# For fitting

source("./khoury_lookup.R")

# Give names to elements of x from priors
named_params = function(priors, constants)
{
    param_names = c(names(priors), names(constants));
    param_defs = c(
        paste0("x[", seq_along(priors) - 1, "]"),
        unlist(as.character(constants))
    )
    paste(mapply(function(nm, def) paste0("double x_", nm, " = ", def, "; (void) x_", nm, ";"), param_names, param_defs), collapse = "\n");
}

# initialiser list with values from x
cpp_vec = function(x) 
{
    paste("{", paste(x, collapse = ", "), "}")
}

asc = function(x, y0, y1, s0, s1)
{
    xx = s0 + x * (s1 - s0);
    h0 = exp(s0) / (1 + exp(s0));
    h1 = exp(s1) / (1 + exp(s1));
    h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);
    return (y0 + (y1 - y0) * h);
}

# Wrap Cpp funcs for multi region fitting
wrap_region = function(model_v2, name, region_names)
{
    ret = 'if (false) {}';
    for (i in seq_along(region_names))
    {
        ret = c(ret,
            glue::glue('else if (P.pop[0].name == "{region_names[i]}")'),
            '{',
            model_v2[[i]][[name]],
            '}'
        )
    }
    return (ret);
}

named_schedules = function()
{
    "enum { CH_CONTACT, CH_TIER_2, CH_TIER_3, CH_CONTACT_CHANGE, CH_XMAS_FUDGE, CH_ADJUST, CH_BA2 };"
}

# create cpp changes
cpp_chgI_voc = function(priors, constants, seasonality, v2, v2_relu, v2_latdur, v2_serial, v2_infdur, v2_immesc, v2_ch_u, v3_relu, v3_severity = 1, v3)
{
    glue::glue(
        named_params(priors, constants),
        named_schedules(),
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'auto interp = [&](double y, vector<double>& curve) {',
        '    if (y < 0) return curve[0];',
        '    unsigned int i = (unsigned int)(y * 100);',
        '    if (i >= curve.size() - 1) return curve.back();',
        '    double f = y * 100 - i;',
        '    return (1 - f) * curve[i] + f * curve[i + 1];',
        '};',
        
        'auto odds = [&](double v, double lo) {',
        '    double a = v / (1 - v);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',
        
        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',
        
        # Do contact adjustment
        # # contact multiplier for contact adjustment
        # adjust_days = as.numeric(ymd(date_fitting) - ymd("2020-01-01"))
        # paramsI$schedule[[6]] = list(
        #     parameter = "contact",
        #     pops = 0,
        #     mode = "multiply",
        #     values = rep(list(rep(1, 8)), adjust_days),
        #     times = 0:(adjust_days - 1)
        # )
        
        # Contact adjustment
        'vector<double> adjustments = { 1.0, 1.0, x_f102, x_f144, x_f186, x_f228, x_f270, x_f312, x_f354, x_f396, x_f438, x_f480, x_f522, x_f564, x_f606, x_f648, x_f690, x_f732, x_f774, x_f816 };',
        'for (unsigned int t = 60; t < P.changes.ch[CH_ADJUST].times.size(); ++t) {',
        '    unsigned int index_left = (t - 39) / 42 + 0;',
        '    unsigned int index_right = (t - 39) / 42 + 1;',
        '    if (index_right >= adjustments.size()) {',
        '        P.changes.ch[CH_ADJUST].values[t] = vector<double>(8, adjustments.back());',
        '    } else {',
        '        double adj1 = adjustments[index_left];',
        '        double adj2 = adjustments[index_right];',
        '        double tx = double((t - 39) % 42) / 42.0;',
        '        P.changes.ch[CH_ADJUST].values[t] = vector<double>(8, asc(tx, adj1, adj2, -10.0, 10.0));',
        '    }',
        '}',
        
        # Set age-specific length of stay
        'P.pop[0].lHosp2.assign(16, Discrete());',
        'P.pop[0].lHosp2[0] = delay_lnorm(5.11, 1.33, 60, 0.25);',
        'P.pop[0].lHosp2[1] = delay_lnorm(5.11, 1.33, 60, 0.25);',
        'P.pop[0].lHosp2[2] = delay_lnorm(5.11, 1.33, 60, 0.25);',
        'P.pop[0].lHosp2[3] = delay_lnorm(5.11, 1.33, 60, 0.25);',
        'P.pop[0].lHosp2[4] = delay_lnorm(5.12, 1.30, 60, 0.25);',
        'P.pop[0].lHosp2[5] = delay_lnorm(5.48, 1.34, 60, 0.25);',
        'P.pop[0].lHosp2[6] = delay_lnorm(6.46, 1.38, 60, 0.25);',
        'P.pop[0].lHosp2[7] = delay_lnorm(7.41, 1.28, 60, 0.25);',
        'P.pop[0].lHosp2[8] = delay_lnorm(8.36, 1.33, 60, 0.25);',
        'P.pop[0].lHosp2[9] = delay_lnorm(9.16, 1.27, 60, 0.25);',
        'P.pop[0].lHosp2[10] = delay_lnorm(9.82, 1.21, 60, 0.25);',
        'P.pop[0].lHosp2[11] = delay_lnorm(11.01, 1.24, 60, 0.25);',
        'P.pop[0].lHosp2[12] = delay_lnorm(11.32, 1.19, 60, 0.25);',
        'P.pop[0].lHosp2[13] = delay_lnorm(12.00, 1.19, 60, 0.25);',
        'P.pop[0].lHosp2[14] = delay_lnorm(12.09, 1.09, 60, 0.25);',
        'P.pop[0].lHosp2[15] = delay_lnorm(11.88, 1.08, 60, 0.25);',
        
        if (seasonality) {
        '    P.pop[0].season_A[0] = x_seasonality;'
        '    P.pop[0].season_T[0] = 365.25;'
        '    P.pop[0].season_phi[0] = 0;'
        } else '',

        'for (unsigned int g = 0; g < P.pop[0].ifr1.size(); ++g) {',
        '    // To death',
        '    P.pop[0].ifr1[g] = odds(P.pop[0].ifr1[g], x_cfr_rlo);',
        '    P.pop[0].ifr2[g] = odds(P.pop[0].ifr2[g], x_cfr_rlo${ if (v2) "+ x_v2_cfr_rlo" else "" });',
        '    P.pop[0].ifr3[g] = min(1.0, odds(P.pop[0].ifr3[g], x_cfr_rlo${ if (v2) "+ x_v2_cfr_rlo" else "" }) * ${ v3_severity });',
        '    // To hospital',
        '    P.pop[0].ihr1[g] = odds(P.pop[0].ihr1[g], x_hosp_rlo);',
        '    P.pop[0].ihr2[g] = odds(P.pop[0].ihr2[g], x_hosp_rlo${ if (v2) "+ x_v2_hosp_rlo" else "" });',
        '    P.pop[0].ihr3[g] = min(1.0, odds(P.pop[0].ihr3[g], x_hosp_rlo${ if (v2) "+ x_v2_hosp_rlo" else "" }) * ${ v3_severity });',
        '    // To ICU',
        '    P.pop[0].iir1[g] = odds(P.pop[0].iir1[g], x_icu_rlo);',
        '    P.pop[0].iir2[g] = odds(P.pop[0].iir2[g], x_icu_rlo${ if (v2) "+ x_v2_icu_rlo" else "" });',
        '    P.pop[0].iir3[g] = min(1.0, odds(P.pop[0].iir3[g], x_icu_rlo${ if (v2) "+ x_v2_icu_rlo" else "" }) * ${ v3_severity });',
        '    // Relative susceptibility',
        '    P.pop[0].u[g]  = P.pop[0].u[g]  * x_u;',
        '    P.pop[0].u3[g] = P.pop[0].u3[g] * x_u * x_v3_relu * x_v2_relu;',
        if (v2_relu) {
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x_u * x_v2_relu;'
        } else {
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x_u;'
        },
        if (v2_immesc) {
        '    P.pop[0].pi2_r[g] = x_v2_immesc;'
        } else '',
        '}',

        if (v2_ch_u) {
        'P.pop[0].u2[0] = P.pop[0].u2[0] * x_v2_ch_u;
        P.pop[0].u2[1] = P.pop[0].u2[1] * x_v2_ch_u;
        P.pop[0].u2[2] = P.pop[0].u2[2] * x_v2_ch_u;
        P.pop[0].u2[3] = P.pop[0].u2[3] * x_v2_ch_u;'
        } else '',

        if (v2_latdur) {
        'P.pop[0].dE2 = delay_gamma(x_v2_latdur * 2.5, 2.5, 15, 0.25);'
        } else '',

        if (v2_serial) {
        'P.pop[0].dE2 = delay_gamma(x_v2_serial * 2.5, 2.5, 15, 0.25);
        P.pop[0].dIp2 = delay_gamma(x_v2_serial * 2.5, 4.0, 30, 0.25);
        P.pop[0].dIs2 = delay_gamma(x_v2_serial * 2.5, 4.0, 30, 0.25);
        P.pop[0].dIa2 = delay_gamma(x_v2_serial * 5.0, 4.0, 30, 0.25);
        for (unsigned int g = 0; g < P.pop[0].u2.size(); ++g) {
            P.pop[0].u2[g] /= x_v2_serial;
        }'
        } else '',

        if (v2_infdur) {
        'P.pop[0].dIp2 = delay_gamma(x_v2_infdur * 2.5, 4.0, 30, 0.25);
        P.pop[0].dIs2 = delay_gamma(x_v2_infdur * 2.5, 4.0, 30, 0.25);
        P.pop[0].dIa2 = delay_gamma(x_v2_infdur * 5.0, 4.0, 30, 0.25);'
        } else '',
        
        '// Delays to death, hospital, and ICU',
        'std::vector<double> dE = delay_gamma(2.5, 2.5, 60, 0.25);',
        'P.pop[0].dDeath = delay_convolve(dE, delay_gamma(x_death_mean, 2.2, 60, 0.25));',
        'P.pop[0].dHosp  = delay_convolve(dE, delay_gamma(x_hosp_admission, 0.71, 60, 0.25));',
        'P.pop[0].dICU   = delay_convolve(dE, delay_gamma(x_icu_admission, 1.91, 60, 0.25));',

        '// Seeding of original and variant strain',
        'P.pop[0].seed_times = seq((int)x_tS, (int)x_tS + 27);',
        if (v2) {
        'P.pop[0].seed_times2 = vector<double>(10, x_v2_when);'
        } else {
        'P.pop[0].seed_times2 = vector<double>(1, 99999);'
        },
        if (v3) {
            'P.pop[0].seed_times3 = vector<double>(10, x_v3_when);'
        } else {
            'P.pop[0].seed_times3 = vector<double>(10, 99999);'
        },
        
        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[CH_CONTACT_CHANGE].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[CH_CONTACT_CHANGE].times.size() - 1.0);',
        '    P.changes.ch[CH_CONTACT_CHANGE].values[i] = vector<double>(8, asc(tx, 1.0, x_contact_final, -x_contact_s0, x_contact_s1));',
        '}',

        # Xmas fudge
        'P.changes.ch[CH_XMAS_FUDGE].values[0] = vector<double>(8, x_xmas_fudge);',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { CH_CONTACT, CH_TIER_2, CH_TIER_3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        double scho = P.changes.ch[k].values[i][5];',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '}',
        .sep = "\n    ", .open = "${", .close = "}"
    )
}

# changes func for extra scenario parameters
cpp_chgI_xparams = function(X_PARAMS)
{
    if (!is.null(X_PARAMS$SDR) & is.null(X_PARAMS$SDD)) {
        stop("Need to specify both SDR and SDD, or neither")
    }
    
    glue::glue(
        if (!is.null(X_PARAMS$VL)) {
        'P.vax_limit = ${ X_PARAMS$VL };'
        } else '',
        
        if (!is.null(X_PARAMS$SDR)) {
        'P.changes.ch[CH_CONTACT_CHANGE].times.resize(731, 0);
        P.changes.ch[CH_CONTACT_CHANGE].values.resize(731, vector<double>());
        
        for (unsigned int i = 366; i < 731; ++i) {
            P.changes.ch[CH_CONTACT_CHANGE].times[i] = i;
            if (i < ${ X_PARAMS$SDD })
                P.changes.ch[CH_CONTACT_CHANGE].values[i] = vector<double>(8, x_contact_final);
            else if (i < ${ X_PARAMS$SDD } + 14)
                P.changes.ch[CH_CONTACT_CHANGE].values[i] = vector<double>(8, x_contact_final + (double(i - ${X_PARAMS$SDD}) / 14.0) * (1.0 - x_contact_final) * ${X_PARAMS$SDR});
            else
                P.changes.ch[CH_CONTACT_CHANGE].values[i] = vector<double>(8, x_contact_final + (1.0 - x_contact_final) * ${X_PARAMS$SDR});
        }'
        } else '',

        # # Contact adjustment
        # 'for (unsigned int i = 0; i < P.changes.ch[CH_CONTACT_CHANGE].times.size(); ++i) {',
        # '    double tx = double(i) / (P.changes.ch[CH_CONTACT_CHANGE].times.size() - 1.0);',
        # '    P.changes.ch[CH_CONTACT_CHANGE].values[i] = vector<double>(8, asc(tx, 1.0, x_contact_final, -x_contact_s0, x_contact_s1));',
        # '}',
        .sep = "\n    ", .open = "${", .close = "}"
    )
}

# create c++ likelihood components - new version for fitting to omicron
cpp_likI_voc_omi = function(params, ld, sitreps, sero, virus, sgtfd, popid, max_date, priors, constants, death_cutoff, use_sgtf, delta, omi, gamdisp)
{
    refdate = ld[, max(date)] + 1;
    ld = copy(ld[date <= max_date]);
    sitreps = copy(sitreps[date <= max_date]);
    sero = copy(sero[End.date <= max_date]);
    virus = copy(virus[End.date <= max_date]);
    sgtfd = copy(sgtfd[date <= max_date]);
    if (!is.null(delta)) {
        delta = copy(delta[date <= max_date]);
    }
    omi = copy(omi[date <= max_date]);
    
    glue::glue(
        named_params(priors, constants),
        named_schedules(),
        'auto restrict = [=](std::vector<double>& t, std::vector<double>& data) {',
        '    size_t start = std::lower_bound(t.begin(), t.end(), t_min) - t.begin();',
        '    size_t end   = std::upper_bound(t.begin(), t.end(), t_max) - t.begin();',
        '    data.assign(data.begin() + start, data.begin() + end);',
        '};',
        '',
        'std::vector<double> death_v =  ${ cpp_vec(ld[!is.na(N), N]) };',
        'std::vector<double> death_t =  ${ cpp_vec(ld[!is.na(N), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> hosp_i_v = ${ cpp_vec(sitreps[!is.na(n_admitted_diagnosed), n_admitted_diagnosed]) };',
        'std::vector<double> hosp_i_t = ${ cpp_vec(sitreps[!is.na(n_admitted_diagnosed), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> hosp_p_v = ${ cpp_vec(sitreps[!is.na(n_in_all_beds), n_in_all_beds]) };',
        'std::vector<double> hosp_p_t = ${ cpp_vec(sitreps[!is.na(n_in_all_beds), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> icu_p_v =  ${ cpp_vec(sitreps[!is.na(n_in_itu), n_in_itu]) };',
        'std::vector<double> icu_p_t =  ${ cpp_vec(sitreps[!is.na(n_in_itu), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> sero_n =   ${ cpp_vec(sero[, n_approx]) };',
        'std::vector<double> sero_p =   ${ cpp_vec(sero[, p]) };',
        'std::vector<double> sero_t0 =   ${ cpp_vec(sero[, as.numeric(ymd(Start.date) - ymd(params$date0))]) };',
        'std::vector<double> sero_t1 =   ${ cpp_vec(sero[, as.numeric(ymd(End.date) - ymd(params$date0))]) };',
        'std::vector<double> virus_n =  ${ cpp_vec(virus[, n_approx]) };',
        'std::vector<double> virus_p =  ${ cpp_vec(virus[, p]) };',
        'std::vector<double> virus_t0 =  ${ cpp_vec(virus[, as.numeric(ymd(Start.date) - ymd(params$date0))]) };',
        'std::vector<double> virus_t1 =  ${ cpp_vec(virus[, as.numeric(ymd(End.date) - ymd(params$date0))]) };',
        'restrict(death_t, death_v);',
        'restrict(death_t, death_t);',
        'restrict(hosp_i_t, hosp_i_v);',
        'restrict(hosp_i_t, hosp_i_t);',
        'restrict(hosp_p_t, hosp_p_v);',
        'restrict(hosp_p_t, hosp_p_t);',
        'restrict(icu_p_t, icu_p_v);',
        'restrict(icu_p_t, icu_p_t);',
        'restrict(sero_t1, sero_n);',
        'restrict(sero_t1, sero_p);',
        'restrict(sero_t1, sero_t0);',
        'restrict(sero_t1, sero_t1);',
        'restrict(virus_t1, virus_n);',
        'restrict(virus_t1, virus_p);',
        'restrict(virus_t1, virus_t0);',
        'restrict(virus_t1, virus_t1);',
        if (use_sgtf) {
            'std::vector<double> sgtf_s = ${ cpp_vec(sgtfd[!is.na(sgtf), sgtf]) };
        std::vector<double> sgtf_o = ${ cpp_vec(sgtfd[!is.na(other), other]) };
        std::vector<double> sgtf_t = ${ cpp_vec(sgtfd[!is.na(sgtf), as.numeric(ymd(date) - ymd(params$date0))]) };
        restrict(sgtf_t, sgtf_s);
        restrict(sgtf_t, sgtf_o);
        restrict(sgtf_t, sgtf_t);'
        } else '',
        if (!is.null(delta)) {
            'std::vector<double> delta_d = ${ cpp_vec(delta[!is.na(delta), delta]) };
        std::vector<double> delta_o = ${ cpp_vec(delta[!is.na(delta), other]) };
        std::vector<double> delta_t = ${ cpp_vec(delta[!is.na(delta), as.numeric(ymd(date) - ymd(params$date0))]) };
        restrict(delta_t, delta_d);
        restrict(delta_t, delta_o);
        restrict(delta_t, delta_t);'
        } else '',
        if (!is.null(omi)) {
            'std::vector<double> omi_s = ${ cpp_vec(omi[!is.na(sgtf), sgtf]) };
        std::vector<double> omi_o = ${ cpp_vec(omi[!is.na(other), other]) };
        std::vector<double> omi_t = ${ cpp_vec(omi[!is.na(sgtf), as.numeric(ymd(date) - ymd(params$date0))]) };
        restrict(omi_t, omi_s);
        restrict(omi_t, omi_o);
        restrict(omi_t, omi_t);'
        } else '',
        'auto size_param = [](double disp) { return 1.0 / (disp * disp); };',
        
        'll += rate_penalty(death_t,  death_v,  dyn, "deaths",            0,                   0.1, 7, ${gamdisp});',
        'll += rate_penalty(hosp_i_t, hosp_i_v, dyn, "hosp_undetected_o", 0,                   0.1, 7, ${gamdisp});',
        'll += rate_penalty(hosp_p_t, hosp_p_v, dyn, "hosp_bed",          "hosp_undetected_p", 0.1, 7, ${gamdisp});',
        'll += rate_penalty(icu_p_t,  icu_p_v,  dyn, "icu_bed",           0,                   0.1, 7, ${gamdisp});',
        'for (unsigned int i = 0; i < virus_t0.size(); ++i) {',
        '    double virus_prev = 0;',
        '    for (double tt = virus_t0[i]; tt <= virus_t1[i]; ++tt)',
        '        virus_prev += dyn("pcr_positive_p", tt, {}, {}) / (virus_t1[i] - virus_t0[i] + 1);',
        '    ll += bbinom(virus_p[i] * virus_n[i], virus_n[i], virus_prev / ${ sum(params$pop[[1]]$size) }, 5000);', # was 5000
        '    //ll += binom(virus_p[i] * virus_n[i], virus_n[i], virus_prev / ${ sum(params$pop[[1]]$size) });',
        '}',
        'for (unsigned int i = 0; i < sero_t0.size(); ++i) {',
        '    double sero_prev = 0;',
        '    for (double tt = sero_t0[i]; tt <= sero_t1[i]; ++tt)',
        '        sero_prev += (',
        '            0.8 * dyn("lfia_positive_p", tt, {}, {3}) + ',
        '            0.8 * dyn("vaccsero_a_p", tt, {}, {3}) + ',
        '            0.8 * dyn("vaccsero_b_p", tt, {}, {3}) + ',
        '            dyn("lfia_positive_p", tt, {}, {4,5,6,7,8,9,10,11,12,13,14,15}) + ',
        '            dyn("vaccsero_a_p", tt, {}, {4,5,6,7,8,9,10,11,12,13,14,15}) + ',
        '            dyn("vaccsero_b_p", tt, {}, {4,5,6,7,8,9,10,11,12,13,14,15})',
        '                ) / (sero_t1[i] - sero_t0[i] + 1);',
        '    ll += bbinom(sero_p[i] * sero_n[i], sero_n[i], sero_prev / ${ 0.8 * params$pop[[1]]$size[4] + sum(params$pop[[1]]$size[5:16]) }, 30);',
        '}',
        if (use_sgtf) {
            'for (unsigned int i = 0; i < sgtf_t.size(); ++i) {
                double s1 = dyn("test_o",  sgtf_t[i], {}, {});
                double s2 = dyn("test2_o", sgtf_t[i], {}, {});
                double p2 = (s1 + s2) > 0 ? s2 / (s1 + s2) : 0;
                double model_sgtf = p2 + (1 - p2) * x_v2_sgtf0;
                ll += bbinom(sgtf_s[i], sgtf_s[i] + sgtf_o[i], model_sgtf, size_param(x_v2_disp));
            }'
        } else '',
        
        if (!is.null(delta)) {
            'for (unsigned int i = 0; i < delta_t.size(); ++i) {
                double d1 = dyn("test3_o",  delta_t[i], {}, {});
                double d0 = dyn("test_o", delta_t[i], {}, {}) + dyn("test2_o", delta_t[i], {}, {});
                double model_delta = (d1 + d0) > 0 ? d1 / (d1 + d0) : 0;
                ll += bbinom(delta_d[i], delta_d[i] + delta_o[i], model_delta, 30.);
            }'
        } else '',
        if (!is.null(omi)) {
            'for (unsigned int i = 0; i < omi_t.size(); ++i) {
                double o1  = dyn("test_o",  omi_t[i], {}, {});
                double o2  = dyn("test2_o", omi_t[i], {}, {});
                double o3  = dyn("test3_o", omi_t[i], {}, {});
                double po2 = (o1 + o2 + o3) > 0 ? o1 / (o1 + o2 + o3) : 0;
                double model_omi = po2 + (1 - po2) * x_v4_sgtf0;
                ll += bbinom(omi_s[i], omi_s[i] + omi_o[i], model_omi, size_param(x_v4_disp));
            }'
        } else '',
        
        .sep = "\n", .open = "${", .close = "}"
        )
}

# create c++ observer function
cpp_obsI_voc = function(concentration, v2, P.death, P.critical, priors, constants, v3_severity = 1)
{
    glue::glue(
        named_params(priors, constants),
        named_schedules(),
        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',
        
        # 'auto odds = [&](double v, double lo) {',
        # '    double a = v / (1 - v);',
        # '    return a * exp(lo) / (a * exp(lo) + 1);',
        # '};',
        
        # 'auto clamp = [&](double v) {',
        # '    return max(0.0, min(1.0, v));',
        # '};',

        'dyn.Obs(t, 0, 0, 0) = estimate_Rt(P, dyn, t, 0, 50);',
        'dyn.Obs(t, 0, 3, 0) = estimate_R0(P, dyn, t, 0, 50);',
        
        # increase in young person mobility
        if (concentration) {
        'if (t == 182) {
            double mode = 0.2;
            double conc = x_concentration1;
            double constant = 0.2;
            for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {
                P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;
                P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;
                P.pop[0].u3[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;
            }
        }
        if (t == 213) {
            double mode = 0.2;
            double conc_prev = x_concentration1;
            double conc = x_concentration2;
            double constant = 0.2;
            for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {
                P.pop[0].u[a]  /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;
                P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;
                P.pop[0].u2[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;
                P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;
                P.pop[0].u3[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;
                P.pop[0].u3[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;
            }
        }
        if (t == 244) {
            double mode = 0.2;
            double conc_prev = x_concentration2;
            double conc = x_concentration3;
            double constant = 0.2;
            for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {
                P.pop[0].u[a]  /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;
                P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;
                P.pop[0].u2[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;
                P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;
                P.pop[0].u3[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;
                P.pop[0].u3[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;
            }
        }'
        } else '',

        'if (t == 419) {',
        '    P.pop[0].dVa1 = delay_gamma(57.8, 1000, 70.0, P.time_step);',
        '    P.pop[0].dVb1 = delay_gamma(57.8, 1000, 70.0, P.time_step);',
        '}',

        # changing CFR and ICU admission
        # 'if ((int)t % 7 == 0) {',
        # '    double ifr_table[] = ${ cpp_vec(P.death) };',
        # '    double icr_table[] = ${ cpp_vec(P.critical) };',
        # '    double adj_f1 = asc(clamp(t / 190.0),          1.0, 0.0, -4.0, 1.0);',
        # '    double adj_f3 = asc(clamp((t - 240.0) / 90.0), 0.0, 1.0, -5.0, 5.0);',
        # '    double adj_f2 = 1.0 - adj_f1 - adj_f3;',
        # '    double adj_c = asc(clamp(t/366.0), 0.0, 1.0, -6, 6);',
        # '    for (unsigned int g = 0; g < P.pop[0].ifr1.size(); ++g) {',
        # '        double ifr = ifr_table[g];',
        # '        double icr = icr_table[g];',
        # '        // To death',
        # '        P.pop[0].ifr1[g] = odds(ifr, adj_f1 * x_cfr_rlo + adj_f2 * x_cfr_rlo2 + adj_f3 * x_cfr_rlo3);',
        # '        P.pop[0].ifr2[g] = odds(ifr, adj_f1 * x_cfr_rlo + adj_f2 * x_cfr_rlo2 + adj_f3 * x_cfr_rlo3${ if (v2) "+ x_v2_cfr_rlo" else "" });',
        # '        P.pop[0].ifr3[g] = min(1.0, odds(ifr, adj_f1 * x_cfr_rlo + adj_f2 * x_cfr_rlo2 + adj_f3 * x_cfr_rlo3${ if (v2) "+ x_v2_cfr_rlo" else "" }) * ${ v3_severity });',
        # 
        # '        // To ICU',
        # '        P.pop[0].iir1[g] = odds(icr, (1 - adj_c) * x_icu_rlo + adj_c * x_icu_rlo2);',
        # '        P.pop[0].iir2[g] = odds(icr, (1 - adj_c) * x_icu_rlo + adj_c * x_icu_rlo2${ if (v2) "+ x_v2_icu_rlo" else "" });',
        # '        P.pop[0].iir3[g] = min(1.0, odds(icr, (1 - adj_c) * x_icu_rlo + adj_c * x_icu_rlo2${ if (v2) "+ x_v2_icu_rlo" else "" }) * ${ v3_severity });',
        # '    }',
        # '}',
        # changing detection rate in patients admitted to hospital
        'double detection = asc(min(t / 365.0, 1.0), 14, 1, -5.86, 33.4);',
        'P.processes[9].delays[0] = delay_gamma(detection, 0.59, 60, 0.25);',
        .sep = "\n    ", .open = "${", .close = "}"
    )
}

# Observer for vaccination programmes
cpp_obsI_vax = function(params, vacc)
{
    glue::glue(
        'if (t == 0) { dyn.scratch["vaxphase"] = 0; }',
        '{',
        'vector<double> vt  = ${ cpp_vec(as.numeric(ymd(vacc$vt) - ymd(params$date0))) };',
        'vector<double> va1   = ${ cpp_vec(unlist(vacc$va1)) };',
        'vector<double> vb1   = ${ cpp_vec(unlist(vacc$vb1)) };',
        'unsigned int phase = dyn.scratch["vaxphase"];',
        'if (vt.size() > phase && t >= vt[phase]) {',
        '    int age_groups = P.pop[0].size.size();',
        '    P.pop[0].va1   = vector<double>(va1.begin()   + age_groups * phase, va1.begin()   + age_groups * (phase + 1));',
        '    P.pop[0].vb1   = vector<double>(vb1.begin()   + age_groups * phase, vb1.begin()   + age_groups * (phase + 1));',
        '    dyn.scratch["vaxphase"] = phase + 1;',
        '}',
        'dyn.Obs(t, 0, 5, 0) = dyn.scratch["vaxphase"];',
        '}',
        .sep = "\n", .open = "${", .close = "}")
}

# Observer for forced seasonality
cpp_obsI_seasonality = function(forced_seasonality, seas_start_t)
{
    glue::glue(
        'if (t == ${seas_start_t}) {',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  /= ${1 + forced_seasonality};',
        '        P.pop[0].u2[a] /= ${1 + forced_seasonality};',
        '        P.pop[0].u3[a] /= ${1 + forced_seasonality};',
        '    }',
        '    P.pop[0].season_A[0] = ${forced_seasonality};',
        '    P.pop[0].season_T[0] = 365.25;',
        '    P.pop[0].season_phi[0] = 0;',
        '}',
        .sep = "\n", .open = "${", .close = "}")
}

# Observer for extra relaxation of R
cpp_obsI_relax = function(relax_t, factor)
{
    glue::glue(
        'if (t == ${relax_t}) {',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  *= ${factor};',
        '        P.pop[0].u2[a] *= ${factor};',
        '        P.pop[0].u3[a] *= ${factor};',
        '    }',
        '}',
        .sep = "\n", .open = "${", .close = "}")
}

# Observer for sensitivity analyses for autumn/winter scenarios
cpp_obsI_aw = function(seasonality_aw, cert, planB_start, planB_end)
{
    glue::glue(
        'if (t == 639) { // October 1st, 2021', 
        '    P.pop[0].season_A[0] = ${seasonality_aw};',
        '    P.pop[0].season_T[0] = 365.25;',
        '    P.pop[0].season_phi[0] = 0;',
        '}',
        'if (${cert} != 0) {',
        '    if (t == ${ifelse(is.na(planB_start), 9999999, as.numeric(ymd(planB_start) - ymd("2020-01-01")))} - 14) {',
        '        double u18_19 = ${cert} * P.pop[0].size[3] * 0.4 / 21.0;',
        '        double u20_24 = ${cert} * P.pop[0].size[4] / 21.0;',
        '        double u25_29 = ${cert} * P.pop[0].size[5] / 21.0;',
        '        P.extra_vb1[3] = u18_19;',
        '        P.extra_vb1[4] = u20_24;',
        '        P.extra_vb1[5] = u25_29;',
        '    }',
        '    if (t == ${ifelse(is.na(planB_start), 9999999, as.numeric(ymd(planB_start) - ymd("2020-01-01")))} + 7) {',
        '        P.extra_vb1[3] = 0;',
        '        P.extra_vb1[4] = 0;',
        '        P.extra_vb1[5] = 0;',
        '    }',
        '}',
        .sep = "\n", .open = "${", .close = "}")

}

# create c++ observer function for fitting to Omicron data
cpp_obsI_voc_omi = function(setup_t, x_protection, vax_factor, vax_assumption, n_seeds_per_day, omi_sev, omi_crit)
{
    if (vax_assumption == "khoury") {
        delta_ei = list()
        
        delta_ei[["az1"]] = 43
        delta_ei[["az2"]] = 63
        delta_ei[["azW"]] = 32
        delta_ei[["pf1"]] = 62
        delta_ei[["pf2"]] = 80
        delta_ei[["pfW"]] = 43
        delta_ei[["nat"]] = 100
        
        omicron_esev = list()
        
        for (imm in c("az1", "az2", "azW", "pf1", "pf2", "pfW")) {
            omicron_esev[[imm]] = 0.01 * khoury_get_3a(delta_ei[[imm]] * vax_factor)$mid
        }
        for (imm in c("nat")) {
            omicron_esev[[imm]] = 0.01 * khoury_get_3a(delta_ei[[imm]] * x_protection)$mid
        }
        
    } else if (vax_assumption != "conditional") {
        stop("vax_assumption must be 'conditional' or 'khoury'")
    }
    
    glue::glue(
        'if (t == ${setup_t}) {',
            'auto VEfunc = [](double x, double given) { return (x - given) / (1.0 - given); };',
            'for (unsigned int day = x_v4_when; day <= x_v4_when + 13; day++){',
                'P.pop[0].seed_times.insert(P.pop[0].seed_times.end(), 10, day);',
            '}',
            'for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
                # Move individuals between compartments
                'Mp.pops[0].R3[a] += Mp.pops[0].R[a];',
                'Mp.pops[0].R[a] = 0;',
                
                # Set u to zero in prep for later emergence
                'P.pop[0].u[a] = 0;',
                
                # cross protection
                # 'P.pop[0].pi_r[a] = ${x_protection};',
                'P.pop[0].pi_r2[a] = ${x_protection};',
                'P.pop[0].pi_r3[a] = ${x_protection};',
                'P.pop[0].pd_r2i[a] = P.pop[0].ed_vb2i3[a];',
                'P.pop[0].pd_r3i[a] = P.pop[0].ed_vb2i3[a];',
                'P.pop[0].pt_r2i[a] = P.pop[0].et_vb2i3[a];',
                'P.pop[0].pt_r3i[a] = P.pop[0].et_vb2i3[a];',
                
                # Reset protection against Delta given Omicron immunity to 1 (was protection against Delta given WT immunity prior to setup_t)
                'P.pop[0].pi3_r[a] = 1.0;',
                
                'if (${vax_assumption == "conditional"}) {',
                '    P.pop[0].ph_r2d[a] = P.pop[0].eh_vb2d3[a];',
                '    P.pop[0].ph_r3d[a] = P.pop[0].eh_vb2d3[a];',
                '    P.pop[0].pm_r2d[a] = P.pop[0].em_vb2d3[a];',
                '    P.pop[0].pm_r3d[a] = P.pop[0].em_vb2d3[a];',
                '} else {',
                '    P.pop[0].ph_r2d[a] = VEfunc(${ omicron_esev[["nat"]] }, 1.0 - (1.0 - P.pop[0].pi_r2[a]) * (1.0 - P.pop[0].pd_r2i[a]));',
                '    P.pop[0].ph_r3d[a] = VEfunc(${ omicron_esev[["nat"]] }, 1.0 - (1.0 - P.pop[0].pi_r3[a]) * (1.0 - P.pop[0].pd_r3i[a]));',
                '    P.pop[0].pm_r2d[a] = VEfunc(${ omicron_esev[["nat"]] }, 1.0 - (1.0 - P.pop[0].pi_r2[a]) * (1.0 - P.pop[0].pd_r2i[a]));',
                '    P.pop[0].pm_r3d[a] = VEfunc(${ omicron_esev[["nat"]] }, 1.0 - (1.0 - P.pop[0].pi_r3[a]) * (1.0 - P.pop[0].pd_r3i[a]));',
                '}',
        
            # do we want to change any of the below cross-protection params?
            # pi2_r, 
            # # pi2_r2, pi2_r3, 
            # pi3_r, 
            # # pi3_r2, pi3_r3
            
            # vaccine protection
            'P.pop[0].ei_va1[a]  = P.pop[0].ei3_va1[a] * ${vax_factor};',
            'P.pop[0].ei_va2[a]  = P.pop[0].ei3_va2[a] * ${vax_factor};',
            'P.pop[0].ei_va3[a]  = P.pop[0].ei3_va3[a] * ${vax_factor};',
            'P.pop[0].ei_vb1[a]  = P.pop[0].ei3_vb1[a] * ${vax_factor};',
            'P.pop[0].ei_vb2[a]  = P.pop[0].ei3_vb2[a] * ${vax_factor};',
            'P.pop[0].ei_vb3[a]  = P.pop[0].ei3_vb3[a] * ${vax_factor};',
            
            'P.pop[0].ed_va1i[a] = P.pop[0].ed_va1i3[a];',
            'P.pop[0].ed_va2i[a] = P.pop[0].ed_va2i3[a];',
            'P.pop[0].ed_va3i[a] = P.pop[0].ed_va3i3[a];',
            'P.pop[0].ed_vb1i[a] = P.pop[0].ed_vb1i3[a];',
            'P.pop[0].ed_vb2i[a] = P.pop[0].ed_vb2i3[a];',
            'P.pop[0].ed_vb3i[a] = P.pop[0].ed_vb3i3[a];',
            
            'P.pop[0].et_va1i[a] = P.pop[0].et_va1i3[a];',
            'P.pop[0].et_va2i[a] = P.pop[0].et_va2i3[a];',
            'P.pop[0].et_va3i[a] = P.pop[0].et_va3i3[a];',
            'P.pop[0].et_vb1i[a] = P.pop[0].et_vb1i3[a];',
            'P.pop[0].et_vb2i[a] = P.pop[0].et_vb2i3[a];',
            'P.pop[0].et_vb3i[a] = P.pop[0].et_vb3i3[a];',
            
            'if (${vax_assumption == "conditional"}) {',
            '    P.pop[0].eh_va1d[a] = P.pop[0].eh_va1d3[a];',
            '    P.pop[0].eh_va2d[a] = P.pop[0].eh_va2d3[a];',
            '    P.pop[0].eh_va3d[a] = P.pop[0].eh_va3d3[a];',
            '    P.pop[0].eh_vb1d[a] = P.pop[0].eh_vb1d3[a];',
            '    P.pop[0].eh_vb2d[a] = P.pop[0].eh_vb2d3[a];',
            '    P.pop[0].eh_vb3d[a] = P.pop[0].eh_vb3d3[a];',
            
            '    P.pop[0].em_va1d[a] = P.pop[0].em_va1d3[a];',
            '    P.pop[0].em_va2d[a] = P.pop[0].em_va2d3[a];',
            '    P.pop[0].em_va3d[a] = P.pop[0].em_va3d3[a];',
            '    P.pop[0].em_vb1d[a] = P.pop[0].em_vb1d3[a];',
            '    P.pop[0].em_vb2d[a] = P.pop[0].em_vb2d3[a];',
            '    P.pop[0].em_vb3d[a] = P.pop[0].em_vb3d3[a];',
            '} else {',
            '    P.pop[0].eh_va1d[a] = VEfunc(${ omicron_esev[["az1"]] }, 1.0 - (1.0 - P.pop[0].ei_va1[a]) * (1.0 - P.pop[0].ed_va1i[a]));',
            '    P.pop[0].eh_va2d[a] = VEfunc(${ omicron_esev[["az2"]] }, 1.0 - (1.0 - P.pop[0].ei_va2[a]) * (1.0 - P.pop[0].ed_va2i[a]));',
            '    P.pop[0].eh_va3d[a] = VEfunc(${ omicron_esev[["azW"]] }, 1.0 - (1.0 - P.pop[0].ei_va3[a]) * (1.0 - P.pop[0].ed_va3i[a]));',
            '    P.pop[0].eh_vb1d[a] = VEfunc(${ omicron_esev[["pf1"]] }, 1.0 - (1.0 - P.pop[0].ei_vb1[a]) * (1.0 - P.pop[0].ed_vb1i[a]));',
            '    P.pop[0].eh_vb2d[a] = VEfunc(${ omicron_esev[["pf2"]] }, 1.0 - (1.0 - P.pop[0].ei_vb2[a]) * (1.0 - P.pop[0].ed_vb2i[a]));',
            '    P.pop[0].eh_vb3d[a] = VEfunc(${ omicron_esev[["pfW"]] }, 1.0 - (1.0 - P.pop[0].ei_vb3[a]) * (1.0 - P.pop[0].ed_vb3i[a]));',
            
            '    P.pop[0].em_va1d[a] = VEfunc(${ omicron_esev[["az1"]] }, 1.0 - (1.0 - P.pop[0].ei_va1[a]) * (1.0 - P.pop[0].ed_va1i[a]));',
            '    P.pop[0].em_va2d[a] = VEfunc(${ omicron_esev[["az2"]] }, 1.0 - (1.0 - P.pop[0].ei_va2[a]) * (1.0 - P.pop[0].ed_va2i[a]));',
            '    P.pop[0].em_va3d[a] = VEfunc(${ omicron_esev[["azW"]] }, 1.0 - (1.0 - P.pop[0].ei_va3[a]) * (1.0 - P.pop[0].ed_va3i[a]));',
            '    P.pop[0].em_vb1d[a] = VEfunc(${ omicron_esev[["pf1"]] }, 1.0 - (1.0 - P.pop[0].ei_vb1[a]) * (1.0 - P.pop[0].ed_vb1i[a]));',
            '    P.pop[0].em_vb2d[a] = VEfunc(${ omicron_esev[["pf2"]] }, 1.0 - (1.0 - P.pop[0].ei_vb2[a]) * (1.0 - P.pop[0].ed_vb2i[a]));',
            '    P.pop[0].em_vb3d[a] = VEfunc(${ omicron_esev[["pfW"]] }, 1.0 - (1.0 - P.pop[0].ei_vb3[a]) * (1.0 - P.pop[0].ed_vb3i[a]));',
            '}',
        
        
            # severity...
            'P.pop[0].ifr1[a] = P.pop[0].ifr3[a] * ${omi_crit};',
            'P.pop[0].ihr1[a] = P.pop[0].ihr3[a] * ${omi_sev};',
            'P.pop[0].iir1[a] = P.pop[0].iir3[a] * ${omi_crit};',
        
            # change any delays??
            # "dDeath"          "dHosp"           "lHosp"           "dICU"            "lICU"
        
            '}',
        
        '}',
        
        'if (t == int(x_v4_when)) {',
            'for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
                'P.pop[0].u[a] = P.pop[0].u3[a] * x_v4_relu;',
            '}',
        '}',
        
        .sep = "\n", .open = "${", .close = "}")
}

# create c++ observer function for nu voc, November 2021
cpp_obsI_voc_nu = function(setup_t, nu_voc_t, tx_factor, x_protection, vax_factor, vax_assumption, sev_factor, n_seeds_per_day)
{
    nu_seed_times = rep(0:19 + nu_voc_t, each = n_seeds_per_day)

    if (vax_assumption == "khoury") {
        delta_ei = list()

        delta_ei[["az1"]] = 43
        delta_ei[["az2"]] = 63
        delta_ei[["azW"]] = 36
        delta_ei[["pf1"]] = 62
        delta_ei[["pf2"]] = 80
        delta_ei[["pfW"]] = 45
        delta_ei[["nat"]] = 100

        omicron_esev = list()

        for (imm in c("az1", "az2", "azW", "pf1", "pf2", "pfW")) {
            omicron_esev[[imm]] = 0.01 * khoury_get_3a(delta_ei[[imm]] * vax_factor)$mid
        }
        for (imm in c("nat")) {
            omicron_esev[[imm]] = 0.01 * khoury_get_3a(delta_ei[[imm]] * x_protection)$mid
        }
        
    } else if (vax_assumption != "conditional") {
        stop("vax_assumption must be 'conditional' or 'khoury'")
    }
    
    glue::glue(
        'if (t == ${setup_t}) {',
        'auto VEfunc = [](double x, double given) { return (x - given) / (1.0 - given); };',
        'P.pop[0].seed_times = ${cpp_vec(nu_seed_times)};',
        'for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        
        # Move individuals between compartments
        'Mp.pops[0].R3[a] += Mp.pops[0].R[a];',
        'Mp.pops[0].R[a] = 0;',
        
        # Set u to zero in prep for later emergence
        'P.pop[0].u[a] = 0;',
        
        # cross protection
        # 'P.pop[0].pi_r[a] = ${x_protection};',
        'P.pop[0].pi_r2[a] = ${x_protection};',
        'P.pop[0].pi_r3[a] = ${x_protection};',
        'P.pop[0].pd_r2i[a] = P.pop[0].ed_vb2i3[a];',
        'P.pop[0].pd_r3i[a] = P.pop[0].ed_vb2i3[a];',
        'P.pop[0].pt_r2i[a] = P.pop[0].et_vb2i3[a];',
        'P.pop[0].pt_r3i[a] = P.pop[0].et_vb2i3[a];',

        # Reset protection against Delta given Omicron immunity to 1 (was protection against Delta given WT immunity prior to setup_t)
        'P.pop[0].pi3_r[a] = 1.0;',

        'if (${vax_assumption == "conditional"}) {',
        '    P.pop[0].ph_r2d[a] = P.pop[0].eh_vb2d3[a];',
        '    P.pop[0].ph_r3d[a] = P.pop[0].eh_vb2d3[a];',
        '    P.pop[0].pm_r2d[a] = P.pop[0].em_vb2d3[a];',
        '    P.pop[0].pm_r3d[a] = P.pop[0].em_vb2d3[a];',
        '} else {',
        '    P.pop[0].ph_r2d[a] = VEfunc(${ omicron_esev[["nat"]] }, 1.0 - (1.0 - P.pop[0].pi_r2[a]) * (1.0 - P.pop[0].pd_r2i[a]));',
        '    P.pop[0].ph_r3d[a] = VEfunc(${ omicron_esev[["nat"]] }, 1.0 - (1.0 - P.pop[0].pi_r3[a]) * (1.0 - P.pop[0].pd_r3i[a]));',
        '    P.pop[0].pm_r2d[a] = VEfunc(${ omicron_esev[["nat"]] }, 1.0 - (1.0 - P.pop[0].pi_r2[a]) * (1.0 - P.pop[0].pd_r2i[a]));',
        '    P.pop[0].pm_r3d[a] = VEfunc(${ omicron_esev[["nat"]] }, 1.0 - (1.0 - P.pop[0].pi_r3[a]) * (1.0 - P.pop[0].pd_r3i[a]));',
        '}',
        
        # do we want to change any of the below cross-protection params?
        # pi2_r, 
        # # pi2_r2, pi2_r3, 
        # pi3_r, 
        # # pi3_r2, pi3_r3
        
        # vaccine protection
        'P.pop[0].ei_va1[a]  = P.pop[0].ei3_va1[a] * ${vax_factor};',
        'P.pop[0].ei_va2[a]  = P.pop[0].ei3_va2[a] * ${vax_factor};',
        'P.pop[0].ei_va3[a]  = P.pop[0].ei3_va3[a] * ${vax_factor};',
        'P.pop[0].ei_vb1[a]  = P.pop[0].ei3_vb1[a] * ${vax_factor};',
        'P.pop[0].ei_vb2[a]  = P.pop[0].ei3_vb2[a] * ${vax_factor};',
        'P.pop[0].ei_vb3[a]  = P.pop[0].ei3_vb3[a] * ${vax_factor};',

        'P.pop[0].ed_va1i[a] = P.pop[0].ed_va1i3[a];',
        'P.pop[0].ed_va2i[a] = P.pop[0].ed_va2i3[a];',
        'P.pop[0].ed_va3i[a] = P.pop[0].ed_va3i3[a];',
        'P.pop[0].ed_vb1i[a] = P.pop[0].ed_vb1i3[a];',
        'P.pop[0].ed_vb2i[a] = P.pop[0].ed_vb2i3[a];',
        'P.pop[0].ed_vb3i[a] = P.pop[0].ed_vb3i3[a];',
        
        'P.pop[0].et_va1i[a] = P.pop[0].et_va1i3[a];',
        'P.pop[0].et_va2i[a] = P.pop[0].et_va2i3[a];',
        'P.pop[0].et_va3i[a] = P.pop[0].et_va3i3[a];',
        'P.pop[0].et_vb1i[a] = P.pop[0].et_vb1i3[a];',
        'P.pop[0].et_vb2i[a] = P.pop[0].et_vb2i3[a];',
        'P.pop[0].et_vb3i[a] = P.pop[0].et_vb3i3[a];',
        
        'if (${vax_assumption == "conditional"}) {',
        '    P.pop[0].eh_va1d[a] = P.pop[0].eh_va1d3[a];',
        '    P.pop[0].eh_va2d[a] = P.pop[0].eh_va2d3[a];',
        '    P.pop[0].eh_va3d[a] = P.pop[0].eh_va3d3[a];',
        '    P.pop[0].eh_vb1d[a] = P.pop[0].eh_vb1d3[a];',
        '    P.pop[0].eh_vb2d[a] = P.pop[0].eh_vb2d3[a];',
        '    P.pop[0].eh_vb3d[a] = P.pop[0].eh_vb3d3[a];',
        
        '    P.pop[0].em_va1d[a] = P.pop[0].em_va1d3[a];',
        '    P.pop[0].em_va2d[a] = P.pop[0].em_va2d3[a];',
        '    P.pop[0].em_va3d[a] = P.pop[0].em_va3d3[a];',
        '    P.pop[0].em_vb1d[a] = P.pop[0].em_vb1d3[a];',
        '    P.pop[0].em_vb2d[a] = P.pop[0].em_vb2d3[a];',
        '    P.pop[0].em_vb3d[a] = P.pop[0].em_vb3d3[a];',
        '} else {',
        '    P.pop[0].eh_va1d[a] = VEfunc(${ omicron_esev[["az1"]] }, 1.0 - (1.0 - P.pop[0].ei_va1[a]) * (1.0 - P.pop[0].ed_va1i[a]));',
        '    P.pop[0].eh_va2d[a] = VEfunc(${ omicron_esev[["az2"]] }, 1.0 - (1.0 - P.pop[0].ei_va2[a]) * (1.0 - P.pop[0].ed_va2i[a]));',
        '    P.pop[0].eh_va3d[a] = VEfunc(${ omicron_esev[["azW"]] }, 1.0 - (1.0 - P.pop[0].ei_va3[a]) * (1.0 - P.pop[0].ed_va3i[a]));',
        '    P.pop[0].eh_vb1d[a] = VEfunc(${ omicron_esev[["pf1"]] }, 1.0 - (1.0 - P.pop[0].ei_vb1[a]) * (1.0 - P.pop[0].ed_vb1i[a]));',
        '    P.pop[0].eh_vb2d[a] = VEfunc(${ omicron_esev[["pf2"]] }, 1.0 - (1.0 - P.pop[0].ei_vb2[a]) * (1.0 - P.pop[0].ed_vb2i[a]));',
        '    P.pop[0].eh_vb3d[a] = VEfunc(${ omicron_esev[["pfW"]] }, 1.0 - (1.0 - P.pop[0].ei_vb3[a]) * (1.0 - P.pop[0].ed_vb3i[a]));',
        
        '    P.pop[0].em_va1d[a] = VEfunc(${ omicron_esev[["az1"]] }, 1.0 - (1.0 - P.pop[0].ei_va1[a]) * (1.0 - P.pop[0].ed_va1i[a]));',
        '    P.pop[0].em_va2d[a] = VEfunc(${ omicron_esev[["az2"]] }, 1.0 - (1.0 - P.pop[0].ei_va2[a]) * (1.0 - P.pop[0].ed_va2i[a]));',
        '    P.pop[0].em_va3d[a] = VEfunc(${ omicron_esev[["azW"]] }, 1.0 - (1.0 - P.pop[0].ei_va3[a]) * (1.0 - P.pop[0].ed_va3i[a]));',
        '    P.pop[0].em_vb1d[a] = VEfunc(${ omicron_esev[["pf1"]] }, 1.0 - (1.0 - P.pop[0].ei_vb1[a]) * (1.0 - P.pop[0].ed_vb1i[a]));',
        '    P.pop[0].em_vb2d[a] = VEfunc(${ omicron_esev[["pf2"]] }, 1.0 - (1.0 - P.pop[0].ei_vb2[a]) * (1.0 - P.pop[0].ed_vb2i[a]));',
        '    P.pop[0].em_vb3d[a] = VEfunc(${ omicron_esev[["pfW"]] }, 1.0 - (1.0 - P.pop[0].ei_vb3[a]) * (1.0 - P.pop[0].ed_vb3i[a]));',
        '}',

        
        # severity...
        'P.pop[0].ifr1[a] = P.pop[0].ifr3[a] * ${sev_factor};',
        'P.pop[0].ihr1[a] = P.pop[0].ihr3[a] * ${sev_factor};',
        'P.pop[0].iir1[a] = P.pop[0].iir3[a] * ${sev_factor};',
        
        # change any delays??
         # "dDeath"          "dHosp"           "lHosp"           "dICU"            "lICU"
        
        '}',
        
        '}',
        
        'if (t == ${nu_voc_t}) {',
        'for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        # relative transmissibility
        'P.pop[0].u[a] = P.pop[0].u3[a] * ${tx_factor};',
        '}',
        '}',
        
        .sep = "\n", .open = "${", .close = "}")
}

# create c++ observer function for BA.2
cpp_obsI_voc_ba2 = function(ba2t, ba2f)
{
    glue::glue(
        'vector <int> ba2t = ${cpp_vec(ba2t)};',
        'vector <double> ba2f = ${cpp_vec(ba2f)};',
        'for (unsigned int i = 0; i < ba2t.size(); i++){',
        'int thist = ba2t[i];',
        'if (t == thist) {',
        'for (unsigned int a = 0; a < P.pop[0].u.size(); a++){',
        'P.pop[0].u[a] *= ba2f[i];',
        '}',
        '}',
        '}',
      
        .sep = "\n", .open = "${", .close = "}")
}

cpp_obsI_booster = function(target_phase1 = 229000, target_phase2 = 229000, 
    proportion_booster = rep(0.9, 16), booster_fold = 1, booster_om_fold = 1,
    booster_duration = 180)
{
    bsched = get_booster_schedule(target_phase1, target_phase2, proportion_booster);
    start_days = bsched$boost_start
    end_days = start_days + booster_duration
    earliest_day = min(start_days);
    earliest_age = bsched$group[which.min(start_days)] - 1;

    glue::glue(
        '{',
        'vector<double> start_days = ${cpp_vec(start_days)};',
        'vector<double> end_days = ${cpp_vec(end_days)};',
        'double boost_fold = ${booster_fold};',
        'double boost_om_fold = ${booster_om_fold};',
        'vector<double> proportion_booster = ${cpp_vec(proportion_booster)};',
        
        'auto VEfunc = [](double x, double given) { return (x - given) / (1.0 - given); };',
        
        'auto boost = [=](double& eff_inf, double eff_dis, double& eff_hosp, double& eff_mort, double extra_fold_boost, double frac_boosted = 1.0)',
        '{',
        '    if (eff_inf == 1.0) return;',
        '    double eff_boosted = efficacy_fold(eff_inf, boost_fold * extra_fold_boost);',
        '    double eff_sev = infection_to_severe(eff_boosted);',
        '    eff_inf = (1.0 - frac_boosted) * eff_inf + frac_boosted * eff_boosted;',
        '    eff_hosp = (1.0 - frac_boosted) * eff_hosp + frac_boosted * VEfunc(eff_sev, 1.0 - (1.0 - eff_boosted) * (1.0 - eff_dis));',
        '    eff_mort = (1.0 - frac_boosted) * eff_mort + frac_boosted * VEfunc(eff_sev, 1.0 - (1.0 - eff_boosted) * (1.0 - eff_dis));',
        '};',

        'for (unsigned int a = 0; a < 16; ++a) {',
        '    double R   = dyn(t, 0, a, 5) ;',
        '    double R2  = dyn(t, 0, a, 10);',
        '    double R3  = dyn(t, 0, a, 15);',
        '    double Va1 = dyn(t, 0, a, 16);',
        '    double Va2 = dyn(t, 0, a, 17);',
        '    double Va3 = dyn(t, 0, a, 18);',
        '    double Vb1 = dyn(t, 0, a, 19);',
        '    double Vb2 = dyn(t, 0, a, 20);',
        '    double Vb3 = dyn(t, 0, a, 21);',
        '    double EV  = dyn(t, 0, a, 36);',
        '    double frac_V_in_R = max(0.0, min(1.0, (EV - (Va1 + Va2 + Va3 + Vb1 + Vb2 + Vb3)) / (R + R2 + R3)));',
        '    dyn.Obs(t, 0, a, 1) = frac_V_in_R;',

        '    double frac_V_in_R_old = min(1.0, (EV - (Va1 + Va2 + Va3 + Vb1 + Vb2 + Vb3)) / (R + R2 + R3));',
        '    if (frac_V_in_R_old < 0) {',
        '        std::cout << "a = " << a << "; t = " << t << "; frac_V_in_R_old = " << frac_V_in_R_old << "\\n";',
        '        Rcpp::stop("frac_V_in_R_old < 0");',
        '    }',
        
        
        '    if (t == start_days[a]) {',
        '        if (t == ${earliest_day} && a == ${earliest_age}) {',
        '            // SAVE',
        '            P.pop[0].X_pi_r =       P.pop[0].pi_r;',
        '            P.pop[0].X_pi_r2 =      P.pop[0].pi_r2;',
        '            P.pop[0].X_pi_r3 =      P.pop[0].pi_r3;',
        '            P.pop[0].X_pi2_r =      P.pop[0].pi2_r;',
        '            P.pop[0].X_pi2_r2 =     P.pop[0].pi2_r2;',
        '            P.pop[0].X_pi2_r3 =     P.pop[0].pi2_r3;',
        '            P.pop[0].X_pi3_r =      P.pop[0].pi3_r;',
        '            P.pop[0].X_pi3_r2 =     P.pop[0].pi3_r2;',
        '            P.pop[0].X_pi3_r3 =     P.pop[0].pi3_r3;',
        '            P.pop[0].X_ei_va1 =     P.pop[0].ei_va1;',
        '            P.pop[0].X_ei2_va1 =    P.pop[0].ei2_va1;',
        '            P.pop[0].X_ei3_va1 =    P.pop[0].ei3_va1;',
        '            P.pop[0].X_ei_va2 =     P.pop[0].ei_va2;',
        '            P.pop[0].X_ei2_va2 =    P.pop[0].ei2_va2;',
        '            P.pop[0].X_ei3_va2 =    P.pop[0].ei3_va2;',
        '            P.pop[0].X_ei_va3 =     P.pop[0].ei_va3;',
        '            P.pop[0].X_ei2_va3 =    P.pop[0].ei2_va3;',
        '            P.pop[0].X_ei3_va3 =    P.pop[0].ei3_va3;',
        '            P.pop[0].X_ei_vb1 =     P.pop[0].ei_vb1;',
        '            P.pop[0].X_ei2_vb1 =    P.pop[0].ei2_vb1;',
        '            P.pop[0].X_ei3_vb1 =    P.pop[0].ei3_vb1;',
        '            P.pop[0].X_ei_vb2 =     P.pop[0].ei_vb2;',
        '            P.pop[0].X_ei2_vb2 =    P.pop[0].ei2_vb2;',
        '            P.pop[0].X_ei3_vb2 =    P.pop[0].ei3_vb2;',
        '            P.pop[0].X_ei_vb3 =     P.pop[0].ei_vb3;',
        '            P.pop[0].X_ei2_vb3 =    P.pop[0].ei2_vb3;',
        '            P.pop[0].X_ei3_vb3 =    P.pop[0].ei3_vb3;',
        '            P.pop[0].X_ph_rd =      P.pop[0].ph_rd;',
        '            P.pop[0].X_ph_rd2 =     P.pop[0].ph_rd2;',
        '            P.pop[0].X_ph_rd3 =     P.pop[0].ph_rd3;',
        '            P.pop[0].X_ph_r2d =     P.pop[0].ph_r2d;',
        '            P.pop[0].X_ph_r2d2 =    P.pop[0].ph_r2d2;',
        '            P.pop[0].X_ph_r2d3 =    P.pop[0].ph_r2d3;',
        '            P.pop[0].X_ph_r3d =     P.pop[0].ph_r3d;',
        '            P.pop[0].X_ph_r3d2 =    P.pop[0].ph_r3d2;',
        '            P.pop[0].X_ph_r3d3 =    P.pop[0].ph_r3d3;',
        '            P.pop[0].X_pm_rd =      P.pop[0].pm_rd;',
        '            P.pop[0].X_pm_rd2 =     P.pop[0].pm_rd2;',
        '            P.pop[0].X_pm_rd3 =     P.pop[0].pm_rd3;',
        '            P.pop[0].X_pm_r2d =     P.pop[0].pm_r2d;',
        '            P.pop[0].X_pm_r2d2 =    P.pop[0].pm_r2d2;',
        '            P.pop[0].X_pm_r2d3 =    P.pop[0].pm_r2d3;',
        '            P.pop[0].X_pm_r3d =     P.pop[0].pm_r3d;',
        '            P.pop[0].X_pm_r3d2 =    P.pop[0].pm_r3d2;',
        '            P.pop[0].X_pm_r3d3 =    P.pop[0].pm_r3d3;',
        '            P.pop[0].X_eh_va1d  = P.pop[0].eh_va1d;',
        '            P.pop[0].X_eh_va1d2 = P.pop[0].eh_va1d2;',
        '            P.pop[0].X_eh_va1d3 = P.pop[0].eh_va1d3;',
        '            P.pop[0].X_eh_va2d  = P.pop[0].eh_va2d;',
        '            P.pop[0].X_eh_va2d2 = P.pop[0].eh_va2d2;',
        '            P.pop[0].X_eh_va2d3 = P.pop[0].eh_va2d3;',
        '            P.pop[0].X_eh_va3d  = P.pop[0].eh_va3d;',
        '            P.pop[0].X_eh_va3d2 = P.pop[0].eh_va3d2;',
        '            P.pop[0].X_eh_va3d3 = P.pop[0].eh_va3d3;',
        '            P.pop[0].X_eh_vb1d  = P.pop[0].eh_vb1d;',
        '            P.pop[0].X_eh_vb1d2 = P.pop[0].eh_vb1d2;',
        '            P.pop[0].X_eh_vb1d3 = P.pop[0].eh_vb1d3;',
        '            P.pop[0].X_eh_vb2d  = P.pop[0].eh_vb2d;',
        '            P.pop[0].X_eh_vb2d2 = P.pop[0].eh_vb2d2;',
        '            P.pop[0].X_eh_vb2d3 = P.pop[0].eh_vb2d3;',
        '            P.pop[0].X_eh_vb3d  = P.pop[0].eh_vb3d;',
        '            P.pop[0].X_eh_vb3d2 = P.pop[0].eh_vb3d2;',
        '            P.pop[0].X_eh_vb3d3 = P.pop[0].eh_vb3d3;',
        '            P.pop[0].X_em_va1d  = P.pop[0].em_va1d;',
        '            P.pop[0].X_em_va1d2 = P.pop[0].em_va1d2;',
        '            P.pop[0].X_em_va1d3 = P.pop[0].em_va1d3;',
        '            P.pop[0].X_em_va2d  = P.pop[0].em_va2d;',
        '            P.pop[0].X_em_va2d2 = P.pop[0].em_va2d2;',
        '            P.pop[0].X_em_va2d3 = P.pop[0].em_va2d3;',
        '            P.pop[0].X_em_va3d  = P.pop[0].em_va3d;',
        '            P.pop[0].X_em_va3d2 = P.pop[0].em_va3d2;',
        '            P.pop[0].X_em_va3d3 = P.pop[0].em_va3d3;',
        '            P.pop[0].X_em_vb1d  = P.pop[0].em_vb1d;',
        '            P.pop[0].X_em_vb1d2 = P.pop[0].em_vb1d2;',
        '            P.pop[0].X_em_vb1d3 = P.pop[0].em_vb1d3;',
        '            P.pop[0].X_em_vb2d  = P.pop[0].em_vb2d;',
        '            P.pop[0].X_em_vb2d2 = P.pop[0].em_vb2d2;',
        '            P.pop[0].X_em_vb2d3 = P.pop[0].em_vb2d3;',
        '            P.pop[0].X_em_vb3d  = P.pop[0].em_vb3d;',
        '            P.pop[0].X_em_vb3d2 = P.pop[0].em_vb3d2;',
        '            P.pop[0].X_em_vb3d3 = P.pop[0].em_vb3d3;',
        '        }',
        '        boost(P.pop[0].ei_va2[a] , P.pop[0].ed_va2i[a],  P.pop[0].eh_va2d[a] , P.pop[0].em_va2d[a] , boost_om_fold, 1.0);',
        '        boost(P.pop[0].ei2_va2[a], P.pop[0].ed_va2i2[a], P.pop[0].eh_va2d2[a], P.pop[0].em_va2d2[a], 1,             1.0);',
        '        boost(P.pop[0].ei3_va2[a], P.pop[0].ed_va2i3[a], P.pop[0].eh_va2d3[a], P.pop[0].em_va2d3[a], 1,             1.0);',
        '        boost(P.pop[0].ei_vb2[a] , P.pop[0].ed_vb2i[a],  P.pop[0].eh_vb2d[a] , P.pop[0].em_vb2d[a] , boost_om_fold, 1.0);',
        '        boost(P.pop[0].ei2_vb2[a], P.pop[0].ed_vb2i2[a], P.pop[0].eh_vb2d2[a], P.pop[0].em_vb2d2[a], 1,             1.0);',
        '        boost(P.pop[0].ei3_vb2[a], P.pop[0].ed_vb2i3[a], P.pop[0].eh_vb2d3[a], P.pop[0].em_vb2d3[a], 1,             1.0);',
        '        boost(P.pop[0].pi_r[a],    P.pop[0].pd_ri[a],    P.pop[0].ph_rd[a],    P.pop[0].pm_rd[a],    boost_om_fold, frac_V_in_R * proportion_booster[a]);',
        '        boost(P.pop[0].pi2_r[a],   P.pop[0].pd_ri2[a],   P.pop[0].ph_rd2[a],   P.pop[0].pm_rd2[a],   1,             frac_V_in_R * proportion_booster[a]);',
        '        boost(P.pop[0].pi3_r[a],   P.pop[0].pd_ri3[a],   P.pop[0].ph_rd3[a],   P.pop[0].pm_rd3[a],   1,             frac_V_in_R * proportion_booster[a]);',
        '        boost(P.pop[0].pi_r2[a],   P.pop[0].pd_r2i[a],   P.pop[0].ph_r2d[a],   P.pop[0].pm_r2d[a],   boost_om_fold, frac_V_in_R * proportion_booster[a]);',
        '        boost(P.pop[0].pi2_r2[a],  P.pop[0].pd_r2i2[a],  P.pop[0].ph_r2d2[a],  P.pop[0].pm_r2d2[a],  1,             frac_V_in_R * proportion_booster[a]);',
        '        boost(P.pop[0].pi3_r2[a],  P.pop[0].pd_r2i3[a],  P.pop[0].ph_r2d3[a],  P.pop[0].pm_r2d3[a],  1,             frac_V_in_R * proportion_booster[a]);',
        '        boost(P.pop[0].pi_r3[a],   P.pop[0].pd_r3i[a],   P.pop[0].ph_r3d[a],   P.pop[0].pm_r3d[a],   boost_om_fold, frac_V_in_R * proportion_booster[a]);',
        '        boost(P.pop[0].pi2_r3[a],  P.pop[0].pd_r3i2[a],  P.pop[0].ph_r3d2[a],  P.pop[0].pm_r3d2[a],  1,             frac_V_in_R * proportion_booster[a]);',
        '        boost(P.pop[0].pi3_r3[a],  P.pop[0].pd_r3i3[a],  P.pop[0].ph_r3d3[a],  P.pop[0].pm_r3d3[a],  1,             frac_V_in_R * proportion_booster[a]);',
        '    }',
        '    if (t == end_days[a]) {',
        '        // RESTORE',
        '        P.pop[0].pi_r[a] =       P.pop[0].X_pi_r[a];',
        '        P.pop[0].pi_r2[a] =      P.pop[0].X_pi_r2[a];',
        '        P.pop[0].pi_r3[a] =      P.pop[0].X_pi_r3[a];',
        '        P.pop[0].pi2_r[a] =      P.pop[0].X_pi2_r[a];',
        '        P.pop[0].pi2_r2[a] =     P.pop[0].X_pi2_r2[a];',
        '        P.pop[0].pi2_r3[a] =     P.pop[0].X_pi2_r3[a];',
        '        P.pop[0].pi3_r[a] =      P.pop[0].X_pi3_r[a];',
        '        P.pop[0].pi3_r2[a] =     P.pop[0].X_pi3_r2[a];',
        '        P.pop[0].pi3_r3[a] =     P.pop[0].X_pi3_r3[a];',
        '        P.pop[0].ei_va1[a] =     P.pop[0].X_ei_va1[a];',
        '        P.pop[0].ei2_va1[a] =    P.pop[0].X_ei2_va1[a];',
        '        P.pop[0].ei3_va1[a] =    P.pop[0].X_ei3_va1[a];',
        '        P.pop[0].ei_va2[a] =     P.pop[0].X_ei_va2[a];',
        '        P.pop[0].ei2_va2[a] =    P.pop[0].X_ei2_va2[a];',
        '        P.pop[0].ei3_va2[a] =    P.pop[0].X_ei3_va2[a];',
        '        P.pop[0].ei_va3[a] =     P.pop[0].X_ei_va3[a];',
        '        P.pop[0].ei2_va3[a] =    P.pop[0].X_ei2_va3[a];',
        '        P.pop[0].ei3_va3[a] =    P.pop[0].X_ei3_va3[a];',
        '        P.pop[0].ei_vb1[a] =     P.pop[0].X_ei_vb1[a];',
        '        P.pop[0].ei2_vb1[a] =    P.pop[0].X_ei2_vb1[a];',
        '        P.pop[0].ei3_vb1[a] =    P.pop[0].X_ei3_vb1[a];',
        '        P.pop[0].ei_vb2[a] =     P.pop[0].X_ei_vb2[a];',
        '        P.pop[0].ei2_vb2[a] =    P.pop[0].X_ei2_vb2[a];',
        '        P.pop[0].ei3_vb2[a] =    P.pop[0].X_ei3_vb2[a];',
        '        P.pop[0].ei_vb3[a] =     P.pop[0].X_ei_vb3[a];',
        '        P.pop[0].ei2_vb3[a] =    P.pop[0].X_ei2_vb3[a];',
        '        P.pop[0].ei3_vb3[a] =    P.pop[0].X_ei3_vb3[a];',
        '        P.pop[0].ph_rd[a] =      P.pop[0].X_ph_rd[a];',
        '        P.pop[0].ph_rd2[a] =     P.pop[0].X_ph_rd2[a];',
        '        P.pop[0].ph_rd3[a] =     P.pop[0].X_ph_rd3[a];',
        '        P.pop[0].ph_r2d[a] =     P.pop[0].X_ph_r2d[a];',
        '        P.pop[0].ph_r2d2[a] =    P.pop[0].X_ph_r2d2[a];',
        '        P.pop[0].ph_r2d3[a] =    P.pop[0].X_ph_r2d3[a];',
        '        P.pop[0].ph_r3d[a] =     P.pop[0].X_ph_r3d[a];',
        '        P.pop[0].ph_r3d2[a] =    P.pop[0].X_ph_r3d2[a];',
        '        P.pop[0].ph_r3d3[a] =    P.pop[0].X_ph_r3d3[a];',
        '        P.pop[0].pm_rd[a] =      P.pop[0].X_pm_rd[a];',
        '        P.pop[0].pm_rd2[a] =     P.pop[0].X_pm_rd2[a];',
        '        P.pop[0].pm_rd3[a] =     P.pop[0].X_pm_rd3[a];',
        '        P.pop[0].pm_r2d[a] =     P.pop[0].X_pm_r2d[a];',
        '        P.pop[0].pm_r2d2[a] =    P.pop[0].X_pm_r2d2[a];',
        '        P.pop[0].pm_r2d3[a] =    P.pop[0].X_pm_r2d3[a];',
        '        P.pop[0].pm_r3d[a] =     P.pop[0].X_pm_r3d[a];',
        '        P.pop[0].pm_r3d2[a] =    P.pop[0].X_pm_r3d2[a];',
        '        P.pop[0].pm_r3d3[a] =    P.pop[0].X_pm_r3d3[a];',
        '        P.pop[0].eh_va1d [a] = P.pop[0].X_eh_va1d[a];',
        '        P.pop[0].eh_va1d2[a] = P.pop[0].X_eh_va1d2[a];',
        '        P.pop[0].eh_va1d3[a] = P.pop[0].X_eh_va1d3[a];',
        '        P.pop[0].eh_va2d [a] = P.pop[0].X_eh_va2d[a];',
        '        P.pop[0].eh_va2d2[a] = P.pop[0].X_eh_va2d2[a];',
        '        P.pop[0].eh_va2d3[a] = P.pop[0].X_eh_va2d3[a];',
        '        P.pop[0].eh_va3d [a] = P.pop[0].X_eh_va3d[a];',
        '        P.pop[0].eh_va3d2[a] = P.pop[0].X_eh_va3d2[a];',
        '        P.pop[0].eh_va3d3[a] = P.pop[0].X_eh_va3d3[a];',
        '        P.pop[0].eh_vb1d [a] = P.pop[0].X_eh_vb1d[a];',
        '        P.pop[0].eh_vb1d2[a] = P.pop[0].X_eh_vb1d2[a];',
        '        P.pop[0].eh_vb1d3[a] = P.pop[0].X_eh_vb1d3[a];',
        '        P.pop[0].eh_vb2d [a] = P.pop[0].X_eh_vb2d[a];',
        '        P.pop[0].eh_vb2d2[a] = P.pop[0].X_eh_vb2d2[a];',
        '        P.pop[0].eh_vb2d3[a] = P.pop[0].X_eh_vb2d3[a];',
        '        P.pop[0].eh_vb3d [a] = P.pop[0].X_eh_vb3d[a];',
        '        P.pop[0].eh_vb3d2[a] = P.pop[0].X_eh_vb3d2[a];',
        '        P.pop[0].eh_vb3d3[a] = P.pop[0].X_eh_vb3d3[a];',
        '        P.pop[0].em_va1d [a] = P.pop[0].X_em_va1d[a];',
        '        P.pop[0].em_va1d2[a] = P.pop[0].X_em_va1d2[a];',
        '        P.pop[0].em_va1d3[a] = P.pop[0].X_em_va1d3[a];',
        '        P.pop[0].em_va2d [a] = P.pop[0].X_em_va2d[a];',
        '        P.pop[0].em_va2d2[a] = P.pop[0].X_em_va2d2[a];',
        '        P.pop[0].em_va2d3[a] = P.pop[0].X_em_va2d3[a];',
        '        P.pop[0].em_va3d [a] = P.pop[0].X_em_va3d[a];',
        '        P.pop[0].em_va3d2[a] = P.pop[0].X_em_va3d2[a];',
        '        P.pop[0].em_va3d3[a] = P.pop[0].X_em_va3d3[a];',
        '        P.pop[0].em_vb1d [a] = P.pop[0].X_em_vb1d[a];',
        '        P.pop[0].em_vb1d2[a] = P.pop[0].X_em_vb1d2[a];',
        '        P.pop[0].em_vb1d3[a] = P.pop[0].X_em_vb1d3[a];',
        '        P.pop[0].em_vb2d [a] = P.pop[0].X_em_vb2d[a];',
        '        P.pop[0].em_vb2d2[a] = P.pop[0].X_em_vb2d2[a];',
        '        P.pop[0].em_vb2d3[a] = P.pop[0].X_em_vb2d3[a];',
        '        P.pop[0].em_vb3d [a] = P.pop[0].X_em_vb3d[a];',
        '        P.pop[0].em_vb3d2[a] = P.pop[0].X_em_vb3d2[a];',
        '        P.pop[0].em_vb3d3[a] = P.pop[0].X_em_vb3d3[a];',
        '    }',
        '}',
        '}',

        .sep = "\n", .open = "${", .close = "}")
}

cpp_obsI_printVE = function(file_basename, p, t_step = 14)
{
    glue::glue(
        '{',
        '    if ((unsigned int)t % ${t_step} == 0) {',
        '        std::string filename = "${file_basename}" + to_string(${p}) + ".csv";',
        '        std::ofstream fout(filename, ios_base::app);',
        '        if (t == 0) {',
        '            P.pop[0].OutputVELine(-1, -1, true, fout);',
        '        }',
        '        P.pop[0].OutputVELine(t, ${p}, false, fout);',
        '    }',
        '}',

        .sep = "\n", .open = "${", .close = "}")
}
