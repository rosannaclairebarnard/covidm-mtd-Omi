library(data.table)

# test = w
# testname = "deaths"
# dispname = "disp_deaths"
# data = ld
# dataname = "N"
# alpha = 0.1
# k = 7
# sigma = 0.3
# extranames = c("deaths_V0", "deaths_V1", "deaths_V2", "deaths_V3")
# adj_file = "./fits/final_adj_2021-12-05.csv"
# dstart = "2021-11-01"
# dend = "2021-11-30"

apply_gamma_multiplier = function(test, testname, dispname, data, dataname, alpha, k, sigma, extranames, 
    adj_file, dstart, dend, adj_ts, ts_start, ts_end)
{
    # Amalgamate data
    mdata = merge(
        test[, .(model = sum(get(testname))), by = .(t, population, run)], 
        data[, .(t = as.numeric(date - ymd("2020-01-01")), data = get(dataname), population = name)], 
        by = c("t", "population"), all.x = TRUE);
    mdata = mdata[order(run, population, t)];
    mdata[, prior_mean := pmax(alpha, model)];
    mdata[, prior_sd := prior_mean * sigma];
    mdata[, prior_shape := prior_mean^2 / prior_sd^2];
    mdata[, prior_rate  := prior_mean / prior_sd^2];
    mdata[!is.na(data), posterior_shape := zoo::rollmean(prior_shape + data, k = k, fill = "extend"), by = .(run, population)];
    mdata[!is.na(data), posterior_rate := prior_rate + 1];
    
    mdata[, gamma_multiplier := (posterior_shape / posterior_rate) / model];
    mdata[, gamma_disp := sqrt(1.0 / posterior_shape)];
    
    mdata[, gamma_multiplier := zoo::na.locf(gamma_multiplier, na.rm = FALSE), by = .(run, population)];
    mdata[, gamma_multiplier := zoo::na.locf(gamma_multiplier, na.rm = FALSE, fromLast = TRUE), by = .(run, population)];

    mdata[, gamma_disp := zoo::na.locf(gamma_disp, na.rm = FALSE), by = .(run, population)];
    mdata[, gamma_disp := zoo::na.locf(gamma_disp, na.rm = FALSE, fromLast = TRUE), by = .(run, population)];
    
    # Bring in fixed final adjustment, if supplied
    if (!is.na(adj_file)) {
        final_adj = fread(adj_file)
        final_adj = final_adj[, .(adj = get(paste0(testname, "_adj"))), by = .(population = NHS.region)]
        mdata = merge(mdata, final_adj, by = "population", all.x = TRUE)
        tstart = as.numeric(ymd(dstart) - ymd("2020-01-01"))
        tend   = as.numeric(ymd(dend)   - ymd("2020-01-01"))
        mdata[, date := ymd("2020-01-01") + t]
        mdata[date >= dstart & date <= dend, 
            gamma_multiplier := ((tend - t) / (tend - tstart)) * gamma_multiplier + ((t - tstart) / (tend - tstart)) * adj]
        mdata[date > dend,
            gamma_multiplier := adj]
        mdata[, date := NULL]
        mdata[, adj := NULL]
            
    }
    
    # Bring in time-series final adjustments, if supplied
    if (!is.na(adj_ts)){
        ts_adj = fread(adj_ts)
        ts_adj = ts_adj[, .(adj = get(paste0(testname, "_adj"))), by = .(t, date, population = NHS.region)]
        mdata = merge(mdata, ts_adj, by = c("t", "population"), all.x = TRUE)
        mdata[, date := ymd("2020-01-01") + t]
        mdata[date >= ts_start & date <= ts_end,
              gamma_multiplier := adj]
        mdata[, date := NULL]
        mdata[, adj := NULL]
    }

    test = merge(test, mdata[, .(t, population, run, gamma_multiplier, gamma_disp)], by = c("t", "population", "run"));
    for (v in c(testname, extranames)) {
        test[!is.na(gamma_multiplier), (v) := get(v) * gamma_multiplier]
    }
    test[!is.na(gamma_multiplier), (dispname) := gamma_disp];

    test[, gamma_multiplier := NULL]
    test[, gamma_disp := NULL]

    return (test)
}

SPIM_output_1var = function(summ, varname, outname, cyear, cmonth, cday, ymd_from, ymd_to, size)
{
    quants = seq(0.05, 0.95, by = 0.05)

    # new code...
    rq = data.table(run = 1:summ[, max(run)], q = runif(summ[, max(run)]))
    summ = merge(summ, rq, by = "run")
    if (is.character(size)) {
        summ[, rv := qnbinom(q, size = 1.0 / get(size)^2, mu = get(varname))]
    } else if (size == 0) { # poisson
        summ[, rv := qpois(q, lambda = get(varname))]
    } else if (size == -1) { # no variation
        summ[, rv := get(varname)]
    } else if (size > 0) { # neg binom
        summ[, rv := qnbinom(q, size = size, mu = get(varname))]
    } else {
        stop("Size must be -1, 0, or a positive number")
    }
    qsumm = summ[, as.list(quantile(rv, quants, na.rm = TRUE)), by = .(t, population, age_group)]

    #qsumm = summ[, as.list(qnbinom(quants, size = 20, mu = quantile(get(varname), quants, na.rm = T))), by = .(t, population, age_group)]
    qsumm[, date := t + ymd("2020-01-01")]
    
    # Add "All" age group if not already present
    if (!"All" %in% qsumm[, unique(age_group)]) {
        if (outname == "sero_prev") {
            nr = qsumm[!age_group %in% c("0-4", "5-14"), .N]
            qsumm = rbind(qsumm,
                          qsumm[!age_group %in% c("0-4", "5-14"), 
                              lapply(.SD, sum), .SDcols = `5%`:`95%`, by = .(t, population, date, age_group = rep("All", nr))],
                          use.names = TRUE, fill = TRUE)
        } else {
            qsumm = rbind(qsumm,
                          qsumm[, 
                              lapply(.SD, sum), .SDcols = `5%`:`95%`, by = .(t, population, date, age_group = rep("All", nrow(qsumm)))],
                          use.names = TRUE, fill = TRUE)
        }
    }

    data_out = qsumm[date %between% c(ymd_from, ymd_to),
        .(Group = "LSHTM", Model = "Transmission", Scenario = "MTP", ModelType = "Multiple", Version = 2,
        `Creation Day` = cday, `Creation Month` = cmonth, `Creation Year` = cyear,
        `Day of Value` = day(date), `Month of Value` = month(date), `Year of Value` = year(date),
        AgeBand = age_group, Geography = population, ValueType = outname, Value = `50%`,
        `5%`, `10%`, `15%`, `20%`, `25%`, `30%`, `35%`, `40%`, `45%`, `50%`, `55%`, `60%`, `65%`, `70%`, `75%`, `80%`, `85%`, `90%`, `95%`)]

    names(data_out)[names(data_out) %like% "\\%$"] = paste("Quantile", quants);
    data_out
}

SPIM_output_full = function(test0, cyear, cmonth, cday, ymd_from, ymd_to)
{
    cat("Restricting time span...\n")
    t0 = as.numeric(ymd(ymd_from) - ymd("2020-01-01"))
    t1 = as.numeric(ymd(ymd_to) - ymd("2020-01-01"))
    test = test0[t %between% c(t0, t1)]

    cat("Adding age groups...\n");
    test[group %between% c(1, 1), age_group := "0-4"]
    test[group %between% c(2, 3), age_group := "5-14"]
    test[group %between% c(4, 5), age_group := "15-24"]
    test[group %between% c(6, 9), age_group := "25-44"]
    test[group %between% c(10, 13), age_group := "45-64"]
    test[group %between% c(14, 15), age_group := "65-74"]
    test[group %between% c(16, 16), age_group := "75+"]
    test[, age_group := factor(age_group, levels = unique(age_group))]

    cat("Summarizing variables...\n");
    summ = test[, .(death_o = sum(deaths), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
    summ3 = test[, .(icu_p = sum(icu_bed), disp = mean(disp_icu_prev)), by = .(run, t, population, age_group)]
    summ4 = test[, .(bed_p = sum(pmax(0, hosp_bed - hosp_undetected_p)), disp = mean(disp_hosp_prev)), by = .(run, t, population, age_group)]
    summ34 = test[, .(admissions = sum(hosp_undetected_o), disp = mean(disp_hosp_inc)), by = .(run, t, population, age_group)]
    summ_inf_i = test[, .(infections_i = sum(pcr_positive_i)), by = .(run, t, population, age_group)]
    summ_inf_p = test[, .(infections_p = sum(pcr_positive_p)), by = .(run, t, population, age_group)];
    summ_sero_p = test[, .(sero_p = sum(lfia_positive_p + vaccsero_a_p + vaccsero_b_p)), by = .(run, t, population, age_group)];
    summ_rt = test[group == 1, .(age_group = "All", Rt = obs0), by = .(run, t, population)]
    summ_r0 = test[group == 4, .(age_group = "All", R0 = obs0), by = .(run, t, population)]
    summ_test = test[, .(test1 = sum(test_o), test2 = sum(test2_o), test3 = sum(test3_o)), by = .(run, t, population, age_group)];
    summ_d_V0 = test[, .(death_o = sum(deaths_V0), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
    summ_d_V1 = test[, .(death_o = sum(deaths_V1), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
    summ_d_V2 = test[, .(death_o = sum(deaths_V2), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
    summ_d_V3 = test[, .(death_o = sum(deaths_V3), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
    summ_a_V0 = test[, .(admissions = sum(hosp_adm_V0), disp = mean(disp_hosp_inc)), by = .(run, t, population, age_group)]
    summ_a_V1 = test[, .(admissions = sum(hosp_adm_V1), disp = mean(disp_hosp_inc)), by = .(run, t, population, age_group)]
    summ_a_V2 = test[, .(admissions = sum(hosp_adm_V2), disp = mean(disp_hosp_inc)), by = .(run, t, population, age_group)]
    summ_a_V3 = test[, .(admissions = sum(hosp_adm_V3), disp = mean(disp_hosp_inc)), by = .(run, t, population, age_group)]
    # adding output to record number of individuals in vaccine compartments
    summ_Va1_i = test[, .(num_vacc = sum(Va1)), by = .(run, t, population, age_group)]
    summ_Va2_i = test[, .(num_vacc = sum(Va2)), by = .(run, t, population, age_group)]
    summ_Va3_i = test[, .(num_vacc = sum(Va3)), by = .(run, t, population, age_group)]
    summ_Vb1_i = test[, .(num_vacc = sum(Vb1)), by = .(run, t, population, age_group)]
    summ_Vb2_i = test[, .(num_vacc = sum(Vb2)), by = .(run, t, population, age_group)]
    summ_Vb3_i = test[, .(num_vacc = sum(Vb3)), by = .(run, t, population, age_group)]


    cat("Running quantiles...\n");
    w = rbind(
        SPIM_output_1var(summ, "death_o", "type28_death_inc_line", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ3, "icu_p", "icu_prev", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ4, "bed_p", "hospital_prev", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ34, "admissions", "hospital_inc", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ_inf_p, "infections_p", "prevalence_mtp", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_inf_i, "infections_i", "infections_inc", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_sero_p, "sero_p", "sero_prev", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_rt, "Rt", "Rt", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_r0, "R0", "R0", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_test, "test1", "inc_v1", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_test, "test2", "inc_v2", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_test, "test3", "inc_v3", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_d_V0, "death_o", "deaths_V0" , cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ_d_V1, "death_o", "deaths_V1", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ_d_V2, "death_o", "deaths_V2", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ_d_V3, "death_o", "deaths_V3", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ_a_V0, "admissions", "admissions_V0", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ_a_V1, "admissions", "admissions_V1", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ_a_V2, "admissions", "admissions_V2", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ_a_V3, "admissions", "admissions_V3", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        # testing getting traces of Va1, Va2, Va3, Vb1, Vb2, Vb3 over time
        SPIM_output_1var(summ_Va1_i, "num_vacc", "Va1", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_Va2_i, "num_vacc", "Va2", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_Va3_i, "num_vacc", "Va3", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_Vb1_i, "num_vacc", "Vb1", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_Vb2_i, "num_vacc", "Vb2", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_Vb3_i, "num_vacc", "Vb3", cyear, cmonth, cday, ymd_from, ymd_to, -1)
    )
    return (w)
}

# next function needs a fix to calculation of sero_p to include vaccine sero
# SPIM_output_full_3strain = function(test0, cyear, cmonth, cday, ymd_from, ymd_to)
# {
#     cat("Restricting time span...\n")
#     t0 = as.numeric(ymd(ymd_from) - ymd("2020-01-01"))
#     t1 = as.numeric(ymd(ymd_to) - ymd("2020-01-01"))
#     test = test0[t %between% c(t0, t1)]
#     
#     cat("Adding age groups...\n");
#     test[group %between% c(1, 1), age_group := "0-4"]
#     test[group %between% c(2, 3), age_group := "5-14"]
#     test[group %between% c(4, 5), age_group := "15-24"]
#     test[group %between% c(6, 9), age_group := "25-44"]
#     test[group %between% c(10, 13), age_group := "45-64"]
#     test[group %between% c(14, 15), age_group := "65-74"]
#     test[group %between% c(16, 16), age_group := "75+"]
#     test[, age_group := factor(age_group, levels = unique(age_group))]
#     
#     cat("Summarizing variables...\n");
#     summd1 = test[, .(death_o = sum(death_o), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
#     summd2 = test[, .(death_o = sum(death2_o), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
#     summd3 = test[, .(death_o = sum(death3_o), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
#     summ = test[, .(death_o = sum(death_o + death2_o + death3_o), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
#     summicu1 = test[, .(icu_p = sum(icu_p), disp = mean(disp_icu_prev)), by = .(run, t, population, age_group)]
#     summicu2 = test[, .(icu_p = sum(icu2_p), disp = mean(disp_icu_prev)), by = .(run, t, population, age_group)]
#     summicu3 = test[, .(icu_p = sum(icu3_p), disp = mean(disp_icu_prev)), by = .(run, t, population, age_group)]
#     summ3 = test[, .(icu_p = sum(icu_p + icu2_p + icu3_p), disp = mean(disp_icu_prev)), by = .(run, t, population, age_group)]
#     summhosp1 = test[, .(bed_p = sum(pmax(0, hosp_p - hosp_undetected_p)), disp = mean(disp_hosp_prev)), by = .(run, t, population, age_group)]
#     summhosp2 = test[, .(bed_p = sum(pmax(0, hosp2_p - hosp_undetected2_p)), disp = mean(disp_hosp_prev)), by = .(run, t, population, age_group)]
#     summhosp3 = test[, .(bed_p = sum(pmax(0, hosp3_p - hosp_undetected3_p)), disp = mean(disp_hosp_prev)), by = .(run, t, population, age_group)]
#     summ4 = test[, .(bed_p = sum(pmax(0, hosp_p + hosp2_p + hosp3_p - hosp_undetected_p - hosp_undetected2_p - hosp_undetected3_p)), disp = mean(disp_hosp_prev)), by = .(run, t, population, age_group)]
#     summhospo1 = test[, .(admissions = sum(hosp_undetected_o), disp = disp_hosp_inc), by = .(run, t, population, age_group)]
#     summhospo2 = test[, .(admissions = sum(hosp_undetected2_o), disp = disp_hosp_inc), by = .(run, t, population, age_group)]
#     summhospo3 = test[, .(admissions = sum(hosp_undetected3_o), disp = disp_hosp_inc), by = .(run, t, population, age_group)]
#     summ34 = test[, .(admissions = sum(hosp_undetected_o + hosp_undetected2_o + hosp_undetected3_o), disp = disp_hosp_inc), by = .(run, t, population, age_group)]
#     summ_inf_i = test[, .(infections_i = sum(pcr_positive_i)), by = .(run, t, population, age_group)]
#     summ_inf_p = test[, .(infections_p = sum(pcr_positive_p)), by = .(run, t, population, age_group)];
#     summ_sero_p = test[, .(sero_p = sum(lfia_positive_p)), by = .(run, t, population, age_group)];
#     summ_rt = test[group == 1, .(age_group = "All", Rt = obs0), by = .(run, t, population)]
#     summ_r0 = test[group == 4, .(age_group = "All", R0 = obs0), by = .(run, t, population)]
#     
#     cat("Running quantiles...\n");
#     w = rbind(
#         SPIM_output_1var(summd1, "death_o", "type28_death_inc_line_s1", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summd2, "death_o", "type28_death_inc_line_s2", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summd3, "death_o", "type28_death_inc_line_s3", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summ, "death_o", "type28_death_inc_line", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summicu1, "icu_p", "icu_prev_s1", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summicu2, "icu_p", "icu_prev_s2", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summicu3, "icu_p", "icu_prev_s3", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summ3, "icu_p", "icu_prev", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summhosp1, "bed_p", "hospital_prev_s1", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summhosp2, "bed_p", "hospital_prev_s2", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summhosp3, "bed_p", "hospital_prev_s3", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summ4, "bed_p", "hospital_prev", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summhospo1, "admissions", "hospital_inc_s1", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summhospo2, "admissions", "hospital_inc_s2", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summhospo3, "admissions", "hospital_inc_s3", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summ34, "admissions", "hospital_inc", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
#         SPIM_output_1var(summ_inf_p, "infections_p", "prevalence_mtp", cyear, cmonth, cday, ymd_from, ymd_to, -1),
#         SPIM_output_1var(summ_inf_i, "infections_i", "infections_inc", cyear, cmonth, cday, ymd_from, ymd_to, -1),
#         SPIM_output_1var(summ_sero_p, "sero_p", "sero_prev", cyear, cmonth, cday, ymd_from, ymd_to, -1),
#         SPIM_output_1var(summ_rt, "Rt", "Rt", cyear, cmonth, cday, ymd_from, ymd_to, -1),
#         SPIM_output_1var(summ_r0, "R0", "R0", cyear, cmonth, cday, ymd_from, ymd_to, -1)
#     )
#     return (w)
# }

make_data = function(ld, sitreps, virus, sero)
{
    rbind(
        ld[, .(ValueType = "type28_death_inc_line", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = N, ymax = NA)],
        sitreps[, .(ValueType = "icu_prev", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_in_itu, ymax = NA)],
        sitreps[, .(ValueType = "hospital_prev", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_in_all_beds, ymax = NA)],
        sitreps[, .(ValueType = "hospital_inc", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_admitted_diagnosed, ymax = NA)],
        virus[Data.source %like% "ONS", .(ValueType = "prevalence_mtp", Geography = NHS.region,
            dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
            ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))],
        sero[!Data.source %like% "NHS BT", .(ValueType = "sero_prev", Geography = NHS.region,
            dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
            ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))]
    )
}
