# v2/cm2_params.R
# parameters for the model

# return translated parameters to work with the backend,
# i.e. fix any times expressed in dates to be expressed in days since date0.
cm_translate_parameters = function(p)
{
    translate_time = function(t) {
        if (is.numeric(t)) {
            return (t)
        } else {
            return (as.numeric(ymd(t) - ymd(p$date0)));
        }
    }

    p$time0 = translate_time(p$time0);
    p$time1 = translate_time(p$time1);

    for (pi in seq_along(p$pop)) {
        p$pop[[pi]]$seed_times = sapply(p$pop[[pi]]$seed_times, translate_time);
        p$pop[[pi]]$seed_times2 = sapply(p$pop[[pi]]$seed_times2, translate_time);
        p$pop[[pi]]$seed_times3 = sapply(p$pop[[pi]]$seed_times3, translate_time);
    }

    for (si in seq_along(p$schedule)) {
        p$schedule[[si]]$times = translate_time(p$schedule[[si]]$times);
    }

    return (p);
}

# check parameters for validity.
cm_check_parameters = function(parameters)
{
    req = function(v, x)
    {
        if (!exists(v, parameters)) {
            stop(paste0("Parameter ", v, " required, but not found."));
        } else if (!eval(parse(text = x), parameters)) {
            stop(paste0("Parameters check: ", x, " failed."));
        }
    }

    reqp = function(i, v, x)
    {
        if (!exists(v, parameters$pop[[i]])) {
            stop(paste0("Population parameter ", v, " required, but not found in population ", i, "."));
        } else if (!eval(parse(text = x), parameters$pop[[i]])) {
            stop(paste0("Parameters check: ", x, " failed in population ", i, "."));
        }
    }

    req("model",            "model == 'SEI3R' || model == 'household'");
    req("time_step",        "is.numeric(time_step) & time_step > 0");
    req("date0",            "is.Date(ymd(date0))");
    req("time0",            "is.numeric(time0)");
    req("time1",            "is.numeric(time1)");
    req("report_every",     "is.numeric(report_every) & report_every == 1. / time_step");
    req("fast_multinomial", "is.logical(fast_multinomial)");
    req("deterministic",    "is.logical(deterministic)");
    req("pop",              "is.list(pop) & length(pop) > 0");
    req("travel",           "is.matrix(travel) & all(travel >= 0) & nrow(travel) == length(parameters$pop) & ncol(travel) == nrow(travel)");

    for (i in 1:length(parameters$pop))
    {
        reqp(i, "dE",   "is.numeric(dE) & all(dE >= 0) & any(dE > 0)");
        reqp(i, "dIp",  "is.numeric(dIp) & all(dIp >= 0) & any(dIp > 0)");
        reqp(i, "dIs",  "is.numeric(dIs) & all(dIs >= 0) & any(dIs > 0)");
        reqp(i, "dIa",  "is.numeric(dIa) & all(dIa >= 0) & any(dIa > 0)");
        reqp(i, "dE2",  "is.numeric(dE2) & all(dE2 >= 0) & any(dE2 > 0)");
        reqp(i, "dIp2", "is.numeric(dIp2) & all(dIp2 >= 0) & any(dIp2 > 0)");
        reqp(i, "dIs2", "is.numeric(dIs2) & all(dIs2 >= 0) & any(dIs2 > 0)");
        reqp(i, "dIa2", "is.numeric(dIa2) & all(dIa2 >= 0) & any(dIa2 > 0)");
        reqp(i, "dE3",  "is.numeric(dE3) & all(dE3 >= 0) & any(dE3 > 0)");
        reqp(i, "dIp3", "is.numeric(dIp3) & all(dIp3 >= 0) & any(dIp3 > 0)");
        reqp(i, "dIs3", "is.numeric(dIs3) & all(dIs3 >= 0) & any(dIs3 > 0)");
        reqp(i, "dIa3", "is.numeric(dIa3) & all(dIa3 >= 0) & any(dIa3 > 0)");
        reqp(i, "dVa1", "is.numeric(dVa1) & all(dVa1 >= 0) & any(dVa1 > 0)");
        reqp(i, "dVb1", "is.numeric(dVb1) & all(dVb1 >= 0) & any(dVb1 > 0)");
        reqp(i, "dVa2", "is.numeric(dVa2) & all(dVa2 >= 0) & any(dVa2 > 0)");
        reqp(i, "dVb2", "is.numeric(dVb2) & all(dVb2 >= 0) & any(dVb2 > 0)");
        reqp(i, "size", "is.numeric(size) & all(size >= 0) & any(size > 0)");
        reqp(i, "imm0", "is.numeric(imm0) & length(imm0) == length(size) & all(imm0 >= 0) & all(imm0 <= 1)");
        reqp(i, "matrices", "is.list(matrices) & length(matrices) > 0");
        for (m in 1:length(parameters$pop[[i]]$matrices))
        {
            reqp(i, "matrices", sprintf("is.matrix(matrices[[%d]]) & all(matrices[[%d]] >= 0) & nrow(matrices[[%d]]) == length(size) & ncol(matrices[[%d]]) == nrow(matrices[[%d]])", m, m, m, m, m));
        }
        reqp(i, "contact",         "is.numeric(contact) & length(contact) == length(matrices) & all(contact >= 0)");
        reqp(i, "contact_mult",    "(is.numeric(contact_mult) & length(contact_mult) == length(matrices) & all(contact_mult >= 0)) | length(contact_mult) == 0");
        reqp(i, "contact_lowerto", "(is.numeric(contact_lowerto) & length(contact_lowerto) == length(matrices) & all(contact_lowerto >= 0)) | length(contact_lowerto) == 0");

        reqp(i, "u",       "is.numeric(u) & length(u) == length(size) & all(u >= 0)");
        reqp(i, "u2",      "is.numeric(u2) & length(u2) == length(size) & all(u2 >= 0)");
        reqp(i, "u3",      "is.numeric(u3) & length(u3) == length(size) & all(u3 >= 0)");
        reqp(i, "y",       "is.numeric(y) & length(y) == length(size) & all(y >= 0)");
        reqp(i, "y2",      "is.numeric(y2) & length(y2) == length(size) & all(y2 >= 0)");
        reqp(i, "y3",      "is.numeric(y3) & length(y3) == length(size) & all(y3 >= 0)");
        reqp(i, "fIp",     "is.numeric(fIp) & length(fIp) == length(size) & all(fIp >= 0)");
        reqp(i, "fIa",     "is.numeric(fIa) & length(fIa) == length(size) & all(fIa >= 0)");
        reqp(i, "fIs",     "is.numeric(fIs) & length(fIs) == length(size) & all(fIs >= 0)");
        reqp(i, "omega",   "is.numeric(omega) & length(omega) == length(size) & all(omega >= 0)");
        reqp(i, "tau",     "is.numeric(tau) & length(tau) == length(size) & all(tau >= 0)");
        reqp(i, "pi_r",    "is.numeric(pi_r) & length(pi_r) == length(size) & all(pi_r >= 0) & all(pi_r <= 1)");
        reqp(i, "pi_r2",   "is.numeric(pi_r2) & length(pi_r2) == length(size) & all(pi_r2 >= 0) & all(pi_r2 <= 1)");
        reqp(i, "pi_r3",   "is.numeric(pi_r3) & length(pi_r3) == length(size) & all(pi_r3 >= 0) & all(pi_r3 <= 1)");
        reqp(i, "pi2_r",   "is.numeric(pi2_r) & length(pi2_r) == length(size) & all(pi2_r >= 0) & all(pi2_r <= 1)");
        reqp(i, "pi2_r2",  "is.numeric(pi2_r2) & length(pi2_r2) == length(size) & all(pi2_r2 >= 0) & all(pi2_r2 <= 1)");
        reqp(i, "pi2_r3",  "is.numeric(pi2_r3) & length(pi2_r3) == length(size) & all(pi2_r3 >= 0) & all(pi2_r3 <= 1)");
        reqp(i, "pi3_r",   "is.numeric(pi3_r) & length(pi3_r) == length(size) & all(pi3_r >= 0) & all(pi3_r <= 1)");
        reqp(i, "pi3_r2",  "is.numeric(pi3_r2) & length(pi3_r2) == length(size) & all(pi3_r2 >= 0) & all(pi3_r2 <= 1)");
        reqp(i, "pi3_r3",  "is.numeric(pi3_r3) & length(pi3_r3) == length(size) & all(pi3_r3 >= 0) & all(pi3_r3 <= 1)");
        reqp(i, "wn",      "is.numeric(wn) & length(wn) == length(size) & all(wn >= 0)");
        reqp(i, "wn2",     "is.numeric(wn2) & length(wn2) == length(size) & all(wn2 >= 0)");
        reqp(i, "wn3",     "is.numeric(wn3) & length(wn3) == length(size) & all(wn3 >= 0)");
        reqp(i, "va1",       "is.numeric(va1) & length(va1) == length(size) & all(va1 >= 0)");
        reqp(i, "wva1",      "is.numeric(wva1) & length(wva1) == length(size) & all(wva1 >= 0)");
        reqp(i, "ei_va1",    "is.numeric(ei_va1) & length(ei_va1) == length(size) & all(ei_va1 >= 0) & all(ei_va1 <= 1)");
        reqp(i, "ei2_va1",   "is.numeric(ei2_va1) & length(ei2_va1) == length(size) & all(ei2_va1 >= 0) & all(ei2_va1 <= 1)");
        reqp(i, "ei3_va1",   "is.numeric(ei3_va1) & length(ei3_va1) == length(size) & all(ei3_va1 >= 0) & all(ei3_va1 <= 1)");
        reqp(i, "ed_va1i",   "is.numeric(ed_va1i) & length(ed_va1i) == length(size) & all(ed_va1i >= 0) & all(ed_va1i <= 1)");
        reqp(i, "ed_va1i2",  "is.numeric(ed_va1i2) & length(ed_va1i2) == length(size) & all(ed_va1i2 >= 0) & all(ed_va1i2 <= 1)");
        reqp(i, "ed_va1i3",  "is.numeric(ed_va1i3) & length(ed_va1i3) == length(size) & all(ed_va1i3 >= 0) & all(ed_va1i3 <= 1)");
        reqp(i, "wva2",      "is.numeric(wva2) & length(wva2) == length(size) & all(wva2 >= 0)");
        reqp(i, "ei_va2",    "is.numeric(ei_va2) & length(ei_va2) == length(size) & all(ei_va2 >= 0) & all(ei_va2 <= 1)");
        reqp(i, "ei2_va2",   "is.numeric(ei2_va2) & length(ei2_va2) == length(size) & all(ei2_va2 >= 0) & all(ei2_va2 <= 1)");
        reqp(i, "ei3_va2",   "is.numeric(ei3_va2) & length(ei3_va2) == length(size) & all(ei3_va2 >= 0) & all(ei3_va2 <= 1)");
        reqp(i, "ed_va2i",   "is.numeric(ed_va2i) & length(ed_va2i) == length(size) & all(ed_va2i >= 0) & all(ed_va2i <= 1)");
        reqp(i, "ed_va2i2",  "is.numeric(ed_va2i2) & length(ed_va2i2) == length(size) & all(ed_va2i2 >= 0) & all(ed_va2i2 <= 1)");
        reqp(i, "ed_va2i3",  "is.numeric(ed_va2i3) & length(ed_va2i3) == length(size) & all(ed_va2i3 >= 0) & all(ed_va2i3 <= 1)");
        reqp(i, "wva3",      "is.numeric(wva3) & length(wva3) == length(size) & all(wva3 >= 0)");
        reqp(i, "ei_va3",    "is.numeric(ei_va3) & length(ei_va3) == length(size) & all(ei_va3 >= 0) & all(ei_va3 <= 1)");
        reqp(i, "ei2_va3",   "is.numeric(ei2_va3) & length(ei2_va3) == length(size) & all(ei2_va3 >= 0) & all(ei2_va3 <= 1)");
        reqp(i, "ei3_va3",   "is.numeric(ei3_va3) & length(ei3_va3) == length(size) & all(ei3_va3 >= 0) & all(ei3_va3 <= 1)");
        reqp(i, "ed_va3i",   "is.numeric(ed_va3i) & length(ed_va3i) == length(size) & all(ed_va3i >= 0) & all(ed_va3i <= 1)");
        reqp(i, "ed_va3i2",  "is.numeric(ed_va3i2) & length(ed_va3i2) == length(size) & all(ed_va3i2 >= 0) & all(ed_va3i2 <= 1)");
        reqp(i, "ed_va3i3",  "is.numeric(ed_va3i3) & length(ed_va3i3) == length(size) & all(ed_va3i3 >= 0) & all(ed_va3i3 <= 1)");
        reqp(i, "vb1",       "is.numeric(vb1) & length(vb1) == length(size) & all(vb1 >= 0)");
        reqp(i, "wvb1",      "is.numeric(wvb1) & length(wvb1) == length(size) & all(wvb1 >= 0)");
        reqp(i, "ei_vb1",    "is.numeric(ei_vb1) & length(ei_vb1) == length(size) & all(ei_vb1 >= 0) & all(ei_vb1 <= 1)");
        reqp(i, "ei2_vb1",   "is.numeric(ei2_vb1) & length(ei2_vb1) == length(size) & all(ei2_vb1 >= 0) & all(ei2_vb1 <= 1)");
        reqp(i, "ei3_vb1",   "is.numeric(ei3_vb1) & length(ei3_vb1) == length(size) & all(ei3_vb1 >= 0) & all(ei3_vb1 <= 1)");
        reqp(i, "ed_vb1i",   "is.numeric(ed_vb1i) & length(ed_vb1i) == length(size) & all(ed_vb1i >= 0) & all(ed_vb1i <= 1)");
        reqp(i, "ed_vb1i2",  "is.numeric(ed_vb1i2) & length(ed_vb1i2) == length(size) & all(ed_vb1i2 >= 0) & all(ed_vb1i2 <= 1)");
        reqp(i, "ed_vb1i3",  "is.numeric(ed_vb1i3) & length(ed_vb1i3) == length(size) & all(ed_vb1i3 >= 0) & all(ed_vb1i3 <= 1)");
        reqp(i, "wvb2",      "is.numeric(wvb2) & length(wvb2) == length(size) & all(wvb2 >= 0)");
        reqp(i, "ei_vb2",    "is.numeric(ei_vb2) & length(ei_vb2) == length(size) & all(ei_vb2 >= 0) & all(ei_vb2 <= 1)");
        reqp(i, "ei2_vb2",   "is.numeric(ei2_vb2) & length(ei2_vb2) == length(size) & all(ei2_vb2 >= 0) & all(ei2_vb2 <= 1)");
        reqp(i, "ei3_vb2",   "is.numeric(ei3_vb2) & length(ei3_vb2) == length(size) & all(ei3_vb2 >= 0) & all(ei3_vb2 <= 1)");
        reqp(i, "ed_vb2i",   "is.numeric(ed_vb2i) & length(ed_vb2i) == length(size) & all(ed_vb2i >= 0) & all(ed_vb2i <= 1)");
        reqp(i, "ed_vb2i2",  "is.numeric(ed_vb2i2) & length(ed_vb2i2) == length(size) & all(ed_vb2i2 >= 0) & all(ed_vb2i2 <= 1)");
        reqp(i, "ed_vb2i3",  "is.numeric(ed_vb2i3) & length(ed_vb2i3) == length(size) & all(ed_vb2i3 >= 0) & all(ed_vb2i3 <= 1)");
        reqp(i, "wvb3",      "is.numeric(wvb3) & length(wvb3) == length(size) & all(wvb3 >= 0)");
        reqp(i, "ei_vb3",    "is.numeric(ei_vb3) & length(ei_vb3) == length(size) & all(ei_vb3 >= 0) & all(ei_vb3 <= 1)");
        reqp(i, "ei2_vb3",   "is.numeric(ei2_vb3) & length(ei2_vb3) == length(size) & all(ei2_vb3 >= 0) & all(ei2_vb3 <= 1)");
        reqp(i, "ei3_vb3",   "is.numeric(ei3_vb3) & length(ei3_vb3) == length(size) & all(ei3_vb3 >= 0) & all(ei3_vb3 <= 1)");
        reqp(i, "ed_vb3i",   "is.numeric(ed_vb3i) & length(ed_vb3i) == length(size) & all(ed_vb3i >= 0) & all(ed_vb3i <= 1)");
        reqp(i, "ed_vb3i2",  "is.numeric(ed_vb3i2) & length(ed_vb3i2) == length(size) & all(ed_vb3i2 >= 0) & all(ed_vb3i2 <= 1)");
        reqp(i, "ed_vb3i3",  "is.numeric(ed_vb3i3) & length(ed_vb3i3) == length(size) & all(ed_vb3i3 >= 0) & all(ed_vb3i3 <= 1)");
        reqp(i, "A",       "is.numeric(A) & length(A) == length(size) & all(A >= 0)");
        reqp(i, "B",       "is.numeric(B) & length(B) == length(size) & all(B >= 0)");
        reqp(i, "D",       "is.numeric(D) & length(D) == length(size) & all(D >= 0)");

        reqp(i, "season_A",   "is.numeric(season_A) & length(season_A) == 1 & all(abs(season_A <= 1))");
        reqp(i, "season_T",   "is.numeric(season_T) & length(season_T) == 1 & all(season_T > 0)");
        reqp(i, "season_phi", "is.numeric(season_phi) & length(season_phi) == 1");

        reqp(i, "seed_times",      "is.numeric(seed_times) & !is.unsorted(seed_times)");
        reqp(i, "seed_times2",     "is.numeric(seed_times2) & !is.unsorted(seed_times2)");
        reqp(i, "seed_times3",     "is.numeric(seed_times3) & !is.unsorted(seed_times3)");
        reqp(i, "dist_seed_ages",  "is.numeric(dist_seed_ages) & length(dist_seed_ages) == length(size)");

        if (!is.null(parameters$pop[[i]]$observer)) {
            if (!(is.function(parameters$pop[[i]]$observer) & length(formals(parameters$pop[[i]]$observer) == 4))) {
                stop(paste0("observer has to be either NULL or a function taking 4 arguments, but is not in population", i));
            }
        }
        reqp(i, "schedule", "is.list(schedule)");
        schedule_times = sapply(parameters$pop[[i]]$schedule, function(x) x$t);
        if (is.unsorted(schedule_times)) {
            stop(paste0("elements t of schedule need to be ordered, but are not in population ", i));
        }
    }
}

# Get demographics for a given location, with error checking.
cm_get_demographics = function(dem_location, n_groups = NULL)
{
    # Get demographics.
    demographics = cm_populations[name == dem_location];
    if (nrow(demographics) == 0) {
        message(paste0("Could not find demographics for dem_location ", dem_location, "."));
        answer = readline(prompt = "View options? (y/n) ");
        if (answer != "n") {
            print(cm_populations[, unique(name)]);
        }
        stop();
    }
    if (nrow(demographics) != demographics[, uniqueN(age)]) {
        stop(paste0("Age not unique in cm_population[name == \"", dem_location, "\"]. This means cm_populations is misspecified for this location."));
    }

    # Adjust number of age groups if needed.
    if (!is.null(n_groups)) {
        if (n_groups > nrow(demographics)) {
            stop(sprintf("Requested %d age groups for model (up to %d+), but demographic data only goes up to %s.",
                         n_groups, (n_groups - 1) * 5, demographics[.N, age]));
        } else if (n_groups < nrow(demographics)) {
            demographics[n_groups]$f = demographics[n_groups:.N, sum(f)];
            demographics[n_groups]$m = demographics[n_groups:.N, sum(m)];
            demographics[n_groups]$age = demographics[n_groups, sub("-[0-9]+", "\\+", age)];
            demographics = demographics[1:n_groups];
        }
    }
    return (demographics);
}

# Get matrices for a given location, with error checking.
cm_get_matrices = function(mat_location, dem_location = NULL)
{
    # If requested, guess mat_location from dem_location.
    guess = F;
    if (mat_location == "guess") {
        if (is.null(dem_location)) {
            stop("cm_get_matrices needs a dem_location to guess the matrices location from.");
        }
        mat_location = dem_location;
        guess = T;
    }

    # Try to find matrices for mat_location.
    if (!mat_location %in% names(cm_matrices)) {
        message(paste0("Could not find matrices for mat_location ", mat_location, "."));
        answer = readline(prompt = "View options? (y/n) ");
        if (answer != "n") {
            print(names(cm_matrices));
        }
        stop();
    }
    if (sum(mat_location %in% names(cm_matrices) > 1)) {
        stop(paste0("Duplicate entries for ", mat_location, " in cm_matrices. This means cm_matrices is misspecified."));
    }
    mat = cm_matrices[[mat_location]];
    if (is.list(mat) & length(mat) > 0 & all(sapply(mat, is.matrix))) {
        return (mat);
    } else {
        if (guess) {
            stop(paste0("Could not guess mat_location for dem_location ", dem_location, "."));
        }
        stop(paste0("No valid entry in cm_matrices for matrix location ", mat_location, "."));
    }
}

# Split matrices
# ex_in: bounds separate into contacts *ex*clusively between lower groups, and *in*clusively between higher groups.
cm_split_matrices_ex_in = function(parameters, bounds)
{
    for (pi in seq_along(parameters$pop))
    {
        ng = nrow(parameters$pop[[pi]]$matrices[[1]]);
        if (any(bounds < 1 | bounds > ng)) {
            stop("Bounds must lie within [1, nrow(mat)] for splitting contact matrices.");
        }
        nmat0 = length(parameters$pop[[pi]]$matrices);
        parameters$pop[[pi]]$matrices = rep(parameters$pop[[pi]]$matrices, length(bounds) + 1);

        for (b in seq_along(bounds))
        {
            lb = floor(bounds[b]);
            fb = bounds[b] %% 1;
            mask1 = matrix(1, nrow = ng, ncol = ng);
            if (lb > 1) {
                mask1[1:(lb - 1), 1:(lb - 1)] = 0;
            }
            mask1[lb, 1:(lb-1)] = 1 - fb;
            mask1[1:(lb-1), lb] = 1 - fb;
            mask1[lb, lb] = (1 - fb)^2;
            mask0 = 1 - mask1;

            for (m in seq_len(nmat0)) {
                names(parameters$pop[[pi]]$matrices)[m + b * nmat0] = paste0(names(parameters$pop[[pi]]$matrices)[m + b * nmat0], b + 1);
                parameters$pop[[pi]]$matrices[[m + (b - 1) * nmat0]] = mask0 * parameters$pop[[pi]]$matrices[[m + (b - 1) * nmat0]];
                parameters$pop[[pi]]$matrices[[m +       b * nmat0]] = mask1 * parameters$pop[[pi]]$matrices[[m +       b * nmat0]];
            }
        }
        parameters$pop[[pi]]$contact = rep_len(parameters$pop[[pi]]$contact, nmat0 * (length(bounds) + 1));
    }

    return (parameters)
}

# TODO
# cm_split_matrices_in_ex
# cm_split_matrices_custom

# Get default population parameters, SEI3R model
cm_base_pop_SEI3R = function(n_groups)
{
    warning("seed_times2 set to 999999 by default as empty numeric vector causes issues. Refactor?")
    list(
        dE  = cm_delay_gamma(4.0, 4.0, t_max = 60, t_step = 0.25)$p, # Derived from Backer et al Eurosurveillance
        dIp = cm_delay_gamma(2.4, 4.0, t_max = 60, t_step = 0.25)$p, # Derived from Backer et al Eurosurveillance
        dIa = cm_delay_gamma(7.0, 4.0, t_max = 60, t_step = 0.25)$p, # Assumed 7 days subclinical shedding
        dIs = cm_delay_gamma(3.2, 3.7, t_max = 60, t_step = 0.25)$p, # Zhang et al 2020
        dE2 = numeric(),
        dIp2 = numeric(),
        dIa2 = numeric(),
        dIs2 = numeric(),
        dE3 = numeric(),
        dIp3 = numeric(),
        dIa3 = numeric(),
        dIs3 = numeric(),
        dVa1 = cm_delay_gamma(21, 100, t_max = 60, t_step = 0.25)$p,
        dVb1 = cm_delay_gamma(21, 100, t_max = 60, t_step = 0.25)$p,
        dVa2 = cm_delay_gamma(222, 100, t_max = 60, t_step = 0.25)$p, # 268 days between 8th December 2020 and 1st September 2021
        dVb2 = cm_delay_gamma(222, 100, t_max = 60, t_step = 0.25)$p, # subtract average first to second dose delay (~60 days), add 14 days for efficacy to kick in

        size = rep(1000, n_groups),
        imm0 = rep(0, n_groups),
        matrices = list(base = diag(n_groups) * 0.5 + 0.5/n_groups),
        contact = 1,
        contact_mult = numeric(),
        contact_lowerto = numeric(),
        u = rep(0.08, n_groups),
        u2 = rep(0.08, n_groups),
        u3 = rep(0.08, n_groups),
        y = rep(0.5, n_groups),
        y2 = rep(0.5, n_groups),
        y3 = rep(0.5, n_groups),
        fIp = rep(1, n_groups),
        fIs = rep(1, n_groups),
        fIa = rep(0.5, n_groups),
        omega = rep(0, n_groups),
        tau = rep(1, n_groups),
        pi_r = rep(1, n_groups),
        pi_r2 = rep(1, n_groups),
        pi_r3 = rep(1, n_groups),
        pi2_r = rep(1, n_groups),
        pi2_r2 = rep(1, n_groups),
        pi2_r3 = rep(1, n_groups),
        pi3_r = rep(1, n_groups),
        pi3_r2 = rep(1, n_groups),
        pi3_r3 = rep(1, n_groups),
        wn = rep(0, n_groups),
        wn2 = rep(0, n_groups),
        wn3 = rep(0, n_groups),
        va1 = rep(0, n_groups),
        wva1 = rep(0, n_groups),
        ei_va1 = rep(1, n_groups),
        ei2_va1 = rep(1, n_groups),
        ei3_va1 = rep(1, n_groups),
        ed_va1i = rep(0, n_groups),
        ed_va1i2 = rep(0, n_groups),
        ed_va1i3 = rep(0, n_groups),
        wva2 = rep(0, n_groups),
        ei_va2 = rep(1, n_groups),
        ei2_va2 = rep(1, n_groups),
        ei3_va2 = rep(1, n_groups),
        ed_va2i = rep(0, n_groups),
        ed_va2i2 = rep(0, n_groups),
        ed_va2i3 = rep(0, n_groups),
        wva3 = rep(0, n_groups),
        ei_va3 = rep(1, n_groups),
        ei2_va3 = rep(1, n_groups),
        ei3_va3 = rep(1, n_groups),
        ed_va3i = rep(0, n_groups),
        ed_va3i2 = rep(0, n_groups),
        ed_va3i3 = rep(0, n_groups),
        vb1 = rep(0, n_groups),
        wvb1 = rep(0, n_groups),
        ei_vb1 = rep(1, n_groups),
        ei2_vb1 = rep(1, n_groups),
        ei3_vb1 = rep(1, n_groups),
        ed_vb1i = rep(0, n_groups),
        ed_vb1i2 = rep(0, n_groups),
        ed_vb1i3 = rep(0, n_groups),
        wvb2 = rep(0, n_groups),
        ei_vb2 = rep(1, n_groups),
        ei2_vb2 = rep(1, n_groups),
        ei3_vb2 = rep(1, n_groups),
        ed_vb2i = rep(0, n_groups),
        ed_vb2i2 = rep(0, n_groups),
        ed_vb2i3 = rep(0, n_groups),
        wvb3 = rep(0, n_groups),
        ei_vb3 = rep(1, n_groups),
        ei2_vb3 = rep(1, n_groups),
        ei3_vb3 = rep(1, n_groups),
        ed_vb3i = rep(0, n_groups),
        ed_vb3i2 = rep(0, n_groups),
        ed_vb3i3 = rep(0, n_groups),
        A = rep(0, n_groups),
        B = rep(0, n_groups),
        D = rep(0, n_groups),
        season_A = 0,
        season_T = 365.25,
        season_phi = 0,

        seed_times = 1,
        seed_times2 = 999999,
        seed_times3 = 999999,
        dist_seed_ages = rep(1, n_groups),

        schedule = list(),
        observer = NULL
    )
}

# Build parameters for a single location, SEI3R model
cm_build_pop_SEI3R = function(dem_location, mat_location = "guess",
    dE = NULL, dIp = NULL, dIs = NULL, dIa = NULL, 
    dE2 = NULL, dIp2 = NULL, dIs2 = NULL, dIa2 = NULL,
    dE3 = NULL, dIp3 = NULL, dIs3 = NULL, dIa3 = NULL, 
    dVa1 = NULL, dVb1 = NULL, dVa2 = NULL, dVb2 = NULL,
    contact = NULL, imm0 = NULL, 
    u = NULL, u2 = NULL, u3 = NULL,
    y = NULL, y2 = NULL, y3 = NULL,
    fIp = NULL, fIa = NULL, fIs = NULL, omega = NULL, tau = NULL,
    pi_r = NULL, pi_r2 = NULL, pi_r3 = NULL,
    pi2_r = NULL, pi2_r2 = NULL, pi2_r3 = NULL,
    pi3_r = NULL, pi3_r2 = NULL, pi3_r3 = NULL,
    wn = NULL, wn2 = NULL, wn3 = NULL,
    va1 = NULL, wva1 = NULL, 
    ei_va1 = NULL, ei2_va1 = NULL, ei3_va1 = NULL, 
    ed_va1i = NULL, ed_va1i2 = NULL, ed_va1i3 = NULL,
                wva2 = NULL, 
    ei_va2 = NULL, ei2_va2 = NULL, ei3_va2 = NULL,
    ed_va2i = NULL, ed_va2i2 = NULL, ed_va2i3 = NULL,
                wva3 = NULL, 
    ei_va3 = NULL, ei2_va3 = NULL, ei3_va3 = NULL,
    ed_va3i = NULL, ed_va3i2 = NULL, ed_va3i3 = NULL,
    vb1 = NULL, wvb1 = NULL, 
    ei_vb1 = NULL, ei2_vb1 = NULL, ei3_vb1 = NULL, 
    ed_vb1i = NULL, ed_vb1i2 = NULL, ed_vb1i3 = NULL,
                wvb2 = NULL, 
    ei_vb2 = NULL, ei2_vb2 = NULL, ei3_vb2 = NULL, 
    ed_vb2i = NULL, ed_vb2i2 = NULL, ed_vb2i3 = NULL,
                wvb3 = NULL, 
    ei_vb3 = NULL, ei2_vb3 = NULL, ei3_vb3 = NULL, 
    ed_vb3i = NULL, ed_vb3i2 = NULL, ed_vb3i3 = NULL,
    A = NULL, B = NULL, D = NULL, 
    season_A = NULL, season_T = NULL, season_phi = NULL,
    seed_times = NULL, seed_times2 = NULL, seed_times3 = NULL,
    dist_seed_ages = NULL, observer = NULL, schedule = NULL)
{
    # Get desired demographics and matrices.
    matrices = cm_get_matrices(mat_location, dem_location);
    n_groups = nrow(matrices[[1]]);
    demographics = cm_get_demographics(dem_location, n_groups);

    # Get base population parameters.
    pop = cm_base_pop_SEI3R(n_groups);

    # Set population parameters.
    assign = function(pop, name, value) {
        if (!is.null(value)) {
            pop[[name]] = value;
        }
        return (pop);
    }

    assign_g = function(pop, name, value, n_groups) {
        if (!is.null(value)) {
            if (length(value) == 1 | length(value) == n_groups) {
                pop[[name]] = rep_len(value, n_groups)
            } else {
                stop(paste0("Parameter ", name, " must be either length 1 or length n_groups = ", n_groups, "."));
            }
        }
        return (pop);
    }

    pop$name = dem_location;
    pop$group_names = colnames(matrices[[1]]);

    pop = assign(pop, "dE", dE);
    pop = assign(pop, "dIp", dIp);
    pop = assign(pop, "dIs", dIs);
    pop = assign(pop, "dIa", dIa);
    pop = assign(pop, "dE2", dE2);
    pop = assign(pop, "dIp2", dIp2);
    pop = assign(pop, "dIs2", dIs2);
    pop = assign(pop, "dIa2", dIa2);
    pop = assign(pop, "dE3", dE3);
    pop = assign(pop, "dIp3", dIp3);
    pop = assign(pop, "dIs3", dIs3);
    pop = assign(pop, "dIa3", dIa3);
    pop = assign(pop, "dVa1", dVa1);
    pop = assign(pop, "dVb1", dVb1);
    pop = assign(pop, "dVa2", dVa2);
    pop = assign(pop, "dVb2", dVb2);
    pop$size = demographics[, round((f + m) * 1000)];
    pop$matrices = matrices;
    pop$contact = rep(1, length(matrices));
    pop = assign(pop, "contact", contact);

    pop = assign_g(pop, "imm0", imm0, n_groups);
    pop = assign_g(pop, "u", u, n_groups);
    pop = assign_g(pop, "u2", u2, n_groups);
    pop = assign_g(pop, "u3", u3, n_groups);
    pop = assign_g(pop, "y", y, n_groups);
    pop = assign_g(pop, "y2", y2, n_groups);
    pop = assign_g(pop, "y3", y3, n_groups);
    pop = assign_g(pop, "fIp", fIp, n_groups);
    pop = assign_g(pop, "fIs", fIs, n_groups);
    pop = assign_g(pop, "fIa", fIa, n_groups);
    pop = assign_g(pop, "omega", omega, n_groups);
    pop = assign_g(pop, "tau", tau, n_groups);
    pop = assign_g(pop, "pi_r", pi_r, n_groups);
    pop = assign_g(pop, "pi_r2", pi_r2, n_groups);
    pop = assign_g(pop, "pi_r3", pi_r3, n_groups);
    pop = assign_g(pop, "pi2_r", pi2_r, n_groups);
    pop = assign_g(pop, "pi2_r2", pi2_r2, n_groups);
    pop = assign_g(pop, "pi2_r3", pi2_r3, n_groups);
    pop = assign_g(pop, "pi3_r", pi3_r, n_groups);
    pop = assign_g(pop, "pi3_r2", pi3_r2, n_groups);
    pop = assign_g(pop, "pi3_r3", pi3_r3, n_groups);
    pop = assign_g(pop, "wn", wn, n_groups);
    pop = assign_g(pop, "wn2", wn2, n_groups);
    pop = assign_g(pop, "wn3", wn3, n_groups);
    pop = assign_g(pop, "va1", va1, n_groups);
    pop = assign_g(pop, "wva1", wva1, n_groups);
    pop = assign_g(pop, "ei_va1", ei_va1, n_groups);
    pop = assign_g(pop, "ei2_va1", ei2_va1, n_groups);
    pop = assign_g(pop, "ei3_va1", ei3_va1, n_groups);
    pop = assign_g(pop, "ed_va1i", ed_va1i, n_groups);
    pop = assign_g(pop, "ed_va1i2", ed_va1i2, n_groups);
    pop = assign_g(pop, "ed_va1i3", ed_va1i3, n_groups);
    pop = assign_g(pop, "wva2", wva2, n_groups);
    pop = assign_g(pop, "ei_va2", ei_va2, n_groups);
    pop = assign_g(pop, "ei2_va2", ei2_va2, n_groups);
    pop = assign_g(pop, "ei3_va2", ei3_va2, n_groups);
    pop = assign_g(pop, "ed_va2i", ed_va2i, n_groups);
    pop = assign_g(pop, "ed_va2i2", ed_va2i2, n_groups);
    pop = assign_g(pop, "ed_va2i3", ed_va2i3, n_groups);
    pop = assign_g(pop, "wva3", wva3, n_groups);
    pop = assign_g(pop, "ei_va3", ei_va3, n_groups);
    pop = assign_g(pop, "ei2_va3", ei2_va3, n_groups);
    pop = assign_g(pop, "ei3_va3", ei3_va3, n_groups);
    pop = assign_g(pop, "ed_va3i", ed_va3i, n_groups);
    pop = assign_g(pop, "ed_va3i2", ed_va3i2, n_groups);
    pop = assign_g(pop, "ed_va3i3", ed_va3i3, n_groups);
    pop = assign_g(pop, "vb1", vb1, n_groups);
    pop = assign_g(pop, "wvb1", wvb1, n_groups);
    pop = assign_g(pop, "ei_vb1", ei_vb1, n_groups);
    pop = assign_g(pop, "ei2_vb1", ei2_vb1, n_groups);
    pop = assign_g(pop, "ei3_vb1", ei3_vb1, n_groups);
    pop = assign_g(pop, "ed_vb1i", ed_vb1i, n_groups);
    pop = assign_g(pop, "ed_vb1i2", ed_vb1i2, n_groups);
    pop = assign_g(pop, "ed_vb1i3", ed_vb1i3, n_groups);
    pop = assign_g(pop, "wvb2", wvb2, n_groups);
    pop = assign_g(pop, "ei_vb2", ei_vb2, n_groups);
    pop = assign_g(pop, "ei2_vb2", ei2_vb2, n_groups);
    pop = assign_g(pop, "ei3_vb2", ei3_vb2, n_groups);
    pop = assign_g(pop, "ed_vb2i", ed_vb2i, n_groups);
    pop = assign_g(pop, "ed_vb2i2", ed_vb2i2, n_groups);
    pop = assign_g(pop, "ed_vb2i3", ed_vb2i3, n_groups);
    pop = assign_g(pop, "wvb3", wvb3, n_groups);
    pop = assign_g(pop, "ei_vb3", ei_vb3, n_groups);
    pop = assign_g(pop, "ei2_vb3", ei2_vb3, n_groups);
    pop = assign_g(pop, "ei3_vb3", ei3_vb3, n_groups);
    pop = assign_g(pop, "ed_vb3i", ed_vb3i, n_groups);
    pop = assign_g(pop, "ed_vb3i2", ed_vb3i2, n_groups);
    pop = assign_g(pop, "ed_vb3i3", ed_vb3i3, n_groups);
    pop = assign_g(pop, "A", A, n_groups);
    pop = assign_g(pop, "B", B, n_groups);
    pop = assign_g(pop, "D", D, n_groups);
    pop = assign(pop, "season_A", season_A);
    pop = assign(pop, "season_T", season_T);
    pop = assign(pop, "season_phi", season_phi);

    pop = assign(pop, "seed_times", seed_times);
    pop = assign(pop, "seed_times2", seed_times2);
    pop = assign(pop, "seed_times3", seed_times3);
    pop = assign(pop, "dist_seed_ages", dist_seed_ages);
    pop = assign(pop, "observer", observer);
    pop = assign(pop, "schedule", schedule);

    return (pop)
}

# Get default simulation parameters, SEI3R model
cm_base_parameters_SEI3R = function(n_groups = 1, pop = cm_base_pop_SEI3R(n_groups))
{
    # If just a single population, rather than a list of populations, has been passed to this function, rectify that.
    if (is.character(pop$type)) {
        pop = list(pop);
    }

    list(
        model = "SEI3R",
        time_step = 0.25,
        date0 = "2020-01-01",
        time0 = 0,
        time1 = 365,
        report_every = 4,
        fast_multinomial = F,
        deterministic = T,
        pop = pop,
        travel = diag(length(pop)),
        processes = NULL
    )
}

# Build parameters for one or several locations, SEI3R model
cm_parameters_SEI3R = function(dem_locations, mat_locations = "guess", 
                               date_start = "2020-03-01", 
                               date_end = "2021-03-01", 
                               deterministic = T, processes = NULL,
    dE = NULL, dIp = NULL, dIs = NULL, dIa = NULL, 
    dE2 = NULL, dIp2 = NULL, dIs2 = NULL, dIa2 = NULL, 
    dE3 = NULL, dIp3 = NULL, dIs3 = NULL, dIa3 = NULL, 
    dVa1 = NULL, dVb1 = NULL, dVa2 = NULL, dVb2 = NULL,
    contact = NULL, imm0 = NULL,
    u = NULL, u2 = NULL, u3 = NULL, 
    y = NULL, y2 = NULL, y3 = NULL,
    fIp = NULL, fIa = NULL, fIs = NULL, omega = NULL, tau = NULL,
    pi_r = NULL, pi_r2 = NULL, pi_r3 = NULL,
    pi2_r = NULL, pi2_r2 = NULL, pi2_r3 = NULL,
    pi3_r = NULL, pi3_r2 = NULL, pi3_r3 = NULL,
    wn = NULL, wn2 = NULL, wn3 = NULL,
    va1 = NULL, wva1 = NULL, 
    ei_va1 = NULL, ei2_va1 = NULL, ei3_va1 = NULL,
    ed_va1i = NULL, ed_va1i2 = NULL, ed_va1i3 = NULL,
                wva2 = NULL, 
    ei_va2 = NULL, ei2_va2 = NULL, ei3_va2 = NULL,
    ed_va2i = NULL, ed_va2i2 = NULL, ed_va2i3 = NULL,
                wva3 = NULL, 
    ei_va3 = NULL, ei2_va3 = NULL, ei3_va3 = NULL,
    ed_va3i = NULL, ed_va3i2 = NULL, ed_va3i3 = NULL,
    vb1 = NULL, wvb1 = NULL, 
    ei_vb1 = NULL, ei2_vb1 = NULL, ei3_vb1 = NULL,
    ed_vb1i = NULL, ed_vb1i2 = NULL, ed_vb1i3 = NULL,
                wvb2 = NULL, 
    ei_vb2 = NULL, ei2_vb2 = NULL, ei3_vb2 = NULL,
    ed_vb2i = NULL, ed_vb2i2 = NULL, ed_vb2i3 = NULL,
                wvb3 = NULL, 
    ei_vb3 = NULL, ei2_vb3 = NULL, ei3_vb3 = NULL,
    ed_vb3i = NULL, ed_vb3i2 = NULL, ed_vb3i3 = NULL,
    A = NULL, B = NULL, D = NULL, 
    season_A = NULL, season_T = NULL, season_phi = NULL,
    seed_times = NULL, seed_times2 = NULL, seed_times3 = NULL,
    dist_seed_ages = NULL, observer = NULL, schedule = NULL)
{
    # Check parameters
    if (length(mat_locations) != length(dem_locations)) {
        if (length(mat_locations) == 1) {
            mat_locations = rep_len(mat_locations, length(dem_locations));
        } else {
            stop("dem_locations and mat_locations must have the same length; or mat_locations can have length 1 and will be used for all dem_locations.");
        }
    }

    # Get population parameters.
    pop = list();
    for (i in seq_along(dem_locations))
    {
        pop[[i]] = cm_build_pop_SEI3R(dem_locations[i], mat_locations[i],
            dE = dE, dIp = dIp, dIs = dIs, dIa = dIa, 
            dE2 = dE2, dIp2 = dIp2, dIs2 = dIs2, dIa2 = dIa2, 
            dE3 = dE3, dIp3 = dIp3, dIs3 = dIs3, dIa3 = dIa3, 
            dVa1 = dVa1, dVb1 = dVb1, dVa2 = dVa2, dVb2 = dVb2,
            contact = contact, imm0 = imm0,
            u = u, u2 = u2, u3 = u3,
            y = y, y2 = y2, y3 = y3,
            fIp = fIp, fIa = fIa, fIs = fIs, omega = omega, tau = tau,
            pi_r = pi_r, pi_r2 = pi_r2, pi_r3 = pi_r3,
            pi2_r = pi2_r, pi2_r2 = pi2_r2, pi2_r3 = pi2_r3,
            pi3_r = pi3_r, pi3_r2 = pi3_r2, pi3_r3 = pi3_r3,
            wn = wn, wn2 = wn2, wn3 = wn3,
            va1 = va1, wva1 = wva1, 
            ei_va1 = ei_va1, ei2_va1 = ei2_va1, ei3_va1 = ei3_va1,
            ed_va1i = ed_va1i, ed_va1i2 = ed_va1i2, ed_va1i3 = ed_va1i3,
                       wva2 = wva2, 
            ei_va2 = ei_va2, ei2_va2 = ei2_va2, ei3_va2 = ei3_va2,
            ed_va2i = ed_va2i, ed_va2i2 = ed_va2i2, ed_va2i3 = ed_va2i3,
                       wva3 = wva3, 
            ei_va3 = ei_va3, ei2_va3 = ei2_va3, ei3_va3 = ei3_va3,
            ed_va3i = ed_va3i, ed_va3i2 = ed_va3i2, ed_va3i3 = ed_va3i3,
            vb1 = vb1, wvb1 = wvb1, 
            ei_vb1 = ei_vb1, ei2_vb1 = ei2_vb1, ei3_vb1 = ei3_vb1,
            ed_vb1i = ed_vb1i, ed_vb1i2 = ed_vb1i2, ed_vb1i3 = ed_vb1i3,
                       wvb2 = wvb2, 
            ei_vb2 = ei_vb2, ei2_vb2 = ei2_vb2, ei3_vb2 = ei3_vb2,
            ed_vb2i = ed_vb2i, ed_vb2i2 = ed_vb2i2, ed_vb2i3 = ed_vb2i3,
                       wvb3 = wvb3, 
            ei_vb3 = ei_vb3, ei2_vb3 = ei2_vb3, ei3_vb3 = ei3_vb3,
            ed_vb3i = ed_vb3i, ed_vb3i2 = ed_vb3i2, ed_vb3i3 = ed_vb3i3,
            A = A, B = B, D = D, 
            season_A = season_A, season_T = season_T, season_phi = season_phi,
            seed_times = seed_times, seed_times2 = seed_times2, seed_times3 = seed_times3,
            dist_seed_ages = dist_seed_ages, observer = observer, schedule = schedule);
    }

    # Build simulation parameters around this population.
    parameters = cm_base_parameters_SEI3R(n_groups, pop);
    parameters$date0 = date_start;
    parameters$time0 = 0;
    parameters$time1 = date_end;
    parameters$deterministic = deterministic;
    parameters$processes = processes;

    return (parameters)
}

# Get regions for the UK.
cm_uk_locations = function(country, level) {
    # Check country code
    country = toupper(country);
    if (country == "UK") {
        country = "EWSN";
    }
    if (!country %like% "^(UK|[EWSN]+)$") {
        stop("country must be UK, or a combination of E, W, S, and/or N.");
    }

    # Interpret level
    level = as.integer(level);
    if (level < 0 | level > 4) {
        stop("level must be 0, 1, 2, 3, or 4");
    }

    if (level == 0) {
        if (country != "EWSN") {
            stop("For level 0, country must be UK.");
        }
        return ("UK | UNITED KINGDOM");
    } else if (level == 1) {
        gE = "Country";
        gW = "Country";
        gS = "Country";
        gN = "Country";
    } else if (level == 2) {
        gE = "Region";
        gW = "Country";
        gS = "Country";
        gN = "Country";
    } else if (level == 3) {
        gE = c("Metropolitan County", "County", "Unitary Authority", "London Borough");
        gW = "Unitary Authority";
        gS = "Council Area";
        gN = "Local Government District";
    } else if (level == 4) {
        gE = c("Metropolitan District", "Non-metropolitan District", "Unitary Authority", "London Borough");
        gW = "Unitary Authority";
        gS = "Council Area";
        gN = "Local Government District";
    }

    # Extract locations
    locs = NULL;
    if (country %like% "E") { locs = c(locs, cm_structure_UK[Code %like% "^E" & Geography1 %in% gE, Name]); }
    if (country %like% "W") { locs = c(locs, cm_structure_UK[Code %like% "^W" & Geography1 %in% gW, Name]); }
    if (country %like% "S") { locs = c(locs, cm_structure_UK[Code %like% "^S" & Geography1 %in% gS, Name]); }
    if (country %like% "N") { locs = c(locs, cm_structure_UK[Code %like% "^N" & Geography1 %in% gN, Name]); }
    return (paste0("UK | ", locs));
}

# get prevalence of morbidities that increase risk of Covid per age-group
# by default, returns prevalence of every included morbidity per age-group
# to return overall prevalence, set aggregate = TRUE
#  to account for mutli-morbidities, set correlation of people with multiple morbidities multimorbidity_corr
cm_high_risk_prevalence <- function(country, aggregate = FALSE, multimorbidity_corr = 0.8){
    cnt <- country
    country_data <- cm_highrisk[country == cnt]
    if(aggregate){
        country_data <- country_data[, .(maxprev = max(prevalence, na.rm = T), cmbprev = 1-prod(1-prevalence, na.rm = T)), by=c("country", "age_from", "age_to")]
        country_data <- country_data[, .(highrisk = maxprev+((cmbprev-maxprev)*(1-multimorbidity_corr))), by=c("country", "age_from", "age_to")]
    }
    return(country_data)
}
