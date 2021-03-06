#
# PLOT FITS.
#

gen_fit = function(test, parametersI, ld, sitreps, virus, sero, populations, sero_cut_off)
{
    test = copy(test);
    test = test[population %in% populations];
    sero = copy(sero);
    sero = sero[NHS.region %in% populations];
    sero = sero[sero$Start.date < sero_cut_off]
    virus = copy(virus);
    virus = virus[NHS.region %in% populations];
    ld = copy(ld);
    ld = ld[name %in% populations];
    sitreps = copy(sitreps);
    sitreps = sitreps[name %in% populations]
    
    # Calculate total population
    popsize = NULL
    for (i in seq_along(parametersI)) {
        if (!is.null(parametersI[[i]])) {
            popsize = rbind(popsize,
                data.table(Geography = parametersI[[i]]$pop[[1]]$name, 
                    population_size = sum(parametersI[[i]]$pop[[1]]$size),
                    population_size_15plus = sum(parametersI[[i]]$pop[[1]]$size[-(1:3)]))
            )
        }
    }
    
    # Create formatted output
    output = SPIM_output_full(test, 2020, 11, 1, "2020-01-01", as.character(ymd("2020-01-01") + max(test$t)))
    output[, d := make_date(`Year of Value`, `Month of Value`, `Day of Value`)]
    output = merge(output, popsize, by = "Geography")
    
    adj_output = function(output, val_type, div, pop = 0, pop15 = 0) {
        output[ValueType == val_type, Value := Value / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.05` := `Quantile 0.05` / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.1`  := `Quantile 0.1`  / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.15` := `Quantile 0.15` / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.2`  := `Quantile 0.2`  / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.25` := `Quantile 0.25` / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.3`  := `Quantile 0.3`  / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.35` := `Quantile 0.35` / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.4`  := `Quantile 0.4`  / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.45` := `Quantile 0.45` / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.5`  := `Quantile 0.5`  / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.55` := `Quantile 0.55` / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.6`  := `Quantile 0.6`  / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.65` := `Quantile 0.65` / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.7`  := `Quantile 0.7`  / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.75` := `Quantile 0.75` / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.8`  := `Quantile 0.8`  / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.85` := `Quantile 0.85` / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.9`  := `Quantile 0.9`  / (div + population_size * pop + population_size_15plus * pop15)]
        output[ValueType == val_type, `Quantile 0.95` := `Quantile 0.95` / (div + population_size * pop + population_size_15plus * pop15)]
    }
    
    adj_output(output, "hospital_inc", 1)
    adj_output(output, "hospital_prev", 1)
    adj_output(output, "icu_prev", 1)
    adj_output(output, "prevalence_mtp", 0, 0.01)
    adj_output(output, "sero_prev", 0, 0, 0.01)
    adj_output(output, "type28_death_inc_line", 1)
    
    # Make data to output
    data = make_data(ld, sitreps, virus, sero)
    data = merge(data, popsize, by = "Geography")
    
    adj_data = function(data, val_type, div, pop = 0) {
        data[ValueType == val_type, ymin := ymin / (div + population_size * pop)]
        data[ValueType == val_type, y    := y    / (div + population_size * pop)]
        data[ValueType == val_type, ymax := ymax / (div + population_size * pop)]
    }
    
    adj_data(data, "hospital_inc", 1)
    adj_data(data, "hospital_prev", 1)
    adj_data(data, "icu_prev", 1)
    adj_data(data, "prevalence_mtp", 0.01)
    adj_data(data, "sero_prev", 0.01)
    adj_data(data, "type28_death_inc_line", 1)
    
    output[ValueType == "hospital_inc", ValueType := "Hospital\nadmissions"]
    output[ValueType == "hospital_prev", ValueType := "Hospital beds\noccupied"]
    output[ValueType == "icu_prev", ValueType := "ICU beds\noccupied"]
    output[ValueType == "infections_inc", ValueType := "Infection\nincidence"]
    output[ValueType == "prevalence_mtp", ValueType := "PCR\nprevalence (%)"]
    output[ValueType == "sero_prev", ValueType := "Seroprevalence\n(%)"]
    output[ValueType == "type28_death_inc_line", ValueType := "Deaths"]

    data[ValueType == "hospital_inc", ValueType := "Hospital\nadmissions"]
    data[ValueType == "hospital_prev", ValueType := "Hospital beds\noccupied"]
    data[ValueType == "icu_prev", ValueType := "ICU beds\noccupied"]
    data[ValueType == "infections_inc", ValueType := "Infection\nincidence"]
    data[ValueType == "prevalence_mtp", ValueType := "PCR\nprevalence (%)"]
    data[ValueType == "sero_prev", ValueType := "Seroprevalence\n(%)"]
    data[ValueType == "type28_death_inc_line", ValueType := "Deaths"]
    
    return (list(data, output))
}

check_fit = function(test, parametersI, ld, sitreps, virus, sero, populations, death_cutoff, max_date, min_date = NULL, sero_cut_off)
{
    # Calculate total population
    popsize = NULL
    for (i in seq_along(parametersI)) {
        if (!is.null(parametersI[[i]])) {
            popsize = rbind(popsize,
                            data.table(Geography = parametersI[[i]]$pop[[1]]$name, 
                                       population_size = sum(parametersI[[i]]$pop[[1]]$size),
                                       population_size_15plus = sum(parametersI[[i]]$pop[[1]]$size[-(1:3)]))
            )
        }
    }
    adj_data = function(data, val_type, div, pop = 0) {
        data[ValueType == val_type, ymin := ymin / (div + population_size * pop)]
        data[ValueType == val_type, y    := y    / (div + population_size * pop)]
        data[ValueType == val_type, ymax := ymax / (div + population_size * pop)]
    }
    extra_sero = sero[sero$Start.date >= sero_cut_off]
    extra_sero = extra_sero[NHS.region %in% populations]
    extra_sero = extra_sero[!Data.source %like% "NHS BT", .(ValueType = "Seroprevalence\n(%)", Geography = NHS.region,
                                                            dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
                                                            ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))]
    extra_sero = merge(extra_sero, popsize, by = "Geography")
    adj_data(extra_sero, "Seroprevalence\n(%)", 0.01)
    
    fit = gen_fit(test, parametersI, ld, sitreps, virus, sero, populations, sero_cut_off)
    data = fit[[1]]
    output = fit[[2]]
    output = output[d <= max_date]
    
    if (!is.null(min_date)) {
        output = output[d >= min_date]
        data = data[d >= min_date]
    }
    
    # Augment deaths data
    cutoff_date = ld[, max(date)] - death_cutoff;
    data = data[ValueType != "Deaths" | (d <= cutoff_date)]

    # Make plot
    theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))
    
    linetypes = c("Deaths", "Hospital\nadmissions", "Hospital beds\noccupied", "ICU beds\noccupied")
    
    plot = ggplot(output[d > "2020-03-01" & AgeBand == "All"]) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = ValueType), alpha = 0.5) +
        geom_line(aes(x = d, y = Value, colour = ValueType)) +
        geom_line(data = data[ValueType %in% linetypes], aes(x = d, y = y), size = 0.2) +
        geom_point(data = data[!ValueType %in% linetypes], aes(x = d, y = y), size = 0.01, shape = 20) +
        geom_linerange(data = data, aes(x = d, ymin = ymin, ymax = ymax), size = 0.2) +
        geom_linerange(data = data, aes(xmin = dmin, xmax = dmax, y = y), size = 0.2) +
        geom_point(data = extra_sero, aes(x = d, y = y), size = 0.01, shape = 20, color = 'dodgerblue2', fill = 'dodgerblue2') +
        geom_linerange(data = extra_sero, aes(x = d, ymin = ymin, ymax = ymax), size = 0.2, color = 'dodgerblue2') +
        geom_linerange(data = extra_sero, aes(xmin = dmin, xmax = dmax, y = y), size = 0.2, color = 'dodgerblue2') +
        facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
        theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank()) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        labs(x = NULL, y = NULL)
    
    return (plot)
}

check_fit_small_output = function(test, parametersI, ld, sitreps, virus, sero, populations, death_cutoff, max_date, min_date = NULL, sero_cut_off)
{
    # Calculate total population
    popsize = NULL
    for (i in seq_along(parametersI)) {
        if (!is.null(parametersI[[i]])) {
            popsize = rbind(popsize,
                            data.table(Geography = parametersI[[i]]$pop[[1]]$name, 
                                       population_size = sum(parametersI[[i]]$pop[[1]]$size),
                                       population_size_15plus = sum(parametersI[[i]]$pop[[1]]$size[-(1:3)]))
            )
        }
    }
    adj_data = function(data, val_type, div, pop = 0) {
        data[ValueType == val_type, ymin := ymin / (div + population_size * pop)]
        data[ValueType == val_type, y    := y    / (div + population_size * pop)]
        data[ValueType == val_type, ymax := ymax / (div + population_size * pop)]
    }
    extra_sero = sero[sero$Start.date >= sero_cut_off]
    extra_sero = extra_sero[NHS.region %in% populations]
    extra_sero = extra_sero[!Data.source %like% "NHS BT", .(ValueType = "Seroprevalence\n(%)", Geography = NHS.region,
                                                            dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
                                                            ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))]
    extra_sero = merge(extra_sero, popsize, by = "Geography")
    adj_data(extra_sero, "Seroprevalence\n(%)", 0.01)
    
    fit = gen_fit(test, parametersI, ld, sitreps, virus, sero, populations, sero_cut_off)
    data = fit[[1]]
    output = fit[[2]]
    output = output[d <= max_date]
    
    if (!is.null(min_date)) {
        output = output[d >= min_date]
        data = data[d >= min_date]
    }
    
    # Augment deaths data
    cutoff_date = ld[, max(date)] - death_cutoff;
    data = data[ValueType != "Deaths" | (d <= cutoff_date)]
    
    # Make plot
    theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))
    
    linetypes = c("Deaths", "Hospital\nadmissions", "Hospital beds\noccupied", "ICU beds\noccupied")
    valuetypes = c("Deaths", "Hospital\nadmissions", "Hospital beds\noccupied", "ICU beds\noccupied", "PCR\nprevalence (%)", "Seroprevalence\n(%)")
    
    plot = ggplot(output[d > "2020-03-01" & AgeBand == "All" & ValueType %in% valuetypes]) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = ValueType), alpha = 0.5) +
        geom_line(aes(x = d, y = Value, colour = ValueType)) +
        geom_line(data = data[ValueType %in% linetypes], aes(x = d, y = y), size = 0.2) +
        geom_point(data = data[!ValueType %in% linetypes], aes(x = d, y = y), size = 0.01, shape = 20) +
        geom_linerange(data = data, aes(x = d, ymin = ymin, ymax = ymax), size = 0.2) +
        geom_linerange(data = data, aes(xmin = dmin, xmax = dmax, y = y), size = 0.2) +
        geom_point(data = extra_sero, aes(x = d, y = y), size = 0.01, shape = 20, color = 'dodgerblue2', fill = 'dodgerblue2') +
        geom_linerange(data = extra_sero, aes(x = d, ymin = ymin, ymax = ymax), size = 0.2, color = 'dodgerblue2') +
        geom_linerange(data = extra_sero, aes(xmin = dmin, xmax = dmax, y = y), size = 0.2, color = 'dodgerblue2') +
        facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
        theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank()) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        labs(x = NULL, y = NULL)
    
    return (plot)
}



compare_fit = function(test, test0, ld, sitreps, virus, sero, populations, populations0, max_date)
{
    fit = gen_fit(test, ld, sitreps, virus, sero, populations)
    output = fit[[2]]

    fit = gen_fit(test0, ld, sitreps, virus, sero, populations0)
    data = fit[[1]]
    output0 = fit[[2]]
    
    output[, kind := "With VOC"]
    output0[, kind := "Without VOC"]
    output = rbind(output, output0)
    output = output[d <= max_date]
    
    # Make plot
    theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))
    
    linetypes = c("Deaths", "Hospital\nadmissions", "Hospital beds\noccupied", "ICU beds\noccupied")
    
    plot = ggplot(output[d > "2020-03-01" & AgeBand == "All"]) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = ValueType, group = kind), alpha = 0.5) +
        geom_line(aes(x = d, y = Value, colour = ValueType, linetype = kind)) +
        geom_line(data = data[ValueType %in% linetypes], aes(x = d, y = y), size = 0.2) +
        geom_point(data = data[!ValueType %in% linetypes], aes(x = d, y = y), size = 0.01, shape = 20) +
        geom_linerange(data = data, aes(x = d, ymin = ymin, ymax = ymax), size = 0.2) +
        geom_linerange(data = data, aes(xmin = dmin, xmax = dmax, y = y), size = 0.2) +
        facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
        theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank()) +
        scale_x_date(date_breaks = "2 months", date_labels = "%b") +
        labs(x = NULL, y = NULL)
    
    return (plot)
}
