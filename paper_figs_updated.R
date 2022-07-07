library(data.table)
library(ggplot2)
library(qs)
library(lubridate)
library(stringr)
library(httr)

options(scipen = 999)

popUK = readRDS("./fitting_data/popNHS.rds");
nhs_regions = popUK[, unique(name)]

###############################################################################

# Figure 1 - aggregated model fit to whole of England

SERO_CUT_OFF = '2020-12-01'

# load data
all_data = qread("./fitting_data/processed-data-20220506122858.qs")
ld = all_data[[1]]
sitreps = all_data[[2]]
virus = all_data[[3]][!Data.source %like% "7a|7b|6a|6b|9a|9b"]
sero = all_data[[4]]
# we only want to fit to seroprevalence data prior to December 2020
sero_to_fit = sero[sero$Start.date < SERO_CUT_OFF]
sgtf = all_data[[5]]
delta = all_data[[6]]
omi = all_data[[7]]
# remove any omicron sgtf data occurring before 1st October 2021
omi = omi[omi$date >= '2021-10-01']

# five colours
cols = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')

# load full output file
w = qread('./output/paper/may22/adjusted1_basecase_220511.qs')

source("./spim_output.R")

w[, sero_positive := lfia_positive_p + vaccsero_a_p + vaccsero_b_p]

popsize = w[population != "England" & t == 0 & run == 1, .(population_size = sum(S), population_size_15plus = sum(S[-(1:3)])), by = .(Geography = population)]
popsize = rbind(popsize, popsize[, .(Geography = "England", population_size = sum(population_size), population_size_15plus = sum(population_size_15plus))])
w = merge(w, popsize[, .(population = Geography, population_size, population_size_15plus)], by = "population")

england_x1 = w[, lapply(.SD, sum), .SDcols = c("deaths", "hosp_adm", "hosp_bed", "icu_bed", "pcr_positive_p", "sero_positive", "hosp_undetected_p", "hosp_undetected_o", "lfia_positive_p", "vaccsero_a_p", "vaccsero_b_p"), by = .(t, run, group)]
england_x2 = w[, lapply(.SD, weighted.mean, population_size), .SDcols = c("disp_deaths", "disp_hosp_inc", "disp_hosp_prev", "disp_icu_prev"), by = .(t, run, group)]
england_x1[, population := "England"]
england = merge(england_x1, england_x2, by = c("t", "run", "group"))
england[, NHS.region := "England"]

w = rbind(w, england, fill = TRUE)

var1 = readline(prompt = "Load existing spim_output? (Y/N): ")
while (var1 != "N" && var1 != "Y" && var1 != "n" && var1 != "y"){
    var1 = readline(prompt = "Load existing spim_output? (Y/N): ")
}
if (var1 %in% c("Y", "y")){
    cat('Loading existing spim_output file')
    spim_output = fread('./output/paper/may22/basecase_220511_spim_output_20220511190051.csv')
} else if (var1 %in% c('N', 'n')) {
    cat('Generating new spim_output')
    spim_output = SPIM_output_full(w, 2022, 05, 11, "2020-01-01", "2022-12-31")
} else {
    cat('Error: spim_output not loaded')
}

spim_output[, d := make_date(`Year of Value`, `Month of Value`, `Day of Value`)]
spim_output = merge(spim_output, popsize, by = "Geography")

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

adj_output(spim_output, "hospital_inc", 1)
adj_output(spim_output, "hospital_prev", 1)
adj_output(spim_output, "icu_prev", 1)
adj_output(spim_output, "prevalence_mtp", 0, 0.01, 0)
adj_output(spim_output, "sero_prev", 0, 0, 0.01)
adj_output(spim_output, "type28_death_inc_line", 1)

# Make data to output
pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100
data = make_data(ld, sitreps, virus, sero_to_fit)
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

spim_output[ValueType == "hospital_inc", ValueType := "Hospital\nadmissions"]
spim_output[ValueType == "hospital_prev", ValueType := "Hospital beds\noccupied"]
spim_output[ValueType == "icu_prev", ValueType := "ICU beds\noccupied"]
spim_output[ValueType == "infections_inc", ValueType := "Infection\nincidence"]
spim_output[ValueType == "prevalence_mtp", ValueType := "PCR\nprevalence (%)"]
spim_output[ValueType == "sero_prev", ValueType := "Sero-\nprevalence (%)"]
spim_output[ValueType == "type28_death_inc_line", ValueType := "Deaths"]

data[ValueType == "hospital_inc", ValueType := "Hospital\nadmissions"]
data[ValueType == "hospital_prev", ValueType := "Hospital beds\noccupied"]
data[ValueType == "icu_prev", ValueType := "ICU beds\noccupied"]
data[ValueType == "infections_inc", ValueType := "Infection\nincidence"]
data[ValueType == "prevalence_mtp", ValueType := "PCR\nprevalence (%)"]
data[ValueType == "sero_prev", ValueType := "Sero-\nprevalence (%)"]
data[ValueType == "type28_death_inc_line", ValueType := "Deaths"]

max_date = '2022-05-06'
date_breaks_opt = "4 months"

F1_deaths = ggplot(spim_output[d >= "2020-03-01" & d <= max_date & AgeBand == "All" & ValueType == "Deaths" & Geography == "England"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
    geom_line(aes(x = d, y = Value, colour = "Modelled")) +
    geom_line(data = data[ValueType == "Deaths" & Geography == "England" & d <= max_date], aes(x = d, y = y, colour = "Data"), size = 0.25) +
    cowplot::theme_cowplot(font_size = 7, font_family = "sans") +
    theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
    labs(x = NULL, y = "Deaths", fill = NULL, colour = NULL) +
    scale_colour_manual(values = c(Data = "black", Modelled = cols[1]), aesthetics = c("fill", "colour")) +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))

F1_hosp_adm = ggplot(spim_output[d >= "2020-03-01" & d <= max_date & AgeBand == "All" & ValueType == "Hospital\nadmissions" & Geography == "England"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
    geom_line(aes(x = d, y = Value, colour = "Modelled")) +
    geom_line(data = data[ValueType == "Hospital\nadmissions" & Geography == "England" & d <= max_date], aes(x = d, y = y, colour = "Data"), size = 0.2) +
    cowplot::theme_cowplot(font_size = 7) +
    theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
    labs(x = NULL, y = "Hospital\nadmissions", fill = NULL, colour = NULL) +
    scale_colour_manual(values = c(Data = "black", Modelled = cols[2]), aesthetics = c("fill", "colour")) +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))

F1_hosp_bed = ggplot(spim_output[d >= "2020-03-01" & d <= max_date & AgeBand == "All" & ValueType == "Hospital beds\noccupied" & Geography == "England"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
    geom_line(aes(x = d, y = Value, colour = "Modelled")) +
    geom_line(data = data[ValueType == "Hospital beds\noccupied" & Geography == "England" & d <= max_date], aes(x = d, y = y, colour = "Data"), size = 0.2) +
    cowplot::theme_cowplot(font_size = 7) +
    theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
    labs(x = NULL, y = "Hospital beds\noccupied", fill = NULL, colour = NULL) +
    scale_colour_manual(values = c(Data = "black", Modelled = cols[3]), aesthetics = c("fill", "colour")) +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))

F1_icu_bed = ggplot(spim_output[d >= "2020-03-01" & d <= max_date & AgeBand == "All" & ValueType == "ICU beds\noccupied" & Geography == "England"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
    geom_line(aes(x = d, y = Value, colour = "Modelled")) +
    geom_line(data = data[ValueType == "ICU beds\noccupied" & Geography == "England" & d <= max_date], aes(x = d, y = y, colour = "Data"), size = 0.2) +
    cowplot::theme_cowplot(font_size = 7) +
    theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
    labs(x = NULL, y = "ICU beds\noccupied", fill = NULL, colour = NULL) +
    scale_colour_manual(values = c(Data = "black", Modelled = cols[4]), aesthetics = c("fill", "colour")) +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))

F1_pcr = ggplot(spim_output[d >= "2020-03-01" & d <= max_date & AgeBand == "All" & ValueType == "PCR\nprevalence (%)" & Geography == "England"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
    geom_line(aes(x = d, y = Value, colour = "Modelled")) +
    geom_ribbon(data = data[ValueType == "PCR\nprevalence (%)" & Geography == "England" & d <= max_date], aes(x = d, ymin = ymin, ymax = ymax, fill = "Data"), alpha = 0.75) +
    cowplot::theme_cowplot(font_size = 7) +
    theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
    labs(x = NULL, y = "PCR\nprevalence (%)", fill = NULL, colour = NULL) +
    scale_colour_manual(values = c(Data = "black", Modelled = cols[5]), aesthetics = c("fill", "colour")) +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))

fig_england_fit = cowplot::plot_grid(
    F1_deaths,
    F1_hosp_adm,
    F1_hosp_bed,
    F1_icu_bed,
    F1_pcr,
    ncol = 1, labels = letters, label_size = 7, align = "v", axis = "bottom")
fig_england_fit

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/july22/fig2_", datetime, ".png"), fig_england_fit, width = 88, height = 130, units = "mm")
ggsave(paste0("./output/paperfigs/july22/fig2_", datetime, ".pdf"), fig_england_fit, width = 88, height = 130, units = "mm", useDingbats = TRUE)

date_breaks_opt = "2 months"
F1_deaths_2col = F1_deaths +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))

F1_hosp_adm_2col = F1_hosp_adm +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))

F1_hosp_bed_2col = F1_hosp_bed +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))

F1_icu_bed_2col = F1_icu_bed +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))

F1_pcr_2col = F1_pcr +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))

fig_england_fit_2col = cowplot::plot_grid(
    F1_deaths_2col,
    F1_hosp_adm_2col,
    F1_hosp_bed_2col,
    F1_icu_bed_2col,
    F1_pcr_2col,
    ncol = 1, labels = letters, label_size = 7, align = "v", axis = "bottom")
fig_england_fit_2col

ggsave(paste0("./output/paperfigs/july22/fig2_", datetime, "_2column.png"), fig_england_fit_2col, width = 180, height = 185, units = "mm")
ggsave(paste0("./output/paperfigs/july22/fig2_", datetime, "_2column.pdf"), fig_england_fit_2col, width = 180, height = 185, units = "mm", useDingbats = TRUE)

if(0){
    # plot same output as above but don't restrict spim_output to max_date
    F12_deaths = ggplot(spim_output[d >= "2020-03-01" & AgeBand == "All" & ValueType == "Deaths" & Geography == "England"]) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
        geom_line(aes(x = d, y = Value, colour = "Modelled")) +
        geom_line(data = data[ValueType == "Deaths" & Geography == "England" & d <= max_date], aes(x = d, y = y, colour = "Data"), size = 0.25) +
        cowplot::theme_cowplot(font_size = 9) +
        theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
        labs(x = NULL, y = "Deaths", fill = NULL, colour = NULL) +
        scale_colour_manual(values = c(Data = "black", Modelled = cols[1]), aesthetics = c("fill", "colour")) +
        scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))
    
    F12_hosp_adm = ggplot(spim_output[d >= "2020-03-01" & AgeBand == "All" & ValueType == "Hospital\nadmissions" & Geography == "England"]) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
        geom_line(aes(x = d, y = Value, colour = "Modelled")) +
        geom_line(data = data[ValueType == "Hospital\nadmissions" & Geography == "England" & d <= max_date], aes(x = d, y = y, colour = "Data"), size = 0.2) +
        cowplot::theme_cowplot(font_size = 9) +
        theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
        labs(x = NULL, y = "Hospital\nadmissions", fill = NULL, colour = NULL) +
        scale_colour_manual(values = c(Data = "black", Modelled = cols[2]), aesthetics = c("fill", "colour")) +
        scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))
    
    F12_hosp_bed = ggplot(spim_output[d >= "2020-03-01" & AgeBand == "All" & ValueType == "Hospital beds\noccupied" & Geography == "England"]) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
        geom_line(aes(x = d, y = Value, colour = "Modelled")) +
        geom_line(data = data[ValueType == "Hospital beds\noccupied" & Geography == "England" & d <= max_date], aes(x = d, y = y, colour = "Data"), size = 0.2) +
        cowplot::theme_cowplot(font_size = 9) +
        theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
        labs(x = NULL, y = "Hospital beds\noccupied", fill = NULL, colour = NULL) +
        scale_colour_manual(values = c(Data = "black", Modelled = cols[3]), aesthetics = c("fill", "colour")) +
        scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))
    
    F12_icu_bed = ggplot(spim_output[d >= "2020-03-01" & AgeBand == "All" & ValueType == "ICU beds\noccupied" & Geography == "England"]) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
        geom_line(aes(x = d, y = Value, colour = "Modelled")) +
        geom_line(data = data[ValueType == "ICU beds\noccupied" & Geography == "England" & d <= max_date], aes(x = d, y = y, colour = "Data"), size = 0.2) +
        cowplot::theme_cowplot(font_size = 9) +
        theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
        labs(x = NULL, y = "ICU beds\noccupied", fill = NULL, colour = NULL) +
        scale_colour_manual(values = c(Data = "black", Modelled = cols[4]), aesthetics = c("fill", "colour")) +
        scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))
    
    F12_pcr = ggplot(spim_output[d >= "2020-03-01" & AgeBand == "All" & ValueType == "PCR\nprevalence (%)" & Geography == "England"]) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
        geom_line(aes(x = d, y = Value, colour = "Modelled")) +
        geom_ribbon(data = data[ValueType == "PCR\nprevalence (%)" & Geography == "England" & d <= max_date], aes(x = d, ymin = ymin, ymax = ymax, fill = "Data"), alpha = 0.75) +
        cowplot::theme_cowplot(font_size = 9) +
        theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
        labs(x = NULL, y = "PCR\nprevalence (%)", fill = NULL, colour = NULL) +
        scale_colour_manual(values = c(Data = "black", Modelled = cols[5]), aesthetics = c("fill", "colour")) +
        scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))
    
    fig2_england_fit = cowplot::plot_grid(
        F12_deaths,
        F12_hosp_adm,
        F12_hosp_bed,
        F12_icu_bed,
        F12_pcr,
        ncol = 1, labels = letters, align = "v", axis = "bottom")
    fig2_england_fit
    datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
    ggsave(paste0("./output/paperfigs/may22/dummyfig_aggregated_fit_to_england_and_basecase_projection", datetime, ".png"), fig_england_fit, width = 15, height = 20, units = "cm")
    ggsave(paste0("./output/paperfigs/may22/dummyfig_aggregated_fit_to_england_and_basecase_projection", datetime, ".pdf"), fig_england_fit, width = 15, height = 20, units = "cm", useDingbats = TRUE)
    
}

###############################################################################

# Figure S1 - model fits to NHS England regions

data$Geography[data$Geography == 'North East and Yorkshire'] = 'NE & Y'
spim_output$Geography[spim_output$Geography == 'North East and Yorkshire'] = 'NE & Y'

valuetypes = c("Hospital\nadmissions", "Hospital beds\noccupied", 
               "ICU beds\noccupied", "PCR\nprevalence (%)", 
               "Seroprevalence (%)\n & seroconversions", "Deaths")
linetypes = c("Deaths", "Hospital\nadmissions", "Hospital beds\noccupied", "ICU beds\noccupied")

six_colours = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02')

extra_sero = sero[sero$Start.date >= SERO_CUT_OFF]
extra_sero = extra_sero[NHS.region %in% nhs_regions]
extra_sero = extra_sero[!Data.source %like% "NHS BT", .(ValueType = "Sero-\nprevalence (%)", Geography = NHS.region,
                                                        dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
                                                        ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))]
extra_sero = merge(extra_sero, popsize, by = "Geography")
adj_data(extra_sero, "Sero-\nprevalence (%)", 0.01)
extra_sero$Geography[extra_sero$Geography == 'North East and Yorkshire'] = 'NE & Y'

spim_output$ValueType[spim_output$ValueType == 'Sero-\nprevalence (%)'] = 'Seroprevalence (%)\n & seroconversions'
data$ValueType[data$ValueType == 'Sero-\nprevalence (%)'] = 'Seroprevalence (%)\n & seroconversions'
extra_sero$ValueType[extra_sero$ValueType == 'Sero-\nprevalence (%)'] = 'Seroprevalence (%)\n & seroconversions'

plot = ggplot(spim_output[d > "2020-03-01" & d < max_date & AgeBand == "All" & ValueType %in% valuetypes & Geography != "England"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = ValueType), alpha = 0.5) +
    geom_line(aes(x = d, y = Value, colour = ValueType)) +
    geom_line(data = data[ValueType %in% linetypes & Geography != "England"], aes(x = d, y = y), size = 0.2) +
    geom_point(data = data[!ValueType %in% linetypes & Geography != "England"], aes(x = d, y = y), size = 0.01, shape = 20) +
    geom_linerange(data = data[Geography != "England"], aes(x = d, ymin = ymin, ymax = ymax), size = 0.2) +
    geom_linerange(data = data[Geography != "England"], aes(xmin = dmin, xmax = dmax, y = y), size = 0.2) +
    geom_point(data = extra_sero[Geography != "England"], aes(x = d, y = y), size = 0.01, shape = 20, color = 'dodgerblue2', fill = 'dodgerblue2') +
    geom_linerange(data = extra_sero[Geography != "England"], aes(x = d, ymin = ymin, ymax = ymax), size = 0.2, color = 'dodgerblue2') +
    geom_linerange(data = extra_sero[Geography != "England"], aes(xmin = dmin, xmax = dmax, y = y), size = 0.2, color = 'dodgerblue2') +
    facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
    cowplot::theme_cowplot(font_size = 9) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank(),
          panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
    scale_x_date(breaks = ymd(c("2020-07-01", "2021-03-01", "2021-11-01")), date_labels = "%m/%y") +
    labs(x = NULL, y = NULL) +
    scale_fill_manual(values=six_colours) +
    scale_color_manual(values=six_colours)
plot
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figS1_modelfit_regional_", datetime, ".pdf"), plot, width = 30, height = 20, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figS1_modelfit_regional_", datetime, ".png"), plot, width = 30, height = 20, units = "cm")

# try splitting regions up into a group of 4 and a group of 3
regionsA = c("East of England", "London", "Midlands", "NE & Y")
regionsB = c("North West", "South East", "South West")

plotA = ggplot(spim_output[Geography %in% regionsA & d > "2020-03-01" & d < max_date & AgeBand == "All" & ValueType %in% valuetypes]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = ValueType), alpha = 0.5) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = ValueType), alpha = 0.75) +
    geom_line(aes(x = d, y = Value, colour = ValueType)) +
    geom_line(data = data[ValueType %in% linetypes & Geography %in% regionsA], aes(x = d, y = y), size = 0.2) +
    geom_linerange(data = data[Geography %in% regionsA], aes(x = d, ymin = ymin, ymax = ymax), size = 0.2, alpha = 0.5) +
    geom_linerange(data = data[Geography %in% regionsA], aes(xmin = dmin, xmax = dmax, y = y), size = 0.2, alpha = 0.5) +
    geom_linerange(data = extra_sero[Geography %in% regionsA], aes(x = d, ymin = ymin, ymax = ymax), size = 0.2, color = 'dodgerblue2', alpha = 0.5) +
    geom_linerange(data = extra_sero[Geography %in% regionsA], aes(xmin = dmin, xmax = dmax, y = y), size = 0.2, color = 'dodgerblue2', alpha = 0.5) +
    facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
    cowplot::theme_cowplot(font_size = 9) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank(),
          panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
    scale_x_date(breaks = ymd(c("2020-07-01", "2021-01-01", "2021-07-01", "2022-01-01")), date_labels = "%m/%y") +
    labs(x = NULL, y = NULL) +
    scale_fill_manual(values=six_colours) +
    scale_color_manual(values=six_colours)
plotA
plotB = ggplot(spim_output[Geography %in% regionsB & d > "2020-03-01" & d < max_date & AgeBand == "All" & ValueType %in% valuetypes]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = ValueType), alpha = 0.5) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = ValueType), alpha = 0.75) +
    geom_line(aes(x = d, y = Value, colour = ValueType)) +
    geom_line(data = data[ValueType %in% linetypes & Geography %in% regionsB], aes(x = d, y = y), size = 0.2) +
    geom_linerange(data = data[Geography %in% regionsB], aes(x = d, ymin = ymin, ymax = ymax), size = 0.2, alpha = 0.5) +
    geom_linerange(data = data[Geography %in% regionsB], aes(xmin = dmin, xmax = dmax, y = y), size = 0.2, alpha = 0.5) +
    geom_linerange(data = extra_sero[Geography %in% regionsB], aes(x = d, ymin = ymin, ymax = ymax), size = 0.2, color = 'dodgerblue2', alpha = 0.5) +
    geom_linerange(data = extra_sero[Geography %in% regionsB], aes(xmin = dmin, xmax = dmax, y = y), size = 0.2, color = 'dodgerblue2', alpha = 0.5) +
    facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
    cowplot::theme_cowplot(font_size = 9) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank(),
          panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5)) +
    scale_x_date(breaks = ymd(c("2020-07-01", "2021-01-01", "2021-07-01", "2022-01-01")), date_labels = "%m/%y") +
    labs(x = NULL, y = NULL) +
    scale_fill_manual(values=six_colours) +
    scale_color_manual(values=six_colours)
plotB
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figS1A_modelfit_regional_", datetime, ".pdf"), plotA, width = 20, height = 20, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figS1A_modelfit_regional_", datetime, ".png"), plotA, width = 20, height = 20, units = "cm")
ggsave(paste0("./output/paperfigs/may22/figS1B_modelfit_regional_", datetime, ".pdf"), plotB, width = 15, height = 20, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figS1B_modelfit_regional_", datetime, ".png"), plotB, width = 15, height = 20, units = "cm")


###############################################################################

# Figure 2 - mobility scenarios and transmission adjustments

# load traces
pf = qread("./fits/pf_relu_yeswane_sev2.0_22050605_20220511163014.qs")
traces = rbindlist(pf$traces, idcol = "population")
posteriors = rbindlist(pf$posteriors, idcol = "population")
traces = merge(traces, posteriors[, .(population, u)], by = "population")

# traces before correspond to traces before the particle filter time limit, 
# which currently corresponds to the maximum PCR prevalence datapoint, to the 
# nearest 5 days -> 30th April 2022 (so t=766 rounded to t=765)
date_fitting = as.Date(max(virus$Start.date))
t_fitting = as.numeric(date_fitting - as.Date('2020-01-01')) 
projection_end = as.Date('2022-12-31')
projection_t = as.numeric(projection_end - as.Date('2020-01-01'))
tr_before = traces[t <= t_fitting, .(adj = weighted.mean(adj, u)), by = .(t, population)]
tr_before_mean = traces[t <= t_fitting, .(adj = weighted.mean(adj, u)), by = t]
tr_all_mean = traces[t <= projection_t & i == 0, .(adj = weighted.mean(adj, u)), by = t]
tr_q = traces[t >= t_fitting & t <= projection_t  & i != 0, .(
    q025 = quantile(adj, 0.025),
    q25  = quantile(adj, 0.25),
    q75  = quantile(adj, 0.75),
    q975 = quantile(adj, 0.975)
), by = .(t)]
tr_sample = traces[t >= t_fitting & t <= projection_t & i == 10, .(adj = weighted.mean(adj, u)), by = .(t, population)]
tr_sample_mean = traces[t >= t_fitting & t <= projection_t & i == 10, .(adj = weighted.mean(adj, u)), by = .(t)]
# set up labels for lockdowns and roadmap steps
marks = fread(
    "label,date
    L1,2020-03-26
    L2,2020-11-05
    L3,2021-01-05
    S1,2021-03-08
    S2,2021-04-12
    S3,2021-05-17
    S4,2021-07-19
    PBA,2021-12-08
    PBE,2022-01-27") # Plan B measures were announced on Wednesday 8th December 2021; 
                     # masks required from Friday 10th, WFH from Monday 13th; COVID passes from Wednesday 15th
                     # Plan B measures will end on 27th January 2022 (announced on 19th January 2022)
marks[, date := ymd(date)]
marks[, y := 3.0]

colours = c(
    "East of England" =  pals::tableau20(20)[1],
    "London" =  pals::tableau20(20)[2],
    "Midlands" =  pals::tableau20(20)[3],
    "NE & Y" =  pals::tableau20(20)[4],
    "North West" =  pals::tableau20(20)[5],
    "South East" =  pals::tableau20(20)[6],
    "South West" = pals::tableau20(20)[7],
    "Projection" = pals::tableau20(20)[14]
)

colfac = names(colours)
nhs_regions[nhs_regions == "North East and Yorkshire"] = "NE & Y"
tr_before[, region := factor(nhs_regions[population], levels = colfac)]
# tr_before_mean[, region := factor(nhs_regions[population], levels = colfac)]
tr_sample[, region := factor(nhs_regions[population], levels = colfac)]
tr_q[, region := factor("Projection", levels = colfac)]
nhs_regions[nhs_regions == "NE & Y"] == "North East and Yorkshire"

pl_tadj = ggplot() + 
    geom_hline(aes(yintercept = 1), colour = "#d0d0d0") +
    geom_ribbon(data = tr_q, aes(t + ymd("2020-01-01") - 5, ymin = q25, ymax = q75, fill = region), alpha = 0.4) +
    geom_line(data = tr_before, aes(t + ymd("2020-01-01") - 5, adj, colour = region), size = 0.3) + #colour = "Fitted", 
    geom_line(data = tr_before_mean, aes(t + ymd("2020-01-01") - 5, adj), size = 0.4) + #colour = "Fitted", 
    geom_line(data = tr_sample[population == 1], aes(t + ymd("2020-01-01") - 5, adj, colour = region), size = 0.3) +
    scale_colour_manual(values = colours[-8], guide = guide_legend(order = 1)) +
    scale_fill_manual(values = colours[8], guide = guide_legend(order = 2)) +
    geom_line(data = tr_sample_mean, aes(t + ymd("2020-01-01") - 5, adj), size = 0.4) +
    geom_vline(data = marks, aes(xintercept = marks$date), linetype = "22", size = 0.25) +
    geom_label(data = marks, aes(x = marks$date, y = 7.5, label = label), size = 2.5, label.padding = unit(0.15, "lines")) +
    annotate("blank", x = ymd(c("2020-03-01"), c("2022-09-30"))) +
    scale_x_date(date_breaks = "2 months", date_labels = "%m/%y", expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0.5, 9, by = 1), limits = c(0.25, 9.25)) +
    labs(x = NULL, y = "Multiplicative factor", colour = "Region", fill = NULL, title = "Transmission adjustment") +
    cowplot::theme_cowplot(font_size = 9) +
    theme(panel.background = element_rect(fill = "#f4f4f4"), strip.background = element_blank(),
          panel.grid = element_line(colour = "#ffffff"), 
          text=element_text(size=7, family="sans"))
pl_tadj

get_overall = function(roadmap, tag = NULL)
{
    overall = roadmap[region_name == "England", .(mob = mean(c(retrec, transit, workplace, school))), by = date]
    overall[, mob := zoo::rollmean(mob, 7, fill = "extend")]
    overall = merge(overall, tr_all_mean[, .(date = t - 5 + ymd("2020-01-01"), adj)], all.x = TRUE)
    overall[, adj := zoo::na.locf(adj, na.rm = FALSE)]
    overall[, adj := zoo::na.locf(adj, na.rm = FALSE, fromLast = TRUE)]
    overall[, adj := zoo::rollmean(adj, 7, fill = "extend")]
    if (!is.null(tag)) {
        overall[, tag := ..tag]
    }
    return (overall)
}

##################################### get mobility schedules ##################

source("build_schedule_google_mobility.R")

xs = make_xs("~/Desktop/Global_Mobility_Report-2022-05-04.csv", TRUE, NA)
# save xs for later
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
write.csv(xs, paste0('./output/paper/may22/', datetime, '_xs_built_from_mobility_2022-05-04.csv'), row.names = FALSE)
schools_data = fread('./fitting_data/schoolattendancedata-England-20220506120901.csv')
schools_data[, date := as.Date(date)]
y_schools = make_y(xs, schools_data);
HOLD_FROM = max(xs$date)

paper_aw_3wk = steps_y(
    xs,
    list(c(1.0, 1.0, 1.0, 1.0)), 
    hold_from = HOLD_FROM, 
    ramp_start= c(HOLD_FROM), 
    ramp_end =  c(HOLD_FROM+21),
    SMOOTH_ADJ = 7) # HOLD_FROM + 21 falls on the same day 3 weeks later
paper_aw_3wk[, school := NULL]
paper_aw_3wk = merge(paper_aw_3wk, y_schools[, .(date, region_name, school)], by = c("date", "region_name"), all.x = TRUE)
setcolorder(paper_aw_3wk, names(y_schools))
setorder(paper_aw_3wk, pop, date)

paper_aw_3mo = steps_y(
    xs,
    list(c(1.0, 1.0, 1.0, 1.0)), 
    hold_from = HOLD_FROM, 
    ramp_start= c(HOLD_FROM), 
    ramp_end =  c(HOLD_FROM+91),
    SMOOTH_ADJ = 7) # HOLD_FROM + 91 falls on the same day 3 months later
paper_aw_3mo[, school := NULL]
paper_aw_3mo = merge(paper_aw_3mo, y_schools[, .(date, region_name, school)], by = c("date", "region_name"), all.x = TRUE)
setcolorder(paper_aw_3mo, names(y_schools))
setorder(paper_aw_3mo, pop, date)

paper_aw_6mo = steps_y(
    xs,
    list(c(1.0, 1.0, 1.0, 1.0)), 
    hold_from = HOLD_FROM, 
    ramp_start= c(HOLD_FROM), 
    ramp_end =  c(HOLD_FROM+183),
    SMOOTH_ADJ = 7) # HOLD_FROM + 183 falls on the same day 6 months later
paper_aw_6mo[, school := NULL]
paper_aw_6mo = merge(paper_aw_6mo, y_schools[, .(date, region_name, school)], by = c("date", "region_name"), all.x = TRUE)
setcolorder(paper_aw_6mo, names(y_schools))
setorder(paper_aw_6mo, pop, date)

paper_aw_flt = steps_y(
    xs,
    list(c(1.0, 1.0, 1.0, 1.0)), 
    hold_from = HOLD_FROM, 
    ramp_start= c('2023-01-01'), 
    ramp_end =  c('2023-01-31'),
    SMOOTH_ADJ = 7)
paper_aw_flt[, school := NULL]
paper_aw_flt = merge(paper_aw_flt, y_schools[, .(date, region_name, school)], by = c("date", "region_name"), all.x = TRUE)
setcolorder(paper_aw_flt, names(y_schools))
setorder(paper_aw_flt, pop, date)


##################################### end get mobility schedules ###############

overall_3wk = get_overall(paper_aw_3wk)
overall_3wk$scenario = rep('3 weeks', dim(overall_3wk)[1])
overall_6mo = get_overall(paper_aw_6mo)
overall_6mo$scenario = rep('6 months', dim(overall_6mo)[1])
overall_3mo = get_overall(paper_aw_3mo)
overall_3mo$scenario = rep('3 months', dim(overall_3mo)[1])
overall_flt = get_overall(paper_aw_flt)
overall_flt$scenario = rep('None', dim(overall_flt)[1])

overall = rbind(overall_3wk, overall_6mo, overall_3mo, overall_flt)
overall$scenario = factor(overall$scenario, levels = c('3 weeks', '3 months', '6 months', 'None'))

col3 = rev(c("#33c5ff", "#8a03fa", "#8f0c0c"))

# update height of lockdown labels
marks[, y2 := 0.1]

trpo_xtra = ggplot(overall[date >= "2020-03-01" & date <= "2022-12-31"]) +
    geom_line(data = overall[date >= '2020-03-01' & date <= max(xs$date)], aes(date, mob * adj, group = scenario)) +
    geom_line(data = overall[date > max(xs$date) & date <= '2022-12-31'], aes(date, mob * adj, group = scenario, colour = scenario)) +
    geom_vline(data = marks, aes(xintercept = marks$date), linetype = "22", size = 0.25) +
    geom_label(data = marks, aes(x = marks$date, y = 3.5, label = label), size = 2.5, label.padding = unit(0.15, "lines")) +
    scale_x_date(date_breaks = "2 months", date_labels = "%m/%y", expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 5.5, by = 0.5), limits = c(0, 5.5)) +
    labs(x = "Date", y = "Transmission potential", colour = 'Return to\nbaseline', title = "Overall") +
    cowplot::theme_cowplot(font_size = 9) +
    theme(panel.background = element_rect(fill = "#f4f4f4"), strip.background = element_blank(),
          panel.grid = element_line(colour = "#ffffff"),
          text=element_text(size=7, family="sans")) +
    scale_colour_manual(values = c('#1b9e77','#d95f02','#7570b3','#e7298a'), aesthetics = c("colour", "fill"))
trpo_xtra

# paw_full
paw_full = plot_y_fancy(list(
    `3 weeks`  = paper_aw_3wk,
    `3 months` = paper_aw_3mo,
    `6 months` = paper_aw_6mo,
    `None`     = paper_aw_flt
),  proj_start_date = max(xs$date), 
xmark_label = c("S4", "PBA", "PBE"), xmark_date = c('2021-07-19', '2021-12-08', '2022-01-27'),
start_date = "2020-03-01", end_date = "2022-12-31", 
dbreaks = "2 months", dlabels = "%m/%y", ybreaks = seq(0.2, 1.4, by = 0.2), ylim = c(0.2, 1.4),
colours_list = c('#1b9e77','#d95f02','#7570b3','#e7298a')) + 
    labs(title = "Mobility", colour = "Return to\nbaseline", x = NULL) +
    theme(text=element_text(size=7, family="sans"))
paw_full


fig1_xtra = cowplot::plot_grid(paw_full, pl_tadj, trpo_xtra, nrow = 3, rel_heights = c(3, 1.2, 1.3), align = "hv", axis = "lr")
fig1_xtra
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/july22/fig3_", datetime, "_2column.pdf"), fig1_xtra, width = 180, height = 195, units = "mm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/july22/fig3_", datetime, "_2column.png"), fig1_xtra, width = 180, height = 195, units = "mm")

###############################################################################

# Figure S2 - regional model fits to Alpha VOC

# load full output file
f = qread('./output/paper/may22/output1_basecase_220511.qs')

# load fit to get posteriors
fit = qread('./fits/pf_relu_yeswane_sev2.0_22050605_20220511163014.qs')
posterior = rbindlist(fit$posteriors, idcol = "population")
posterior = posterior[, .(v2_sgtf0 = mean(v2_sgtf0), v2_disp = mean(v2_disp),
                          v4_sgtf0 = mean(v4_sgtf0), v4_disp = mean(v4_disp)), by = population]
f = merge(f, posterior, by = "population")

################################################################################

# plot fit to SGTF data
sgtf[, qlo := qbeta(0.025, sgtf + 1, other + 1)]
sgtf[, qhi := qbeta(0.975, sgtf + 1, other + 1)]
vmodel = f[, .(I1 = sum(test_o + test3_o), I2 = sum(test2_o), sgtf0 = mean(v2_sgtf0), disp = mean(v2_disp)), by = .(t, population, run)]
vmodel[, conc := 1/(disp^2)]
vmodel[, p2 := I2 / (I1 + I2)]
vmodel[is.nan(p2), p2 := 0]
vmodel[, sgtf := (1 - p2) * sgtf0 + p2];
vmodel[, alpha := sgtf * (conc - 2) + 1]
vmodel[, beta := (1 - sgtf) * (conc - 2) + 1]
vmodel[, q025 := qbeta(0.025, alpha, beta)]
vmodel[, q500 := qbeta(0.500, alpha, beta)]
vmodel[, q975 := qbeta(0.975, alpha, beta)]
vmodel = vmodel[, lapply(.SD, mean), .SDcols = c("q025", "q500", "q975"), by = .(nhs_name = population, t)]
vmodel$nhs_name[vmodel$nhs_name == 1] = 'East of England'
vmodel$nhs_name[vmodel$nhs_name == 3] = 'London'
vmodel$nhs_name[vmodel$nhs_name == 4] = 'Midlands'
vmodel$nhs_name[vmodel$nhs_name == 5] = 'North East and Yorkshire'
vmodel$nhs_name[vmodel$nhs_name == 6] = 'North West'
vmodel$nhs_name[vmodel$nhs_name == 9] = 'South East'
vmodel$nhs_name[vmodel$nhs_name == 10] = 'South West'

sgtf[date >= "2021-08-15", qlo := NA]

sgtf_regions = ggplot(sgtf[(pid + 1) %in% c(1,3,4,5,6,9,10)]) +
    geom_ribbon(aes(x = date, ymin = qlo, ymax = qhi), fill = "black", alpha = 0.1) +
    geom_ribbon(data = vmodel[t + ymd("2020-01-01") >= "2020-10-01"], 
                aes(x = ymd("2020-01-01") + t, ymin = q025, ymax = q975), fill = "darkorchid", alpha = 0.5) +
    geom_line(data = vmodel[t + ymd("2020-01-01") >= "2020-10-01"], 
              aes(x = ymd("2020-01-01") + t, y = q500), colour = "darkorchid") +
    geom_line(aes(x = date, y = sgtf / (sgtf + other)), size = 0.25) +
    facet_wrap(~nhs_name, ncol = 3) +
    labs(x = NULL, y = "Relative frequency of S gene target failure") +
    scale_x_date(date_breaks = "3 months", date_labels = "%m/%y", limits = c(min(sgtf$date), '2021-09-30')) +
    cowplot::theme_cowplot(font_size = 12) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank())
sgtf_regions

sgtf_regions2 = sgtf_regions + 
    labs(x = NULL, y = "Relative frequency of\nS gene target failure")
sgtf_regions2
# Save final plot
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figS2_sgtf_fit_regional_", datetime, ".pdf"), sgtf_regions, width = 26, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figS2_sgtf_fit_regional_", datetime, ".png"), sgtf_regions, width = 26, height = 12, units = "cm")

################################################################################

# Figure S3 - regional model fits to Delta VOC
delta[, qlo := qbeta(0.025, delta + 1, other + 1)]
delta[, qhi := qbeta(0.975, delta + 1, other + 1)]

dmodel = f[, .(D1 = sum(test3_o), D2 = sum(test_o), D3 = sum(test2_o), delta0 = 0, conc = 30), by = .(t, population, run)]
dmodel[, dprop := D1 / (D1+D2+D3)]
dmodel[is.nan(dprop), dprop := 0]
dmodel[, delta := (1-dprop) * delta0 + dprop]
dmodel[, alpha := delta * (conc - 2) + 1]
dmodel[, beta := (1 - delta) * (conc - 2) + 1]
dmodel[, q025 := qbeta(0.025, alpha, beta)]
dmodel[, q500 := qbeta(0.500, alpha, beta)]
dmodel[, q975 := qbeta(0.975, alpha, beta)]
dmodel = dmodel[, lapply(.SD, mean), .SDcols = c("q025", "q500", "q975"), by = .(nhs_name = population, t)]
dmodel$nhs_name[dmodel$nhs_name == 1] = 'East of England'
dmodel$nhs_name[dmodel$nhs_name == 3] = 'London'
dmodel$nhs_name[dmodel$nhs_name == 4] = 'Midlands'
dmodel$nhs_name[dmodel$nhs_name == 5] = 'North East and Yorkshire'
dmodel$nhs_name[dmodel$nhs_name == 6] = 'North West'
dmodel$nhs_name[dmodel$nhs_name == 9] = 'South East'
dmodel$nhs_name[dmodel$nhs_name == 10] = 'South West'

delta[date >= "2021-09-14", qlo := NA]

delta_regions = ggplot(delta[(pid + 1) %in% c(1,3,4,5,6,9,10)]) +
    geom_ribbon(aes(x = date, ymin = qlo, ymax = qhi), fill = "black", alpha = 0.1) +
    geom_ribbon(data = dmodel[t + ymd("2020-01-01") >= "2021-01-01"], 
                aes(x = ymd("2020-01-01") + t, ymin = q025, ymax = q975), fill = "darkorchid", alpha = 0.5) +
    geom_line(data = dmodel[t + ymd("2020-01-01") >= "2021-01-01"], 
              aes(x = ymd("2020-01-01") + t, y = q500), colour = "darkorchid") +
    geom_line(aes(x = date, y = delta / (delta + other)), size = 0.25) +
    #scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99), limits = c(0.01, 0.99)) +
    facet_wrap(~nhs_name, ncol = 3) +
    cowplot::theme_cowplot(font_size = 12) +
    labs(x = NULL, y = "Relative frequency of Delta B.1.617.2 VOC") +
    scale_x_date(date_breaks = "2 months", date_labels = "%m/%y", limits = c(as.Date('2021-02-01'), '2021-09-30')) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank())
delta_regions

delta_regions2 = delta_regions + 
    labs(x = NULL, y = "Relative frequency\nof Delta B.1.617.2 VOC")
delta_regions2
# Save final plot
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figS3_delta_fit_regional_", datetime, ".pdf"), delta_regions, width = 26, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figS3_delta_fit_regional_", datetime, ".png"), delta_regions, width = 26, height = 12, units = "cm")

################################################################################

# Figure S4v2 - regional model fits to Omicron VOC

omi[, qlo := qbeta(0.025, sgtf + 1, other + 1)]
omi[, qhi := qbeta(0.975, sgtf + 1, other + 1)]

omodel = f[, .(IO = sum(test_o), INO = sum(test2_o + test3_o), sgtfomi0 = v4_sgtf0[1], conc = 1/(v4_disp[1]*v4_disp[1])), by = .(t, population, run)]
omodel[, o2 := IO / (IO + INO)]
omodel[is.nan(o2), o2 := 0]
omodel[, sgtf := (1 - o2) * sgtfomi0 + o2]
omodel[, alpha := sgtf * (conc - 2) + 1]
omodel[, beta := (1 - sgtf) * (conc - 2) + 1]
omodel[, q025 := qbeta(0.025, alpha, beta)]
omodel[, q500 := qbeta(0.500, alpha, beta)]
omodel[, q975 := qbeta(0.975, alpha, beta)]
omodel = omodel[, lapply(.SD, mean), .SDcols = c("q025", "q500", "q975"), by = .(nhs_name = population, t)]
omodel$nhs_name[omodel$nhs_name == 1] = 'East of England'
omodel$nhs_name[omodel$nhs_name == 3] = 'London'
omodel$nhs_name[omodel$nhs_name == 4] = 'Midlands'
omodel$nhs_name[omodel$nhs_name == 5] = 'North East and Yorkshire'
omodel$nhs_name[omodel$nhs_name == 6] = 'North West'
omodel$nhs_name[omodel$nhs_name == 9] = 'South East'
omodel$nhs_name[omodel$nhs_name == 10] = 'South West'

plotO = ggplot(omi[(pid + 1) %in% c(1,3,4,5,6,9,10)]) +
    geom_ribbon(aes(x = date, ymin = qlo, ymax = qhi), fill = "black", alpha = 0.1) +
    geom_ribbon(data = omodel[t + ymd("2020-01-01") >= "2021-10-01" & t + ymd("2020-01-01") < "2022-02-01"],
                aes(x = ymd("2020-01-01") + t, ymin = q025, ymax = q975), fill = "darkorchid", alpha = 0.5) +
    geom_line(data = omodel[t + ymd("2020-01-01") >= "2021-10-01" & t + ymd("2020-01-01") < "2022-02-01"],
              aes(x = ymd("2020-01-01") + t, y = q500), colour = "darkorchid") +
    geom_line(aes(x = date, y = sgtf / (sgtf + other)), size = 0.25) +
    #scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99), limits = c(0.01, 0.99)) +
    facet_wrap(~nhs_name) +
    cowplot::theme_cowplot(font_size = 12) +
    labs(x = NULL, y = "Relative frequency of S gene target failure") +
    scale_x_date(date_breaks = "1 month", date_labels = "%m/%y", limits = as.Date(c("2021-11-01", "2022-01-15"))) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank())
plotO

plotO2 = plotO +
    labs(x = NULL, y = "Relative frequency\nof S gene target failure")
plotO2
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figS4_omicron_fit_regional_", datetime, ".pdf"), plotO, width = 26, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figS4_omicron_fit_regional_", datetime, ".png"), plotO, width = 26, height = 12, units = "cm")

# join alpha, delta and omicron plots together

VOCplot = cowplot::plot_grid(sgtf_regions2, delta_regions2, plotO2, ncol = 1, labels = c('a','b','c'), align = 'hv')
VOCplot

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figVOCs_fit_regional_", datetime, ".pdf"), VOCplot, width = 20, height = 25, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figVOCs_fit_regional_", datetime, ".png"), VOCplot, width = 20, height = 25, units = "cm")


################################################################################

# Figure 3 - story figure showing key messages of paper

# see plot_story_fig.R / plot_story_fig_updated.R

################################################################################

# Functions to plot modelled burdens output and generate tables

burden_plots = function(bnlist, dstart = NULL, dend = NULL, pstart = '2021-05-01', pend = '2022-09-30', vline = 10000, dbreaks = '3 months', peak_lines = FALSE)
{
    ww = list()
    for (i in seq_along(bnlist))
    {
        basename = names(bnlist)[i]
        title = unname(bnlist[i])
        
        w = fread(paste0("./output/paper/may22/spim_", basename, ".csv"))
        w[, tag := title]
        ww[[i]] = w
    }
    ww = rbindlist(ww)
    ww = ww[variable %in% c("pcr_positive_i", "hosp_adm", "deaths")]
    ww[variable == "pcr_positive_i", name := "Infections"]
    ww[variable == "hosp_adm", name := "Admissions"]
    ww[variable == "deaths", name := "Deaths"]
    ww = ww[date >= "2021-05-01"]
    
    ww[, tag := factor(tag, bnlist)]
    
    if (peak_lines == FALSE){
        plot1 = function(ww, nm, col, thou)
        {
            if (thou == TRUE) {
                unit = 1000
                header = paste(nm, "(thousands)")
            } else {
                unit = 1
                header = nm
            }
            ggplot(ww[name == nm & date >= pstart & date <= pend]) + 
                geom_ribbon(aes(date, ymin = q05/unit, ymax = q95/unit), fill = col, alpha = 0.2) +
                geom_ribbon(aes(date, ymin = q25/unit, ymax = q75/unit), fill = col, alpha = 0.4) +
                geom_line(aes(date, q50/unit), colour = col) +
                geom_line(aes(date, va/unit), colour = col, linetype = "dashed") +
                geom_hline(aes(yintercept = 0)) +
                geom_vline(aes(xintercept = as.Date('2020-01-01')+vline), size = .5, linetype = "dotted") +
                scale_x_date(date_breaks = dbreaks, date_labels = "%m/%y", expand = c(0, 0)) +
                scale_y_continuous(expand = c(0.0, 0.0)) +
                labs(x = NULL, y = NULL, title = header) +
                facet_wrap(~tag, ncol = 1) +
                cowplot::theme_cowplot(font_size = 10) +
                theme(strip.background = element_blank())
        }
    } else if (peak_lines == TRUE) {
        plot1 = function(ww, nm, col, thou)
        {
            if (thou == TRUE) {
                unit = 1000
                header = paste(nm, "(thousands)")
            } else {
                unit = 1
                header = nm
            }
            ggplot(ww[name == nm & date >= pstart & date <= pend]) + 
                geom_ribbon(aes(date, ymin = q05/unit, ymax = q95/unit), fill = col, alpha = 0.2) +
                geom_ribbon(aes(date, ymin = q25/unit, ymax = q75/unit), fill = col, alpha = 0.4) +
                geom_line(aes(date, q50/unit), colour = col) +
                geom_line(aes(date, va/unit), colour = col, linetype = "dashed") +
                geom_hline(data = peaks[metric == nm], aes(yintercept = peak), linetype = "dashed", alpha = 0.5) +
                geom_hline(aes(yintercept = 0)) +
                geom_vline(aes(xintercept = as.Date('2020-01-01')+vline), size = .5, linetype = "dotted") +
                scale_x_date(date_breaks = dbreaks, date_labels = "%m/%y", expand = c(0, 0)) +
                scale_y_continuous(expand = c(0.0, 0.0)) +
                labs(x = NULL, y = NULL, title = header) +
                facet_wrap(~tag, ncol = 1) +
                cowplot::theme_cowplot(font_size = 10) +
                theme(strip.background = element_blank())
        }
    }
    
    
    intervention = function(gg, dstart, dend, height)
    {
        if (!is.null(dstart))
        {
            gg +
                annotate("ribbon", x = ymd(c(dstart, dend)), ymin = c(0, 0), ymax = c(height, height),
                         fill = "#000000", alpha = 0.1)
        }
        else
        {
            gg
        }
    }
    
    cowplot::plot_grid(
        intervention(plot1(ww, "Infections", "#bb88ff", TRUE), dstart, dend, 200),
        intervention(plot1(ww, "Admissions", "#3399bb", FALSE), dstart, dend, 2000),
        intervention(plot1(ww, "Deaths", "#aa2222", FALSE), dstart, dend, 400),
        nrow = 1,
        align = "hv"
    )
}

burden_plots_alt = function(bnlist, dstart = NULL, dend = NULL, pstart = '2021-05-01', pend = '2022-09-30', vline = 10000, dbreaks = '3 months')
{
    ww = list()
    for (i in seq_along(bnlist))
    {
        basename = names(bnlist)[i]
        title = unname(bnlist[i])
        
        w = fread(paste0("./output/paper/may22/spim_", basename, ".csv"))
        w[, tag := title]
        ww[[i]] = w
    }
    ww = rbindlist(ww)
    ww = ww[variable %in% c("pcr_positive_i", "hosp_adm", "deaths")]
    ww[variable == "pcr_positive_i", name := "Infections"]
    ww[variable == "hosp_adm", name := "Admissions"]
    ww[variable == "deaths", name := "Deaths"]
    ww = ww[date >= "2021-05-01" & date <= "2022-09-30"]
    
    ww[, tag := factor(tag, bnlist)]
    
    plot1 = function(ww, nm, col, thou)
    {
        if (thou == TRUE) {
            unit = 1000
            header = paste(nm, "(thousands)")
        } else {
            unit = 1
            header = nm
        }
        ggplot(ww[name == nm & date >= pstart & date <= pend]) + 
            geom_ribbon(aes(date, ymin = q05/unit, ymax = q95/unit), fill = col, alpha = 0.2) +
            geom_ribbon(aes(date, ymin = q25/unit, ymax = q75/unit), fill = col, alpha = 0.4) +
            geom_line(aes(date, q50/unit), colour = col) +
            geom_line(aes(date, va/unit), colour = col, linetype = "dashed") +
            geom_hline(aes(yintercept = 0)) +
            geom_vline(aes(xintercept = as.Date('2020-01-01')+vline), size = .5, linetype = "dotted") +
            scale_x_date(date_breaks = dbreaks, date_labels = "%m/%y", expand = c(0, 0)) +
            scale_y_continuous(expand = c(0.0, 0.0)) +
            labs(x = NULL, y = NULL, title = header) +
            facet_wrap(~tag, nrow = 1) +
            cowplot::theme_cowplot(font_size = 10) +
            theme(strip.background = element_blank())
    }
    
    intervention = function(gg, dstart, dend, height)
    {
        if (!is.null(dstart))
        {
            gg +
                annotate("ribbon", x = ymd(c(dstart, dend)), ymin = c(0, 0), ymax = c(height, height),
                         fill = "#000000", alpha = 0.1)
        }
        else
        {
            gg
        }
    }
    
    cowplot::plot_grid(
        intervention(plot1(ww, "Infections", "#bb88ff", TRUE), dstart, dend, 200),
        intervention(plot1(ww, "Admissions", "#3399bb", FALSE), dstart, dend, 2000),
        intervention(plot1(ww, "Deaths", "#aa2222", FALSE), dstart, dend, 400),
        ncol = 1,
        align = "hv"
    )
}

tables = function(bnlist, wh = "from_oct_21")
{
    ww = list()
    for (i in seq_along(bnlist))
    {
        basename = names(bnlist)[i]
        title = unname(bnlist[i])
        
        w = fread(paste0("./output/paper/may22/totals_", basename, ".csv"))
        w[, tag := title]
        pew = function(x) prettyNum(signif(x, 3), ",")
        w[, result := paste0(pew(median), " (", pew(q05), " - ", pew(q95), ")")]
        ww[[i]] = w
    }
    ww = rbindlist(ww)
    ww = ww[variable %in% c("pcr_positive_i", "hosp_adm", "deaths")]
    ww[variable == "pcr_positive_i", name := "Infections"]
    ww[variable == "hosp_adm", name := "Admissions"]
    ww[variable == "deaths", name := "Deaths"]
    ww = ww[, .(tag, name, when, result)]
    ww = ww[when == wh]
    ww[, tag := factor(tag, unique(tag))]
    ww[, name := factor(name, unique(name))]
    
    return (dcast(ww, tag ~ name, value.var = "result"))
}

################BEHAVIOUR PLOTS#################################################

# Figure showing different assumptions for future behaviour

# Strip of transmission potential
ov = rbind(
    get_overall(paper_aw_3wk, "3 weeks"),
    get_overall(paper_aw_3mo, "3 months"),
    get_overall(paper_aw_6mo, "6 months"),
    get_overall(paper_aw_flt, "None")
)

ov[, tag := factor(tag, c("3 weeks", "3 months", "6 months", "None"))]
ov[, tp := mob * adj]

mobility_strip = ggplot(ov[date >= "2021-10-01" & date <= "2022-12-31"]) +
    geom_line(aes(date, tp)) +
    geom_vline(aes(xintercept = as.Date('2020-01-01')+t_fitting), size = .5, linetype = "dotted") +
    facet_wrap(~tag, ncol = 1) +
    scale_x_date(date_breaks = "3 months", date_labels = "%m/%y", expand = c(0, 0)) +
    ylim(NA, 5.5) +
    labs(x = NULL, y = "Transmission potential", title = "Return to baseline") +
    cowplot::theme_cowplot(font_size = 10) +
    theme(panel.background = element_rect(fill = "#f4f4f4"), strip.background = element_blank(),
          panel.grid = element_line(colour = "#ffffff"))
mobility_strip

fig2 = cowplot::plot_grid(
    mobility_strip,
    burden_plots(c(mob_3wk_220511 = "3 weeks", mob_3mo_220511 = "3 months", basecase_220511 = "6 months", mob_flt_220511 = "None"),
                 pstart = '2021-10-01', pend = '2022-12-31', vline = t_fitting, dbreaks = '3 months'),
    nrow = 1,
    rel_widths = c(1, 3))
fig2
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/fig2_Oct21onwards_", datetime, ".pdf"), fig2, width = 24, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/fig2_Oct21onwards_", datetime, ".png"), fig2, width = 24, height = 12, units = "cm")

# now reproduce same plot but for January 2022 to December 2022 only

fig2B = cowplot::plot_grid(
    mobility_strip + scale_x_date(limits = c(as.Date('2022-01-01'), as.Date('2022-12-31')), date_breaks = "2 months", date_labels = "%m/%y", expand = c(0, 0)),
    burden_plots(c(mob_3wk_220511 = "3 weeks", mob_3mo_220511 = "3 months", basecase_220511 = "6 months", mob_flt_220511 = "None"),
                 pstart = '2022-01-01', pend = '2022-12-31', vline = t_fitting, dbreaks = '2 months'),
    nrow = 1,
    rel_widths = c(1, 3))
fig2B
ggsave(paste0("./output/paperfigs/may22/fig2_Jan22onwards_", datetime, ".pdf"), fig2B, width = 24, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/fig2_Jan22onwards_", datetime, ".png"), fig2B, width = 24, height = 12, units = "cm")

# now reproduce same plot but for March 2022 to December 2022 only
fig2C = cowplot::plot_grid(
    mobility_strip + scale_x_date(limits = c(as.Date('2022-03-01'), as.Date('2022-12-31')), date_breaks = "2 months", date_labels = "%m/%y", expand = c(0, 0)),
    burden_plots(c(mob_3wk_220511 = "3 weeks", mob_3mo_220511 = "3 months", basecase_220511 = "6 months", mob_flt_220511 = "None"),
                 pstart = '2022-03-01', pend = '2022-12-31', vline = t_fitting, dbreaks = '2 months'),
    nrow = 1,
    rel_widths = c(1, 3))
fig2C
ggsave(paste0("./output/paperfigs/may22/fig2_Mar22onwards_", datetime, ".pdf"), fig2C, width = 24, height = 16, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/fig2_Mar22onwards_", datetime, ".png"), fig2C, width = 24, height = 16, units = "cm")


# generate tables summarising total outcomes
bhv_tab = tables(c(mob_3wk_220511 = "3 weeks", mob_3mo_220511 = "3 months", basecase_220511 = "6 months", mob_flt_220511 = "None"))
bhv_tab
bhv_tabjn = tables(c(mob_3wk_220511 = "3 weeks", mob_3mo_220511 = "3 months", basecase_220511 = "6 months", mob_flt_220511 = "None"), wh = "from_jan_22")
bhv_tabjn

################BOOSTER ROLLOUT PLOTS###########################################

# Figure showing different assumptions for booster vaccination rollout

# add horizontal lines showing previous measured peak daily hospital admissions
# and deaths for England (from coronavirus.data.gov.uk dashboard)
dshbrd_data = as.data.table(content(GET("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newDeaths28DaysByDeathDate&metric=newAdmissions&format=csv")))
# 4134 COVID-19 hospital admissions on 12th January 2021 in England
# 1251 COVID-19 deaths (within 28 days of test)
peaks = data.table(peak = c(max(dshbrd_data$newAdmissions, na.rm=T), 
                            max(dshbrd_data$newDeaths28DaysByDeathDate, na.rm=T)),
                   metric = c("Admissions", "Deaths"))
peaks
figB = burden_plots(c(boostNO_220511 = "No boosters",
                      boost50_220511 = "50+ boosters",
                      basecase_220511 = "Actual boosters",
                      boostHI_220511 = "Higher uptake"), pstart = '2021-10-01', pend = '2022-12-31',
                    vline = t_fitting, dbreaks = '2 months', peak_lines = TRUE)
figB

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figBOOST_Oct21onwards_", datetime, ".pdf"), figB, width = 24, height = 16, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figBOOST_Oct21onwards_", datetime, ".png"), figB, width = 24, height = 16, units = "cm")

# do the same but without plotting peak for deaths
peaks = peaks[!(metric == 'Deaths'),]
figB = burden_plots(c(boostNO_220511 = "No boosters",
                      boost50_220511 = "50+ boosters",
                      basecase_220511 = "Actual boosters",
                      boostHI_220511 = "Higher uptake"), pstart = '2021-10-01', pend = '2022-12-31',
                    vline = t_fitting, dbreaks = '2 months', peak_lines = TRUE)
figB
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figBOOST_Oct21onwards_nodeathspeakline_", datetime, ".pdf"), figB, width = 24, height = 16, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figBOOST_Oct21onwards_nodeathspeakline_", datetime, ".png"), figB, width = 24, height = 16, units = "cm")

# generate tables summarising total outcomes
bst_tab = tables(c(boostNO_220511 = "No boosters",
                   boost50_220511 = "50+ boosters",
                   basecase_220511 = "Actual boosters",
                   boostHI_220511 = "Higher uptake"))
bst_tab
bst_tabjn = tables(c(boostNO_220511 = "No boosters",
                     boost50_220511 = "50+ boosters",
                     basecase_220511 = "Actual boosters",
                     boostHI_220511 = "Higher uptake"), wh = "from_jan_22")
bst_tabjn

################################################################################

# Figure showing different assumptions for waning immunity

figW = burden_plots(c(basecase_220511 = "Basecase waning",
                      waneHI_220511   = "Higher waning",
                      waneVHI_220511  = "Very high waning"), pstart = '2021-10-01', pend = '2022-12-31',
                    vline = t_fitting, dbreaks = '2 months')
figW

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figWANE_Oct21onwards_", datetime, ".pdf"), figW, width = 24, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figWANE_Oct21onwards_", datetime, ".png"), figW, width = 24, height = 12, units = "cm")

figWB = burden_plots(c(basecase_220511 = "Basecase waning",
                      waneHI_220511   = "Higher waning",
                      waneVHI_220511  = "Very high waning"), pstart = '2022-01-01', pend = '2022-12-31',
                    vline = t_fitting, dbreaks = '2 months')
figWB

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figWANE_Jan22onwards_", datetime, ".pdf"), figWB, width = 24, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figWANE_Jan22onwards_", datetime, ".png"), figWB, width = 24, height = 12, units = "cm")

figWC = burden_plots(c(basecase_220511 = "Basecase waning",
                       waneHI_220511   = "Higher waning",
                       waneVHI_220511  = "Very high waning"), pstart = '2022-03-01', pend = '2022-12-31',
                     vline = t_fitting, dbreaks = '2 months')
figWC

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figWANE_Mar22onwards_", datetime, ".pdf"), figWC, width = 24, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figWANE_Mar22onwards_", datetime, ".png"), figWC, width = 24, height = 12, units = "cm")


# generate tables summarising total outcomes
wan_tab = tables(c(basecase_220511 = "Basecase waning",
                   waneHI_220511   = "Higher waning",
                   waneVHI_220511  = "Very high waning"))
wan_tab
wan_tabjn = tables(c(basecase_220511 = "Basecase waning",
                     waneHI_220511   = "Higher waning",
                     waneVHI_220511  = "Very high waning"), wh = "from_jan_22")
wan_tabjn

wan_tabst = tables(c(basecase_220511 = "Basecase waning",
                   waneHI_220511   = "Higher waning",
                   waneVHI_220511  = "Very high waning"), wh = "from_oct_21_to_dec_21")
wan_tabst

################################################################################

# Figure showing different assumptions for seasonality

figS = burden_plots(c(seas_10_220511  = "10% seasonality",
                      basecase_220511 = "20% seasonality*",
                      seas_30_220511  = "30% seasonality",
                      seas_40_220511  = "40% seasonality"), pstart = '2021-10-01', pend = '2022-12-31',
                    vline = t_fitting, dbreaks = '2 months')
figS

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figSEAS_Oct21onwards_", datetime, ".pdf"), figS, width = 24, height = 16, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figSEAS_Oct21onwards_", datetime, ".png"), figS, width = 24, height = 16, units = "cm")

figSB = burden_plots(c(seas_10_220511  = "10% seasonality",
                      basecase_220511 = "20% seasonality*",
                      seas_30_220511  = "30% seasonality",
                      seas_40_220511  = "40% seasonality"), pstart = '2022-01-01', pend = '2022-12-31',
                    vline = t_fitting, dbreaks = '1 months')
figSB

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figSEAS_Jan22onwards_", datetime, ".pdf"), figSB, width = 24, height = 16, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figSEAS_Jan22onwards_", datetime, ".png"), figSB, width = 24, height = 16, units = "cm")

figSC = burden_plots(c(seas_10_220511  = "10% seasonality",
                       basecase_220511 = "20% seasonality*",
                       seas_30_220511  = "30% seasonality",
                       seas_40_220511  = "40% seasonality"), pstart = '2022-03-01', pend = '2022-12-31',
                     vline = t_fitting, dbreaks = '2 months')
figSC

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figSEAS_Mar22onwards_", datetime, ".pdf"), figSC, width = 24, height = 16, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figSEAS_Mar22onwards_", datetime, ".png"), figSC, width = 24, height = 16, units = "cm")


# generate tables summarising total outcomes
sea_tab = tables(c(seas_10_220511  = "10% seasonality",
                    basecase_220511 = "20% seasonality*",
                    seas_30_220511  = "30% seasonality",
                    seas_40_220511  = "40% seasonality"))
sea_tab
sea_tabjn = tables(c(seas_10_220511  = "10% seasonality",
                     basecase_220511 = "20% seasonality*",
                     seas_30_220511  = "30% seasonality",
                     seas_40_220511  = "40% seasonality"), wh = "from_jan_22")
sea_tabjn

sea_tabst = tables(c(seas_10_220511  = "10% seasonality",
                     basecase_220511 = "20% seasonality*",
                     seas_30_220511  = "30% seasonality",
                     seas_40_220511  = "40% seasonality"), wh = "from_oct_21_to_dec_21")
sea_tabst

################################################################################

# Figure showing different assumptions for vaccinating children

figVK = burden_plots(c(basecase_220511  = "5+, 80% uptake",
                       vax0550_220511   = "5+, 50% uptake"), pstart = '2021-10-01', pend = '2022-12-31',
                     vline = t_fitting, dbreaks = '2 months')
figVK

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figVAX_Oct21onwards_", datetime, ".pdf"), figVK, width = 24, height = 8, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figVAX_Oct21onwards_", datetime, ".png"), figVK, width = 24, height = 8, units = "cm")


figVKB = burden_plots(c(basecase_220511  = "5+, 80% uptake",
                       vax0550_220511    = "5+, 50% uptake"), pstart = '2022-01-01', pend = '2022-12-31',
                     vline = t_fitting, dbreaks = '2 months')
figVKB

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figVAX_Jan22onwards_", datetime, ".pdf"), figVKB, width = 24, height = 8, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figVAX_Jan22onwards_", datetime, ".png"), figVKB, width = 24, height = 8, units = "cm")

figVKC = burden_plots(c(basecase_220511  = "5+, 80% uptake",
                        vax0550_220511    = "5+, 50% uptake"), pstart = '2022-03-01', pend = '2022-12-31',
                      vline = t_fitting, dbreaks = '2 months')
figVKC

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figVAX_Mar22onwards_", datetime, ".pdf"), figVKC, width = 24, height = 8, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figVAX_Mar22onwards_", datetime, ".png"), figVKC, width = 24, height = 8, units = "cm")



# generate tables summarising total outcomes
vax_tab = tables(c(basecase_220511  = "5+, 80% uptake",
                   vax0550_220511   = "5+, 50% uptake"))
vax_tab
vax_tabjn = tables(c(basecase_220511  = "5+, 80% uptake",
                     vax0550_220511   = "5+, 50% uptake"), wh = "from_jan_22")
vax_tabjn

# add extra column to plots showing different vaccine schedules

two_scheds = read.csv('./output/paper/may22/vaccine_schedules_20220506094327.csv')
setDT(two_scheds)
# calculate sums over all NHS England regions

eng_sched = two_scheds[, .(num = sum(num)), by = .(date, age.group, dose, title)]

eng_sched$age[eng_sched$age.group == 1]  = '0-4'
eng_sched$age[eng_sched$age.group == 2]  = '5-9' 
eng_sched$age[eng_sched$age.group == 3]  = '10-14'
eng_sched$age[eng_sched$age.group == 4]  = '15-19'
eng_sched$age[eng_sched$age.group == 5]  = '20-24'
eng_sched$age[eng_sched$age.group == 6]  = '25-29'
eng_sched$age[eng_sched$age.group == 7]  = '30-34'
eng_sched$age[eng_sched$age.group == 8]  = '35-39'
eng_sched$age[eng_sched$age.group == 9]  = '40-44'
eng_sched$age[eng_sched$age.group == 10] = '45-49'
eng_sched$age[eng_sched$age.group == 11] = '50-54'
eng_sched$age[eng_sched$age.group == 12] = '55-59'
eng_sched$age[eng_sched$age.group == 13] = '60-64'
eng_sched$age[eng_sched$age.group == 14] = '65-69'
eng_sched$age[eng_sched$age.group == 15] = '70-74'
eng_sched$age[eng_sched$age.group == 16] = '75+' 
eng_sched$age = factor(eng_sched$age, levels = c('0-4'   ,
                                                 '5-9' ,
                                                 '10-14',
                                                 '15-19',
                                                 '20-24',
                                                 '25-29',
                                                 '30-34',
                                                 '35-39',
                                                 '40-44',
                                                 '45-49',
                                                 '50-54',
                                                 '55-59',
                                                 '60-64',
                                                 '65-69',
                                                 '70-74',
                                                 '75+' ))


path_to_newcovidvax  <- "."
uk_covid_data_path   <- paste0(path_to_newcovidvax, "/fitting_data/")
datapath = function(x) paste0(uk_covid_data_path, x)
popUK = readRDS(datapath("popNHS.rds"))
engpop = popUK[popUK$name == 'England',]
engpop$total = (engpop$f + engpop$m)*1000
over75s = sum(engpop$total[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
over75sf = sum(engpop$f[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
over75sm = sum(engpop$m[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
new_row = data.table(name = 'England',
                     age = '75+',
                     f = over75sf,
                     m = over75sm,
                     location_type = 5,
                     country_code = 2222,
                     total = over75s)
engpop = rbind(engpop, new_row)
engpop = engpop[!(engpop$age %in% c('75-79', '80-84', '85-89', '90+')),]
engpop$age.group = NULL
engpop$age.group[engpop$age == '0-4'] = 1
engpop$age.group[engpop$age == '5-9'] = 2
engpop$age.group[engpop$age == '10-14'] = 3
engpop$age.group[engpop$age == '15-19'] = 4
engpop$age.group[engpop$age == '20-24'] = 5
engpop$age.group[engpop$age == '25-29'] = 6
engpop$age.group[engpop$age == '30-34'] = 7
engpop$age.group[engpop$age == '35-39'] = 8
engpop$age.group[engpop$age == '40-44'] = 9
engpop$age.group[engpop$age == '45-49'] = 10
engpop$age.group[engpop$age == '50-54'] = 11
engpop$age.group[engpop$age == '55-59'] = 12
engpop$age.group[engpop$age == '60-64'] = 13
engpop$age.group[engpop$age == '65-69'] = 14
engpop$age.group[engpop$age == '70-74'] = 15
engpop$age.group[engpop$age == '75+'] = 16

all = merge(eng_sched[, .(date, age.group, dose, title, num, age)], engpop[, .(total, age.group)], by = "age.group")

all[,cumulative := cumsum(num), by = .(age.group, dose, title)]
all$coverage = (all$cumulative / all$total)*100


all$title[all$title == '5plus_80%'] = '5+, 80% uptake'
all$title[all$title == '5plus_50%'] = '5+, 50% uptake'
all$title = factor(all$title, levels = c( '5+, 80% uptake',
                                          '5+, 50% uptake'))

vax_schedule_col = ggplot(all[all$dose == 'First' & all$date > '2022-03-01' & all$date < '2022-12-31' & all$age.group %in% c(2,3,4),]) + 
    geom_line(aes(as.Date(date), coverage, colour=factor(age), group = age)) +
    facet_wrap(~title, ncol = 1) +
    scale_x_date(date_breaks = '2 months', date_labels = "%m/%y", expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = 'Vaccination coverage (%)', colour = 'Age group') +
    cowplot::theme_cowplot(font_size = 10) +
    theme(strip.background = element_blank(), 
          legend.position = c(0.7, 0.7),
          legend.background = element_rect(fill = "white", # Background
                                           colour = 1),
          legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
          legend.title = element_text(size=8), #change legend title font size
          legend.text = element_text(size=8)) +
    geom_vline(aes(xintercept = as.Date('2020-01-01')+t_fitting), size = .5, linetype = "dotted")
vax_schedule_col

figVKC_and_schedule = cowplot::plot_grid(vax_schedule_col, figVKC, nrow = 1, rel_widths = c(1,3))
figVKC_and_schedule

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figVAX_Mar22onwards_plusschedule_", datetime, ".pdf"), figVKC_and_schedule, width =24, height = 8, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figVAX_Mar22onwards_plusschedule_", datetime, ".png"), figVKC_and_schedule, width = 24, height = 8, units = "cm")

################################################################################

# Figure showing different assumptions for booster duration

figBD = burden_plots(c(shrtbst_220511   = "90 days",
                       basecase_220511  = "180 days",
                       longbst_220511   = "270 days"), pstart = '2021-10-01', pend = '2022-12-31',
                     vline = t_fitting, dbreaks = '2 months')
figBD

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figBOOSTDUR_Oct21onwards_", datetime, ".pdf"), figBD, width = 24, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figBOOSTDUR_Oct21onwards_", datetime, ".png"), figBD, width = 24, height = 12, units = "cm")


figBDB = burden_plots(c(shrtbst_220511   = "90 days",
                       basecase_220511  = "180 days",
                       longbst_220511   = "270 days"), pstart = '2022-01-01', pend = '2022-12-31',
                     vline = t_fitting, dbreaks = '2 months')
figBDB

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figBOOSTDUR_Jan22onwards_", datetime, ".pdf"), figBDB, width = 24, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figBOOSTDUR_Jan22onwards_", datetime, ".png"), figBDB, width = 24, height = 12, units = "cm")

figBDC = burden_plots(c(shrtbst_220511   = "90 days",
                        basecase_220511  = "180 days",
                        longbst_220511   = "270 days"), pstart = '2022-03-01', pend = '2022-12-31',
                      vline = t_fitting, dbreaks = '2 months')
figBDC

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figBOOSTDUR_Mar22onwards_", datetime, ".pdf"), figBDC, width = 24, height = 12, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figBOOSTDUR_Mar22onwards_", datetime, ".png"), figBDC, width = 24, height = 12, units = "cm")


# generate tables summarising total outcomes
bstdur_tab = tables(c(shrtbst_220511   = "90 days",
                      basecase_220511  = "180 days",
                      longbst_220511   = "270 days"))
bstdur_tab
bstdur_tabjn = tables(c(shrtbst_220511   = "90 days",
                        basecase_220511  = "180 days",
                        longbst_220511   = "270 days"), wh = "from_jan_22")
bstdur_tabjn

fwrite(bstdur_tabjn, "./output/paper/may22/table_booster_duration_Jan22onwards.csv")

bstdur_tabst = tables(c(shrtbst_220511   = "90 days",
                        basecase_220511  = "180 days",
                        longbst_220511   = "270 days"), wh = "from_oct_21_to_dec_21")
bstdur_tabst

# we also want to check the total outcomes between January and April/May 2022

# use the burdens files for each and generate a new_totals ... .csv file
list = c("shrtbst_220511", 
         "basecase_220511",
         "longbst_220511")

for (i in list){
    
    print(i)
    w0 = qread(paste0('./output/paper/may22/burdens_', i, '.qs'))
    
    tots = rbind(
        w0[,                                      .(when = "all", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2021-10-01", .(when = "from_oct_21", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2021-10-01" & ymd("2020-01-01") + t <= "2021-12-31",
           .(when = "from_oct_21_to_dec_21", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2021-12-01" & ymd("2020-01-01") + t <= "2022-04-30",
           .(when = "from_dec_21_to_apr_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-01-01", .(when = "from_jan_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-01-01" & ymd("2020-01-01") + t <= "2022-05-31", .(when = "from_jan_22_to_may_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-06-01" & ymd("2020-01-01") + t <= "2022-12-31", .(when = "from_jun_22_to_dec_22", tot = sum(value)), by = .(variable, run)]
    )
    
    tots = tots[, .(
        median = median(tot),
        mean = mean(tot),
        q05 = quantile(tot, 0.05),
        q25 = quantile(tot, 0.25),
        q75 = quantile(tot, 0.75),
        q95 = quantile(tot, 0.95)
    ), by = .(variable, when)]
    
    fwrite(tots, paste0("./output/paper/may22/newtotals_", i, ".csv"))
}

newtables = function(bnlist, wh = "from_oct_21")
{
    ww = list()
    for (i in seq_along(bnlist))
    {
        basename = names(bnlist)[i]
        title = unname(bnlist[i])
        
        w = fread(paste0("./output/paper/may22/newtotals_", basename, ".csv"))
        w[, tag := title]
        pew = function(x) prettyNum(signif(x, 3), ",")
        w[, result := paste0(pew(median), " (", pew(q05), " - ", pew(q95), ")")]
        ww[[i]] = w
    }
    ww = rbindlist(ww)
    ww = ww[variable %in% c("pcr_positive_i", "hosp_adm", "deaths")]
    ww[variable == "pcr_positive_i", name := "Infections"]
    ww[variable == "hosp_adm", name := "Admissions"]
    ww[variable == "deaths", name := "Deaths"]
    ww = ww[, .(tag, name, when, result)]
    ww = ww[when == wh]
    ww[, tag := factor(tag, unique(tag))]
    ww[, name := factor(name, unique(name))]
    
    return (dcast(ww, tag ~ name, value.var = "result"))
}

bstdur_tab_JanMay22 = newtables(c(shrtbst_220511   = "90 days",
                        basecase_220511  = "180 days",
                        longbst_220511   = "270 days"), wh = "from_jan_22_to_may_22")
bstdur_tab_JanMay22

bstdur_tab_JunDec22 = newtables(c(shrtbst_220511   = "90 days",
                                  basecase_220511  = "180 days",
                                  longbst_220511   = "270 days"), wh = "from_jun_22_to_dec_22")
bstdur_tab_JunDec22
################################################################################

# Function to generate immunity traces 

immunity_trace = function(basename, tag)
{
    w = qread(paste0("./output/paper/may22/output1_", basename, ".qs"))
    w = w[run == 1]
    w = w[, .(
        sus = sum(S),
        inf = sum(E+Ip+Is+Ia+E2+Ip2+Is2+Ia2+E3+Ip3+Is3+Ia3),
        rec = sum(R + R2 + R3),
        vd1 = sum(Va1 + Vb1),
        vd2 = sum(Va2 + Vb2),
        vwn = sum(Va3 + Vb3),
        ev  = sum(EV)
    ), by = .(t)]
    w[, date := ymd("2020-01-01") + t]
    w[, tag := ..tag]
    w = melt(w[], id.vars = c("t", "date", "tag"))
    w[variable == "sus", name := "Susceptible"]
    w[variable == "inf", name := "Infected"]
    w[variable == "rec", name := "Natural protection"]
    w[variable == "vd1", name := "Vax 1-dose"]
    w[variable == "vd2", name := "Vax 2-dose"]
    w[variable == "vwn", name := "Vax waned"]
    w[variable == "ev",  name := "Ever vaccinated"]
    w[, name := factor(name, unique(name))]
    return (w[])
}

immunity_trace_some_ages = function(basename, tag, age_groups)
{
    w = qread(paste0("./output/paper/may22/output1_", basename, ".qs"))
    w = w[run == 1]
    w = w[group %in% age_groups]
    w = w[, .(
        sus = sum(S),
        inf = sum(E+Ip+Is+Ia+E2+Ip2+Is2+Ia2+E3+Ip3+Is3+Ia3),
        rec = sum(R + R2 + R3),
        vd1 = sum(Va1 + Vb1),
        vd2 = sum(Va2 + Vb2),
        vwn = sum(Va3 + Vb3),
        ev  = sum(EV)
    ), by = .(t)]
    w[, date := ymd("2020-01-01") + t]
    w[, tag := ..tag]
    w = melt(w[], id.vars = c("t", "date", "tag"))
    w[variable == "sus", name := "Susceptible"]
    w[variable == "inf", name := "Infected"]
    w[variable == "rec", name := "Natural protection"]
    w[variable == "vd1", name := "Vax 1-dose"]
    w[variable == "vd2", name := "Vax 2-dose"]
    w[variable == "vwn", name := "Vax waned"]
    w[variable == "ev",  name := "Ever vaccinated"]
    w[, name := factor(name, unique(name))]
    return (w[])
}

immunity_trace_byage = function(basename, tag)
{
    w = qread(paste0("./output/paper/may22/output1_", basename, ".qs"))
    w = w[run == 1]
    # w = w[group %in% age_groups]
    w = w[, .(
        sus = sum(S),
        inf = sum(E+Ip+Is+Ia+E2+Ip2+Is2+Ia2+E3+Ip3+Is3+Ia3),
        rec = sum(R + R2 + R3),
        vd1 = sum(Va1 + Vb1),
        vd2 = sum(Va2 + Vb2),
        vwn = sum(Va3 + Vb3)
    ), by = .(t, group)]
    w[, date := ymd("2020-01-01") + t]
    w[, tag := ..tag]
    w = melt(w[], id.vars = c("t", "date", "group", "tag"))
    w[variable == "sus", name := "Susceptible"]
    w[variable == "inf", name := "Infected"]
    w[variable == "rec", name := "Natural protection"]
    w[variable == "vd1", name := "Vax 1-dose"]
    w[variable == "vd2", name := "Vax 2-dose"]
    w[variable == "vwn", name := "Vax waned"]
    w[, name := factor(name, unique(name))]
    return (w[])
}

################################################################################

# Figure showing immunity profiles by age group for the basecase scenario

it_bc = immunity_trace_byage("basecase_220511", "Basecase scenario")

it_bc[group == 1, age_group := '0-4']
it_bc[group == 2, age_group := '5-9']
it_bc[group == 3, age_group := '10-14']
it_bc[group == 4, age_group := '15-19']
it_bc[group == 5, age_group := '20-24']
it_bc[group == 6, age_group := '25-29']
it_bc[group == 7, age_group := '30-34']
it_bc[group == 8, age_group := '35-39']
it_bc[group == 9, age_group := '40-44']
it_bc[group == 10, age_group := '45-49']
it_bc[group == 11, age_group := '50-54']
it_bc[group == 12, age_group := '55-59']
it_bc[group == 13, age_group := '60-64']
it_bc[group == 14, age_group := '65-69']
it_bc[group == 15, age_group := '70-74']
it_bc[group == 16, age_group := '75+']

it_bc$age_group = factor(it_bc$age_group, levels = c('0-4','5-9','10-14','15-19',
                                                     '20-24','25-29','30-34','35-39',
                                                     '40-44','45-49','50-54','55-59',
                                                     '60-64','65-69','70-74','75+'))

it_full = ggplot(it_bc) +
    geom_area(aes(date, value, fill = name), position = position_fill()) +
    scale_x_date(date_breaks = "8 months", date_labels = "%m/%y", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c(
        "Susceptible" = "#eeeeee",
        "Infected" = "#ff0000",
        "Natural protection" = "#eedd99",
        "Vax 1-dose" = "#aaaaff",
        "Vax 2-dose" = "#6666ff",
        "Vax waned" = "#ffaaff",
        "Ever vaccinated" = "#aaffaa")) +
    facet_wrap(~age_group, ncol = 4) +
    labs(x = NULL, y = "Proportion of population", fill = NULL, title = "Basecase") +
    cowplot::theme_cowplot(font_size = 9) +
    theme(strip.background = element_blank(), legend.position = c(0.03, 0.92), 
          legend.background = element_rect(fill = "#ffffff", colour = "#000000"), legend.margin = margin(1, 3, 3, 3),
          legend.text = element_text(size = 7), legend.key.size = unit(5, units = "pt"))
it_full

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/figIMMTRACE_full_", datetime, ".pdf"), it_full, width = 24, height = 18, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/may22/figIMMTRACE_full_", datetime, ".png"), it_full, width = 24, height = 18, units = "cm")

# it_jan22 = it_full + scale_x_date(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0), limits = as.Date(c('2020-01-01', '2022-01-31')))
# it_jan22
# 
# datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
# ggsave(paste0("./output/paperfigs/jan22/figIMMTRACE_tojan22_", datetime, ".pdf"), it_jan22, width = 24, height = 18, units = "cm", useDingbats = FALSE)
# ggsave(paste0("./output/paperfigs/jan22/figIMMTRACE_tojan22_", datetime, ".png"), it_jan22, width = 24, height = 18, units = "cm")

################################################################################

# Figure showing gamma multiplier adjustments for the basecase scenario

adjustment_fig = function(codenm, vlinedate)
{
    outfile = paste0("./output/paper/may22/output1_", codenm, ".qs");
    adjfile = paste0("./output/paper/may22/adjusted1_", codenm, ".qs");
    wo = qread(outfile)
    wa = qread(adjfile)
    
    ir = merge(
        wo[, .(t, NHS.region, run, group, deaths_o = deaths, hosp_adm_o = hosp_adm, hosp_bed_o = hosp_bed - hosp_undetected_p, icu_bed_o = icu_bed, pcr_incidence_o = pcr_positive_i)], 
        wa[, .(t, NHS.region, run, group, deaths_a = deaths, hosp_adm_a = hosp_adm, hosp_bed_a = hosp_bed - hosp_undetected_p, icu_bed_a = icu_bed, pcr_incidence_a = pcr_positive_i)],
        by = c("t", "NHS.region", "run", "group"))
    
    ir = ir[, .(
        incidence_o = sum(pcr_incidence_o),
        incidence_a = sum(pcr_incidence_a),
        deaths_o = sum(deaths_o),
        deaths_a = sum(deaths_a),
        hosp_adm_o = sum(hosp_adm_o),
        hosp_adm_a = sum(hosp_adm_a),
        hosp_bed_o = sum(hosp_bed_o),
        hosp_bed_a = sum(hosp_bed_a),
        icu_bed_o = sum(icu_bed_o),
        icu_bed_a = sum(icu_bed_a)),
        by = .(t, NHS.region, run)]
    
    ir[, date := t + ymd("2020-01-01")]
    vline = data.table(date = ymd(vlinedate))
    
    
    adj1 = ggplot(ir[date >= "2020-03-01"]) +
        geom_vline(data = vline, aes(xintercept = date), linetype = "dashed") +
        geom_line(aes(date, deaths_a / deaths_o, colour = NHS.region, group = paste(NHS.region, run))) +
        cowplot::theme_cowplot(font_size = 10) +
        theme(legend.position = c(0.75, 0.6), legend.text = element_text(size = 9)) +
        ylim(0, 5) +
        scale_x_date(date_breaks = "3 months", date_labels = "%m/%y") +
        labs(x = NULL, y = "Fatality rate\nadjustment", colour = "NHS England\nregion")
    
    adj2 = ggplot(ir[date >= "2020-03-01"]) +
        geom_vline(data = vline, aes(xintercept = date), linetype = "dashed") +
        geom_line(aes(date, hosp_adm_a / hosp_adm_o, colour = NHS.region, group = paste(NHS.region, run))) +
        cowplot::theme_cowplot(font_size = 10) +
        theme(legend.position = "none") +
        ylim(0, 5) +
        scale_x_date(date_breaks = "3 months", date_labels = "%m/%y") +
        labs(x = NULL, y = "Admissions rate\nadjustment", colour = "NHS England\nregion")
    
    adj3 = ggplot(ir[date >= "2020-03-01"]) +
        geom_vline(data = vline, aes(xintercept = date), linetype = "dashed") +
        geom_line(aes(date, hosp_bed_a / hosp_bed_o, colour = NHS.region, group = paste(NHS.region, run))) +
        cowplot::theme_cowplot(font_size = 10) +
        theme(legend.position = "none") +
        ylim(0, 5) +
        scale_x_date(date_breaks = "3 months", date_labels = "%m/%y") +
        labs(x = NULL, y = "Hospital bed occupancy\nadjustment", colour = "NHS England\nregion")
    
    adj4 = ggplot(ir[date >= "2020-03-01"]) +
        geom_vline(data = vline, aes(xintercept = date), linetype = "dashed") +
        geom_line(aes(date, icu_bed_a / icu_bed_o, colour = NHS.region, group = paste(NHS.region, run))) +
        cowplot::theme_cowplot(font_size = 10) +
        theme(legend.position = "none") +
        ylim(0, 5) +
        scale_x_date(date_breaks = "3 months", date_labels = "%m/%y") +
        labs(x = NULL, y = "ICU bed occupancy\nadjustment", colour = "NHS England\nregion")
    
    adj_plot = cowplot::plot_grid(adj1, adj2, adj3, adj4, nrow = 4, labels = letters[1:4], label_size = 12)
    datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
    ggsave(paste0("./output/paperfigs/may22/supp_adj_", codenm, datetime, ".png"), adj_plot, width = 25, height = 25, units = "cm")
    ggsave(paste0("./output/paperfigs/may22/supp_adj_", codenm, datetime, ".pdf"), adj_plot, width = 25, height = 25, units = "cm", useDingbats = FALSE)
    return(adj_plot)
}

adj_plot = adjustment_fig("basecase_220511", date_fitting)

################################################################################

# Figure showing immunity traces for booster vaccination scenarios

itbstbc = immunity_trace("basecase_220511", "Actual booster uptake*")
itbstHI = immunity_trace("boostHI_220511", "High booster uptake")
itbstNO = immunity_trace("boostNO_220511", "No booster uptake")
itbst50 = immunity_trace("boost50_220511", "Boosters for 50+")

itbst_all = rbind(itbstbc, itbstHI, itbstNO, itbst50)
itbst_all[, tag := factor(tag, rev(unique(tag)))]

ggplot(itbst_all[variable != 'ev' & date > '2021-10-01']) +
    geom_area(aes(date, value, fill = name), position = position_fill()) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n'%y", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c(
        "Susceptible" = "#eeeeee",
        "Infected" = "#ff0000",
        "Natural protection" = "#eedd99",
        "Vax 1-dose" = "#aaaaff",
        "Vax 2-dose" = "#6666ff",
        "Vax waned" = "#ffaaff",
        "Ever vaccinated" = "#aaffaa")) +
    facet_wrap(~tag, ncol = 1) +
    labs(x = NULL, y = "Proportion of population", fill = NULL, title = "Booster uptake") +
    cowplot::theme_cowplot(font_size = 9) +
    theme(strip.background = element_blank(), legend.position = c(0.8, 0.92),
          legend.background = element_rect(fill = "#ffffff", colour = "#000000"), legend.margin = margin(1, 3, 3, 3),
          legend.text = element_text(size = 7), legend.key.size = unit(5, units = "pt"))

# plot lines showing proportion in each state over time

itbst_all$popsize = itbst_all$value[itbst_all$t==0 & 
                    itbst_all$tag == 'Actual booster uptake*' &
                    itbst_all$variable == 'sus']

itbst_all$prop = itbst_all$value / itbst_all$popsize

itbst_all$tag = factor(itbst_all$tag, levels = c('No booster uptake',
                                                   'Boosters for 50+',
                                                   'Actual booster uptake*',
                                                   'High booster uptake'))

itbst_scenarios = ggplot(itbst_all[date > '2021-10-01' & variable %in% c('sus', 'inf', 'rec', 'vd2')]) +
    geom_line(aes(date, prop, colour = tag, group = tag)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%m/%y", expand = c(0, 0)) +
    facet_wrap(~name, ncol = 1, scales = 'free_y') +
    labs(x = NULL, y = "Proportion of population", colour = 'Scenario', title = "Booster vaccination scenarios") +
    cowplot::theme_cowplot(font_size = 9) +
    theme(strip.background = element_blank(), legend.position = c(0.7, 0.3),
          legend.background = element_rect(fill = "#ffffff", colour = "#000000"), legend.margin = margin(1, 3, 3, 3),
          legend.text = element_text(size = 7), legend.key.size = unit(5, units = "pt"))
itbst_scenarios
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/may22/supp_bst_scenarios_immunity_trace_", datetime, ".png"), itbst_scenarios, width = 12, height = 15, units = "cm")
ggsave(paste0("./output/paperfigs/may22/supp_bst_scenarios_immunity_trace_", datetime, ".pdf"), itbst_scenarios, width = 12, height = 15, units = "cm", useDingbats = FALSE)

################################################################################

# Figure showing immunity traces for vaccination of 5+, 12+ children scenarios

itbc =      immunity_trace("basecase_220214", "12+; 80%*")
it12pl_50 = immunity_trace("vax1250_220214",  "12+; 50%")
it5pl_80 =  immunity_trace("vax0580_220214",  "5+; 80%")
it5pl_50 =  immunity_trace("vax0550_220214",  "5+; 50%")

it_all = rbind(itbc, it12pl_50, it5pl_80, it5pl_50)
it_all[, tag := factor(tag, rev(unique(tag)))]

ggplot(it_all[variable != 'ev' & date > '2022-01-01']) +
    geom_area(aes(date, value, fill = name), position = position_fill()) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n'%y", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c(
        "Susceptible" = "#eeeeee",
        "Infected" = "#ff0000",
        "Natural protection" = "#eedd99",
        "Vax 1-dose" = "#aaaaff",
        "Vax 2-dose" = "#6666ff",
        "Vax waned" = "#ffaaff",
        "Ever vaccinated" = "#aaffaa")) +
    facet_wrap(~tag, ncol = 1) +
    labs(x = NULL, y = "Proportion of population", fill = NULL, title = "Booster uptake") +
    cowplot::theme_cowplot(font_size = 9) +
    theme(strip.background = element_blank(), legend.position = c(0.03, 0.92),
          legend.background = element_rect(fill = "#ffffff", colour = "#000000"), legend.margin = margin(1, 3, 3, 3),
          legend.text = element_text(size = 7), legend.key.size = unit(5, units = "pt"))

# try looking at younger age groups only (ages 0-19)

Yitbc =      immunity_trace_some_ages("basecase_220214", "12+; 80%*", c(1,2,3,4))
Yit12pl_50 = immunity_trace_some_ages("vax1250_220214",  "12+; 50%", c(1,2,3,4))
Yit5pl_80 =  immunity_trace_some_ages("vax0580_220214",  "5+; 80%", c(1,2,3,4))
Yit5pl_50 =  immunity_trace_some_ages("vax0550_220214",  "5+; 50%", c(1,2,3,4))

Yit_all = rbind(Yitbc, Yit12pl_50, Yit5pl_80, Yit5pl_50)
Yit_all[, tag := factor(tag, rev(unique(tag)))]

ggplot(Yit_all[variable != 'ev' & date > '2022-01-01']) +
    geom_area(aes(date, value, fill = name), position = position_fill()) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n'%y", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c(
        "Susceptible" = "#eeeeee",
        "Infected" = "#ff0000",
        "Natural protection" = "#eedd99",
        "Vax 1-dose" = "#aaaaff",
        "Vax 2-dose" = "#6666ff",
        "Vax waned" = "#ffaaff",
        "Ever vaccinated" = "#aaffaa")) +
    facet_wrap(~tag, ncol = 1) +
    labs(x = NULL, y = "Proportion of population", fill = NULL, title = "Booster uptake") +
    cowplot::theme_cowplot(font_size = 9) +
    theme(strip.background = element_blank(), legend.position = c(0.03, 0.92),
          legend.background = element_rect(fill = "#ffffff", colour = "#000000"), legend.margin = margin(1, 3, 3, 3),
          legend.text = element_text(size = 7), legend.key.size = unit(5, units = "pt"))

# plot lines showing proportion in each state over time

Yit_all_noEV = Yit_all[Yit_all$variable != 'ev']

pop_size = Yit_all_noEV$value[Yit_all_noEV$t==1 & 
                                  Yit_all_noEV$tag == '12+; 80%*' & 
                                  Yit_all_noEV$variable == 'sus']
Yit_all_noEV$popsize = pop_size

ggplot(Yit_all_noEV[date > '2022-01-01' & variable %in% c('sus', 'rec', 'vd1', 'vd2')]) +
    geom_line(aes(date, value/pop_size, colour = tag, group = tag)) +
    # geom_area(aes(date, value, fill = name), position = position_fill()) +
    scale_x_date(date_breaks = "1 month", date_labels = "%m/%y", expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0)) +
    # scale_fill_manual(values = c(
    #     "Susceptible" = "#eeeeee",
    #     "Infected" = "#ff0000",
    #     "Natural protection" = "#eedd99",
    #     "Vax 1-dose" = "#aaaaff",
    #     "Vax 2-dose" = "#6666ff",
    #     "Vax waned" = "#ffaaff",
    #     "Ever vaccinated" = "#aaffaa")) +
    facet_wrap(~name, ncol = 1, scales = 'free_y') +
    labs(x = NULL, y = "Proportion of population aged 0-19", colour = 'Scenario', title = "Vaccination scenarios") +
    # cowplot::theme_cowplot(font_size = 9) +
    theme(strip.background = element_blank(), legend.position = c(0.03, 0.92),
          legend.background = element_rect(fill = "#ffffff", colour = "#000000"), legend.margin = margin(1, 3, 3, 3),
          legend.text = element_text(size = 7), legend.key.size = unit(5, units = "pt"))









######################################### extra totals calculation 13th May 2022

list = c("basecase_220511",
         "mob_3wk_220511",
         "mob_3mo_220511",
         "mob_flt_220511",
         "waneHI_220511", 
         "waneVHI_220511",
         "seas_10_220511",
         "seas_30_220511",
         "seas_40_220511",
         "vax0550_220511",
         "shrtbst_220511",
         "longbst_220511",
         "boostHI_220511",
         "boostNO_220511",
         "boost50_220511")

for (i in list){
    
    print(i)
    w0 = qread(paste0('./output/paper/may22/burdens_', i, '.qs'))
    
    tots = rbind(
        w0[,                                      .(when = "all", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2021-10-01", .(when = "from_oct_21", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2021-10-01" & ymd("2020-01-01") + t <= "2021-12-31",
           .(when = "from_oct_21_to_dec_21", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2021-12-01" & ymd("2020-01-01") + t <= "2022-04-30",
           .(when = "from_dec_21_to_apr_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-01-01", .(when = "from_jan_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-01-01" & ymd("2020-01-01") + t <= "2022-05-31", .(when = "from_jan_22_to_may_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-06-01" & ymd("2020-01-01") + t <= "2022-12-31", .(when = "from_jun_22_to_dec_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-03-01" & ymd("2020-01-01") + t <= "2022-12-31", .(when = "from_mar_22_to_dec_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-04-01" & ymd("2020-01-01") + t <= "2022-12-31", .(when = "from_apr_22_to_dec_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-05-01" & ymd("2020-01-01") + t <= "2022-12-31", .(when = "from_may_22_to_dec_22", tot = sum(value)), by = .(variable, run)]
    )
    
    tots = tots[, .(
        median = median(tot),
        mean = mean(tot),
        q05 = quantile(tot, 0.05),
        q25 = quantile(tot, 0.25),
        q75 = quantile(tot, 0.75),
        q95 = quantile(tot, 0.95)
    ), by = .(variable, when)]
    
    fwrite(tots, paste0("./output/paper/may22/extratotals_", i, ".csv"))
}

extratables = function(bnlist, wh = "from_oct_21")
{
    ww = list()
    for (i in seq_along(bnlist))
    {
        basename = names(bnlist)[i]
        title = unname(bnlist[i])
        
        w = fread(paste0("./output/paper/may22/extratotals_", basename, ".csv"))
        w[, tag := title]
        pew = function(x) prettyNum(signif(x, 3), ",")
        w[, result := paste0(pew(median), " (", pew(q05), " - ", pew(q95), ")")]
        ww[[i]] = w
    }
    ww = rbindlist(ww)
    ww = ww[variable %in% c("pcr_positive_i", "hosp_adm", "deaths")]
    ww[variable == "pcr_positive_i", name := "Infections"]
    ww[variable == "hosp_adm", name := "Admissions"]
    ww[variable == "deaths", name := "Deaths"]
    ww = ww[, .(tag, name, when, result)]
    ww = ww[when == wh]
    ww[, tag := factor(tag, unique(tag))]
    ww[, name := factor(name, unique(name))]
    
    return (dcast(ww, tag ~ name, value.var = "result"))
}

behavMay22 = extratables(c(mob_3wk_220511 = "3 weeks",
                           mob_3mo_220511 = "3 months",
                           basecase_220511 = "6 months*",
                           mob_flt_220511 = "No change"), wh = "from_may_22_to_dec_22")
behavMay22
fwrite(behavMay22, "./output/paper/may22/tab_behav.csv")

behavMar22 = extratables(c(mob_3wk_220511 = "3 weeks",
                           mob_3mo_220511 = "3 months",
                           basecase_220511 = "6 months*",
                           mob_flt_220511 = "No change"), wh = "from_mar_22_to_dec_22")
behavMar22

waningMay22 = extratables(c(basecase_220511 = "Basecase waning*",
                            waneHI_220511 = 'High', 
                            waneVHI_220511 = 'Very high'), wh = "from_may_22_to_dec_22")
waningMay22

waningMar22 = extratables(c(basecase_220511 = "Basecase waning*",
                            waneHI_220511 = 'High', 
                            waneVHI_220511 = 'Very high'), wh = "from_mar_22_to_dec_22")
waningMar22

seasMay22 = extratables(c(seas_10_220511 = '10%',
                          basecase_220511 = "20%*",
                          seas_30_220511 = '30%',
                          seas_40_220511 = '40%'), wh = "from_may_22_to_dec_22")
seasMay22

seasMar22 = extratables(c(seas_10_220511 = '10%',
                          basecase_220511 = "20%*",
                          seas_30_220511 = '30%',
                          seas_40_220511 = '40%'), wh = "from_mar_22_to_dec_22")
seasMar22

vaxMay22 = extratables(c(basecase_220511 = "5+, 80% uptake",
                          vax0550_220511 = '5+, 50% uptake'), wh = "from_may_22_to_dec_22")
vaxMay22

