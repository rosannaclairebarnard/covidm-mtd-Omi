# testing script to produce story figure (Figure 3 in main manuscript)

library(ggplot2)
library(data.table)
library(stringr)
library(qs)
library(colorspace)
library(lubridate)
library(cowplot)
theme_set(theme_classic())
theme_set(theme_cowplot(font_size=7))
options(scipen = 999)

# variables to choose from for plotting are "pcr_positive_i" "pcr_positive_p" "hosp_adm"       "hosp_bed"       "icu_bed"        "deaths"

# # list of scenarios that were run for resubmission, May 2022

# first set of scenarios for resubmission following Omicron updates (11th May 2022)
# scenario_sheet = fread(
#     "codenm,          mobility,                    wane_yn, vac_eff, wane_scen, bdur, bage, pby,  pbo,  pbst_all, seas, cert, mask,  wfh_date,   title,                                                                           vax_file,                                      adj_ts
# basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Base case: 6-month relaxation; central waning; actual boosters; 20% seasonality, vax-covidm20220505205235.rds,                  NA
# mob_3wk_220511,       paper_aw_3wk_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         3-week relaxation,                                                               vax-covidm20220505205235.rds,                  NA
# mob_3mo_220511,       paper_aw_3mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         3-month relaxation,                                                              vax-covidm20220505205235.rds,                  NA
# mob_flt_220511,       paper_aw_flt_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         No relaxation,                                                                   vax-covidm20220505205235.rds,                  NA
# waneHI_220511,        paper_aw_6mo_20220509105131, yeswane, central, hiwane,    180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         High waning,                                                                     vax-covidm20220505205235.rds,                  NA
# waneVHI_220511,       paper_aw_6mo_20220509105131, yeswane, central, vhiwane,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Very high waning,                                                                vax-covidm20220505205235.rds,                  NA
# seas_10_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.05,  0,    0,     NA,         10% seasonality,                                                                 vax-covidm20220505205235.rds,                  NA
# seas_30_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.15,  0,    0,     NA,         30% seasonality,                                                                 vax-covidm20220505205235.rds,                  NA
# seas_40_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.20,  0,    0,     NA,         40% seasonality,                                                                 vax-covidm20220505205235.rds,                  NA
# vax0550_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Vax 5+ 50% uptake,                                                               vax-covidm202205052222505plus_50percent.rds,   NA
# shrtbst_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,    90,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Short booster duration (90 days),                                                vax-covidm20220505205235.rds,                  NA
# longbst_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   270,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Long booster duration (270 days),                                                vax-covidm20220505205235.rds,                  NA
# ")

# # second set of scenarios for resubmission following Omicron updates (11th May 2022)
# scenario_sheet = fread(
#     "codenm,          mobility,                    wane_yn, vac_eff, wane_scen, bdur, bage, pby,  pbo,  pbst_all, seas, cert, mask,  wfh_date,   title,                                                                           vax_file,                                      adj_ts
# boostHI_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.90, 0.98, NULL,     0.1,  0,    0,     NA,         High booster uptake: 90% in under 50s and 98% in 50+,                            vax-covidm20220505205235.rds,                  ./fits/final_adj_tseries_basecase_220511.csv
# boostNO_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,    0.0,  0.0, NULL,     0.1,  0,    0,     NA,         No booster uptake,                                                               vax-covidm20220505205235.rds,                  ./fits/final_adj_tseries_basecase_220511.csv
# boost50_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,    0.0, 0.95, NULL,     0.1,  0,    0,     NA,         Boosters for 50+ only: 95% uptake,                                               vax-covidm20220505205235.rds,                  ./fits/final_adj_tseries_basecase_220511.csv
# ")

###############################################################################

# basecase model fit (aggregated to England) and projection panel

popUK = readRDS("./fitting_data/popNHS.rds");
nhs_regions = popUK[, unique(name)]

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
# not sure if ba2 is needed here.... but including it just in case
ba2 = all_data[[8]]
# process ba2 data so that multiplicative effect is not applied recursively
BA2_RELU = 1.5
ba2v = 1.0 * (1 - ba2$ba2) + BA2_RELU * ba2$ba2
ba2f = ba2v
for (j in length(ba2v):2){
    ba2f[j] = ba2v[j] / ba2v[j-1]
}
ba2t = ba2$t

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

# instead of calculating quantiles across the means only, we need to calculate 
# quantiles taking into account the mean and sizes for deaths, hospital admissions,
# hospital bed and ICU bed occupancy, e.g.:
# spim_output = SPIM_output_full(w, 2022, 01, 31, "2020-01-01", "2022-09-30")
spim_output = read.csv('./output/paper/may22/basecase_220511_spim_output_20220511190051.csv', check.names = FALSE)
# save output for later
# datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
# write.csv(spim_output, paste0('./output/paper/jan22/basecase_220214_spim_output_', datetime, '.csv'))
setDT(spim_output)
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

# plot options

max_date = '2022-12-31'
date_breaks_opt = "2 months"
cols = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')

# hospital admissions
bc_ha = ggplot(spim_output[d >= "2020-03-01" & d <= max_date & AgeBand == "All" & ValueType == "Hospital\nadmissions" & Geography == "England"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
    geom_line(aes(x = d, y = Value, colour = "Modelled")) +
    geom_line(data = data[ValueType == "Hospital\nadmissions" & Geography == "England" & d <= max_date], aes(x = d, y = y, colour = "Data"), size = 0.2) +
    cowplot::theme_cowplot(font_size = 7) +
    theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5),
          text = element_text(size = 7, family = "sans")) +
    labs(x = NULL, y = "Hospital\nadmissions", fill = NULL, colour = NULL, tag = 'a', title = 'Basecase') +
    scale_colour_manual(values = c(Data = "black", Modelled = cols[2]), aesthetics = c("fill", "colour")) +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))
bc_ha

# infection prevalence
bc_ip = ggplot(spim_output[d >= "2020-03-01" & d <= max_date & AgeBand == "All" & ValueType == "PCR\nprevalence (%)" & Geography == "England"]) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = "Modelled"), alpha = 0.5) +
    geom_ribbon(aes(x = d, ymin = `Quantile 0.25`, ymax = `Quantile 0.75`, fill = "Modelled"), alpha = 0.75) +
    geom_line(aes(x = d, y = Value, colour = "Modelled")) +
    geom_ribbon(data = data[ValueType == "PCR\nprevalence (%)" & Geography == "England" & d <= max_date], aes(x = d, ymin = ymin, ymax = ymax, fill = "Data"), alpha = 0.5) +
    cowplot::theme_cowplot(font_size = 7) +
    theme(legend.position = 'none', panel.background = element_rect(fill = "#f4f4f4"), panel.grid = element_line(colour = "#ffffff", size = 0.5),
          text = element_text(size = 7, family = "sans")) +
    labs(x = NULL, y = "PCR\nprevalence (%)", fill = NULL, colour = NULL, tag = 'b', title = 'Basecase') +
    scale_colour_manual(values = c(Data = "black", Modelled = cols[5]), aesthetics = c("fill", "colour")) +
    scale_x_date(date_breaks = date_breaks_opt, date_labels = "%m/%y", expand = c(0.01, 0.01))
bc_ip

top_panel = cowplot::plot_grid(bc_ha, bc_ip, nrow = 2, align = 'hv')
top_panel
################################################################################

################################################################################

# panel A - behaviour plots

panelA_sheet = fread(
    "codenm,         mobility,                    wane_yn, vac_eff, wane_scen, bdur, bage, pby,  pbo,  pbst_all, seas, cert, mask,  wfh_date,   title,                                                                           vax_file,                         result_type, basecaseYN, basecase2YN
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Base case: 6-month relaxation; central waning; actual boosters; 20% seasonality, vax-covidm20220505205235.rds,    Behaviour,   Y,          N
mob_3wk_220511,       paper_aw_3wk_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         3-week relaxation,                                                               vax-covidm20220505205235.rds,    Behaviour,   N,          N
mob_3mo_220511,       paper_aw_3mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         3-month relaxation,                                                              vax-covidm20220505205235.rds,    Behaviour,   N,          N
mob_flt_220511,       paper_aw_flt_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         No relaxation,                                                                   vax-covidm20220505205235.rds,    Behaviour,   N,          N
")

panelA_scenarios <- function(variables_to_plot, scenario_sheet, list){
    all_d = data.table(variable = NULL, median = NULL, q05 = NULL, q25 = NULL, q75 = NULL, q95 = NULL, scenario = NULL, result_type = NULL, basecaseYN = NULL, basecase2YN = NULL, title = NULL)
    for (i in 1:length(list)){
        this_d = list[[i]]
        this_d = this_d[this_d$variable %in% variables_to_plot,]
        this_d$scenario = rep(scenario_sheet[i, codenm], dim(this_d)[1])
        this_d$result_type = rep(scenario_sheet[i, result_type], dim(this_d)[1])
        this_d$basecaseYN = rep(scenario_sheet[i, basecaseYN], dim(this_d)[1])
        this_d$basecase2YN = rep(scenario_sheet[i, basecase2YN], dim(this_d)[1])
        this_d$title = rep(scenario_sheet[i, title], dim(this_d)[1])
        print(scenario_sheet[i, title])
        all_d = rbind(all_d, this_d)
    }
    return(all_d)
}

panelA1_list = list()
for (ROW in 1:nrow(panelA_sheet)) {
    print(ROW)
    cat(panelA_sheet[ROW, codenm], "\n");
    panelA1_list[[ROW]] = fread(paste0('./output/paper/may22/cumul_', panelA_sheet[ROW, codenm], '.csv'))
}

d1 = panelA_scenarios(c('deaths','pcr_positive_i', 'hosp_adm'), panelA_sheet, panelA1_list)
d1[variable == "pcr_positive_i", name := "Cumulative infections"]
d1[variable == "hosp_adm", name := "Cumulative admissions"]
d1[variable == "deaths", name := "Cumulative deaths"]

panelA2_list = list()
for (ROW in 1:nrow(panelA_sheet)) {
    print(ROW)
    cat(panelA_sheet[ROW, codenm], "\n");
    panelA2_list[[ROW]] = fread(paste0('./output/paper/may22/spim_', panelA_sheet[ROW, codenm], '.csv'))
}
d2 = panelA_scenarios(c('deaths','pcr_positive_i', 'hosp_adm'), panelA_sheet, panelA2_list)
d2[variable == "pcr_positive_i", name := "Infections"]
d2[variable == "hosp_adm", name := "Admissions"]
d2[variable == "deaths", name := "Deaths"]

d1$type = rep('Cumulative deaths', dim(d1)[1])
d1$unit = rep(1, dim(d1)[1])
d1[variable == "pcr_positive_i", unit := 1000000]
d1[variable == "deaths", unit := 1000]
d1[variable == "hosp_adm", unit := 1000]
d1copy = d1
d1copy2 = d1
d1copy3 = d1
# value_Oct1_2021 = d1$q50[d1$date == '2021-10-01']
d2$type = rep('Deaths', dim(d2)[1])

cum_deaths_Oct1_2021 = mean(d1$q50[d1$date == '2021-10-01' & d1$variable == 'deaths'])
d1copy$q25[d1copy$variable == 'deaths' & d1copy$date > '2021-10-01'] = d1copy$q25[d1copy$variable == 'deaths' & d1copy$date > '2021-10-01'] - cum_deaths_Oct1_2021
d1copy$q50[d1copy$variable == 'deaths' & d1copy$date > '2021-10-01'] = d1copy$q50[d1copy$variable == 'deaths' & d1copy$date > '2021-10-01'] - cum_deaths_Oct1_2021
d1copy$q75[d1copy$variable == 'deaths' & d1copy$date > '2021-10-01'] = d1copy$q75[d1copy$variable == 'deaths' & d1copy$date > '2021-10-01'] - cum_deaths_Oct1_2021
d1copy$type = rep('Cumulative deaths\n(thousands)', dim(d1)[1])

cum_deaths_Jan1_2022 = mean(d1$q50[d1$date == '2022-01-01' & d1$variable == 'deaths'])
d1copy2$q25[d1copy2$variable == 'deaths' & d1copy2$date > '2022-01-01'] = d1copy2$q25[d1copy2$variable == 'deaths' & d1copy2$date > '2022-01-01'] - cum_deaths_Jan1_2022
d1copy2$q50[d1copy2$variable == 'deaths' & d1copy2$date > '2022-01-01'] = d1copy2$q50[d1copy2$variable == 'deaths' & d1copy2$date > '2022-01-01'] - cum_deaths_Jan1_2022
d1copy2$q75[d1copy2$variable == 'deaths' & d1copy2$date > '2022-01-01'] = d1copy2$q75[d1copy2$variable == 'deaths' & d1copy2$date > '2022-01-01'] - cum_deaths_Jan1_2022
d1copy2$type = rep('Cumulative deaths\n(thousands)', dim(d1)[1])

cum_deaths_Mar1_2022 = mean(d1$q50[d1$date == '2022-03-01' & d1$variable == 'deaths'])
d1copy3$q25[d1copy2$variable == 'deaths' & d1copy3$date > '2022-03-01'] = d1copy3$q25[d1copy3$variable == 'deaths' & d1copy3$date > '2022-03-01'] - cum_deaths_Mar1_2022
d1copy3$q50[d1copy2$variable == 'deaths' & d1copy3$date > '2022-03-01'] = d1copy3$q50[d1copy3$variable == 'deaths' & d1copy3$date > '2022-03-01'] - cum_deaths_Mar1_2022
d1copy3$q75[d1copy2$variable == 'deaths' & d1copy3$date > '2022-03-01'] = d1copy3$q75[d1copy3$variable == 'deaths' & d1copy3$date > '2022-03-01'] - cum_deaths_Mar1_2022
d1copy3$type = rep('Cumulative deaths\n(thousands)', dim(d1)[1])

d2$unit = rep(1, dim(d2)[1])
all = rbind(d1,d2)
all2 = rbind(d1copy, d2)
all3 = rbind(d1copy2, d2)
all4 = rbind(d1copy3, d2)
all[title == 'Base case: 6-month relaxation; central waning; actual boosters; 20% seasonality', title := '6-month relaxation*']
all2[title == 'Base case: 6-month relaxation; central waning; actual boosters; 20% seasonality', title := '6-month relaxation*']
all3[title == 'Base case: 6-month relaxation; central waning; actual boosters; 20% seasonality', title := '6-month relaxation*']
all4[title == 'Base case: 6-month relaxation; central waning; actual boosters; 20% seasonality', title := '6-month relaxation*']

colours3 = c('#1b9e77','#d95f02','#7570b3')
colours4 = c('#1b9e77','#d95f02','#7570b3','#e7298a')
colours3b = c('#66c2a5', '#fc8d62', '#8da0cb')

panelAalt2 = ggplot(all3[date > '2022-03-01' & variable == 'deaths'], aes(x=date)) +
    geom_line(aes(y = q50/unit, group = title, colour = title)) +
    geom_ribbon(aes(ymin = q25/unit, ymax = q75/unit, fill = title), alpha = 0.4) +
    labs(title='Behaviour', x = NULL, y = NULL, color = 'Scenario', fill = 'Scenario', tag = 'd') + 
    facet_wrap(~ type, scales = 'free_y', strip.position = 'left', ncol = 1) +
    #theme(strip.placement = 'outside', strip.background = element_blank(), axis.line.x = element_line(), text = element_text(size = 11), legend.position = c(0.2, 0.85)) +
    theme(strip.placement = 'outside', strip.background = element_blank(), axis.line.x = element_line(), legend.position = c(0.6, 0.5),
          text = element_text(size = 7, family = "sans")) +
    scale_color_manual(values = colours4) +
    scale_fill_manual(values = colours4) +
    #scale_x_date(breaks = c('2021-10-01', '2022-01-01', '2022-04-01', '2022-07-01'), labels = c('10/21', '01/22', '04/22', '07/22'))
    #scale_x_date(breaks = four_breaks, labels = four_labels) 
    scale_x_date(date_breaks = '1 month', date_labels = '%m/%y')
panelAalt2

################################################################################

# panel B - uncertainty plots

u = qread('./output/paper/may22/burdens_waneVHI_220511.qs')
u$date = u$t + as.Date('2020-01-01')
u = u[date > '2022-01-01']
u$run = factor(u$run, levels = (1:100))

# plot showing epidemic traces
colours100 = diverge_hcl(100, h = c(0, 300), c = 100, l = 75)
panelBa = ggplot(u[date > '2022-03-01' & variable == 'deaths'], aes(x = date)) +
    geom_line(aes(x = date, y = value, group = run, colour = run)) +
    labs(x = NULL, y = 'Deaths', tag = 'c', title = 'Very high waning') + 
    scale_color_manual(values = colours100) +
    theme(legend.position = 'none', axis.text.x = NULL,
          text = element_text(size = 7, family = "sans")) +
    #theme(text = element_text(size = 11), legend.position = 'none', axis.text.x = NULL) +
    scale_x_date(date_breaks = '1 month', date_labels = '%m/%y', limits = c(as.Date('2022-03-01'), as.Date('2022-12-31')))
panelBa

# we want to superimpose the median trajectory on top of panelBa
s = fread('./output/paper/may22/spim_waneVHI_220511.csv')
panelB = panelBa + geom_line(data = s[s$variable == 'deaths' & date > '2022-03-01'], aes(x = date, y = q50), size = 0.8) +
    geom_ribbon(data = s[s$variable == 'deaths' & date > '2021-10-01'], aes(x = date, ymin = q25, ymax = q75), alpha = 0.6) +
    geom_ribbon(data = s[s$variable == 'deaths' & date > '2021-10-01'], aes(x = date, ymin = q05, ymax = q95), alpha = 0.4) 
panelB

# list of objects so far
# 1. top_panel
# 1. panelAalt2
# 2. panelB (very high waning showing different epidemic traces)

row1 = cowplot::plot_grid(panelB, panelAalt2, nrow = 1, rel_widths = c(1.25,1))
row1
top_half = cowplot::plot_grid(top_panel, row1, nrow = 2, rel_heights = c(1, 0.75))
top_half

################################################################################

# plots showing booster and waning scenarios

panelC_sheet = fread(
"codenm,             mobility,                    wane_yn, vac_eff, wane_scen, bdur, bage, pby,  pbo,  pbst_all, seas, cert, mask,  wfh_date,   title,                     vax_file,                      result_type,      basecaseYN, basecase2YN
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,        Central waning*,           vax-covidm20220505205235.rds,  Waning,           Y,          N
waneHI_220511,        paper_aw_6mo_20220509105131, yeswane, central, hiwane,    180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,        High waning,               vax-covidm20220505205235.rds,  Waning,           N,          N
waneVHI_220511,       paper_aw_6mo_20220509105131, yeswane, central, vhiwane,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,        Very high waning,          vax-covidm20220505205235.rds,  Waning,           N,          N
boostNO_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,    0.0,  0.0, NULL,     0.1,  0,    0,     NA,        No booster uptake,         vax-covidm20220505205235.rds,  Boosters,         N,          N
boost50_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,    0.0, 0.95, NULL,     0.1,  0,    0,     NA,        Boosters for 50+,          vax-covidm20220505205235.rds,  Boosters,         N,          N
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,        Actual booster uptake*,    vax-covidm20220505205235.rds,  Boosters,         Y,          N
boostHI_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.90, 0.98, NULL,     0.1,  0,    0,     NA,        High booster uptake,       vax-covidm20220505205235.rds,  Boosters,         N,          N
shrtbst_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,    90,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,        Short booster duration,    vax-covidm20220505205235.rds,  Booster duration, N,          N
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,        Central booster duration*, vax-covidm20220505205235.rds,  Booster duration, Y,          N
longbst_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   270,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,        Long booster duration,     vax-covidm20220505205235.rds,  Booster duration, N,          N
seas_10_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.05,  0,    0,     NA,        10% seasonality,           vax-covidm20220505205235.rds,  Seasonality,      N,          N
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,        20% seasonality*,          vax-covidm20220505205235.rds,  Seasonality,      Y,          N
seas_30_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.15,  0,    0,     NA,        30% seasonality,           vax-covidm20220505205235.rds,  Seasonality,      N,          N
seas_40_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.20,  0,    0,     NA,        40% seasonality,           vax-covidm20220505205235.rds,  Seasonality,      N,          N
")

panelC_scenarios <- function(variables_to_plot, scenario_sheet, list){
    all_d = data.table(variable = NULL, median = NULL, q05 = NULL, q25 = NULL, q75 = NULL, q95 = NULL, scenario = NULL, result_type = NULL, basecaseYN = NULL, basecase2YN = NULL, title = NULL)
    for (i in 1:length(list)){
        this_d = list[[i]]
        this_d = this_d[this_d$variable %in% variables_to_plot,]
        this_d$scenario = rep(scenario_sheet[i, codenm], dim(this_d)[1])
        this_d$result_type = rep(scenario_sheet[i, result_type], dim(this_d)[1])
        this_d$basecaseYN = rep(scenario_sheet[i, basecaseYN], dim(this_d)[1])
        this_d$basecase2YN = rep(scenario_sheet[i, basecase2YN], dim(this_d)[1])
        this_d$title = rep(scenario_sheet[i, title], dim(this_d)[1])
        print(scenario_sheet[i, title])
        all_d = rbind(all_d, this_d)
    }
    return(all_d)
}

panelC_list = list()
for (ROW in 1:nrow(panelC_sheet)) {
    print(ROW)
    cat(panelC_sheet[ROW, codenm], "\n");
    panelC_list[[ROW]] = fread(paste0('./output/paper/may22/spim_', panelC_sheet[ROW, codenm], '.csv'))
}

c = panelC_scenarios(c('deaths','pcr_positive_i', 'hosp_adm'), panelC_sheet, panelC_list)
c[variable == "pcr_positive_i", name := "Infections"]
c[variable == "hosp_adm", name := "Admissions"]
c[variable == "deaths", name := "Deaths"]
c$name = factor(c$name, levels = c('Infections', 'Admissions', 'Deaths'))
c$date = c$t + as.Date('2020-01-01')
c$title = factor(c$title, levels = c("Central waning*",          
                                     "High waning",              
                                     "Very high waning",
                                     "No booster uptake",        
                                     "Boosters for 50+",         
                                     "Actual booster uptake*",   
                                     "High booster uptake",      
                                     "Short booster duration",   
                                     "Central booster duration*",
                                     "Long booster duration",    
                                     "10% seasonality",          
                                     "20% seasonality*",         
                                     "30% seasonality",          
                                     "40% seasonality"))

# let's calculate cumulative burdens
panelC_list_cumul = list()
for (ROW in 1:nrow(panelC_sheet)) {
    print(ROW)
    cat(panelC_sheet[ROW, codenm], "\n");
    panelC_list_cumul[[ROW]] = fread(paste0('./output/paper/may22/cumul_', panelC_sheet[ROW, codenm], '.csv'))
}

c2 = panelC_scenarios(c('deaths','pcr_positive_i', 'hosp_adm'), panelC_sheet, panelC_list_cumul)
c2[variable == "pcr_positive_i", name := "Cumulative infections"]
c2[variable == "hosp_adm", name := "Cumulative admissions"]
c2[variable == "deaths", name := "Cumulative deaths"]
c2$date = c2$t + as.Date('2020-01-01')
c2$title = factor(c2$title, levels = c("Central waning*",          
                                       "High waning",              
                                       "Very high waning",
                                       "No booster uptake",        
                                       "Boosters for 50+",         
                                       "Actual booster uptake*",   
                                       "High booster uptake",      
                                       "Short booster duration",   
                                       "Central booster duration*",
                                       "Long booster duration",    
                                       "10% seasonality",          
                                       "20% seasonality*",         
                                       "30% seasonality",          
                                       "40% seasonality"))


# for retrospective booster policy scenarios, we want cumulative burdens since 1st October 2021

# for other scenarios we want cumulative burdens since 1st January 2022

# so we need to calculate both :)

# firstly, cumulative burdens since 1st October 2021

# try to subtract cumulative infections, admissions and deaths from 1st October 2021 before plotting
c2copy = c2
cum_vals_Oct1_2021 = c2[c2$date == '2021-10-01']
variables = unique(c2$variable)
titles = unique(c2$title)

for (i in variables){
    for (j in titles){
        cum_vals = cum_vals_Oct1_2021[cum_vals_Oct1_2021$variable == i & cum_vals_Oct1_2021$title == j]
        c2copy$q05[c2copy$variable == i & c2copy$title == j & c2copy$date > '2021-10-01'] = c2copy$q05[c2copy$variable == i & c2copy$title == j & c2copy$date > '2021-10-01'] - cum_vals$q05
        c2copy$q25[c2copy$variable == i & c2copy$title == j & c2copy$date > '2021-10-01'] = c2copy$q25[c2copy$variable == i & c2copy$title == j & c2copy$date > '2021-10-01'] - cum_vals$q25
        c2copy$q50[c2copy$variable == i & c2copy$title == j & c2copy$date > '2021-10-01'] = c2copy$q50[c2copy$variable == i & c2copy$title == j & c2copy$date > '2021-10-01'] - cum_vals$q50
        c2copy$q75[c2copy$variable == i & c2copy$title == j & c2copy$date > '2021-10-01'] = c2copy$q75[c2copy$variable == i & c2copy$title == j & c2copy$date > '2021-10-01'] - cum_vals$q75
        c2copy$q95[c2copy$variable == i & c2copy$title == j & c2copy$date > '2021-10-01'] = c2copy$q95[c2copy$variable == i & c2copy$title == j & c2copy$date > '2021-10-01'] - cum_vals$q95
        
    }
}

# split up into result types
be = c2copy[c2copy$result_type == 'Boosters']
w = c2copy[c2copy$result_type == 'Waning']
bd = c2copy[c2copy$result_type == 'Booster duration']
se = c2copy[c2copy$result_type == 'Seasonality']

be$newname = be$name
bd$newname = bd$name
se$newname = se$name
w$newname = w$name

be[name == 'Cumulative infections', newname := 'Cumulative infections\n(millions)']
be[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']
bd[name == 'Cumulative infections', newname := 'Cumulative infections\n(millions)']
bd[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']
se[name == 'Cumulative infections', newname := 'Cumulative infections\n(millions)']
se[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']
w[name =='Cumulative infections', newname := 'Cumulative infections\n(millions)']
w[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']

be$newname = factor(be$newname, levels = c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))
w$newname = factor(w$newname, levels = c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))
se$newname = factor(se$newname, levels = c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))
bd$newname = factor(bd$newname, levels = c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))

be$unit = rep(1, dim(be)[1])
be[name == 'Cumulative infections', unit := 1000000]
be[name == 'Cumulative deaths', unit := 1000]
w$unit = rep(1, dim(w)[1])
w[name == 'Cumulative infections', unit := 1000000]
w[name == 'Cumulative deaths', unit := 1000]

bd$unit = rep(1, dim(bd)[1])
bd[name == 'Cumulative infections', unit := 1000000]
bd[name == 'Cumulative deaths', unit := 1000]
se$unit = rep(1, dim(se)[1])
se[name == 'Cumulative infections', unit := 1000000]
se[name == 'Cumulative deaths', unit := 1000]


be$newtitle = be$title
be[title == "No booster uptake", newtitle := "No boosters"]
be[title == "Boosters for 50+", newtitle := "Boost 50+ only"]
be[title == "Actual booster uptake*", newtitle := "Actual boosters*"]
be[title == "High booster uptake", newtitle := "Higher uptake"]
be$newtitle = factor(be$newtitle, levels = c("No boosters",
                                             "Boost 50+ only",
                                             "Actual boosters*",
                                             "Higher uptake"))

plotboost = ggplot(be[be$date > '2021-10-01' & be$name != 'Cumulative admissions']) + 
    geom_ribbon(aes(date, ymin = q05/unit, ymax = q95/unit, fill = newtitle), alpha = 0.2) +
    geom_ribbon(aes(date, ymin = q25/unit, ymax = q75/unit, fill = newtitle), alpha = 0.4) +
    geom_line(aes(date, q50/unit, colour = newtitle)) +
    #geom_line(aes(date, va/unit), colour = col, linetype = "dashed") +
    geom_hline(aes(yintercept = 0), size = .5) +
    scale_x_date(date_breaks = "2 months", date_labels = "%m/%y", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0.0, 0.0)) +
    labs(x = NULL, y = NULL, tag = 'e', title = 'Booster vaccinations') +
    facet_wrap(~newname, scales = 'free', ncol = 1, strip.position = 'left') +
    #facet_grid(name ~ result_type, scales = 'free', switch = 'y') +
    #cowplot::theme_cowplot(font_size = 9) +
    #theme(strip.background = element_blank(), strip.placement = 'outside', legend.position = 'bottom', text = element_text(size = 11)) +
    theme(strip.background = element_blank(), strip.placement = 'outside', legend.position = 'bottom',
          text = element_text(size = 7, family = "sans")) +
    labs(fill = 'Scenario', colour = 'Scenario')
plotboost

# now do the same 4 plots, but show cumulative numbers from January 2022 onwards
c2copy2 = c2
cum_vals_Jan1_2022 = c2[c2$date == '2022-01-01']

for (i in variables){
    for (j in titles){
        cum_vals = cum_vals_Jan1_2022[cum_vals_Jan1_2022$variable == i & cum_vals_Jan1_2022$title == j]
        c2copy2$q05[c2copy2$variable == i & c2copy2$title == j & c2copy2$date > '2022-01-01'] = c2copy2$q05[c2copy2$variable == i & c2copy2$title == j & c2copy2$date > '2022-01-01'] - cum_vals$q05
        c2copy2$q25[c2copy2$variable == i & c2copy2$title == j & c2copy2$date > '2022-01-01'] = c2copy2$q25[c2copy2$variable == i & c2copy2$title == j & c2copy2$date > '2022-01-01'] - cum_vals$q25
        c2copy2$q50[c2copy2$variable == i & c2copy2$title == j & c2copy2$date > '2022-01-01'] = c2copy2$q50[c2copy2$variable == i & c2copy2$title == j & c2copy2$date > '2022-01-01'] - cum_vals$q50
        c2copy2$q75[c2copy2$variable == i & c2copy2$title == j & c2copy2$date > '2022-01-01'] = c2copy2$q75[c2copy2$variable == i & c2copy2$title == j & c2copy2$date > '2022-01-01'] - cum_vals$q75
        c2copy2$q95[c2copy2$variable == i & c2copy2$title == j & c2copy2$date > '2022-01-01'] = c2copy2$q95[c2copy2$variable == i & c2copy2$title == j & c2copy2$date > '2022-01-01'] - cum_vals$q95
        
    }
}

be2 = c2copy2[c2copy2$result_type == 'Boosters']
w2 =  c2copy2[c2copy2$result_type == 'Waning']
bd2 = c2copy2[c2copy2$result_type == 'Booster duration']
se2 = c2copy2[c2copy2$result_type == 'Seasonality']

be2$newname = be2$name
bd2$newname = bd2$name
se2$newname = se2$name
w2$newname = w2$name

be2[name == 'Cumulative infections', newname := 'Cumulative infections\n(millions)']
be2[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']
bd2[name == 'Cumulative infections', newname := 'Cumulative infections\n(millions)']
bd2[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']
se2[name == 'Cumulative infections', newname := 'Cumulative infections\n(millions)']
se2[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']
w2[name =='Cumulative infections', newname := 'Cumulative infections\n(millions)']
w2[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']

be2$newname = factor(be2$newname, levels = c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))
w2$newname = factor(w2$newname, levels =   c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))
se2$newname = factor(se2$newname, levels = c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))
bd2$newname = factor(bd2$newname, levels = c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))

be2$unit = rep(1, dim(be2)[1])
be2[name == 'Cumulative infections', unit := 1000000]
be2[name == 'Cumulative deaths', unit := 1000]
w2$unit = rep(1, dim(w2)[1])
w2[name == 'Cumulative infections', unit := 1000000]
w2[name == 'Cumulative deaths', unit := 1000]

bd2$unit = rep(1, dim(bd2)[1])
bd2[name == 'Cumulative infections', unit := 1000000]
bd2[name == 'Cumulative deaths', unit := 1000]
se2$unit = rep(1, dim(se2)[1])
se2[name == 'Cumulative infections', unit := 1000000]
se2[name == 'Cumulative deaths', unit := 1000]


be2$newtitle = be2$title
be2[title == "No booster uptake", newtitle := "No boosters"]
be2[title == "Boosters for 50+", newtitle := "Boost 50+ only"]
be2[title == "Actual booster uptake*", newtitle := "Actual boosters*"]
be2[title == "High booster uptake", newtitle := "Higher uptake"]
be2$newtitle = factor(be2$newtitle, levels = c("No boosters", 
                                               "Boost 50+ only",
                                               "Actual boosters*",
                                               "Higher uptake"))

plotwane2 = ggplot(w2[w2$date > '2022-01-01' & w2$name != 'Cumulative admissions']) + 
    geom_ribbon(aes(date, ymin = q05/unit, ymax = q95/unit, fill = title), alpha = 0.2) +
    geom_ribbon(aes(date, ymin = q25/unit, ymax = q75/unit, fill = title), alpha = 0.4) +
    geom_line(aes(date, q50/unit, colour = title)) +
    #geom_line(aes(date, va/unit), colour = col, linetype = "dashed") +
    geom_hline(aes(yintercept = 0), size = .5) +
    scale_x_date(date_breaks = "1 month", date_labels = "%m/%y", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0.0, 0.0)) +
    labs(x = NULL, y = NULL, tag = 'f', title = 'Waning immunity') +
    facet_wrap(~newname, scales = 'free', ncol = 1, strip.position = 'left') +
    #facet_grid(name ~ result_type, scales = 'free', switch = 'y') +
    #cowplot::theme_cowplot(font_size = 9) +
    #theme(strip.background = element_blank(), strip.placement = 'outside', legend.position = 'bottom', text = element_text(size = 11)) +
    theme(strip.background = element_blank(), strip.placement = 'outside', legend.position = 'bottom',
          text = element_text(size = 7, family = "sans")) +
    labs(fill = 'Scenario', colour = 'Scenario')
plotwane2

################################################################################

# now do the same 4 plots, but show cumulative numbers from March 2022 onwards
c2copy3 = c2
cum_vals_Mar1_2022 = c2[c2$date == '2022-03-01']

for (i in variables){
    for (j in titles){
        cum_vals = cum_vals_Mar1_2022[cum_vals_Mar1_2022$variable == i & cum_vals_Mar1_2022$title == j]
        c2copy3$q05[c2copy3$variable == i & c2copy3$title == j & c2copy3$date > '2022-03-01'] = c2copy3$q05[c2copy3$variable == i & c2copy3$title == j & c2copy3$date > '2022-03-01'] - cum_vals$q05
        c2copy3$q25[c2copy3$variable == i & c2copy3$title == j & c2copy3$date > '2022-03-01'] = c2copy3$q25[c2copy3$variable == i & c2copy3$title == j & c2copy3$date > '2022-03-01'] - cum_vals$q25
        c2copy3$q50[c2copy3$variable == i & c2copy3$title == j & c2copy3$date > '2022-03-01'] = c2copy3$q50[c2copy3$variable == i & c2copy3$title == j & c2copy3$date > '2022-03-01'] - cum_vals$q50
        c2copy3$q75[c2copy3$variable == i & c2copy3$title == j & c2copy3$date > '2022-03-01'] = c2copy3$q75[c2copy3$variable == i & c2copy3$title == j & c2copy3$date > '2022-03-01'] - cum_vals$q75
        c2copy3$q95[c2copy3$variable == i & c2copy3$title == j & c2copy3$date > '2022-03-01'] = c2copy3$q95[c2copy3$variable == i & c2copy3$title == j & c2copy3$date > '2022-03-01'] - cum_vals$q95
        
    }
}

be3 = c2copy3[c2copy3$result_type == 'Boosters']
w3 =  c2copy3[c2copy3$result_type == 'Waning']
bd3 = c2copy3[c2copy3$result_type == 'Booster duration']
se3 = c2copy3[c2copy3$result_type == 'Seasonality']

be3$newname = be3$name
bd3$newname = bd3$name
se3$newname = se3$name
w3$newname = w3$name

be3[name == 'Cumulative infections', newname := 'Cumulative infections\n(millions)']
be3[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']
bd3[name == 'Cumulative infections', newname := 'Cumulative infections\n(millions)']
bd3[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']
se3[name == 'Cumulative infections', newname := 'Cumulative infections\n(millions)']
se3[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']
w3[name =='Cumulative infections', newname := 'Cumulative infections\n(millions)']
w3[name == 'Cumulative deaths', newname := 'Cumulative deaths\n(thousands)']

be3$newname = factor(be3$newname, levels = c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))
w3$newname = factor(w3$newname, levels =   c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))
se3$newname = factor(se3$newname, levels = c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))
bd3$newname = factor(bd3$newname, levels = c('Cumulative infections\n(millions)', 'Cumulative admissions', 'Cumulative deaths\n(thousands)'))

be3$unit = rep(1, dim(be3)[1])
be3[name == 'Cumulative infections', unit := 1000000]
be3[name == 'Cumulative deaths', unit := 1000]
w3$unit = rep(1, dim(w3)[1])
w3[name == 'Cumulative infections', unit := 1000000]
w3[name == 'Cumulative deaths', unit := 1000]

bd3$unit = rep(1, dim(bd3)[1])
bd3[name == 'Cumulative infections', unit := 1000000]
bd3[name == 'Cumulative deaths', unit := 1000]
se3$unit = rep(1, dim(se3)[1])
se3[name == 'Cumulative infections', unit := 1000000]
se3[name == 'Cumulative deaths', unit := 1000]


be3$newtitle = be3$title
be3[title == "No booster uptake", newtitle := "No boosters"]
be3[title == "Boosters for 50+", newtitle := "Boost 50+ only"]
be3[title == "Actual booster uptake*", newtitle := "Actual boosters*"]
be3[title == "High booster uptake", newtitle := "Higher uptake"]
be3$newtitle = factor(be3$newtitle, levels = c("No boosters", 
                                               "Boost 50+ only",
                                               "Actual boosters*",
                                               "Higher uptake"))

plotwane3 = ggplot(w3[w3$date > '2022-03-01' & w3$name != 'Cumulative admissions']) + 
    geom_ribbon(aes(date, ymin = q05/unit, ymax = q95/unit, fill = title), alpha = 0.2) +
    geom_ribbon(aes(date, ymin = q25/unit, ymax = q75/unit, fill = title), alpha = 0.4) +
    geom_line(aes(date, q50/unit, colour = title)) +
    #geom_line(aes(date, va/unit), colour = col, linetype = "dashed") +
    geom_hline(aes(yintercept = 0), size = .5) +
    scale_x_date(date_breaks = "1 month", date_labels = "%m/%y", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0.0, 0.0)) +
    labs(x = NULL, y = NULL, tag = 'f', title = 'Waning immunity') +
    facet_wrap(~newname, scales = 'free', ncol = 1, strip.position = 'left') +
    #facet_grid(name ~ result_type, scales = 'free', switch = 'y') +
    #cowplot::theme_cowplot(font_size = 9) +
    #theme(strip.background = element_blank(), strip.placement = 'outside', legend.position = 'bottom', text = element_text(size = 11)) +
    theme(strip.background = element_blank(), strip.placement = 'outside', legend.position = 'bottom',
          text = element_text(size = 7, family = "sans")) +
    labs(fill = 'Scenario', colour = 'Scenario')
plotwane3

row2_mix = cowplot::plot_grid(plotboost, plotwane3, nrow = 1)
row2_mix
first_two_rows_3 = cowplot::plot_grid(row1, row2_mix, nrow = 2, rel_heights = c(1,2))
first_two_rows_3

first_three_rows_3 = cowplot::plot_grid(top_half, row2_mix, nrow = 2, rel_heights = c(1.75, 1.25))
first_three_rows_3

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0('./output/paperfigs/july22/fig4_', datetime, '_2column.pdf'), first_three_rows_3, width = 180, height = 190, units = "mm", useDingbats = FALSE)
ggsave(paste0('./output/paperfigs/july22/fig4_', datetime, '_2column.png'), first_three_rows_3, width = 180, height = 190, units = "mm")

