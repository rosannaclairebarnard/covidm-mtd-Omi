# this script runs N covidm (3-strain) fits with vaccines on AWS, from the command line
#
# e.g., call
#
# Rscript fit.R mechanism set REP_START REP_END WANE_YN VAC_EFF V3_SEVERITY DFR SEAS_YN SEAS_AMP START_FILE BFOLD OMI_PROTECTION OMI_SEV OMI_CRIT BA2_RELU GAM_DISP
#
# currently this script is only setup to fit mechanism = relu for the second and third strains
#
# where
#
# mechanism : (relu novoc)
# set : (all else mnsw test)
# REP_START : first new fit number
# REP_END : last new fit number
# WANE_YN: determines whether natural and vaccine induced immunity are allowed to wane ("yeswane", "nowane")
# VAC_EFF: "hi" or "lo"
# V3_SEVERITY: relative severity of third variant (compared to second)
# DFR: Fold reduction in neutralisation for Delta
# SEAS_YN: "seasyes" or "seasno" or "seaslate" (currently implemented on day 1, 20% amplitude peak to trough; seaslate is 1st April 2021)
# SEAS_AMP: amplitude of seasonality
# START_FILE: either "previous", "scratch", or the base name of a fit file to start from in ./fits/

# currently, the script is set up to run on English regions only (TODO: add devolved administrations)

library(data.table)
library(ggplot2)
library(lubridate)
library(here)
library(cowplot)
library(readxl)
library(sn)
library(qs)
library(stringr)
library(mgcv)
library(binom)

theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

######################## CHANGE HERE FOR EACH RUN ########################

creation_year <- 2022
creation_month <- 05
creation_day <- 06
forecast_start <- paste0(creation_year,"-",creation_month,"-",creation_day)
forecast_end <- "2022-06-29" # >= 6 weeks following forecast_start for SPI-M
analysis_end <- "2022-06-29" # for other analyses, make sure larger than forecast_end
if (forecast_end > analysis_end){
    stop("Error: end of SPI-M forecast exceeds analysis_end date")
}

BURN_IN = 2500
ITER = 250
BURN_IN_FINAL = 2500
ITER_FINAL = 250

data_file = "processed-data-20220506122858.qs"
mobility_file = "schedule3-MTPs-20220506121302.rds"
date_fitting = "2022-05-06"
vax_file = "vax-covidm20220505205235.rds"
sgtf_stop = "2021-02-15"

# SERO_CUT_OFF determines the date at which seroprevalence data is no longer used to fit to (we use the start date of estimates for cut off)
SERO_CUT_OFF = '2020-12-01'

# Our assumptions on booster uptake rely on NHS England statistics on vaccines 
# delivered - we calculate booster dose uptake relative to second dose uptake by
# dividing the total number of booster vaccination doses delivered by the total 
# number of second vaccine doses delivered by age group. This calculation was 
# last updated using data from 10th March 2022, available at:
# https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/
# The numbers below are to 3 decimal places, workings are in the .xlsx file in
# newcovid3/data/COVID-19-monthly-announced-vaccinations-10-March-2022_EDIT.xlsx

# OLD calculations from 13th January 2022 data (kept here for comparison only)
# booster_uptake_measured = fread(
#     "under_18, 18_24, 25_29, 30_34, 35_39, 40_44, 45_49, 50_54, 55_59, 60_64, 65_69, 70_74, 75_79, 80plus
#          0.06, 0.396, 0.444, 0.510, 0.573, 0.655, 0.722, 0.811, 0.851, 0.890, 0.929, 0.953, 0.959, 0.945
# ")
# OLD calculations from 10th February 2022 data
# booster_uptake_measured = fread(
#     "under_18, 18_24, 25_29, 30_34, 35_39, 40_44, 45_49, 50_54, 55_59, 60_64, 65_69, 70_74, 75_79, 80plus
#         0.128, 0.516, 0.561, 0.616, 0.677, 0.742, 0.797, 0.854, 0.885, 0.913, 0.944, 0.962, 0.967, 0.956
# ")
# OLD calculations from 10th March 2022 data
# booster_uptake_measured = fread(
#     "under_18, 18_24, 25_29, 30_34, 35_39, 40_44, 45_49, 50_54, 55_59, 60_64, 65_69, 70_74, 75_79, 80plus
#         0.179, 0.544, 0.588, 0.638, 0.695, 0.756, 0.810, 0.862, 0.891, 0.918, 0.947, 0.965, 0.970, 0.960
# ")
# NEW calculations from 14th April 2022 data
# booster_uptake_measured = fread(
#        "16_17, 18_24, 25_29, 30_34, 35_39, 40_44, 45_49, 50_54, 55_59, 60_64, 65_69, 70_74, 75_79, 80plus
#         0.183, 0.544, 0.597, 0.643, 0.698, 0.757, 0.809, 0.861, 0.892, 0.917, 0.946, 0.964, 0.971, 0.963  
# ")

# We base our assumptions on booster uptake on the measured booster dose uptake
# from NHS data above to 14th April 2022, rounding to 3 decimal places.

# For the 75+ covidm age group, we take the average of the measured uptake for 
# age groups 75-79 and 80+: mean(c(0.971,0.963))
# [1] 0.967

# For the 15-19 covidm age group, we assume an uptake of 40%, based on that age
# group comprising 18-19 year olds (who have a measured uptake of 54.4%) and 
# 16-17 year olds (who have a measured uptake of 18.3%). In England, people
# aged 16 and above are currently being offered booster vaccines at least 3 
# months following their second dose.

# values from April 2022
# booster_uptake_assumed_cm_age_groups = fread(
#     "0_4, 5_9, 10_14, 15_19, 20_24, 25_29, 30_34, 35_39, 40_44, 45_49, 50_54, 55_59, 60_64, 65_69, 70_74, 75plus
#        0,   0,     0,   0.4, 0.544, 0.597, 0.643, 0.698, 0.757, 0.809, 0.861, 0.892, 0.917, 0.946, 0.964,  0.967
# ")
PBOOST=c(0,   0,     0,   0.4, 0.544, 0.597, 0.643, 0.698, 0.757, 0.809, 0.861, 0.892, 0.917, 0.946, 0.964,  0.967)
# BOOST_UP = 16
# BOOST_AGE_SPLIT = 50
# P_BOOST_YOUNG = 0.85
# P_BOOST_OLD = 0.95

# set these up properly later
OMI_SETUP_T = 630 # switch variants on 22nd September 2021 (on this day variant 1 is reparameterised from wildtype to become Omicron and any R1 individuals move to R3)
# OMI_X_PROTECTION = 0.551
# OMI_VAX_FACTOR = 0.551
OMI_VAX_ASSUMPTION = 'khoury'
OMI_SEEDS_PER_DAY = 10

############################## END CHANGES ##############################

which_pops = c(1, 3, 4, 5, 6, 9, 10)
set_id = ""
# c_args = c('relu', 'all', '1', '1', 'yeswane', 'central', '2.0', '1.0', 'seaslate', '0.1', 'scratch', '2.5', '0.551', '0.5', '0.5', '1.3', '0.1')
# Command line
c_args = commandArgs(trailingOnly = TRUE);
if (length(c_args) != 17) {
    stop("17 arguments required.")
}
FIT_TYPE = c_args[[1]];
POP_SET = c_args[[2]];
REP_START = as.numeric(c_args[[3]]);
REP_END = as.numeric(c_args[[4]]);
WANE_YN = c_args[[5]]
VAC_EFF = c_args[[6]]
V3_SEVERITY = c_args[[7]]
DFR = as.numeric(c_args[[8]])
SEAS_YN = c_args[[9]]
if (!SEAS_YN %in% c("seasyes", "seasno", "seaslate")){
    stop("Seasonality option should be seasyes, seasno, seaslate")
}
SEAS_AMP = as.numeric(c_args[[10]])
START_FILE = c_args[[11]]
BFOLD = as.numeric(c_args[[12]])
OMI_PROTECTION = as.numeric(c_args[[13]])
OMI_X_PROTECTION = OMI_PROTECTION
OMI_VAX_FACTOR = OMI_PROTECTION
OMI_SEV = as.numeric(c_args[[14]])
OMI_CRIT = as.numeric(c_args[[15]])
BA2_RELU = as.numeric(c_args[[16]])
GAMDISP = as.numeric(c_args[[17]])

opt_conc = TRUE;
opt_seas = FALSE;

opt_v2 = TRUE;
opt_relu = FALSE;
opt_latdur = FALSE;
opt_serial = FALSE;
opt_infdur = FALSE;
opt_immesc = FALSE;
opt_ch_u = FALSE;
extra_priors = list();

if (FIT_TYPE == "relu") {
    extra_priors = list(v2_relu = "L 0.0 0.4 T 0.25 4",
                        v3_relu = "L 0.0 0.4 T 0.25 4",
                        v4_relu = "L 0.4 0.1 T 0.25 4");
    opt_relu = TRUE;
    opt_relu3 = TRUE;
    opt_relu4 = TRUE;
    opt_v3 = TRUE;
    opt_v4 = TRUE;
} else if (FIT_TYPE == "novoc") {
    opt_v2 = FALSE;
    opt_v3 = FALSE;
    opt_v4 = FALSE;
} else {
    stop("Need to specify fit type at command line.");
}

if (POP_SET == "else") {
    which_pops = c(1, 3, 9)
    pop_letter = "_ELSE"
} else if (POP_SET == "mnsw") {
    which_pops = c(4, 5, 6, 10)
    pop_letter = "_MNSW"
} else if (POP_SET == "test"){
    which_pops = c(3)
    pop_letter = "_test"
} else if (POP_SET == "all") {
    which_pops = c(1, 3, 4, 5, 6, 9, 10)
    pop_letter = ""
} else {
    stop("POP_SET must be else or all.");
}

set_id = paste0(FIT_TYPE, pop_letter, "_", WANE_YN, "_sev", V3_SEVERITY, "_");

uk_covid_data_path = "./fitting_data/";
datapath = function(x) paste0(uk_covid_data_path, x)

#
# SETUP
#

# set up covidm
cm_path = "./covidm_for_fitting/";
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 3;
source(paste0(cm_path, "/R/covidm.R"))
popUK = readRDS(datapath("popNHS.rds"));
matricesUK = readRDS(datapath("matricesNHS.rds"));

cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)

source("./spim_output.R");
source("./check_fit.R")
source("./params.R")
source("./booster_schedule.R")

#
# DATA
#

# PBOOST = c(rep(P_BOOST_YOUNG, BOOST_AGE_SPLIT %/% 5), 
#            rep(P_BOOST_OLD, 16 - BOOST_AGE_SPLIT %/% 5))
# PBOOST = PBOOST * cm_age_coefficients(BOOST_UP, 80, seq(0, 80, by = 5))


vacc <- readRDS(datapath(vax_file))

nhs_regions = popUK[, unique(name)]
pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100

all_data = qread(datapath(data_file))
ld = all_data[[1]]
sitreps = all_data[[2]]
virus = all_data[[3]][!Data.source %like% "7a|7b|6a|6b|9a|9b"]
sero = all_data[[4]]
# we only want to fit to seroprevalence data prior to December 2020
sero_to_fit = sero[sero$Start.date < SERO_CUT_OFF]
sgtf = all_data[[5]]
delta = all_data[[6]]
omi = all_data[[7]]
ba2 = all_data[[8]]
# remove any omicron sgtf data occurring before 1st October 2021
omi = omi[omi$date >= '2021-10-01']

# process ba2 data so that multiplicative effect is not applied recursively
ba2v = 1.0 * (1 - ba2$ba2) + BA2_RELU * ba2$ba2
ba2f = ba2v
for (j in length(ba2v):2){
    ba2f[j] = ba2v[j] / ba2v[j-1]
}
ba2t = ba2$t
# ba2 = fread(
# "t,ba2
# 600,0.01
# 601,0.02
# 602,0.03")

#
# FITTING
#

# NUMBER OF REGIONS TO FIT
N_REG = 12;

# Build parameters for NHS regions
params = cm_parameters_SEI3R(nhs_regions[1:N_REG], deterministic = T, 
                             date_start = "2020-01-01", 
                             date_end = date_fitting,
                             dE  = cm_delay_gamma(2.5, 2.5, 
                                                  t_max = 15, t_step = 0.25)$p,
                             dIp = cm_delay_gamma(2.5, 4.0, 
                                                  t_max = 15, t_step = 0.25)$p,
                             dIs = cm_delay_gamma(2.5, 4.0, 
                                                  t_max = 15, t_step = 0.25)$p,
                             dIa = cm_delay_gamma(5.0, 4.0, 
                                                  t_max = 15, t_step = 0.25)$p)
params = cm_split_matrices_ex_in(params, 15)

# Load age-varying symptomatic rate
covid_scenario = qread(datapath("2-linelist_both_fit_fIa0.5-rbzvih.qs"));
covu = unname(rep(colMeans(covid_scenario[,  5:12]), each = 2));
covy = unname(rep(colMeans(covid_scenario[, 13:20]), each = 2));

# Health burden processes
source("./processes.R")
params$processes = burden_processes


for (i in seq_along(params$pop)) {
    params$pop[[i]]$u = covu / mean(covu);
    params$pop[[i]]$u2 = covu / mean(covu);
    params$pop[[i]]$u3 = covu / mean(covu);
    params$pop[[i]]$y = covy;
    params$pop[[i]]$y2 = covy;
    params$pop[[i]]$y3 = covy;
    
    params = reinfection_defaults(params, i)  #######
    
    # set up waning parameters
    params = waning_scenario(WANE_YN, params, i)
    
    # Assign vaccine efficacy parameters
    params = VE_scenario(VAC_EFF, PBOOST, 
                         params, i, delta_fold_reduction = DFR)
    
    params$pop[[i]]$ifr1 = P.death
    params$pop[[i]]$ihr1 = P.hosp
    params$pop[[i]]$iir1 = P.critical
    params$pop[[i]]$ifr2 = P.death
    params$pop[[i]]$ihr2 = P.hosp
    params$pop[[i]]$iir2 = P.critical
    params$pop[[i]]$ifr3 = P.death
    params$pop[[i]]$ihr3 = P.hosp
    params$pop[[i]]$iir3 = P.critical
    
    params$pop[[i]]$dDeath = cm_delay_lnorm(15, 0.9, 60, 0.25)$p;
    params$pop[[i]]$dHosp  = cm_delay_gamma(6.0 + 2.5, 0.71, 60, 0.25)$p;
    params$pop[[i]]$lHosp  = cm_delay_lnorm(11.08, 1.202, 60, 0.25)$p;
    params$pop[[i]]$dICU   = cm_delay_gamma(9.6 + 2.5, 1.91, 60, 0.25)$p;
    params$pop[[i]]$lICU   = cm_delay_lnorm(13.33, 1.25, 60, 0.25)$p;
}


# changes
schedule_all = readRDS(datapath(mobility_file));
schedule = list();
for (i in seq_along(schedule_all)) {
    if (schedule_all[[i]]$pops < N_REG) {
        schedule[[length(schedule) + 1]] = schedule_all[[i]]
    }
}

# Remove NAs
for (i in seq_along(schedule)) {
    for (j in seq_along(schedule[[i]]$values)) {
        if (any(is.na(schedule[[i]]$values[[j]]))) {
            schedule[[i]]$values[[j]] = ifelse(is.na(schedule[[i]]$values[[j]]), 
                                               prev, schedule[[i]]$values[[j]])
        }
        prev = schedule[[i]]$values[[j]];
    }
}
params$schedule = schedule

#
# Individual fits
#

source("./cpp_funcs.R")

# Fitting
priorsI = list(
    tS = "U 0 60",
    u = "N 0.09 0.02 T 0.04 0.2",
    death_mean = "N 15 2 T 5 30",    # <<< co-cin
    hosp_admission = "N 8 1 T 4 20", # <<< co-cin
    icu_admission = "N 12.5 1 T 8 14", # <<< co-cin
    cfr_rlo = "N 0 0.1 T -2 2",
    cfr_rlo2 = "N 0 0.1 T -2 2",
    cfr_rlo3 = "N 0 0.1 T -2 2",
    hosp_rlo = "N 0 0.1 T -2 2", 
    icu_rlo = "N 0 0.1 T -2 2",
    icu_rlo2 = "N 0 0.1 T -2 2",
    contact_final = "N 1 0.1 T 0 1",
    contact_s0 = "E 0.1 0.1",
    contact_s1 = "E 0.1 0.1",
    disp_deaths = "E 10 10",
    disp_hosp_inc = "E 10 10",
    disp_hosp_prev = "E 10 10",
    disp_icu_prev = "E 10 10",
    concentration1 = "N 2 .3 T 2 10",
    concentration2 = "N 2 .2 T 2 10",
    concentration3 = "N 2 .1 T 2 10",
    xmas_fudge = "N 1 0.25 T 0 2",
    
    f102 = "L 0 0.1 T 0.5 2",
    f144 = "L 0 0.1 T 0.5 2",
    f186 = "L 0 0.1 T 0.5 2",
    f228 = "L 0 0.1 T 0.5 2",
    f270 = "L 0 0.1 T 0.5 2",
    f312 = "L 0 0.1 T 0.5 2",
    f354 = "L 0 0.1 T 0.5 2",
    f396 = "L 0 0.1 T 0.5 2",
    f438 = "L 0 0.1 T 0.5 2",
    f480 = "L 0 0.1 T 0.5 2",
    f522 = "L 0 0.1 T 0.5 2",
    f564 = "L 0 0.1 T 0.5 2",
    f606 = "L 0 0.1 T 0.5 2",
    f648 = "L 0 0.1 T 0.5 2",
    f690 = "L 0 0.1 T 0.5 2",
    f732 = "L 0 0.1 T 0.5 2", # 6-week period beginning 2nd January 2022
    f774 = "L 0 0.1 T 0.5 2", # 6-week period beginning 13th February 2022
    f816 = "L 0 0.1 T 0.5 2"  # 6-week period beginning 27th March 2022 (until 8th May 2022)
);
constants = list();

if (opt_v2) {
    priorsI = c(priorsI, list(
        v2_when = "U 144 365",
        v2_sgtf0 = "B 1.5 15",
        v2_disp = "E 10 10 T 0 0.25",
        v2_hosp_rlo = "N 0 0.1 T -4 4",
        v2_icu_rlo = "N 0 0.1 T -4 4",
        v2_cfr_rlo = "N 0 0.1 T -4 4"
    ));
}

if (opt_v3) {
    priorsI = c(priorsI, list(
        v3_when = "U 366 486" # for now, we just fit relative tx and start date for delta between 1st Jan and 1st May 2021
    ));
}

if (opt_v4) {
    priorsI = c(priorsI, list(
        v4_when = "N 685 7 T 670 700", # fit relative tx and start date for Omicron between 1st November 2021 and 1st December 2021
        v4_sgtf0 = "B 1.5 15", # for now, just copying the same priors we used originally for Alpha's SGTF data
        v4_disp = "E 10 10 T 0 0.25"
    ));
    
}

priorsI = c(priorsI, extra_priors);

posteriorsI = list()
dynamicsI = list()
parametersI = list()

init_previous = TRUE
init_previous_amount = 1

if (START_FILE == "scratch") {
    cat("Starting fit from scratch.\n");
} else if (START_FILE == "previous") {
    saved = qread(paste0("./fits/", set_id, REP_START - 1, ".qs"))
    posteriorsI = saved[[1]]
    parametersI = saved[[2]]
    rm(saved)
} else {
    saved = qread(paste0("./fits/", START_FILE, ".qs"))
    posteriorsI = saved[[1]]
    parametersI = saved[[2]]
    rm(saved)
}

for (i in seq_along(posteriorsI)) {
    if (!is.null(posteriorsI[[i]]) && "v2_conc" %in% names(posteriorsI[[i]])) {
        posteriorsI[[i]][, v2_disp := 1 / sqrt(v2_conc)];
    }
}

# notify Rosie
notify_command = paste0(
    'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
    '--form-string "user=ubdi7mpz6bfiy3a485qkgpt5kzk9mb" ',
    '--form-string "message=Commencing fitting." https://api.pushover.net/1/messages.json')
system(notify_command)

# Elimination of unneeded burden-related parameters, 1st Oct 2021
constants = list(
    # tS
    # u
    # death_mean = 17.72438,
    # hosp_admission = 13.16867,
    # icu_admission = 13.496,
    cfr_rlo = 0,
    cfr_rlo2 = 0,
    cfr_rlo3 = 0,
    hosp_rlo = 0,
    icu_rlo = 0,
    icu_rlo2 = 0,
    contact_final = 1,
    contact_s0 = 1,
    contact_s1 = 1,
    disp_deaths = 0,
    disp_hosp_inc = 0,
    disp_hosp_prev = 0,
    disp_icu_prev = 0,
    concentration1 = 2,
    concentration2 = 2,
    concentration3 = 2,
    xmas_fudge = 1
    # v2_when = 255.1093,
    # v2_sgtf0 = 0.03198553,
    # v2_disp = 0.0603328,
    # v2_hosp_rlo = 0.192403,
    # v2_icu_rlo = 0.4294857,
    # v2_cfr_rlo = 0.02273579,
    # v3_when = 414.4822,
    # v2_relu = 1.619153,
    # v3_relu = 1.928168
)

priorsI = priorsI[setdiff(names(priorsI), names(constants))];

# Define number of threads to use
N_THREADS = length(priorsI) * 2;

for (replic in REP_START:REP_END)
{
    # Loop through regions
    for (pn in which_pops) {
        paramsI = rlang::duplicate(params);
        paramsI$pop = list(rlang::duplicate(params$pop[[pn]]));
        paramsI$travel = matrix(1, nrow = 1, ncol = 1);
        paramsI$schedule = list();
        j = 1;
        for (i in seq_along(params$schedule)) {
            if (pn - 1 == params$schedule[[i]]$pops) {
                paramsI$schedule[[j]] = rlang::duplicate(params$schedule[[i]]);
                paramsI$schedule[[j]]$pops = 0;
                j = j + 1;
            }
        }
        
        # contact placeholder for tier 2
        paramsI$schedule[[2]] = rlang::duplicate(paramsI$schedule[[1]]);
        for (i in seq_along(paramsI$schedule[[2]]$values)) {
            paramsI$schedule[[2]]$values[[i]][1] = paramsI$schedule[[1]]$values[[i]][1] +  0.2497655 / 100;
            paramsI$schedule[[2]]$values[[i]][2] = paramsI$schedule[[1]]$values[[i]][2] + -0.2307939 / 100;
            paramsI$schedule[[2]]$values[[i]][3] = paramsI$schedule[[1]]$values[[i]][3] + -1.5907698 / 100;
            paramsI$schedule[[2]]$values[[i]][4] = paramsI$schedule[[1]]$values[[i]][4] + -3.4866544 / 100;
            paramsI$schedule[[2]]$values[[i]][5] = paramsI$schedule[[1]]$values[[i]][5] + -3.4524518 / 100;
        }
        paramsI$schedule[[2]]$mode = "bypass";
        
        # contact placeholder for tier 3
        paramsI$schedule[[3]] = rlang::duplicate(paramsI$schedule[[1]]);
        for (i in seq_along(paramsI$schedule[[3]]$values)) {
            paramsI$schedule[[3]]$values[[i]][1] = paramsI$schedule[[1]]$values[[i]][1] +  2.080457 / 100;
            paramsI$schedule[[3]]$values[[i]][2] = paramsI$schedule[[1]]$values[[i]][2] + -8.045226 / 100;
            paramsI$schedule[[3]]$values[[i]][3] = paramsI$schedule[[1]]$values[[i]][3] + -2.476266 / 100;
            paramsI$schedule[[3]]$values[[i]][4] = paramsI$schedule[[1]]$values[[i]][4] + -10.144043 / 100;
            paramsI$schedule[[3]]$values[[i]][5] = paramsI$schedule[[1]]$values[[i]][5] + -7.681244 / 100;
        }
        paramsI$schedule[[3]]$mode = "bypass";
        
        # contact multiplier for gradual contact change
        paramsI$schedule[[4]] = list(
            parameter = "contact",
            pops = 0,
            mode = "multiply",
            values = rep(list(rep(1, 8)), 366),
            times = 0:365
        )
        
        # contact multiplier for xmas fudgery
        # between Friday 18th December 2020 and Friday 1st January 2021
        paramsI$schedule[[5]] = list(
            parameter = "contact",
            pops = 0,
            mode = "multiply",
            values = list(rep(1, 8), rep(1, 8)),
            times = c(352, 366)
        )
        
        # contact multiplier for contact adjustment
        adjust_days = as.numeric(ymd(date_fitting) - ymd("2020-01-01"))
        paramsI$schedule[[6]] = list(
            parameter = "contact",
            pops = 0,
            mode = "multiply",
            values = rep(list(rep(1, 8)), adjust_days),
            times = 0:(adjust_days - 1)
        )
        
        # # Transmission rate multiplier for BA.2
        # paramsI$schedule[[7]] = list(
        #     parameter = "u",
        #     pops = 0,
        #     mode = "multiply",
        #     values = lapply(ba2f, rep, 16),
        #     times = ba2$t
        # )
        
        ldI = rlang::duplicate(ld);
        ldI = ldI[pid == pn - 1];
        sitrepsI = rlang::duplicate(sitreps);
        sitrepsI = sitrepsI[pid == pn - 1];
        seroI = rlang::duplicate(sero_to_fit);
        seroI = seroI[pid == pn - 1 & Data.source != "NHSBT"];   # sero: all but NHSBT
        virusI = rlang::duplicate(virus);
        virusI = virusI[pid == pn - 1 & Data.source %like% "ONS"]; # virus: ONS-CIS only
        sgtfI = copy(sgtf);
        sgtfI = sgtfI[pid == pn - 1 & date <= sgtf_stop];
        deltaI = copy(delta)
        deltaI = deltaI[pid == pn - 1 & date > as.Date("2021-04-01")] # first delta sequence was recorded on 2nd April 2021 (checked 21st June 2021)
        omiI = copy(omi)
        omiI = omiI[pid == pn - 1 & date > as.Date("2021-10-01")] # we allow model to fit introduction time of Omicron from 1st October 2021 to January 2022
        
        # load user defined functions
        model_v3_contents = list(
            cpp_changes = cpp_chgI_voc(priorsI, constants, seasonality = opt_seas,
                                       v2 = opt_v2, v2_relu = opt_relu, v2_latdur = opt_latdur, 
                                       v2_serial = opt_serial, v2_infdur = opt_infdur, 
                                       v2_immesc = opt_immesc, v2_ch_u = opt_ch_u,
                                       v3_relu = NULL, v3_severity = V3_SEVERITY, v3 = opt_v3),
            cpp_loglikelihood = cpp_likI_voc_omi(paramsI, ldI, sitrepsI, seroI, virusI, sgtfI, 
                                                 pn, date_fitting, priorsI, constants, death_cutoff = 0,
                                                 use_sgtf = opt_v2, delta = deltaI, omi = omiI, gamdisp = GAMDISP),
            cpp_observer = c(
                cpp_obsI_voc(concentration = opt_conc, v2 = opt_v2, 
                             P.death, P.critical, priorsI, constants, v3_severity = V3_SEVERITY),
                cpp_obsI_voc_omi(OMI_SETUP_T, OMI_X_PROTECTION, OMI_VAX_FACTOR, OMI_VAX_ASSUMPTION, OMI_SEEDS_PER_DAY,
                    omi_sev = OMI_SEV, omi_crit = OMI_CRIT),
                cpp_obsI_vax(paramsI, vacc[[pn]]),
                if (SEAS_YN == "seasyes") cpp_obsI_seasonality(SEAS_AMP, 1) else if (SEAS_YN == "seaslate") cpp_obsI_seasonality(SEAS_AMP, 456) else "",
                # cpp_obsI_aw(seasonality_aw = SEAS_AMP, 0.0, NA, NA),
                # cpp_obsI_voc_nu(
                #     setup_t = 630,
                #     nu_voc_t = 685,
                #     tx_factor = 2.05,
                #     x_protection = 0.551,
                #     vax_factor = 0.551,
                #     vax_assumption = "khoury",
                #     sev_factor = 1.0,
                #     n_seeds_per_day = 10
                # ),
                cpp_obsI_booster(
                    target_phase1 = 229000,
                    target_phase2 = 1000000,
                    proportion_booster = PBOOST,
                    booster_fold = BFOLD,
                    booster_om_fold = 1,
                    booster_duration = 180
                ),
                cpp_obsI_voc_ba2(ba2t, ba2f)
            )
        )
        
        cm_source_backend(
            user_defined = list(
                model_v3 = model_v3_contents
            )
        )
        
        priorsI2 = rlang::duplicate(priorsI)
        if (init_previous) {
            for (k in seq_along(priorsI2)) {
                pname = names(priorsI2)[k];
                if (length(posteriorsI) >= pn && pname %in% names(posteriorsI[[pn]])) {
                    init_values = quantile(posteriorsI[[pn]][[pname]], c(0.025, 0.975));
                    cat(paste0("Using 95% CI ", init_values[1], " - ", init_values[2], " for initial values of parameter ", pname, 
                               " with probability ", init_previous_amount, "\n"));
                    priorsI2[[pname]] = paste0(priorsI2[[pname]], " I ", init_values[1], " ", init_values[2], " ", init_previous_amount);
                    cat(paste0(priorsI2[[pname]], "\n"));
                } else {
                    cat(paste0("Could not find init values for parameter ", pname, "\n"));
                    cat(paste0(priorsI2[[pname]], "\n"));
                }
            }
        }
        
        # Run MCMC, with error handling
        qsave(model_v3_contents, "./last_run_model_v3_contents.qs");
        qsave(paramsI, "./last_run_paramsI.qs");
        qsave(priorsI2, "./last_run_priorsI2.qs");
        
        tryCatch({
            postI = cm_backend_mcmc_test(cm_translate_parameters(paramsI), priorsI2,
                                         seed = 0, 
                                         burn_in = ifelse(replic == REP_END, BURN_IN_FINAL, BURN_IN), 
                                         iterations = ifelse(replic == REP_END, ITER_FINAL, ITER), 
                                         n_threads = N_THREADS, classic_gamma = T);
        },
        error = function(e) {
            
            # Notify Rosie
            notify_command = paste0(
                'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
                '--form-string "user=ubdi7mpz6bfiy3a485qkgpt5kzk9mb" ',
                '--form-string "message=Error caught in rep ', replic, ' region ', pn, '. Stopping." https://api.pushover.net/1/messages.json')
            system(notify_command)
            
            print(e)
            
            stop("Error caught in rep ", replic, " region ", pn, ". Stopping.")
        })
        
        
        setDT(postI)
        # Add constants to posteriors
        postI = cbind(postI, as.data.table(constants))
        posteriorsI[[pn]] = postI
        
        # Notify Rosie
        notify_command = paste0(
            'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
            '--form-string "user=ubdi7mpz6bfiy3a485qkgpt5kzk9mb" ',
            '--form-string "message=Rep ', replic, ' region ', pn, ' done fitting." https://api.pushover.net/1/messages.json')
        system(notify_command)
        
        parametersI[[pn]] = rlang::duplicate(paramsI)
        qsave(rlang::duplicate(list(posteriorsI, parametersI)), paste0("./fits/", set_id, replic, "-progress.qs"))
        
        print(pn)
    }
    
    qsave(rlang::duplicate(list(posteriorsI, parametersI)), paste0("./fits/", set_id, replic, ".qs"))
    
    # saved = qread("./fits/relu_yeswane_sev2_2202110020005.qs")
    # posteriorsI = saved[[1]]
    # parametersI = saved[[2]]
    
    # Generate SPI-M output
    # Sample dynamics from fit
    dynamicsI = list()
    for (pn in which_pops)  {
        cat(paste0("Sampling fit for population ", pn, "...\n"))
        
        # Source backend
        cm_source_backend(
            user_defined = list(
                model_v3 = list(
                    cpp_changes = cpp_chgI_voc(priorsI, constants, seasonality = opt_seas, 
                                               v2 = opt_v2, v2_relu = opt_relu, 
                                               v2_latdur = opt_latdur, 
                                               v2_serial = opt_serial, 
                                               v2_infdur = opt_infdur, 
                                               v2_immesc = opt_immesc, 
                                               v2_ch_u = opt_ch_u,
                                               v3_relu = NULL, 
                                               v3_severity = V3_SEVERITY, 
                                               v3 = opt_v3),
                    cpp_loglikelihood = "",
                    cpp_observer = c(cpp_obsI_voc(concentration = opt_conc, 
                                                  v2 = opt_v2, 
                                                  P.death, 
                                                  P.critical, 
                                                  priorsI, 
                                                  constants,
                                                  v3_severity = V3_SEVERITY),
                                     cpp_obsI_voc_omi(OMI_SETUP_T, OMI_X_PROTECTION, OMI_VAX_FACTOR, OMI_VAX_ASSUMPTION, OMI_SEEDS_PER_DAY,
                                         omi_sev = OMI_SEV, omi_crit = OMI_CRIT),
                                     cpp_obsI_vax(parametersI[[pn]], vacc[[pn]]),
                                     if (SEAS_YN == "seasyes") cpp_obsI_seasonality(SEAS_AMP, 1) else if (SEAS_YN == "seaslate") cpp_obsI_seasonality(SEAS_AMP, 456) else "",
                                     cpp_obsI_booster(
                                         target_phase1 = 229000,
                                         target_phase2 = 1000000,
                                         proportion_booster = PBOOST,
                                         booster_fold = BFOLD,
                                         booster_om_fold = 1,
                                         booster_duration = 180
                                     ),
                                     cpp_obsI_voc_ba2(ba2t, ba2f)
                                     )
                )
            )
        )
        
        # Sampling fits
        paramsI2 = rlang::duplicate(parametersI[[pn]])
        paramsI2$time1 = as.character(analysis_end);
        test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriorsI[[pn]], 100, seed = 0, n_threads = 6); # do we want to change n_threads to something like N_THREADS = 54 + length(extra_priors) * 2; ??? 
        rows = cm_backend_sample_fit_rows(cm_translate_parameters(paramsI2), posteriorsI[[pn]], 100, seed = 0);
        
        test = rbindlist(test)
        test[, population := pn]
        
        # Add unvaccinated outputs
        test[, deaths_V0 := deaths - deaths_V1 - deaths_V2 - deaths_V3]
        test[, hosp_adm_V0 := hosp_adm - hosp_adm_V1 - hosp_adm_V2 - hosp_adm_V3]
        
        # Add dispersion parameters
        disp = posteriorsI[[pn]][rows, .SD, .SDcols = patterns("^disp|v2_conc|v2_disp|v2_sgtf0|v4_sgtf0|v4_disp")]
        disp[, run := .I]
        test = merge(test, disp, by = "run")
        
        dynamicsI[[pn]] = test
    }
    
    # Concatenate dynamics for SPI-M output
    test = rbindlist(dynamicsI, fill = TRUE)
    test[, population := nhs_regions[population]]
    
    # Plot outputs which are fitted to data PRIOR to gamma multiplier adjustment
    plot_pa = check_fit_small_output(test, parametersI, ld, sitreps, virus, sero, nhs_regions[which_pops], death_cutoff = 0, date_fitting, min_date = NULL, sero_cut_off = SERO_CUT_OFF)
    ggsave(paste0("./output/fitsm_noadjust_", set_id, replic, ".pdf"), plot_pa, width = 80 * length(which_pops) / 7, height = 50, units = "cm", useDingbats = FALSE)

    # Adjust for variable infection burden ratios
    test = apply_gamma_multiplier(test, "deaths", "disp_deaths", ld, "N", 0.1, 7, GAMDISP, c("deaths_V1", "deaths_V2", "deaths_V3"), adj_file = NA, "2021-11-01", "2021-11-30", adj_ts = NA, '2020-01-01', '2022-09-30')
    test = apply_gamma_multiplier(test, "hosp_adm", "disp_hosp_inc", sitreps, "n_admitted_diagnosed", 0.1, 7, GAMDISP, c("hosp_adm_V1", "hosp_adm_V2", "hosp_adm_V3", "hosp_undetected_o"), adj_file = NA, "2021-11-01", "2021-11-30", adj_ts = NA, '2020-01-01', '2022-09-30')
    test[, known_hosp_beds := hosp_bed - hosp_undetected_p]
    test = apply_gamma_multiplier(test, "known_hosp_beds", "disp_hosp_prev", sitreps, "n_in_all_beds", 0.1, 7, GAMDISP, c("hosp_bed", "hosp_undetected_p"), adj_file = NA, "2021-11-01", "2021-11-30", adj_ts = NA, '2020-01-01', '2022-09-30')
    test[, known_hosp_beds := NULL]
    test = apply_gamma_multiplier(test, "icu_bed", "disp_icu_prev", sitreps, "n_in_itu", 0.1, 7, GAMDISP, "icu_adm", adj_file = NA, "2021-11-01", "2021-11-30", adj_ts = NA, '2020-01-01', '2022-09-30')
    
    # Visually inspect fit
    
    # Plot outputs which are fitted to data
    plot2 = check_fit_small_output(test, parametersI, ld, sitreps, virus, sero, nhs_regions[which_pops], death_cutoff = 0, date_fitting, min_date = NULL, sero_cut_off = SERO_CUT_OFF)
    ggsave(paste0("./output/fitsm_", set_id, replic, ".pdf"), plot2, width = 80 * length(which_pops) / 7, height = 50, units = "cm", useDingbats = FALSE)
    
    # Plot all outputs
    plot = check_fit(test, parametersI, ld, sitreps, virus, sero, nhs_regions[which_pops], death_cutoff = 0, date_fitting, min_date = NULL, sero_cut_off = SERO_CUT_OFF)
    #plot = plot + geom_vline(aes(xintercept = ymd("2020-12-24")), size = 0.25, linetype = "42")
    ggsave(paste0("./output/fit_", set_id, replic, ".pdf"), plot, width = 80 * length(which_pops) / 7, height = 65, units = "cm", useDingbats = FALSE)
    
    # Posteriors
    post = rbindlist(posteriorsI, idcol = "population", fill = TRUE)
    post[, pop := nhs_regions[population]]
    melted = melt(post, id.vars = c(1:5, ncol(post)))
    plot = ggplot(melted) + geom_density(aes(x = value, colour = pop)) + facet_wrap(~variable, scales = "free")
    ggsave(paste0("./output/post_", set_id, replic, ".pdf"), plot, width = 40, height = 30, units = "cm", useDingbats = FALSE)
    
    # Fit to SGTF data
    if (opt_v2) {
        sgtf[, qlo := qbeta(0.025, sgtf + 1, other + 1)]
        sgtf[, qhi := qbeta(0.975, sgtf + 1, other + 1)]
        vmodel = test[, .(I1 = sum(test_o + test3_o), I2 = sum(test2_o), sgtf0 = v2_sgtf0[1], conc = 1/(v2_disp[1]*v2_disp[1])), by = .(t, population, run)]
        vmodel[, p2 := I2 / (I1 + I2)]
        vmodel[is.nan(p2), p2 := 0]
        vmodel[, sgtf := (1 - p2) * sgtf0 + p2];
        vmodel[, alpha := sgtf * (conc - 2) + 1]
        vmodel[, beta := (1 - sgtf) * (conc - 2) + 1]
        vmodel[, q025 := qbeta(0.025, alpha, beta)]
        vmodel[, q500 := qbeta(0.500, alpha, beta)]
        vmodel[, q975 := qbeta(0.975, alpha, beta)]
        vmodel = vmodel[, lapply(.SD, mean), .SDcols = c("q025", "q500", "q975"), by = .(nhs_name = population, t)]
        plotS = ggplot(sgtf[(pid + 1) %in% which_pops]) +
            geom_ribbon(aes(x = date, ymin = qlo, ymax = qhi), fill = "black", alpha = 0.1) +
            geom_ribbon(data = vmodel[t + ymd("2020-01-01") >= "2020-10-01"],
                        aes(x = ymd("2020-01-01") + t, ymin = q025, ymax = q975), fill = "darkorchid", alpha = 0.5) +
            geom_line(data = vmodel[t + ymd("2020-01-01") >= "2020-10-01"],
                      aes(x = ymd("2020-01-01") + t, y = q500), colour = "darkorchid") +
            geom_line(aes(x = date, y = sgtf / (sgtf + other)), size = 0.25) +
            #scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99), limits = c(0.01, 0.99)) +
            facet_wrap(~nhs_name) +
            labs(x = NULL, y = "Relative frequency of\nS gene target failure") +
            scale_x_date(date_breaks = "1 month", date_labels = "%b")
        ggsave(paste0("./output/sgtf_check_", set_id, replic, ".pdf"), plotS, width = 20, height = 6, units = "cm", useDingbats = FALSE)
    }
    
    # Fit to Delta data
    if (opt_v3) {
        delta[, qlo := qbeta(0.025, delta + 1, other + 1)]
        delta[, qhi := qbeta(0.975, delta + 1, other + 1)]
        
        dmodel = test[, .(D1 = sum(test3_o), D2 = sum(test_o), D3 = sum(test2_o), delta0 = 0, conc = 30), by = .(t, population, run)]
        dmodel[, dprop := D1 / (D1+D2+D3)]
        dmodel[is.nan(dprop), dprop := 0]
        dmodel[, delta := (1-dprop) * delta0 + dprop]
        dmodel[, alpha := delta * (conc - 2) + 1]
        dmodel[, beta := (1 - delta) * (conc - 2) + 1]
        dmodel[, q025 := qbeta(0.025, alpha, beta)]
        dmodel[, q500 := qbeta(0.500, alpha, beta)]
        dmodel[, q975 := qbeta(0.975, alpha, beta)]
        dmodel = dmodel[, lapply(.SD, mean), .SDcols = c("q025", "q500", "q975"), by = .(nhs_name = population, t)]
        
        
        plotS = ggplot(delta[(pid + 1) %in% which_pops]) +
            geom_ribbon(aes(x = date, ymin = qlo, ymax = qhi), fill = "black", alpha = 0.1) +
            geom_ribbon(data = dmodel[t + ymd("2020-01-01") >= "2021-01-01"],
                        aes(x = ymd("2020-01-01") + t, ymin = q025, ymax = q975), fill = "darkorchid", alpha = 0.5) +
            geom_line(data = dmodel[t + ymd("2020-01-01") >= "2021-01-01"],
                      aes(x = ymd("2020-01-01") + t, y = q500), colour = "darkorchid") +
            geom_line(aes(x = date, y = delta / (delta + other)), size = 0.25) +
            #scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99), limits = c(0.01, 0.99)) +
            facet_wrap(~nhs_name) +
            labs(x = NULL, y = "Relative frequency of\nDelta B.1.617.2 VOC") +
            scale_x_date(date_breaks = "1 month", date_labels = "%b")
        ggsave(paste0("./output/delta_check_", set_id, replic, ".pdf"), plotS, width = 20, height = 6, units = "cm", useDingbats = FALSE)
    }
    
    # Fit to SGTF Omicron data
    if (opt_v4) {
        omi[, qlo := qbeta(0.025, sgtf + 1, other + 1)]
        omi[, qhi := qbeta(0.975, sgtf + 1, other + 1)]
       
        omodel = test[, .(IO = sum(test_o), INO = sum(test2_o + test3_o), sgtfomi0 = v4_sgtf0[1], conc = 1/(v4_disp[1]*v4_disp[1])), by = .(t, population, run)]
        omodel[, o2 := IO / (IO + INO)]
        omodel[is.nan(o2), o2 := 0]
        omodel[, sgtf := (1 - o2) * sgtfomi0 + o2]
        omodel[, alpha := sgtf * (conc - 2) + 1]
        omodel[, beta := (1 - sgtf) * (conc - 2) + 1]
        omodel[, q025 := qbeta(0.025, alpha, beta)]
        omodel[, q500 := qbeta(0.500, alpha, beta)]
        omodel[, q975 := qbeta(0.975, alpha, beta)]
        omodel = omodel[, lapply(.SD, mean), .SDcols = c("q025", "q500", "q975"), by = .(nhs_name = population, t)]
        
        plotO = ggplot(omi[(pid + 1) %in% which_pops]) +
            geom_ribbon(aes(x = date, ymin = qlo, ymax = qhi), fill = "black", alpha = 0.1) +
            geom_ribbon(data = omodel[t + ymd("2020-01-01") >= "2021-10-01"],
                        aes(x = ymd("2020-01-01") + t, ymin = q025, ymax = q975), fill = "darkorchid", alpha = 0.5) +
            geom_line(data = omodel[t + ymd("2020-01-01") >= "2021-10-01"],
                      aes(x = ymd("2020-01-01") + t, y = q500), colour = "darkorchid") +
            geom_line(aes(x = date, y = sgtf / (sgtf + other)), size = 0.25) +
            #scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99), limits = c(0.01, 0.99)) +
            facet_wrap(~nhs_name) +
            labs(x = NULL, y = "Relative frequency of\nS gene target failure") +
            scale_x_date(date_breaks = "1 month", date_labels = "%b %Y")
        ggsave(paste0("./output/omi_sgtf_check_", set_id, replic, ".pdf"), plotO, width = 20, height = 6, units = "cm", useDingbats = FALSE)
    }
    
    # save output to inspect later (saveRDS() call is last thing in this script)
    # restrict time span here
    ymd_from <- "2020-01-01"
    ymd_to <- paramsI2$time1
    t0 = as.numeric(ymd(ymd_from) - ymd("2020-01-01"))
    t1 = as.numeric(ymd(ymd_to) - ymd("2020-01-01"))
    test0 = test[t %between% c(t0, t1)]
    cat(paste0("Max number in Va1 compartment is ", max(test0$Va1), "\n"));
    cat(paste0("Max number in Va2 compartment is ", max(test0$Va2), "\n"));
    cat(paste0("Max number in Va3 compartment is ", max(test0$Va3), "\n"));
    cat(paste0("Max number in Vb1 compartment is ", max(test0$Vb1), "\n"));
    cat(paste0("Max number in Vb2 compartment is ", max(test0$Vb2), "\n"));
    cat(paste0("Max number in Vb3 compartment is ", max(test0$Vb3), "\n"));
    
    # save test0 output
    # saveRDS(test0, file = paste0("./output/test_output_", set_id, replic, ".rds"))
    
    # create SPI-M output
    mtp_output = SPIM_output_full(test, creation_year, creation_month, creation_day, forecast_start, forecast_end)
    
    plot = ggplot(mtp_output[AgeBand == "All"], aes(x = make_date(`Year of Value`, `Month of Value`, `Day of Value`))) +
        geom_ribbon(aes(ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = `Geography`)) +
        facet_grid(ValueType ~ Geography, scales = "free") +
        theme(legend.position = "none")
    
    ggsave(paste0("./output/mtp_check_spim_", set_id, replic, ".pdf"), plot, width = 40 * length(which_pops) / 7, height = 25, units = "cm", useDingbats = FALSE)
    
    # Save output
    fwrite(mtp_output, paste0("./output/1SPIM_mtp_", set_id, replic, ".csv"))
    
    # # Save log-scale plots
    # ggsave(paste0("./output/fitsm_log_", set_id, replic, ".pdf"), plot2 + scale_y_log10(limits = c(1,NA)), width = 80 * length(which_pops) / 7, height = 50, units = "cm", useDingbats = FALSE)
    # ggsave(paste0("./output/fit_log_", set_id, replic, ".pdf"), plot + scale_y_log10(limits = c(1,NA)), width = 80 * length(which_pops) / 7, height = 65, units = "cm", useDingbats = FALSE)
    # 
    # Notify Rosie
    notify_command = paste0(
        'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
        '--form-string "user=ubdi7mpz6bfiy3a485qkgpt5kzk9mb" ',
        '--form-string "message=Figures for rep ', replic, ' are saved." https://api.pushover.net/1/messages.json')
    system(notify_command)
}


# Notify Rosie
notify_command = paste0(
    'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
    '--form-string "user=ubdi7mpz6bfiy3a485qkgpt5kzk9mb" ',
    '--form-string "message=All jobs finished." https://api.pushover.net/1/messages.json')
system(notify_command)
