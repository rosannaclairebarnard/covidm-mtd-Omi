# Construction of traces via particle filter
#
# at the command line, use
#
# Rscript pfilter.R FIT_FILE
#
# e.g. Rscript pfilter.R relu_yeswane_sev2.0_22012904

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
library(forecast)
library(zoo)

theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

######################## CHANGE HERE FOR EACH RUN ########################

# c_args = c('relu_test_yeswane_sev2.0_22041306')
c_args = commandArgs(trailingOnly = TRUE)
FIT_FILE = c_args[[1]]

datetime <- str_replace_all(Sys.time(), "[- :BSTGMT]", "")

FIT_FILE_PATH = paste0('./fits/', FIT_FILE, '.qs')
PF_FILE_TO_SAVE = paste0('pf_', FIT_FILE, '_', datetime)

data_file = "processed-data-20220506122858.qs"
mobility_file = "schedule3-MTPs-20220506121302.rds"
date_fitting = "2022-05-06"
vax_file = "vax-covidm20220505205235.rds"
sgtf_stop = "2021-02-15"

# SERO_CUT_OFF determines the date at which seroprevalence data is no longer 
# used to fit to (we use the start date of estimates for cut off)
SERO_CUT_OFF = '2020-12-01'

# default number of particles is 500; set PT_MULTIPLIER=1 for no change
PT_MULTIPLIER = 1

# particle filter time limit is now calculated further down corresponding to the
# last data point in virus prevalence data

# set PT_REDUCE to check if stopping particle filtering earlier achieves better fit
PT_REDUCE = 0 # set PT_REDUCE=0 to have PT limit equivalent to last PCR datapt

UPDATE_PARAMS = TRUE

BFOLD = 2.5 # fold increase in neutralisation titres after boosters (2.5 / 4.9)

DFR = 1.0  # Delta fold reduction (e.g. 1.0 or 3.9)

# probability of receiving a booster dose by covidm 5-year age groups (this 
# should match the PBOOST assumptions used in fit-omi.R)
PBOOST=c(0,   0,     0,   0.4, 0.544, 0.597, 0.643, 0.698, 0.757, 0.809, 0.861, 0.892, 0.917, 0.946, 0.964,  0.967)

OMI_SETUP_T = 630 # switch variants on 22nd September 2021 (on this day variant 
# 1 is reparameterised from wildtype to Omicron and R1 individuals move to R3)
OMI_VAX_ASSUMPTION = 'khoury'
OMI_SEEDS_PER_DAY = 10
OMI_X_PROTECTION = 0.551
OMI_VAX_FACTOR = 0.551
OMI_SEV = 0.5
OMI_CRIT = 0.5
BA2_RELU = 1.5
GAMDISP = 0.3

WHICH_POPS = c(1, 3, 4, 5, 6, 9, 10)

############################## END CHANGES ##############################

which_pops = WHICH_POPS

# Command line
c_args = list("relu", "yeswane", "central", "2.0", "seaslate", "0.1", FIT_FILE_PATH, 
              PF_FILE_TO_SAVE)
if (length(c_args) != 8) {
    stop("8 arguments required.")
}
FIT_TYPE = c_args[[1]];
POP_SET = "all";
WANE_YN = c_args[[2]]
VAC_EFF = c_args[[3]]
V3_SEVERITY = c_args[[4]]
SEAS_YN = c_args[[5]]
if (!SEAS_YN %in% c("seasyes", "seasno", "seaslate")){
    stop("Seasonality option should be seasyes, seasno, seaslate")
}
SEAS_AMP = as.numeric(c_args[[6]])
FIT_FILE = c_args[[7]]
PF_FILE_TO_SAVE = c_args[[8]]

# particle filter file name to be saved at end
PF_FILE = paste0("./fits/", PF_FILE_TO_SAVE, ".qs")

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

# calculate time limit for particle filter based on last virus prevalence datapt
PTLIMIT = as.numeric(as.Date(max(virus$Start.date)) - as.Date('2020-01-01')) - PT_REDUCE

#
# FITTING
#

# NUMBER OF REGIONS TO FIT
N_REG = 12;

# Build parameters for NHS regions
params = cm_parameters_SEI3R(nhs_regions[1:N_REG], deterministic = T, 
                             date_start = "2020-01-01", 
                             date_end = date_fitting,
                             dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
                             dIp = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
                             dIs = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
                             dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 0.25)$p)
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
    f816 = "L 0 0.1 T 0.5 2"
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

for (i in seq_along(posteriorsI)) {
    if (!is.null(posteriorsI[[i]]) && "v2_conc" %in% names(posteriorsI[[i]])) {
        posteriorsI[[i]][, v2_disp := 1 / sqrt(v2_conc)];
    }
}

# notify Rosie
notify_command = paste0(
    'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
    '--form-string "user=ubdi7mpz6bfiy3a485qkgpt5kzk9mb" ',
    '--form-string "message=Commencing particle filtering." https://api.pushover.net/1/messages.json')
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

saved = qread(FIT_FILE)
posteriorsI = saved[[1]]
parametersI = saved[[2]]

which_pops = WHICH_POPS

########
########
########

sched = readRDS(paste0("./fitting_data/", mobility_file))

for (pn in seq_along(parametersI)) {
    if (!is.null(parametersI[[pn]])) {
        parametersI[[pn]]$time1 = date_fitting
        if (UPDATE_PARAMS) {
            parametersI[[pn]] = reinfection_defaults(parametersI[[pn]], 1)
            parametersI[[pn]] = waning_scenario(WANE_YN, parametersI[[pn]], 1)
            parametersI[[pn]] = VE_scenario    (VAC_EFF, PBOOST, 
                                                parametersI[[pn]], 1, 
                                                delta_fold_reduction = DFR)
            
            # SET MOBILITY
            
            # Set mobility scenario for no tier
            parametersI[[pn]]$schedule[[1]] = rlang::duplicate(sched[[pn]])
            parametersI[[pn]]$schedule[[1]]$pops = 0

            # ... for tier 2
            parametersI[[pn]]$schedule[[2]] = rlang::duplicate(parametersI[[pn]]$schedule[[1]]);
            for (i in seq_along(parametersI[[pn]]$schedule[[2]]$values)) {
                parametersI[[pn]]$schedule[[2]]$values[[i]][1] = parametersI[[pn]]$schedule[[1]]$values[[i]][1] +  0.2497655 / 100;
                parametersI[[pn]]$schedule[[2]]$values[[i]][2] = parametersI[[pn]]$schedule[[1]]$values[[i]][2] + -0.2307939 / 100;
                parametersI[[pn]]$schedule[[2]]$values[[i]][3] = parametersI[[pn]]$schedule[[1]]$values[[i]][3] + -1.5907698 / 100;
                parametersI[[pn]]$schedule[[2]]$values[[i]][4] = parametersI[[pn]]$schedule[[1]]$values[[i]][4] + -3.4866544 / 100;
                parametersI[[pn]]$schedule[[2]]$values[[i]][5] = parametersI[[pn]]$schedule[[1]]$values[[i]][5] + -3.4524518 / 100;
            }
            parametersI[[pn]]$schedule[[2]]$mode = "bypass";

            # ... for tier 3
            parametersI[[pn]]$schedule[[3]] = rlang::duplicate(parametersI[[pn]]$schedule[[1]]);
            for (i in seq_along(parametersI[[pn]]$schedule[[3]]$values)) {
                parametersI[[pn]]$schedule[[3]]$values[[i]][1] = parametersI[[pn]]$schedule[[1]]$values[[i]][1] +  2.080457 / 100;
                parametersI[[pn]]$schedule[[3]]$values[[i]][2] = parametersI[[pn]]$schedule[[1]]$values[[i]][2] + -8.045226 / 100;
                parametersI[[pn]]$schedule[[3]]$values[[i]][3] = parametersI[[pn]]$schedule[[1]]$values[[i]][3] + -2.476266 / 100;
                parametersI[[pn]]$schedule[[3]]$values[[i]][4] = parametersI[[pn]]$schedule[[1]]$values[[i]][4] + -10.144043 / 100;
                parametersI[[pn]]$schedule[[3]]$values[[i]][5] = parametersI[[pn]]$schedule[[1]]$values[[i]][5] + -7.681244 / 100;
            }
            parametersI[[pn]]$schedule[[3]]$mode = "bypass";
            
        } else {
            # Probability of having a booster dose
            parametersI[[pn]]$pop[[1]]$p_boost_va2 = rep(1, 16)
            parametersI[[pn]]$pop[[1]]$p_boost_vb2 = rep(1, 16)
            parametersI[[pn]]$pop[[1]]$p_boost_va3 = rep(1, 16)
            parametersI[[pn]]$pop[[1]]$p_boost_vb3 = rep(1, 16)
            
            # Adjust waning
            parametersI[[pn]]$pop[[1]]$dVa2 = cm_delay_erlang(30000, 10)$p
            parametersI[[pn]]$pop[[1]]$dVb2 = cm_delay_erlang(30000, 10)$p
        }
    }
}

dynamics = list()
pfilter = list()
posteriors = list()
parameters = list()
pcrp = list()

for (pn in which_pops) {
    # Get data and parameters
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
    
    paramsI = parametersI[[pn]]
    
    # load user defined functions
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
                cpp_loglikelihood = cpp_likI_voc_omi(paramsI, ldI, sitrepsI, seroI, virusI, sgtfI, 
                                                     pn, date_fitting, priorsI, constants, death_cutoff = 0,
                                                     use_sgtf = opt_v2, delta = deltaI, omi = omiI, gamdisp = GAMDISP),
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

    # Delete fudges
    posteriorsI[[pn]][, f102 := 1.0]
    posteriorsI[[pn]][, f144 := 1.0]
    posteriorsI[[pn]][, f186 := 1.0]
    posteriorsI[[pn]][, f228 := 1.0]
    posteriorsI[[pn]][, f270 := 1.0]
    posteriorsI[[pn]][, f312 := 1.0]
    posteriorsI[[pn]][, f354 := 1.0]
    posteriorsI[[pn]][, f396 := 1.0]
    posteriorsI[[pn]][, f438 := 1.0]
    posteriorsI[[pn]][, f480 := 1.0]
    posteriorsI[[pn]][, f522 := 1.0]
    posteriorsI[[pn]][, f564 := 1.0]
    posteriorsI[[pn]][, f606 := 1.0]
    posteriorsI[[pn]][, f648 := 1.0]
    posteriorsI[[pn]][, f690 := 1.0]
    posteriorsI[[pn]][, f732 := 1.0]
    posteriorsI[[pn]][, f774 := 1.0]
    posteriorsI[[pn]][, f816 := 1.0]

    THE_ROW = posteriorsI[[pn]][, which.max(lp)]
    
    # params_json = jsonlite::toJSON(paramsI, pretty = TRUE)
    # writeLines(params_json, "./pfilter_pop0_params.txt")
    # stop()
    
    print(system.time(
        w <- cm_backend_pfilter(cm_translate_parameters(paramsI), posteriorsI[[pn]],
            row = THE_ROW, seed = 0, n_particles = 500*PT_MULTIPLIER, times = seq(110, PTLIMIT, by = 5), 
            adjust_nu = 3.0, adjust_scale = 0.05, adjust_b = 0.8, n_threads = N_THREADS)
    ))

    pcrp[[pn]] = copy(as.data.table(rbindlist(w[[1]])[, .(run, t, population, group, pcr_positive_p)]))
    # dynamics[[pn]] = copy(as.data.table(rbindlist(w[[1]])))
    pfilter[[pn]] = copy(as.data.table(w[[2]]))
    posteriors[[pn]] = copy(posteriorsI[[pn]][THE_ROW])
    parameters[[pn]] = rlang::duplicate(paramsI)
    
    # Notify Rosie
    notify_command = paste0(
        'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
        '--form-string "user=ubdi7mpz6bfiy3a485qkgpt5kzk9mb" ',
        '--form-string "message=Region ', pn, ' done particle filtering." https://api.pushover.net/1/messages.json')
    system(notify_command)
    
}

# plot = ggplot(pcrp[[1]][, sum(pcr_positive_p), by = .(run, t)]) + 
#     geom_line(aes(t, V1, group = run, colour = run)) + 
#     geom_vline(aes(xintercept = PTLIMIT))
# ggsave("./pfilter_pcr.png", plot, width = 10, height = 5)

# stop()

notify_command = paste0(
    'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
    '--form-string "user=uqmiw2w9jho4hpmc2uqsfdk7hwyccb" ',
    '--form-string "message=Commencing trace making." https://api.pushover.net/1/messages.json')
system(notify_command)

# Make traces
traces = list()
NTRACES = 250
arima_method = "detrend"
detrend_window = 35 # higher values smooth more

set.seed(12345)
for (pn in which_pops)
{
    # Get best-fitting trace
    ptrace = pfilter[[pn]][t == max(t)][which.min(ll_sum), particle]
    for (tt in pfilter[[pn]][, rev(unique(t)[-1])]) {
        ptrace = c(ptrace, pfilter[[pn]][t == tt & particle == tail(ptrace, 1), parent])
    }

    adj_match = data.table(t = pfilter[[pn]][, unique(t)], particle = rev(ptrace));
    adj_match = adj_match[t <= PTLIMIT];
    adj_trace = merge(adj_match, pfilter[[pn]], by = c("t", "particle"));
    adj_trace[, log_adj := log(adj)]
    
    if (arima_method == "detrend") {
        # Get trend
        adj_trace[, log_adj_trend := zoo::rollapply(log_adj, detrend_window, mean, fill = "extend", partial = TRUE)];
        adj_trace[, log_adj_detrend := log_adj - log_adj_trend];
        final_trend = adj_trace[which.max(t), log_adj_trend];
        
        # Get logged time series
        myts = ts(adj_trace[, log_adj_detrend], start = adj_trace[, min(t)], end = adj_trace[, max(t)], deltat = adj_trace[, unique(diff(t))])

        # Fit ARIMA
        myarima = Arima(myts, c(2, 0, 0))
    
        # Simulate
        new.t = seq(from = tsp(myts)[1], by = 1 / tsp(myts)[3], length.out = 2 * ((tsp(myts)[2] - tsp(myts)[1]) * tsp(myts)[3] + 1))
        new.y = exp(c(adj_trace$log_adj, predict(myarima, n.ahead = length(new.t) - length(myts))$pred + final_trend))
        sims = data.table(i = 0, t = new.t, adj = new.y)
        for (i in 1:NTRACES)
        {
            new.y = exp(c(adj_trace$log_adj, simulate(myarima) + final_trend))
            sims = rbind(sims, data.table(i = i, t = new.t, adj = new.y))
        }
    } else {
        # Previous method
        # Get logged time series
        myts = ts(adj_trace[, log_adj], start = adj_trace[, min(t)], end = adj_trace[, max(t)], deltat = adj_trace[, unique(diff(t))])
    
        # Fit ARIMA
        myarima = Arima(myts, c(2, 0, 0))
    
        # Simulate
        new.t = seq(from = tsp(myts)[1], by = 1 / tsp(myts)[3], length.out = 2 * ((tsp(myts)[2] - tsp(myts)[1]) * tsp(myts)[3] + 1))
        new.y = exp(c(myts, predict(myarima, n.ahead = length(new.t) - length(myts))$pred))
        sims = data.table(i = 0, t = new.t, adj = new.y)
        for (i in 1:NTRACES)
        {
            new.y = exp(c(myts, simulate(myarima)))
            sims = rbind(sims, data.table(i = i, t = new.t, adj = new.y))
        }
    }
    
    # ggplot(sims[i==0]) + geom_line(aes(x = t, y = adj, colour = i, group = i))
    
    traces[[pn]] = sims
    
    # Notify Rosie
    notify_command = paste0(
        'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
        '--form-string "user=ubdi7mpz6bfiy3a485qkgpt5kzk9mb" ',
        '--form-string "message=Region ', pn, ' traces finished." https://api.pushover.net/1/messages.json')
    system(notify_command)
}

# Save traces and parameters and posteriors
qsave(list(traces = traces, posteriors = posteriors, parameters = parameters), PF_FILE)

notify_command = paste0(
    'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
    '--form-string "user=uqmiw2w9jho4hpmc2uqsfdk7hwyccb" ',
    '--form-string "message=Traces, parameters and posteriors are saved." https://api.pushover.net/1/messages.json')
system(notify_command)

# Some output
traces2 = rbindlist(traces, idcol = "population")
traces2[, population := nhs_regions[population]]
plot1 = ggplot(traces2) +
    geom_line(aes(t, adj, group = i, colour = ifelse(t < PTLIMIT, 0, i))) +
    facet_wrap(~population) +
    theme(legend.position = "none") + 
    # scale_y_log10() + 
    geom_vline(aes(xintercept=550))
plot1
plot2 = ggplot(traces2[population == "London" & i < 10]) +
    geom_line(aes(t, adj, group = i, colour = ifelse(t < PTLIMIT, 0, i))) +
    theme(legend.position = "none") + 
    facet_wrap(~i) 
# +
#     scale_y_log10()
plot2
ggsave(paste0("./output/", PF_FILE_TO_SAVE, "_pfilter_1_xpro.png"), plot1, width = 20, height = 20, units = "cm")
ggsave(paste0("./output/", PF_FILE_TO_SAVE, "_pfilter_2_xpro.png"), plot2, width = 20, height = 20, units = "cm")

notify_command = paste0(
    'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
    '--form-string "user=uqmiw2w9jho4hpmc2uqsfdk7hwyccb" ',
    '--form-string "message=Output plots are saved. pfilter.R script finished." https://api.pushover.net/1/messages.json')
system(notify_command)


if (0){
    # plot output from >1 particle filter file
    o1 = qread('./fits/pf_relu_yeswane_sev2.0_22012204.qs')
    o2 = qread('./fits/pf_relu_yeswane_sev2.0_22012904_2.qs')
    o3 = qread('./fits/pf_relu_yeswane_sev2.0_22012904_3.qs')
    o4 = qread('./fits/pf_relu_yeswane_sev2.0_22012904_4.qs')
    
    t1 = o1$traces
    t2 = o2$traces
    t3 = o3$traces
    t4 = o4$traces
    
    ts1 = rbindlist(t1, idcol = "population")
    ts2 = rbindlist(t2, idcol = "population")
    ts3 = rbindlist(t3, idcol = "population")
    ts4 = rbindlist(t4, idcol = "population")
    
    ts1[, population := nhs_regions[population]]
    ts2[, population := nhs_regions[population]]
    ts3[, population := nhs_regions[population]]
    ts4[, population := nhs_regions[population]]
    
    ts1$PTLIM = 720
    ts2$PTLIM = 750
    ts3$PTLIM = 700
    ts4$PTLIM = 750
    
    ts1$NPTCL = 500
    ts2$NPTCL = 500
    ts3$NPTCL = 500
    ts4$NPTCL = 500*4
    
    allts = rbindlist(list(ts1, ts2, ts3, ts4), idcol = "fit")
    
    # plot two files with PTLIM = 750 with different numbers of particle
    ptcmpr = ggplot(allts[PTLIM == 750,]) +
        geom_line(aes(t, adj, group = i, colour = ifelse(t < PTLIM, 0, i))) +
        facet_grid(NPTCL~population) +
        theme(legend.position = "none") + 
        scale_y_log10() + geom_vline(aes(xintercept=PTLIM), linetype = "dotted")
    ptcmpr
    
    # same as above but for two regions only
    ptcmpr2 = ggplot(allts[PTLIM == 750 & population %in% c('London', 'North West')]) +
        geom_line(aes(t, adj, group = i, colour = ifelse(t < PTLIM, 0, i))) +
        facet_grid(NPTCL~population) +
        theme(legend.position = "none") + 
        scale_y_log10() + geom_vline(aes(xintercept=PTLIM), linetype = "dotted")
    ptcmpr2
    
    # all regions for all files
    allr = ggplot(allts) +
        geom_line(aes(t, adj, group = i, colour = ifelse(t < PTLIM, 0, i))) +
        facet_grid(fit~population) +
        theme(legend.position = "none") + 
        scale_y_log10() + geom_vline(aes(xintercept=PTLIM), linetype = "dotted")
    allr
    # two regions for all files
    twor = ggplot(allts[population %in% c('London', 'North West')]) +
        geom_line(aes(t, adj, group = i, colour = ifelse(t < PTLIM, 0, i))) +
        facet_grid(PTLIM~population) +
        theme(legend.position = "none") + 
        scale_y_log10() + geom_vline(aes(xintercept=PTLIM), linetype = "dotted")
    twor
    
    twor + scale_x_continuous(limits = c(600, 800))
}
