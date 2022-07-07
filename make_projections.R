# generate output for paper based on existing particle filtered fit

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

theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

######################## CHANGE HERE FOR EACH RUN ########################

PF_FILE = "./fits/pf_relu_yeswane_sev2.0_22050605_20220511163014.qs"

data_file = "processed-data-20220506122858.qs"
mobility_file = "schedule3-MTPs-20220506121302.rds"
date_fitting = "2022-05-06"
vax_file = "vax-covidm20220505205235.rds"
sgtf_stop = "2021-02-15"

# these arguments are not used as we specify adj_file = NA
adj_dstart = "2021-11-01"
adj_dend = "2021-11-30"

# SERO_CUT_OFF determines the date at which sero-prevalence data is no longer
# used to fit to (we use the start date of estimates for cut off)
SERO_CUT_OFF = '2020-12-01'

UPDATE_PARAMS = TRUE

# probability of receiving a booster dose by covidm 5-year age groups (this
# should match the PBOOST assumptions used in fit-omi.R)
PBOOST=c(0,   0,     0,   0.4, 0.544, 0.597, 0.643, 0.698, 0.757, 0.809, 0.861, 0.892, 0.917, 0.946, 0.964,  0.967)

BFOLD = 2.5 # fold increase in neutralisation titres after boosters (2.5 / 4.9)

DFR = 1.0  # Delta fold reduction (e.g. 1.0 or 3.9)

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

############################## END CHANGES ##############################

# Command line
c_args = list("relu", "yeswane", "central", "2.0", "seaslate", "0.1")

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

# these have to go after covidm.R has been sourced
# PBOOST = c(rep(PBY, BAGE %/% 5), rep(PBO, 16 - BAGE %/% 5)) ## Added this
# PBOOST = PBOOST * cm_age_coefficients(BUP, 80, seq(0, 80, by = 5)) ## And this

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

#
# SETUP
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

priorsI = c(priorsI, extra_priors)

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

priorsI = priorsI[setdiff(names(priorsI), names(constants))]
N_THREADS = length(priorsI) * 2;

# Load particle filter fit from saved file
saved = qread(PF_FILE)
traces = saved$traces
posteriors = saved$posteriors
parameters = saved$parameters

# Function to run projections
project = function(
    which_pops = c(1, 3, 4, 5, 6, 9, 10),
    n_traces = 100,
    mobility_file = "./fitting_data/schedule3-MTPs-20220506121302.rds",
    wane_yn = "yeswane",
    vac_eff = "central",
    wane_scen = "central",
    bup = 16,
    boost_age_split = 50,
    p_boost_young = 0.85,
    p_boost_old = 0.95,
    seasonality_aw = 0.1,
    gen_reduce = 0,
    gen_start = NA,
    gen_end = NA,
    gen2_reduce = 0,
    gen2_start = NA,
    gen2_end = NA,
    cert = 0,
    mask = 0,
    wfh_reference_date = NA,
    planB_start = "2022-04-01",
    planB_end = "2022-06-01",
    lockdown_reference_date = NA,
    lockdown_start = "2022-04-01",
    lockdown_end = "2022-06-01",
    n_threads = N_THREADS,
    VAX_SCHEDULE = vax_file,
    bst_rate1 = 229000,
    bst_rate2 = 1000000,
    bst_fold = 1,
    bst_om_fold = 1,
    bst_duration = 180,
    delta_fold_escape = 1,
    VE_debug_file = FALSE,
    p_boost_for_obs = NULL
    )
{

    if (VE_debug_file != FALSE) {
        n_threads = 1
        n_traces = 0
        VE_debug_file = paste0(VE_debug_file, "_", str_replace_all(Sys.time(), "[ :GMTBST-]", ""), "_")
    }

    # VACCINE SCHEDULE
    vacc = readRDS(paste0('./fitting_data/', VAX_SCHEDULE))

    # Load particle filter from saved
    # saved = qread(paste0('./fits/', pf_file));
    # traces = saved$traces;
    # posteriors = saved$posteriors;
    # parameters = saved$parameters;

    planB_t0 = as.numeric(ymd(planB_start) - ymd("2020-01-01"));
    planB_t1 = as.numeric(ymd(planB_end) - ymd("2020-01-01"));
    lockdown_t0 = as.numeric(ymd(lockdown_start) - ymd("2020-01-01"));
    lockdown_t1 = as.numeric(ymd(lockdown_end) - ymd("2020-01-01"));

    gen_reduce_t0  = as.numeric(ymd(gen_start)  - ymd("2020-01-01"));
    gen_reduce_t1  = as.numeric(ymd(gen_end)    - ymd("2020-01-01"));
    gen2_reduce_t0 = as.numeric(ymd(gen2_start) - ymd("2020-01-01"));
    gen2_reduce_t1 = as.numeric(ymd(gen2_end)   - ymd("2020-01-01"));

    if (p_boost_for_obs == "NULL"){
        p_boost_for_obs = c(rep(p_boost_young, boost_age_split %/% 5), rep(p_boost_old, 16 - boost_age_split %/% 5)) ## Added this
        p_boost_for_obs = p_boost_for_obs * cm_age_coefficients(bup, 80, seq(0, 80, by = 5)) ## And this
    } else {
        p_boost_for_obs = PBOOST
    }

    # MOBILITY SCENARIO
    sched = readRDS(mobility_file)

    dyns = list()

    for (pn in which_pops) {
        # Get parameters
        paramsI = rlang::duplicate(parameters[[pn]])
        paramsI$time1 = "2022-12-31"

        # SET WANING / VE
        paramsI = reinfection_defaults(paramsI, 1)
        paramsI = waning_scenario(wane_yn, paramsI, 1)
        paramsI = VE_scenario    (vac_eff, PBOOST, paramsI, 1, delta_fold_reduction = delta_fold_escape)

        if (wane_scen != "central"){
            paramsI = extra_waning_scenario(wane_scen, paramsI, 1)
        }

        # SET MOBILITY

        # Set mobility scenario for no tier
        paramsI$schedule[[1]] = rlang::duplicate(sched[[pn]])
        paramsI$schedule[[1]]$pops = 0

        if (!is.na(lockdown_reference_date)) {
            # Schedule 1, indices 0 and 1
            s1i0 = which(paramsI$schedule[[1]]$times > lockdown_t0)[1] - 1;
            s1i1 = which(paramsI$schedule[[1]]$times >= lockdown_t1)[1];
            # Lockdown lookback time
            ld_lt = as.numeric(ymd(lockdown_reference_date) - ymd("2020-01-01"));
            # Lockdown lookback indices
            lb0 = which(paramsI$schedule[[1]]$times > ld_lt)[1] - 1;
            lb1 = which(paramsI$schedule[[1]]$times >= ld_lt)[1] - 1;
            # Contact rates
            ld_val1 = mean(sapply(paramsI$schedule[[1]]$values[lb0:lb1], `[`, 1)); # wplc
            ld_val2 = mean(sapply(paramsI$schedule[[1]]$values[lb0:lb1], `[`, 2)); # groc
            ld_val3 = mean(sapply(paramsI$schedule[[1]]$values[lb0:lb1], `[`, 3)); # rtrc
            ld_val4 = mean(sapply(paramsI$schedule[[1]]$values[lb0:lb1], `[`, 4)); # trns
            ld_val5 = mean(sapply(paramsI$schedule[[1]]$values[lb0:lb1], `[`, 5)); # scho
            # Overwrite contact rates
            for (s1i in s1i0:s1i1) {
                paramsI$schedule[[1]]$values[[s1i]][1] = ld_val1;
                paramsI$schedule[[1]]$values[[s1i]][2] = ld_val2;
                paramsI$schedule[[1]]$values[[s1i]][3] = ld_val3;
                paramsI$schedule[[1]]$values[[s1i]][4] = ld_val4;
                paramsI$schedule[[1]]$values[[s1i]][5] = ld_val5;
            }
        }
        if (!is.na(wfh_reference_date)) {
            s1i0 = which(paramsI$schedule[[1]]$times > planB_t0)[1] - 1;
            s1i1 = which(paramsI$schedule[[1]]$times >= planB_t1)[1];
            wfh_lt = as.numeric(ymd(wfh_reference_date) - ymd("2020-01-01"));
            lb0 = which(paramsI$schedule[[1]]$times > wfh_lt)[1] - 1;
            lb1 = which(paramsI$schedule[[1]]$times >= wfh_lt)[1] - 1;
            wfh_val = mean(sapply(paramsI$schedule[[1]]$values[lb0:lb1], `[`, 2));
            for (s1i in s1i0:s1i1) {
                paramsI$schedule[[1]]$values[[s1i]][2] = wfh_val;
            }

        }

        # ... for tier 2
        paramsI$schedule[[2]] = rlang::duplicate(paramsI$schedule[[1]]);
        for (i in seq_along(paramsI$schedule[[2]]$values)) {
            paramsI$schedule[[2]]$values[[i]][1] = paramsI$schedule[[1]]$values[[i]][1] +  0.2497655 / 100;
            paramsI$schedule[[2]]$values[[i]][2] = paramsI$schedule[[1]]$values[[i]][2] + -0.2307939 / 100;
            paramsI$schedule[[2]]$values[[i]][3] = paramsI$schedule[[1]]$values[[i]][3] + -1.5907698 / 100;
            paramsI$schedule[[2]]$values[[i]][4] = paramsI$schedule[[1]]$values[[i]][4] + -3.4866544 / 100;
            paramsI$schedule[[2]]$values[[i]][5] = paramsI$schedule[[1]]$values[[i]][5] + -3.4524518 / 100;
        }
        paramsI$schedule[[2]]$mode = "bypass";

        # ... for tier 3
        paramsI$schedule[[3]] = rlang::duplicate(paramsI$schedule[[1]]);
        for (i in seq_along(paramsI$schedule[[3]]$values)) {
            paramsI$schedule[[3]]$values[[i]][1] = paramsI$schedule[[1]]$values[[i]][1] +  2.080457 / 100;
            paramsI$schedule[[3]]$values[[i]][2] = paramsI$schedule[[1]]$values[[i]][2] + -8.045226 / 100;
            paramsI$schedule[[3]]$values[[i]][3] = paramsI$schedule[[1]]$values[[i]][3] + -2.476266 / 100;
            paramsI$schedule[[3]]$values[[i]][4] = paramsI$schedule[[1]]$values[[i]][4] + -10.144043 / 100;
            paramsI$schedule[[3]]$values[[i]][5] = paramsI$schedule[[1]]$values[[i]][5] + -7.681244 / 100;
        }
        paramsI$schedule[[3]]$mode = "bypass";

        # Load user defined functions
        cat(paste0(nhs_regions[pn], "...\n"))
        cm_source_backend(
            user_defined = list(
                model_v3 = list(
                    cpp_changes = cpp_chgI_voc(priorsI, constants, seasonality = opt_seas,
                                               v2 = opt_v2, v2_relu = opt_relu, v2_latdur = opt_latdur,
                                               v2_serial = opt_serial, v2_infdur = opt_infdur,
                                               v2_immesc = opt_immesc, v2_ch_u = opt_ch_u,
                                               v3_relu = NULL, v3_severity = V3_SEVERITY, v3 = opt_v3),
                    cpp_loglikelihood = "",
                    cpp_observer = c(
                        cpp_obsI_voc(concentration = opt_conc, v2 = opt_v2,
                                     P.death, P.critical, priorsI, constants, v3_severity = V3_SEVERITY),
                        cpp_obsI_voc_omi(OMI_SETUP_T, OMI_X_PROTECTION, OMI_VAX_FACTOR, OMI_VAX_ASSUMPTION, OMI_SEEDS_PER_DAY,
                                         omi_sev = OMI_SEV, omi_crit = OMI_CRIT),
                        cpp_obsI_vax(paramsI, vacc[[pn]]),
                        if (SEAS_YN == "seasyes") cpp_obsI_seasonality(seasonality_aw, 1) else if (SEAS_YN == "seaslate") cpp_obsI_seasonality(seasonality_aw, 456) else "",
                        cpp_obsI_booster(
                            target_phase1 = bst_rate1,
                            target_phase2 = bst_rate2,
                            proportion_booster = p_boost_for_obs,
                            booster_fold = bst_fold,
                            booster_om_fold = bst_om_fold,
                            booster_duration = bst_duration
                        ),
                        cpp_obsI_voc_ba2(ba2t, ba2f),
                        # need the following if we want plan B measures??
                        # cpp_obsI_aw(seasonality_aw, cert, planB_start, planB_end),
                        if (VE_debug_file != FALSE) cpp_obsI_printVE(VE_debug_file, pn, t_step = 14) else ""
                    )
                )
            ),
            verbose = FALSE
        )
        cat(paste0(nhs_regions[pn], " sourced. Starting traces.\n"))
        # Do traces
        all_parameters = list()
        for (TRACE in 0:n_traces) {
            # NS = length(all_parameters[[TRACE + 1]]$schedule) + 1;
            all_parameters[[TRACE + 1]] = rlang::duplicate(paramsI);
            all_parameters[[TRACE + 1]]$schedule[[7]] = list(
                parameter = 'contact',
                pops = 0,
                mode = 'multiply',
                values = asplit(matrix(rep(traces[[pn]][i == TRACE, adj], 8), ncol = 8), 1),
                times = traces[[pn]][i == TRACE, t - 5]
            )

            if (mask > 0) {
                s7i0 = which(all_parameters[[TRACE + 1]]$schedule[[7]]$times > planB_t0)[1] - 1;
                s7i1 = which(all_parameters[[TRACE + 1]]$schedule[[7]]$times >= planB_t1)[1];
                for (s7i in s7i0:s7i1) {
                    all_parameters[[TRACE + 1]]$schedule[[7]]$values[[s7i]] =
                        all_parameters[[TRACE + 1]]$schedule[[7]]$values[[s7i]] * (1.0 - mask);
                }
            }

            if (gen_reduce > 0) {
                s7i0 = which(all_parameters[[TRACE + 1]]$schedule[[7]]$times > gen_reduce_t0)[1] - 1;
                s7i1 = which(all_parameters[[TRACE + 1]]$schedule[[7]]$times >= gen_reduce_t1)[1];
                for (s7i in s7i0:s7i1) {
                    all_parameters[[TRACE + 1]]$schedule[[7]]$values[[s7i]] =
                        all_parameters[[TRACE + 1]]$schedule[[7]]$values[[s7i]] * (1.0 - gen_reduce);
                }
            }

            if (gen2_reduce > 0) {
                s7i0 = which(all_parameters[[TRACE + 1]]$schedule[[7]]$times > gen2_reduce_t0)[1] - 1;
                s7i1 = which(all_parameters[[TRACE + 1]]$schedule[[7]]$times >= gen2_reduce_t1)[1];
                for (s7i in s7i0:s7i1) {
                    all_parameters[[TRACE + 1]]$schedule[[7]]$values[[s7i]] =
                        all_parameters[[TRACE + 1]]$schedule[[7]]$values[[s7i]] * (1.0 - gen2_reduce);
                }
            }

            all_parameters[[TRACE + 1]] = cm_translate_parameters(all_parameters[[TRACE + 1]]);
        }

        params_json = jsonlite::toJSON(all_parameters[[1]], pretty = TRUE)
        writeLines(params_json, "./results_pop0_params.txt")
        # stop()

        dyn = cm_backend_sample_fit_2(all_parameters, posteriors[[pn]],
                                      n = n_traces + 1, seed = 0, n_threads = n_threads);
        dyn = rbindlist(dyn);
        dyn[, population := pn];
        dyn[, NHS.region := nhs_regions[population]];
        dyns[[length(dyns) + 1]] = dyn;
    }

    return (rbindlist(dyns))
}

source("./spim_output.R")

do_output = function(w2, titl, filename, adj_file){
    w = copy(w2)
    w[, population := NULL]
    w[, population := NHS.region]

    # Add unvaccinated outputs
    w[, deaths_V0 := deaths - deaths_V1 - deaths_V2 - deaths_V3]
    w[, hosp_adm_V0 := hosp_adm - hosp_adm_V1 - hosp_adm_V2 - hosp_adm_V3]

    # Adjust for variable infection burden ratios
    w = apply_gamma_multiplier(w, "deaths", "disp_deaths", ld, "N", 0.1, 7, GAMDISP, c("deaths_V0", "deaths_V1", "deaths_V2", "deaths_V3"), adj_file = NA, adj_dstart, adj_dend)
    w = apply_gamma_multiplier(w, "hosp_adm", "disp_hosp_inc", sitreps, "n_admitted_diagnosed", 0.1, 7, GAMDISP, c("hosp_adm_V0", "hosp_adm_V1", "hosp_adm_V2", "hosp_adm_V3", "hosp_undetected_o"), adj_file = NA, adj_dstart, adj_dend)
    w[, known_hosp_beds := hosp_bed - hosp_undetected_p]
    w = apply_gamma_multiplier(w, "known_hosp_beds", "disp_hosp_prev", sitreps, "n_in_all_beds", 0.1, 7, GAMDISP, c("hosp_bed", "hosp_undetected_p"), adj_file = NA, adj_dstart, adj_dend)
    w[, known_hosp_beds := NULL]
    w = apply_gamma_multiplier(w, "icu_bed", "disp_icu_prev", sitreps, "n_in_itu", 0.1, 7, GAMDISP, "icu_adm", adj_file = NA, adj_dstart, adj_dend)

    # Add all-England population
    england = w[, lapply(.SD, sum), .SDcols = S:hosp_undetected_o, by = .(run, t, group)]
    england[, population := "England"]
    england[, NHS.region := "England"]
    w0 = rbind(w, england, use.names = TRUE, fill = TRUE)
    w0[, NHS.region := factor(NHS.region, levels = c("England", nhs_regions[c(1, 3, 4, 5, 6, 9, 10)]))]

    w0 = w0[, .(run, NHS.region, t, group, pcr_positive_i, pcr_positive_p, test_o, hosp_adm, hosp_bed, icu_bed, deaths)]
    w0 = w0[, lapply(.SD, sum), .SDcols = pcr_positive_i:deaths, by = .(run, NHS.region, t)]
    w0 = w0[NHS.region == "England"]

    w0 = melt(w0, id.vars = 1:3)
    w0[variable == "pcr_positive_p", value := 100 * value / 56550138]
    w0[, cumulative_value := cumsum(value), by = .(NHS.region, run, variable)]

    summ = w0[, .(
        q05 = quantile(value[-1], 0.05),
        q25 = quantile(value[-1], 0.25),
        q50 = quantile(value[-1], 0.5),
        q75 = quantile(value[-1], 0.75),
        q95 = quantile(value[-1], 0.95),
        va = value[1]), by = .(NHS.region, t, variable)]
    summ[, date := ymd("2020-01-01") + t]

    cumul_summ = w0[, .(
        q05 = quantile(cumulative_value[-1], 0.05),
        q25 = quantile(cumulative_value[-1], 0.25),
        q50 = quantile(cumulative_value[-1], 0.5),
        q75 = quantile(cumulative_value[-1], 0.75),
        q95 = quantile(cumulative_value[-1], 0.95),
        va = cumulative_value[1]), by = .(NHS.region, t, variable)]
    cumul_summ[, date := ymd("2020-01-01") + t]

    tots = rbind(
        w0[,                                      .(when = "all", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-01-01", .(when = "from_jan_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-01-01" & ymd("2020-01-01") + t <= "2022-04-30",
           .(when = "from_jan_22_to_apr_22", tot = sum(value)), by = .(variable, run)]
    )

    tots = tots[, .(
        median = median(tot),
        mean = mean(tot),
        q05 = quantile(tot, 0.05),
        q25 = quantile(tot, 0.25),
        q75 = quantile(tot, 0.75),
        q95 = quantile(tot, 0.95)
    ), by = .(variable, when)]

    return (list(w, w0, summ, tots, cumul_summ))
}

burdens = function(w, adj_file, adj_ts)
{
    w[, population := NULL]
    w[, population := NHS.region]

    # Add unvaccinated outputs
    w[, deaths_V0 := deaths - deaths_V1 - deaths_V2 - deaths_V3]
    w[, hosp_adm_V0 := hosp_adm - hosp_adm_V1 - hosp_adm_V2 - hosp_adm_V3]

    # Adjust for variable infection burden ratios
    w = apply_gamma_multiplier(w, "deaths", "disp_deaths", ld, "N", 0.1, 7, GAMDISP, c("deaths_V0", "deaths_V1", "deaths_V2", "deaths_V3"), adj_file, adj_dstart, adj_dend, adj_ts, '2020-01-01', '2022-09-30')
    w = apply_gamma_multiplier(w, "hosp_adm", "disp_hosp_inc", sitreps, "n_admitted_diagnosed", 0.1, 7, GAMDISP, c("hosp_adm_V0", "hosp_adm_V1", "hosp_adm_V2", "hosp_adm_V3", "hosp_undetected_o"), adj_file, adj_dstart, adj_dend, adj_ts, '2020-01-01', '2022-09-30')
    w[, known_hosp_beds := hosp_bed - hosp_undetected_p]
    w = apply_gamma_multiplier(w, "known_hosp_beds", "disp_hosp_prev", sitreps, "n_in_all_beds", 0.1, 7, GAMDISP, c("hosp_bed", "hosp_undetected_p"), adj_file, adj_dstart, adj_dend, adj_ts, '2020-01-01', '2022-09-30')
    w[, known_hosp_beds := NULL]
    w = apply_gamma_multiplier(w, "icu_bed", "disp_icu_prev", sitreps, "n_in_itu", 0.1, 7, GAMDISP, "icu_adm", adj_file, adj_dstart, adj_dend, adj_ts, '2020-01-01', '2022-09-30')

    # Add all-England population
    england = w[, lapply(.SD, sum), .SDcols = S:hosp_undetected_o, by = .(run, t, group)]
    england[, population := "England"]
    england[, NHS.region := "England"]
    w0 = rbind(w, england, use.names = TRUE, fill = TRUE)
    w0[, NHS.region := factor(NHS.region, levels = c("England", nhs_regions[c(1, 3, 4, 5, 6, 9, 10)]))]

    w0 = w0[, .(run, NHS.region, t, group, pcr_positive_i, pcr_positive_p, test_o, hosp_adm, hosp_bed, icu_bed, deaths)]
    w0 = w0[, lapply(.SD, sum), .SDcols = pcr_positive_i:deaths, by = .(run, NHS.region, t)]
    w0 = w0[NHS.region == "England"]

    w0 = melt(w0, id.vars = 1:3)
    w0[variable == "pcr_positive_p", value := 100 * value / 56550138]
    w0[, cumulative_value := cumsum(value), by = .(NHS.region, run, variable)];

    summ = w0[, .(
        q05 = quantile(value[-1], 0.05),
        q25 = quantile(value[-1], 0.25),
        q50 = quantile(value[-1], 0.5),
        q75 = quantile(value[-1], 0.75),
        q95 = quantile(value[-1], 0.95),
        va = value[1]), by = .(NHS.region, t, variable)]
    summ[, date := ymd("2020-01-01") + t]

    cumul_summ = w0[, .(
        q05 = quantile(cumulative_value[-1], 0.05),
        q25 = quantile(cumulative_value[-1], 0.25),
        q50 = quantile(cumulative_value[-1], 0.5),
        q75 = quantile(cumulative_value[-1], 0.75),
        q95 = quantile(cumulative_value[-1], 0.95),
        va = cumulative_value[1]), by = .(NHS.region, t, variable)]
    cumul_summ[, date := ymd("2020-01-01") + t]

    tots = rbind(
        w0[,                                      .(when = "all", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2021-10-01", .(when = "from_oct_21", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2021-10-01" & ymd("2020-01-01") + t <= "2021-12-31",
           .(when = "from_oct_21_to_dec_21", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2021-12-01" & ymd("2020-01-01") + t <= "2022-04-30",
           .(when = "from_dec_21_to_apr_22", tot = sum(value)), by = .(variable, run)],
        w0[ymd("2020-01-01") + t >= "2022-01-01", .(when = "from_jan_22", tot = sum(value)), by = .(variable, run)]
    )

    tots = tots[, .(
        median = median(tot),
        mean = mean(tot),
        q05 = quantile(tot, 0.05),
        q25 = quantile(tot, 0.25),
        q75 = quantile(tot, 0.75),
        q95 = quantile(tot, 0.95)
    ), by = .(variable, when)]

    return (list(w, w0, summ, tots, cumul_summ))
}

process_burdens = function(basename, adj_file, adj_ts, only_region = NULL){

    if (is.null(only_region)) {
        w = burdens(qread(paste0("./output/paper/may22/output1_", basename, ".qs")), adj_file, adj_ts)
    } else {
        w = burdens(qread(paste0("./output/paper/may22/output1_", basename, ".qs"))[NHS.region == only_region], adj_file, adj_ts)
        basename = paste0(basename, "_", only_region)

        # correct NHS.region
        w[[2]][, NHS.region := only_region]
        w[[3]][, NHS.region := only_region]
        w[[5]][, NHS.region := only_region]
    }
    qsave(w[[1]], paste0("./output/paper/may22/adjusted1_", basename, ".qs"))
    qsave(w[[2]], paste0("./output/paper/may22/burdens_", basename, ".qs"))
    fwrite(w[[3]], paste0("./output/paper/may22/spim_", basename, ".csv"))
    fwrite(w[[4]], paste0("./output/paper/may22/totals_", basename, ".csv"))
    fwrite(w[[5]], paste0("./output/paper/may22/cumul_", basename, ".csv"))
}

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

# second set of scenarios for resubmission following Omicron updates (11th May 2022)
scenario_sheet = fread(
    "codenm,          mobility,                    wane_yn, vac_eff, wane_scen, bdur, bage, pby,  pbo,  pbst_all, seas, cert, mask,  wfh_date,   title,                                                                           vax_file,                                      adj_ts
boostHI_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.90, 0.98, NULL,     0.1,  0,    0,     NA,         High booster uptake: 90% in under 50s and 98% in 50+,                            vax-covidm20220505205235.rds,                  ./fits/final_adj_tseries_basecase_220511.csv
boostNO_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,    0.0,  0.0, NULL,     0.1,  0,    0,     NA,         No booster uptake,                                                               vax-covidm20220505205235.rds,                  ./fits/final_adj_tseries_basecase_220511.csv
boost50_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,    0.0, 0.95, NULL,     0.1,  0,    0,     NA,         Boosters for 50+ only: 95% uptake,                                               vax-covidm20220505205235.rds,                  ./fits/final_adj_tseries_basecase_220511.csv
")

# test basecase scenario for resubmission following BA.2 updates (May 2022)
# scenario_sheet = fread(
#     "codenm,          mobility,                    wane_yn, vac_eff, wane_scen, bdur, bage, pby,  pbo,  pbst_all, seas, cert, mask,  wfh_date,   title,                                                                           vax_file,                                      adj_ts
# basecase_220510,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Base case: 6-month relaxation; central waning; actual boosters; 20% seasonality, vax-covidm20220505205235.rds,                  NA
# ")

# list of scenarios to run for resubmission following BA.2 updates (April 2022)
# scenario_sheet = fread(
#     "codenm,          mobility,                    wane_yn, vac_eff, wane_scen, bdur, bage, pby,  pbo,  pbst_all, seas, cert, mask,  wfh_date,   title,                      vax_file,                      adj_ts
# test_220504,    schedule3-MTPs-20220429122830,     yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         TEST: 4th May 2022,      vax-covidm20220330102111.rds,  NA
# ")

# list of scenarios to run for resubmission following BA.2 updates (5th May 2022)
# scenario_sheet = fread(
#     "codenm,          mobility,                    wane_yn, vac_eff, wane_scen, bdur, bage, pby,  pbo,  pbst_all, seas, cert, mask,  wfh_date,   title,                                                                           vax_file,                                      adj_ts
# basecase_220505,      paper_aw_6mo_20220504164209, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Base case: 6-month relaxation; central waning; actual boosters; 20% seasonality, vax-covidm20220330102111.rds,                  NA
# mob_3wk_220505,       paper_aw_3wk_20220504164209, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         3-week relaxation,                                                               vax-covidm20220330102111.rds,                  NA
# mob_3mo_220505,       paper_aw_3mo_20220504164209, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         3-month relaxation,                                                              vax-covidm20220330102111.rds,                  NA
# mob_flt_220505,     schedule3-MTPs-20220429122830, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         No relaxation,                                                                   vax-covidm20220330102111.rds,                  NA
# 
# ZboostHI_220505,       paper_aw_6mo_20220504164209, yeswane, central, central,   180,  50,   0.90, 0.96, NULL,     0.1,  0,    0,     NA,         High booster uptake: 90% in under 50s and 96% in 50+,                            vax-covidm20220330102111.rds,                  ADD ADJUSTED FILE HERE
# ZboostNO_220505,
# 
# boostHI_220219,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.90, 0.96, NULL,     0.1,  0,    0,     NA,         High booster uptake: 90% in under 50s and 96% in 50+,                            vax-covidm20220117203900.rds,                  ./fits/final_adj_tseries_basecase_220214.csv
# boostNO_220219,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,    0.0,  0.0, NULL,     0.1,  0,    0,     NA,         No booster uptake,                                                               vax-covidm20220117203900.rds,                  ./fits/final_adj_tseries_basecase_220214.csv
# boost50_220219,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,    0.0, 0.96, NULL,     0.1,  0,    0,     NA,         Boosters for 50+ only: 96% uptake,                                               vax-covidm20220117203900.rds,                  ./fits/final_adj_tseries_basecase_220214.csv
# ZwaneHI_220214,       paper_aw_6mo_20220214160607, yeswane, central, hiwane,    180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         High waning,                                                                     vax-covidm20220117203900.rds,                  NA
# ZwaneVHI_220214,      paper_aw_6mo_20220214160607, yeswane, central, vhiwane,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Very high waning,                                                                vax-covidm20220117203900.rds,                  NA
# Zseas_10_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.05,  0,    0,     NA,         10% seasonality,                                                                 vax-covidm20220117203900.rds,                  NA
# Zseas_30_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.15,  0,    0,     NA,         30% seasonality,                                                                 vax-covidm20220117203900.rds,                  NA
# Zseas_40_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.20,  0,    0,     NA,         40% seasonality,                                                                 vax-covidm20220117203900.rds,                  NA
# Zvax1250_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Vax 12+ 50% uptake,                                                              vax-covidm2022011720390012plus_50percent.rds,  NA
# Zvax0580_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Vax 5+  80% uptake,                                                              vax-covidm202201172039005plus_80percent.rds,   NA
# Zvax0550_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Vax 5+  50% uptake,                                                              vax-covidm202201172039005plus_50percent.rds,   NA
# Zshrtbst_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,    90,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Short booster duration (90 days),                                                vax-covidm20220117203900.rds,                  NA
# Zlongbst_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   270,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Long booster duration (270 days),                                                vax-covidm20220117203900.rds,                  NA
# ")

# list of scenarios to run for resubmission
# scenario_sheet = fread(
#     "codenm,          mobility,                    wane_yn, vac_eff, wane_scen, bdur, bage, pby,  pbo,  pbst_all, seas, cert, mask,  wfh_date,   title,                                                                           vax_file,                                      adj_ts
# Zbasecase_220214,     paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Base case: 6-month relaxation; central waning; actual boosters; 20% seasonality, vax-covidm20220117203900.rds,                  NA
# Zmob_3mo_220214,       paper_aw_3mo_20220215184431, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         3-month relaxation,                                                              vax-covidm20220117203900.rds,                  NA
# Zmob_flt_220214,      schedule3-MTPs-2022-02-11,   yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         No relaxation,                                                                   vax-covidm20220117203900.rds,                  NA
# boostHI_220219,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.90, 0.96, NULL,     0.1,  0,    0,     NA,         High booster uptake: 90% in under 50s and 96% in 50+,                            vax-covidm20220117203900.rds,                  ./fits/final_adj_tseries_basecase_220214.csv
# boostNO_220219,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,    0.0,  0.0, NULL,     0.1,  0,    0,     NA,         No booster uptake,                                                               vax-covidm20220117203900.rds,                  ./fits/final_adj_tseries_basecase_220214.csv
# boost50_220219,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,    0.0, 0.96, NULL,     0.1,  0,    0,     NA,         Boosters for 50+ only: 96% uptake,                                               vax-covidm20220117203900.rds,                  ./fits/final_adj_tseries_basecase_220214.csv
# ZwaneHI_220214,       paper_aw_6mo_20220214160607, yeswane, central, hiwane,    180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         High waning,                                                                     vax-covidm20220117203900.rds,                  NA
# ZwaneVHI_220214,      paper_aw_6mo_20220214160607, yeswane, central, vhiwane,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Very high waning,                                                                vax-covidm20220117203900.rds,                  NA
# Zseas_10_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.05,  0,    0,     NA,         10% seasonality,                                                                 vax-covidm20220117203900.rds,                  NA
# Zseas_30_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.15,  0,    0,     NA,         30% seasonality,                                                                 vax-covidm20220117203900.rds,                  NA
# Zseas_40_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.20,  0,    0,     NA,         40% seasonality,                                                                 vax-covidm20220117203900.rds,                  NA
# Zvax1250_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Vax 12+ 50% uptake,                                                              vax-covidm2022011720390012plus_50percent.rds,  NA
# Zvax0580_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Vax 5+  80% uptake,                                                              vax-covidm202201172039005plus_80percent.rds,   NA
# Zvax0550_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Vax 5+  50% uptake,                                                              vax-covidm202201172039005plus_50percent.rds,   NA
# Zshrtbst_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,    90,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Short booster duration (90 days),                                                vax-covidm20220117203900.rds,                  NA
# Zlongbst_220214,      paper_aw_6mo_20220214160607, yeswane, central, central,   270,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Long booster duration (270 days),                                                vax-covidm20220117203900.rds,                  NA
# ")


# Notify Rosie
notify_command = paste0(
    'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
    '--form-string "user=ubdi7mpz6bfiy3a485qkgpt5kzk9mb" ',
    '--form-string "message=Starting projections." https://api.pushover.net/1/messages.json')
system(notify_command)

for (ROW in 1:nrow(scenario_sheet)) {

    if (scenario_sheet[ROW, codenm] %like% "^Z") {
        next;
    }

    cat(scenario_sheet[ROW, codenm], "\n");
    w = project(
        which_pops = which_pops,
        n_traces = 100,
        mobility_file = paste0("./fitting_data/", scenario_sheet[ROW, mobility], ".rds"),
        wane_yn = scenario_sheet[ROW, wane_yn],
        vac_eff = scenario_sheet[ROW, vac_eff],
        wane_scen = scenario_sheet[ROW, wane_scen],
        bup = 16,
        boost_age_split = scenario_sheet[ROW, bage],
        p_boost_young = scenario_sheet[ROW, pby],
        p_boost_old   = scenario_sheet[ROW, pbo],
        seasonality_aw = scenario_sheet[ROW, seas],
        gen_reduce = 0,
        gen_start = NA,
        gen_end = NA,
        gen2_reduce = 0,
        gen2_start = NA,
        gen2_end = NA,
        cert = scenario_sheet[ROW, cert],
        mask = scenario_sheet[ROW, mask],
        wfh_reference_date = ymd(scenario_sheet[ROW, wfh_date]),
        planB_start = "2022-04-01",
        planB_end = "2022-06-01",
        lockdown_reference_date = NA,
        lockdown_start = "2022-04-01",
        lockdown_end = "2022-06-01",
        n_threads = N_THREADS,
        VAX_SCHEDULE = scenario_sheet[ROW, vax_file],
        bst_rate1 = 229000,
        bst_rate2 = 1000000,
        bst_fold = BFOLD,
        bst_om_fold = 1,
        bst_duration = scenario_sheet[ROW, bdur],
        delta_fold_escape = 1,
        VE_debug_file = FALSE,
        p_boost_for_obs = scenario_sheet[ROW, pbst_all])
    #do_output(w, scenario_sheet[ROW, title], paste0("./output/paper/may22/output1_", scenario_sheet[ROW, codenm], ".pdf"))
    outfile = paste0("./output/paper/may22/output1_", scenario_sheet[ROW, codenm], ".qs");
    qsave(w, outfile)
    process_burdens(scenario_sheet[ROW, codenm], adj_file = NA, adj_ts = scenario_sheet[ROW, adj_ts])
    gc()
    # Notify Rosie
    notify_command = paste0(
        'curl -s --form-string "token=arckmxd33fp57d3ze6dujxeg1in48s" ',
        '--form-string "user=ubdi7mpz6bfiy3a485qkgpt5kzk9mb" ',
        '--form-string "message=Scenario ', ROW, ' of ', nrow(scenario_sheet), ' done." https://api.pushover.net/1/messages.json')
    system(notify_command)
}

# extra code to run locally for booster policy scenarios (for these scenarios we
# need to use the measured gamma multiplier from the basecase scenario to adjust
# the model outputs for each scenario). The following code snippet inside the if
# statement saves the gamma adjustments from the basecase scenario and can be 
# run locally
if (0){
  
    measure_final_adj = function(codenm, dstart = "2021-11-01", dend = "2021-11-30"){
        outfile = paste0("./output/paper/may22/output1_", codenm, ".qs");
        adjfile = paste0("./output/paper/may22/adjusted1_", codenm, ".qs");
        wo = qread(outfile)
        wa = qread(adjfile)
        
        ir = merge(
            wo[, .(t, NHS.region, run, group, deaths_o = deaths, hosp_adm_o = hosp_adm, hosp_bed_o = hosp_bed - hosp_undetected_p, icu_bed_o = icu_bed, pcr_incidence_o = pcr_positive_i, population)], 
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
            by = .(t, NHS.region, population, run)]
        
        ir[, date := t + ymd("2020-01-01")]
        
        adj = ir[date >= dstart & date <= dend, .(
            deaths_adj =   mean(deaths_a / deaths_o),
            hosp_adm_adj = mean(hosp_adm_a / hosp_adm_o),
            known_hosp_beds_adj = mean(hosp_bed_a / hosp_bed_o),
            icu_bed_adj =  mean(icu_bed_a / icu_bed_o)
        ), by = .(population, NHS.region)]
        
        return(adj)
    }
    
    adj = measure_final_adj("basecase_220511", dstart = "2021-11-01", dend = "2021-11-30")
    fwrite(adj, "./fits/final_adj_basecase_220214.csv")
    
    measure_final_adj_tseries = function(codenm, dstart = "2020-01-01", dend = "2022-12-31"){
        outfile = paste0("./output/paper/may22/output1_", codenm, ".qs");
        adjfile = paste0("./output/paper/may22/adjusted1_", codenm, ".qs");
        wo = qread(outfile)
        wa = qread(adjfile)
        
        ir = merge(
            wo[, .(t, NHS.region, run, group, deaths_o = deaths, hosp_adm_o = hosp_adm, hosp_bed_o = hosp_bed - hosp_undetected_p, icu_bed_o = icu_bed, pcr_incidence_o = pcr_positive_i, population)], 
            wa[, .(t, NHS.region, run, group, deaths_a = deaths, hosp_adm_a = hosp_adm, hosp_bed_a = hosp_bed - hosp_undetected_p, icu_bed_a = icu_bed, pcr_incidence_a = pcr_positive_i)],
            by = c("t", "NHS.region", "run", "group"))
        
        ira = ir[, .(
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
            by = .(t, NHS.region, population, run)]
        
        ira[, date := t + ymd("2020-01-01")]
        
        adj = ira[date >= dstart & date <= dend, .(
            deaths_adj =   mean(deaths_a / deaths_o),
            hosp_adm_adj = mean(hosp_adm_a / hosp_adm_o),
            known_hosp_beds_adj = mean(hosp_bed_a / hosp_bed_o),
            icu_bed_adj =  mean(icu_bed_a / icu_bed_o)
        ), by = .(t, NHS.region, date)]
        
        return(adj)
    }
    
    adj_tseries = measure_final_adj_tseries("basecase_220511", dstart = "2020-01-01", dend = "2022-12-31")
    fwrite(adj_tseries, "./fits/final_adj_tseries_basecase_220511.csv")
    
}

# old scenario sheet below from first submission of paper

# old_scenario_sheet = fread(
#     "codenm,     mobility, wane_yn,   vac_eff,   bage, pby, pbo, seas, cert,  mask,   wfh_date, title,                  vax_file
# basecase,     6mo,     yeswane,   central,     50, 0.0, 0.9,  0.1,    0,     0,         NA, Base case: 6-month relaxation; medium waning; boosters in 50+; 20% seasonality, vax-covidm20211019180925.rds
# vacc_hi,      6mo,     yeswane, centralHI,     50, 0.0, 0.9,  0.1,    0,     0,         NA, high waning,            vax-covidm20211019180925.rds
# vacc_lo,      6mo,     yeswane, centralLO,     50, 0.0, 0.9,  0.1,    0,     0,         NA, low waning,             vax-covidm20211019180925.rds
# mob_3wk,      3wk,     yeswane,   central,     50, 0.0, 0.9,  0.1,    0,     0,         NA, 3-week relaxation,      vax-covidm20211019180925.rds
# mob_3mo,      3mo,     yeswane,   central,     50, 0.0, 0.9,  0.1,    0,     0,         NA, 3-month relaxation,     vax-covidm20211019180925.rds
# mob_flt,      flt,     yeswane,   central,     50, 0.0, 0.9,  0.1,    0,     0,         NA, no relaxation,          vax-covidm20211019180925.rds
# bst_none,     6mo,     yeswane,   central,     50, 0.0, 0.0,  0.1,    0,     0,         NA, no boosters,            vax-covidm20211019180925.rds
# bst_50,       6mo,     yeswane,   central,     50, 0.0, 0.5,  0.1,    0,     0,         NA, boosters in 50+ (50%),  vax-covidm20211019180925.rds
# bst_5090,     6mo,     yeswane,   central,     50, 0.5, 0.9,  0.1,    0,     0,         NA, 50% in <50; 90% in 50+, vax-covidm20211019180925.rds
# bst_all,      6mo,     yeswane,   central,     50, 0.9, 0.9,  0.1,    0,     0,         NA, boosters for all,       vax-covidm20211019180925.rds
# Gbst_70plus,  6mo,     yeswane,   central,     70, 0.0, 0.9,  0.1,    0,     0,         NA, boosters in 70+ (90%),  vax-covidm20211019180925.rds
# Gbst_70plus,  6mo,     yeswane,   central,     70, 0.5, 0.9,  0.1,    0,     0,         NA, boosters in 70+ (90%),  vax-covidm20211019180925.rds
# bst_40plus,   6mo,     yeswane,   central,     40, 0.0, 0.9,  0.1,    0,     0,         NA, boosters in 40+ (90%),  vax-covidm20211019180925.rds
# bst_30plus,   6mo,     yeswane,   central,     30, 0.0, 0.9,  0.1,    0,     0,         NA, boosters in 30+ (90%),  vax-covidm20211019180925.rds
# B_cert,       6mo,     yeswane, centralLO,     50, 0.0, 0.9,  0.1, 0.17,     0,         NA, Plan B: Certification,  vax-covidm20211019180925.rds
# B_mask,       6mo,     yeswane, centralLO,     50, 0.0, 0.9,  0.1,    0, 0.075,         NA, Plan B: Face coverings, vax-covidm20211019180925.rds
# B_wfh,        6mo,     yeswane, centralLO,     50, 0.0, 0.9,  0.1,    0,     0, 2021-03-15, Plan B: Work from home, vax-covidm20211019180925.rds
# B_all,        6mo,     yeswane, centralLO,     50, 0.0, 0.9,  0.1, 0.17, 0.075, 2021-03-15, Plan B: All measures,   vax-covidm20211019180925.rds
# seas_10,      6mo,     yeswane,   central,     50, 0.0, 0.9, 0.05,    0,     0,         NA, 10% seasonality,        vax-covidm20211019180925.rds
# seas_30,      6mo,     yeswane,   central,     50, 0.0, 0.9, 0.15,    0,     0,         NA, 30% seasonality,        vax-covidm20211019180925.rds
# seas_40,      6mo,     yeswane,   central,     50, 0.0, 0.9, 0.20,    0,     0,         NA, 40% seasonality,        vax-covidm20211019180925.rds
# 12plus_70pc,  6mo,     yeswane,   central,     50, 0.0, 0.9,  0.1,    0,     0,         NA, vax 12+ 70% uptake,     vax-covidm20211019180925.rds
# 12plus_50pc,  6mo,     yeswane,   central,     50, 0.0, 0.9,  0.1,    0,     0,         NA, vax 12+ 50% uptake,     vax-covidm2021101918092512plus_50percent.rds
# 5plus_70pc,   6mo,     yeswane,   central,     50, 0.0, 0.9,  0.1,    0,     0,         NA, vax 5+  70% uptake,     vax-covidm202110191809255plus_70percent.rds
# 5plus_50pc,   6mo,     yeswane,   central,     50, 0.0, 0.9,  0.1,    0,     0,         NA, vax 5+  50% uptake,     vax-covidm202110191809255plus_50percent.rds
# ")
