library(data.table)
library(ggplot2)
library(lubridate)
library(here)
library(cowplot)
library(readxl)
library(sn)
library(qs)
library(stringr)
 
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
# source("./distribution_fit.R");
source("./spim_output.R");

#
# FOR PROPORTION DATA
#
approx_n = function(central, lo95, hi95)
{
    if (length(central) != length(lo95) || length(lo95) != length(hi95)) {
        stop("central, lo95, and hi95 must all be the same length.")
    }
    if (any(c(central, lo95, hi95) < 0) || any(c(central, lo95, hi95) > 1) ||
            any(lo95 > central) || any(central > hi95) || any(lo95 >= hi95)) {
        stop("central, lo95 and hi95 must be between 0 and 1 and be correctly ordered.")
    }
    dist = function(n_est_log10, central, lo95, hi95)
    {
        n_est = 10 ^ n_est_log10;
        alpha = n_est * central;
        beta = n_est * (1 - central);
        w1 = hi95 - lo95;
        w2 = qbeta(0.975, alpha, beta) - qbeta(0.025, alpha, beta);
        (w2 - w1) ^ 2
    }
    
    n_est = rep(0, length(central));
    for (i in seq_along(central))
    {
        n_est[i] = 10 ^ optimize(dist, c(0, 9), central = central[i], lo95 = lo95[i], hi95 = hi95[i])$minimum;
    }
    n_est
}


#
# DATA
#

nhs_regions = popUK[, unique(name)]
pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100

# Process new data
new_data = fread("fitting_data/data-20220506115838.csv")

# removing problematic data with under reporting of hospital incidence
days_to_remove_scot <- 2  ## VISUALLY INSPECT DATA TO DETERMINE THESE VALUES
days_to_remove_wales <- 1 ## VISUALLY INSPECT DATA TO DETERMINE THESE VALUES
scottish_data <- new_data[new_data$location == "Scotland"]
scottish_hospital_inc <- scottish_data[scottish_data$indicator == "hospital_inc"]
max_date_scotland <- max(new_data[new_data$location=="Scotland" & 
                                    new_data$indicator=="hospital_inc"]$date)
dates_to_remove_scotland <- seq(scottish_hospital_inc$date[length(scottish_hospital_inc$date)
                                                           -(days_to_remove_scot-1)], 
                                by = "day", 
                                length.out = (max_date_scotland-
                                                scottish_hospital_inc$date[length(scottish_hospital_inc$date)
                                                                           -(days_to_remove_scot)]))
welsh_data <- new_data[new_data$location == "Wales"]
welsh_hospital_inc <- welsh_data[welsh_data$indicator == "hospital_inc"]
max_date_wales <- max(new_data[new_data$location=="Wales" & 
                                 new_data$indicator=="hospital_inc"]$date)
dates_to_remove_wales <- seq(welsh_hospital_inc$date[length(welsh_hospital_inc$date)
                                                     -(days_to_remove_wales-1)],
                            by = "day",
                             length.out = (max_date_wales-
                                             welsh_hospital_inc$date[length(welsh_hospital_inc$date)
                                                                     -(days_to_remove_wales)]))
new_data <- new_data[!(new_data$location == "Scotland" & 
                         new_data$indicator == "hospital_inc" & 
                         new_data$date %in% dates_to_remove_scotland)]
new_data <- new_data[!(new_data$location == "Wales" & 
                         new_data$indicator == "hospital_inc" & 
                         new_data$date %in% dates_to_remove_wales)]

new_data[, pid := match(location, nhs_regions) - 1]
ld = new_data[indicator == "death_inc_line", .(date, N = value, name = location, pid)];
sitreps = new_data[indicator != "death_inc_line", .(date = ymd(date), value_type = indicator, value, name = location, pid)]
sitreps = dcast(sitreps, date + name + pid ~ value_type, value.var = "value", fill = NA, fun.aggregate = function(x) x[1])
sitreps = sitreps[, .(date, n_in_itu = icu_prev, n_in_all_beds = hospital_prev, n_admitted_diagnosed = hospital_inc, name, pid)];


sero = fread("fitting_data/seroprev_nhs_regions_20220506122452.csv")
sero[, p := pct(Central.estimate)]
sero[, n_approx := approx_n(pct(Central.estimate), pct(Lower.bound), pct(Upper.bound))]
#sero[, Start.date := dmy(Start.date)]
#sero[, End.date := dmy(End.date)]

virus = fread("fitting_data/virusprev_nhs_regions_20220506122021.csv")
#virus[, Start.date := dmy(Start.date)]
#virus[, End.date := dmy(End.date)]
virus[, p := pct(Central.estimate)]
# remove any missing values
virus = virus[!(is.na(virus$Central.estimate)==TRUE),]
virus = virus[!(is.na(virus$Lower.bound)==TRUE),]
virus = virus[!(is.na(virus$Upper.bound)==TRUE),]
virus[, n_approx := approx_n(pct(Central.estimate), pct(Lower.bound), pct(Upper.bound))]
sero[, pid := match(NHS.region, nhs_regions) - 1]
virus[, pid := match(NHS.region, nhs_regions) - 1]


# Add England to deaths series
ld = rbind(ld,
    ld[!name %in% c("Northern Ireland", "Scotland", "Wales"), .(N = sum(N), name = "England", pid = 1), by = date]
)
ld[, date := as.Date(date)]

# Add England to sitrep series
sitreps = rbind(sitreps,
    sitreps[!name %in% c("Northern Ireland", "Scotland", "Wales"),
        .(n_in_itu = sum(n_in_itu, na.rm = T), n_in_all_beds = sum(n_in_all_beds, na.rm = T), n_admitted_diagnosed = sum(n_admitted_diagnosed, na.rm = T),
            name = "England", pid = 1),
        by = date]
)

# SGTF data, add England (this data corresponds to the introduction and decline of Alpha B.1.1.7)
sgtf = fread(datapath("sgtf-2021-12-15.csv"))

sgtf = rbind(sgtf, 
    sgtf[!nhs_name %in% c("Northern Ireland", "Scotland", "Wales"),
        .(sgtf = sum(sgtf, na.rm = T), other = sum(other, na.rm = T), nhs_name = "England"),
        by = date], fill = TRUE
)
sgtf[, pid := match(nhs_name, nhs_regions) - 1]
sgtf[, date := as.Date(date)]

# Variant data from PHE
delta = qread("./fitting_data/delta-20211205151636.qs")
# delta = delta[pillar == "Pillar 2" & v_overall_travel == 0 & NHSER_name != "", 
#     .(delta = sum(v_variant == "VOC-21APR-02"), other = sum(v_variant != "VOC-21APR-02")), keyby = .(date = v_specimen_date_sk, nhs_name = NHSER_name)]
# delta = rbind(delta, 
#     delta[!nhs_name %in% c("Northern Ireland", "Scotland", "Wales"),
#         .(delta = sum(delta, na.rm = T), other = sum(other, na.rm = T), nhs_name = "England"),
#         by = date], fill = TRUE
# )
delta[, pid := match(nhs_name, nhs_regions) - 1]

# SGTF (Omicron only) data, add England
omi = fread(datapath("sgtf-2022-01-10-omicron-only.csv"))

omi = rbind(omi, 
             omi[!nhs_name %in% c("Northern Ireland", "Scotland", "Wales"),
                  .(sgtf = sum(sgtf, na.rm = T), other = sum(other, na.rm = T), nhs_name = "England"),
                  by = date], fill = TRUE
)
omi[, pid := match(nhs_name, nhs_regions) - 1]
omi[, date := as.Date(date)]


# BA.2 proportion from Sanger data
sanger = fread("fitting_data/Proportion_in_England-2022-05-06.csv")
ba2 = merge(
    sanger[lineage %like% "^BA\\.1(\\.|$)", .(ba1 = sum(value)), by = date],
    sanger[lineage %like% "^BA\\.2(\\.|$)", .(ba2 = sum(value)), by = date],
    by = "date")[, .(t = as.numeric(as.Date(date) - ymd("2020-01-01")), ba2 = ba2 / (ba1 + ba2))]
# Remove earlier entries with NaN as ba2 proportion
last_nan = max(ba2[!is.finite(ba2), t])
ba2 = ba2[t > last_nan]

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
qsave(list(ld, sitreps, virus, sero, sgtf, delta, omi, ba2), 
      datapath(paste0("processed-data-", datetime, ".qs")))
