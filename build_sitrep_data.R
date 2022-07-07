library(data.table)
library(readxl)
library(ggplot2)
library(lubridate)
library(ogwrangler)

source("./commit.R")

# Note: this file contains sensitive NHS data, so it is not included wth the repo.
spim_data = read_excel("~/Documents/uk_covid_data_sensitive/ukhsa/20220503_All_SPIM_trust__002.xlsx", "Extracted Data", col_type = "text")
setDT(spim_data)
spim_data[, DateVal := as.Date(DateVal)]
spim_data = cbind(
    spim_data[, .SD, .SDcols = c(1, 5, 6, 7, 8)],
    spim_data[, lapply(.SD, as.numeric), .SDcols = c(2, 3, 4, 9:ncol(spim_data))]
)

regions = c("East of England", "London", "Midlands", "North East and Yorkshire",
    "North West", "South West",  "South East", "Northern Ireland", "Scotland", "Wales")
rlevels = c("Region", "National")

# hospital_prev and icu_prev
hpEN = melt_annotate(spim_data[Geography %in% regions & !(Geography %in% c("Scotland", "Wales")) & ReportLevel %in% rlevels & !is.na(hospital_prev), hospital_prev,
  keyby = .(location = Geography, date = ymd(DateVal))])

hpW = melt_annotate(spim_data[Geography %in% regions & Geography == "Wales" & ReportLevel %in% rlevels & !is.na(hospital_prev), .(hospital_prev = hospital_prev + hospital_prev_recovering),
  keyby = .(location = Geography, date = ymd(DateVal))])

hpS = melt_annotate(spim_data[Geography == "Scotland" & ReportLevel %in% rlevels & !is.na(`hospital_prev_<28days`), .(hospital_prev = `hospital_prev_<28days`),
  keyby = .(location = Geography, date = ymd(DateVal))])

hp = rbind(hpEN, hpW, hpS)

ipENW = melt_annotate(spim_data[Geography %in% regions & Geography != "Scotland" & ReportLevel %in% rlevels & !is.na(icu_prev), icu_prev,
  keyby = .(location = Geography, date = ymd(DateVal))])

ipS = melt_annotate(spim_data[Geography == "Scotland" & ReportLevel %in% rlevels & !is.na(`icu_prev_<28days`), .(icu_prev = `icu_prev_<28days`),
  keyby = .(location = Geography, date = ymd(DateVal))])

ip = rbind(ipENW, ipS)

# hospital_inc
sumNA = function(x, y)
{
    ifelse(is.na(x) & is.na(y), NA, sum(c(x, y), na.rm = T))
}

hiENS = melt_annotate(spim_data[Geography %in% regions & Geography != "Wales" & ReportLevel %in% rlevels & !is.na(sumNA(hospital_inc, hospital_inc_new)),
  .(hospital_inc = sumNA(hospital_inc, hospital_inc_new)),
  keyby = .(location = Geography, date = ymd(DateVal))])

hiW = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(SPIM_hosp_inc_new),
  .(hospital_inc = SPIM_hosp_inc_new),
  keyby = .(location = Geography, date = ymd(DateVal))])

hi = rbindlist(list(hiENS, hiW[date >= "2020-02-01"]))

# death_inc_line
diE = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(PHE_type28_death_inc_line),
  .(death_inc_line = PHE_type28_death_inc_line),
  keyby = .(location = Geography, date = ymd(DateVal))])

diN = melt_annotate(spim_data[Geography == "Northern Ireland" & ReportLevel %in% rlevels & !is.na(SitRep_death_inc_line),
  .(death_inc_line = SitRep_death_inc_line),
  keyby = .(location = Geography, date = ymd(DateVal))])

diS = spim_data[Geography == "Scotland" & ReportLevel %in% rlevels & !is.na(SitRep_death_inc_line),
  .(death_inc_line = SitRep_death_inc_line),
  keyby = .(location = Geography, date = ymd(DateVal))]
zero_filler = data.table(location = "Scotland", date = min(diS$date) + 0:as.numeric(max(diS$date) - min(diS$date)))
diS = merge(zero_filler, diS, by = c("location", "date"), all = TRUE)
diS[is.na(death_inc_line), death_inc_line := 0]
diS = melt_annotate(diS)

diW = melt_annotate(spim_data[Geography == "Wales" & ReportLevel %in% rlevels & !is.na(PHW_Death_inc_line),
  .(death_inc_line = PHW_Death_inc_line),
  keyby = .(location = Geography, date = ymd(DateVal))])

di = rbindlist(list(diE, diN, diS, diW))

# Amalgamate
data = rbindlist(list(hp, ip, hi, di))

# Compare previous build to new build
existing = fread("~/Dropbox/uk_covid_data/data-20220429123339.csv")
existing[, date := as.Date(date)]
committed = commit(existing, data)

if(1){ # 21st April - removing spurious datapoint for Midlands icu_prev
    # inspect following subset to determine which data needs removing
    committed[location == 'Midlands' & indicator == 'icu_prev' & date > '2022-03-28' & date < '2022-04-10']
    # datapoint is 257 ICU beds occupied on 3rd April 2022 
    # (up from 42 the previous day and 46 the day after)
    # removing that day's data for model fitting purposes
    committed = committed[!(location == 'Midlands' & indicator == 'icu_prev' & date == '2022-04-03')]
}

# Save final output
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
fwrite(committed, paste0("~/Dropbox/uk_covid_data/data-", datetime, ".csv"))
fwrite(committed, paste0("./fitting_data/data-", datetime, ".csv"))
