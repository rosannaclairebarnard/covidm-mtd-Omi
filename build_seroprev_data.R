# Process sero(logical) prevalence data @ NHS England regional scales for covidm

# We use a combination of different data sources for England and English regions
# and add new data to existing data on an ad-hoc basis (unlike virus prevalence)
#
# The following data sources are used for regional sero prevalence estimates:
#
# 1. Sero prevalence estimates from the UK Biobank, see:
# https://www.ukbiobank.ac.uk/learn-more-about-uk-biobank/covid-19-hub 
#
# 2. Sero prevalence estimates from Imperial College London's REACT-2 study, see:
# https://www.imperial.ac.uk/medicine/research-and-impact/groups/react-study/
#
# 3. Sero prevalence estimates from the Office for National Statistics, see:
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk

library(readxl)
library(lubridate)
library(stringr)     #
library(data.table)  # these libraries must be loaded for ogwrangler to work
library(glue)        # 
library(ogwrangler)  # see https://github.com/nicholasdavies/ogwrangler
library(ggplot2)
library(ggthemes)
library(cowplot)

# remove existing objects in environment
rm(list=ls())

# Specify existing .csv file containing seroprevalence estimates
file = './fitting_data/seroprev_nhs_regions_20220421090540.csv'

# Load existing file
sero = read.csv(file)

# Make sure dates in sero are consistent formats
if (is.Date(sero$Start.date) == FALSE){
    sero$Start.date = as.Date(sero$Start.date)
}
if (is.Date(sero$End.date) == FALSE){
    sero$End.date   = as.Date(sero$End.date)
}

# Specify path to file, filename and the publication date of new ONS data
newdata = './fitting_data/covid19infectionsurveyantibodies20220504.xlsx'

# We usually find the data required within Excel sheets 1b, but this may change 
# if the format of the published .xlsx file changes, so please check! Note that 
# equivalent data for England as a whole is (currently) found in sheet 1a
s1b = read_excel(newdata, sheet = "1b")
s1a = read_excel(newdata, sheet = "1a")

# find last row containing data - this assumes that specific text will continue 
# to be placed in the row below the last row of data, so please check in future
last_idx_eng = which(s1a$Contents == 'Source: Office for National Statistics – Coronavirus (COVID-19) Infection Survey')
last_idx = which(s1b$Contents == 'Source: Office for National Statistics – Coronavirus (COVID-19) Infection Survey')

# this also assumes row 4 contains appropriate headers, please check in future
s1ad = read_excel(newdata, sheet = "1a", range = cell_rows(6:last_idx_eng))
s1bd = read_excel(newdata, sheet = "1b", range = cell_rows(5:last_idx))

# get correct date ranges for each row entry in s1ad
start_dates_eng = as.Date(NA)
end_dates_eng   = as.Date(NA)
for (i in 1:dim(s1ad)[1]){
    daterange      = s1ad$`Weekly period`[i]
    s1ad_dates     = strsplit(daterange, split = " to ")[[1]]
    start_dates_eng[i] = as.Date(s1ad_dates[1], format = '%d %b %Y')
    end_dates_eng[i]   = as.Date(s1ad_dates[2], format = '%d %b %Y')
}
s1ad$Start.date = start_dates_eng
s1ad$End.date   = end_dates_eng

# get correct date ranges for each row entry in s1bd
start_dates     = as.Date(NA)
end_dates       = as.Date(NA)
for (i in 2:dim(s1bd)[1]){
    daterange      = s1bd$`...1`[i]
    s1bd_dates     = strsplit(daterange, split = " to ")[[1]]
    start_dates[i] = as.Date(s1bd_dates[1], format = '%d %b %Y')
    end_dates[i]   = as.Date(s1bd_dates[2], format = '%d %b %Y')
}
s1bd$Start.date = start_dates
s1bd$End.date   = end_dates

# get PHE regions and codes (used by ONS and REACT)
regions = read.csv("./fitting_data/PHE_England_regions.csv")

# initialise data frames to store resized central, lower, upper bound estimates
rszd_swab_prev = data.frame()
rszd_swab_prev_lb = data.frame()
rszd_swab_prev_ub = data.frame()

# for each date range, convert seroprevalence estimates for all regions
for (i in 2:dim(s1bd)[1]){
    
    # print progress of loop
    print(paste0('Processing date ', i-1, ' of ', dim(s1bd)[1]-1))
    
    # get list of central estimates to convert
    estimates = as.numeric(c(s1bd$`North East`[i],
                             s1bd$`North West`[i],
                             s1bd$`Yorkshire and The Humber`[i],
                             s1bd$`East Midlands`[i],
                             s1bd$`West Midlands`[i],
                             s1bd$`East of England`[i],
                             s1bd$London[i],
                             s1bd$`South East`[i],
                             s1bd$`South West`[i]))
    
    # convert estimates to fit nhs regions
    resized_estimates = ogwrangle(estimates, regions$RGN19CD, "e.reg", 
                                  "e.nhser20", "pop2018", "proportion")
    
    # store converted central estimates
    rszd_swab_prev = rbind(rszd_swab_prev, resized_estimates$estimates)
    
    # get list of lower bound estimates
    lb_estimates = as.numeric(c(s1bd$...3[i],   # North East
                                s1bd$...7[i],   # North West
                                s1bd$...11[i],  # Yorkshire and The Humber
                                s1bd$...15[i],  # East Midlands
                                s1bd$...19[i],  # West Midlands
                                s1bd$...23[i],  # East of England
                                s1bd$...27[i],  # London
                                s1bd$...31[i],  # South East
                                s1bd$...35[i])) # South West
    
    # convert lower bound estimates to fit nhs regions
    resized_lb_estimates = ogwrangle(lb_estimates, regions$RGN19CD, "e.reg", 
                                     "e.nhser20", "pop2018", "proportion")
    
    # store converted lower bound estimates
    rszd_swab_prev_lb = rbind(rszd_swab_prev_lb, 
                              resized_lb_estimates$lb_estimates)
    
    # get list of upper bound estimates (skipping values for England)
    ub_estimates = as.numeric(c(s1bd$...4[i],   # North East
                                s1bd$...8[i],   # North West
                                s1bd$...12[i],  # Yorkshire and The Humber
                                s1bd$...16[i],  # East Midlands
                                s1bd$...20[i],  # West Midlands
                                s1bd$...24[i],  # East of England
                                s1bd$...28[i],  # London
                                s1bd$...32[i],  # South East
                                s1bd$...36[i])) # South West
    
    # convert estimates to fit nhs regions
    resized_ub_estimates = ogwrangle(ub_estimates, regions$RGN19CD, "e.reg", 
                                     "e.nhser20", "pop2018", "proportion")
    
    # store converted upper bound estimates
    rszd_swab_prev_ub = rbind(rszd_swab_prev_ub, 
                              resized_ub_estimates$ub_estimates)
}

# add column headings to all three data frames
names(rszd_swab_prev)    = ogname(resized_estimates$e.nhser20)
names(rszd_swab_prev_lb) = ogname(resized_lb_estimates$e.nhser20)
names(rszd_swab_prev_ub) = ogname(resized_ub_estimates$e.nhser20)

# create list of NHS England regional codes
nhs_regions <- c("E40000009", "E40000010", "E40000008", "E40000007", 
                 "E40000003", "E40000005", "E40000006")
names(nhs_regions) <- ogname(nhs_regions)

# initialise final data frame to store sero prevalence data
sero_prev <- data.frame(NHS.region = NULL, Start.date = NULL, End.date = NULL, 
                         Central.estimate = NULL, Lower.bound = NULL, 
                         Upper.bound = NULL, Test = NULL, Median.or.mean = NULL, 
                         Min.age = NULL, Max.age = NULL, Data.source = NULL, 
                         N.tests = NULL, React2.round = NULL)

# count number of entries of sero prevalence per region in rszd_swab_prev
num_entries <- dim(rszd_swab_prev)[1]

for (i in 1:length(nhs_regions)){
    
    # populate entries for this section of the dataframe
    this_region = rep(names(nhs_regions)[i], num_entries)
    start_date  = as.character(start_dates[1:num_entries+1])
    end_date    = as.character(end_dates[1:num_entries+1])
    central_estimate = rszd_swab_prev[,i]
    lower_bound      = rszd_swab_prev_lb[,i]
    upper_bound      = rszd_swab_prev_ub[,i]
    test <- rep("", num_entries)
    median_mean <- rep("", num_entries)
    min_age <- rep("", num_entries)
    max_age <- rep("", num_entries)
    data_source <- rep("ONS-CIS", num_entries)
    N_tests <- rep("", num_entries)
    react_round <- rep("", num_entries)
    
    # make temporary dataframe to bind to master dataframe
    temp_df <- data.frame(NHS.region       = this_region, 
                          Start.date       = start_date, 
                          End.date         = end_date, 
                          Central.estimate = as.numeric(central_estimate), 
                          Lower.bound      = as.numeric(lower_bound), 
                          Upper.bound      = as.numeric(upper_bound), 
                          Test             = test, 
                          Median.or.mean   = median_mean, 
                          Min.age          = min_age, 
                          Max.age          = max_age, 
                          Data.source      = data_source, 
                          N.tests          = N_tests,
                          React2.round     = react_round)
    
    # bind temporary dataframe to master dataframe
    sero_prev <- rbind(sero_prev, temp_df)
}

# add England-wide estimates to sero_prev
if (1){
    cat("Adding England entries...")
    # populate entries for this section of the dataframe
    num_entries = dim(s1ad)[1]
    this_region = rep('England', num_entries)
    start_date  = as.character(start_dates_eng)
    end_date    = as.character(end_dates_eng)
    central_estimate = s1ad$`Modelled % testing positive for antibodies against SARS-CoV-2`
    lower_bound      = s1ad$`95% Lower credible interval`
    upper_bound      = s1ad$`95% Upper credible interval`
    test <- rep("", num_entries)
    median_mean <- rep("", num_entries)
    min_age <- rep("", num_entries)
    max_age <- rep("", num_entries)
    data_source <- rep("ONS-CIS", num_entries)
    N_tests <- rep("", num_entries)
    react_round <- rep("", num_entries)

    # make temporary dataframe to bind to master dataframe
    temp_df <- data.frame(NHS.region       = this_region,
                          Start.date       = start_date,
                          End.date         = end_date,
                          Central.estimate = as.numeric(central_estimate),
                          Lower.bound      = as.numeric(lower_bound),
                          Upper.bound      = as.numeric(upper_bound),
                          Test             = test,
                          Median.or.mean   = median_mean,
                          Min.age          = min_age,
                          Max.age          = max_age,
                          Data.source      = data_source,
                          N.tests          = N_tests,
                          React2.round     = react_round)

    # bind temporary dataframe to master dataframe
    sero_prev <- rbind(sero_prev, temp_df)
    
}

# append % to estimates (existing estimates have percentage signs which get 
# removed later)
sero_prev$Central.estimate = paste0(sero_prev$Central.estimate, "%")
sero_prev$Lower.bound = paste0(sero_prev$Lower.bound, "%")
sero_prev$Upper.bound = paste0(sero_prev$Upper.bound, "%")

# before binding new ONS-CIS data in sero_prev to old ONS-CIS data in sero, 
# remove any old data from ONS-CIS in sero which already exist within the date
# ranges in the new data sero_prev 
earliest_date = min(start_dates, na.rm = TRUE)
latest_date   = max(end_dates,   na.rm = TRUE)
series_dates  = seq(earliest_date, latest_date, by = 1)

# get row elements to remove - sero$Data.source is ONS-CIS, sero$Start.date 
# occurs in series_dates and sero$End.date occurs in series_dates
rtr = which(sero$Data.source == 'ONS-CIS' & sero$Start.date %in% series_dates & 
                sero$End.date %in% series_dates)

# remove those rows from sero
sero = sero[-(rtr),]

# bind new data to existing data
final_data = rbind(sero, sero_prev)

# save final data 
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
write.csv(final_data, file=paste0("./fitting_data/seroprev_nhs_regions_", 
                                  datetime, ".csv"), row.names = F)

# optional: plots to check data (change to 'if (1)' to execute)
if (0){
    
    # function required to remove % signs
    pct = function(x) as.numeric(str_replace_all(x, "%", ""))
    
    # plot all regions together as single data points
    ggplot(data = final_data[!(final_data$NHS.region %in% c('Wales','Scotland')),]) +
        geom_point(aes(x = Start.date + (End.date-Start.date)/2, y = pct(Central.estimate), color = NHS.region, shape = Data.source), size = 1.5) +
        labs(x = "Date", y = "Seroprevalence (%)", colour = "NHS England region", shape = 'Data source') +
        scale_x_date(date_breaks = "2 months", date_labels = "%b %Y")  +
        scale_color_colorblind()
    
    # plot all regions together with uncertainty included
    ggplot(data = final_data[!(final_data$NHS.region %in% c('Wales','Scotland','England')),]) +
        geom_point(aes(x = Start.date + (End.date-Start.date)/2, y = pct(Central.estimate), color = NHS.region, shape = Data.source), size = 1.25) +
        geom_linerange(aes(xmin = Start.date, xmax = End.date, y = pct(Central.estimate), color = NHS.region, linetype = Data.source), size = 0.2) +
        geom_linerange(aes(x = Start.date + (End.date-Start.date)/2, ymin = pct(Lower.bound), ymax = pct(Upper.bound), color = NHS.region, linetype = Data.source), size = 0.2) + 
        labs(x = "Date", y = "Seroprevalence (%)", colour = "NHS England region", shape = 'Data source') +
        scale_x_date(date_breaks = "2 months", date_labels = "%b %Y")  +
        scale_color_colorblind()
    
    # plot regions separately for comparison
    ggplot(data = final_data[!(final_data$NHS.region %in% c('Wales','Scotland','England')),]) +
        geom_point(aes(x = Start.date + (End.date-Start.date)/2, y = pct(Central.estimate), shape = Data.source, color = Data.source), size = 1.5) +
        geom_ribbon(aes(x = Start.date + (End.date-Start.date)/2, ymin = pct(Lower.bound), 
                        ymax = pct(Upper.bound)), alpha = 0.3) +
        labs(x = "Date", y = "Seroprevalence (%)", shape = 'Data source', colour = 'Data source') +
        scale_x_date(date_breaks = "2 months", date_labels = "%b %Y")  +
    facet_grid(NHS.region ~ ., scales = "free", switch = "y")
    
    # plot regions separately for comparison with uncertainty
    ggplot(data = final_data[!(final_data$NHS.region %in% c('Wales','Scotland','England')),]) +
        geom_point(aes(x = Start.date + (End.date-Start.date)/2, y = pct(Central.estimate), shape = Data.source, color = Data.source), size = 1.5) +
        geom_linerange(aes(xmin = Start.date, xmax = End.date, y = pct(Central.estimate), linetype = Data.source), size = 0.2) +
        geom_linerange(aes(x = Start.date + (End.date-Start.date)/2, ymin = pct(Lower.bound), ymax = pct(Upper.bound), linetype = Data.source), size = 0.2) +
        labs(x = "Date", y = "Seroprevalence (%)", linetype = 'Data source', colour = 'Data source') +
        scale_x_date(date_breaks = "2 months", date_labels = "%b %Y")  +
        facet_grid(NHS.region ~ ., scales = "free", switch = "y")
    
    # plot one region on it's own
    ggplot(data = final_data[(final_data$NHS.region == 'South West'),]) +
        geom_point(aes(x = Start.date + (End.date-Start.date)/2, y = pct(Central.estimate), shape = Data.source, color = Data.source), size = 1.5) +
        geom_linerange(aes(xmin = Start.date, xmax = End.date, y = pct(Central.estimate), linetype = Data.source), size = 0.2) +
        geom_linerange(aes(x = Start.date + (End.date-Start.date)/2, ymin = pct(Lower.bound), ymax = pct(Upper.bound), linetype = Data.source), size = 0.2) +
        labs(x = "Date", y = "Seroprevalence (%)", linetype = 'Data source', colour = 'Data source') +
        scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") 
    
    # plot each NHS England region + England as a whole on separate panels
    p = list()
    for (i in 1:length(unique(final_data$NHS.region))){
        this_region = unique(final_data$NHS.region)[i]
        p[[i]] = ggplot(data = final_data[(final_data$NHS.region == this_region),]) +
            geom_point(aes(x = Start.date + (End.date-Start.date)/2, 
                            y = pct(Central.estimate), shape = Data.source, color = Data.source), size = 1.5) +
            geom_linerange(aes(xmin = Start.date, xmax = End.date, y = pct(Central.estimate), linetype = Data.source), size = 0.2) +
            geom_linerange(aes(x = Start.date + (End.date-Start.date)/2, ymin = pct(Lower.bound), ymax = pct(Upper.bound), linetype = Data.source), size = 0.2) +
            labs(x = "Date", y = "Seroprevalence (%)", title = this_region) +
            scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
            theme(legend.position = 'none')
    }
    c1 = plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], nrow = 5)
    c2 = plot_grid(p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], nrow = 5)
    plot_grid(c1, c2, nrow = 1)
    
    }
    