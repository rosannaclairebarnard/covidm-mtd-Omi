# Process virus prevalence data at NHS England regional scales for covidm

# We use virus prevalence data for England and for English regions from the
# Office for National Statistics' Coronavirus Infection Survey (ONS-CIS):
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata

library(readxl)
library(lubridate)
library(stringr)     #
library(data.table)  # these libraries must be loaded for ogwrangler to work
library(glue)        # 
library(ogwrangler)  # see https://github.com/nicholasdavies/ogwrangler
library(ggplot2)
library(ggthemes)
library(cowplot)

# Specify path to file and filename as well as the publication date of the file
file = './fitting_data/20220819covid19infectionsurveydatasetsengland.xlsx'
publicationdate = as.Date('2022-08-19')

# We usually find the data required within Excel sheets 1b, 1f and 1q, but this 
# may change if the format of the published .xlsx file changes, so please check!
# (Last checked and updated on 1st April 2022)
sheet_name_1b = "1b"
sheet_name_1f = "1f"
sheet_name_1p = '1m'
s1b = read_excel(file, sheet = sheet_name_1b)
s1f = read_excel(file, sheet = sheet_name_1f)
s1p = read_excel(file, sheet = sheet_name_1p)

# As of 3rd September 2021, sheet 1b contains: 'Modelled daily rates of the 
# percentage of the population testing positive for COVID-19, England', sheet 1f 
# contains: 'Modelled daily rates of the percentage of the population testing 
# positive for COVID-19 by region, England' and sheet 1p contains: 'Unrounded 
# modelled daily rates of the percentage of the population testing positive for 
# COVID-19 by region (historic series), England'

# For sheet 1p, make sure each time series is marked with a publication date
# in the cell above the data and that the cell below the last datapoint in the 
# last time series is marked 'Source: ...' (see line 57 in this script)

# sheet 1p contains the historic time series of virus prevalence both for 
# England followed by English regions, split into 6-week publication batches

num_timeseries = 0
pub_dates = NULL
row_elmts = NULL
row_start = NULL
first_ts_last_row = NULL

for (i in 1:dim(s1p)[1]){
    
    if (s1p$Contents[i] %like% 'Publication'){
        
        num_timeseries = num_timeseries + 1
        print(s1p$Contents[i])
        pub_date = s1p$Contents[i]
        row_elmts[num_timeseries] = i
        row_to_check = i+2
        
        while(is.na(s1p$Contents[row_to_check]) == TRUE){
            row_to_check = row_to_check + 1
        }
        
        row_start[num_timeseries] = row_to_check
        
        if (pub_date %like% 'Publication Date:'){
            pub_dates[num_timeseries] = substring(pub_date, 19, nchar(pub_date))
        } else if (pub_date %like% 'Publication date'){
            pub_dates[num_timeseries] = substring(pub_date, 18, nchar(pub_date))
        } else {
            stop('Error: publication date not recorded')
        }
        
    } else if (s1p$Contents[i] %like% 'Source:'){
        first_ts_last_row = i-1
    }
}

pub_dates = as.Date(pub_dates, format = '%d %b %Y')
colnames = as.character(s1p[6,], na.rm = TRUE)
s1p_df = data.frame(matrix(ncol = length(colnames), nrow = 0))
colnames(s1p_df) = colnames

for (i in num_timeseries:1){
    
    pub_date    = pub_dates[i]
    first_row   = row_start[i]
    
    if (i == num_timeseries){
        
        # for the earliest time series, use first_ts_last_row as the end of data
        this_data = s1p[(first_row:first_ts_last_row),]
        this_data$publication_date = rep(pub_date, dim(this_data)[1])
        
        if (nchar(this_data$Contents[1]) > 5){
            this_data$date = as.Date(this_data$Contents, format = '%d %b %Y')
        } else {
            this_data$date = as.Date(as.numeric(this_data$Contents), 
                                     origin = '1899-12-30')
        }
    
        # bind data to output dataframe
        s1p_df    = rbind(s1p_df, this_data)
        
    } else {
        
        # for the remaining time series, use previous pub_date as end of data
        last_row  = row_elmts[i+1]-2
        this_data = s1p[(first_row:last_row),]
        this_data$publication_date = rep(pub_date, dim(this_data)[1])
        
        if (nchar(this_data$Contents[1]) > 5){
            this_data$date = as.Date(this_data$Contents, format = '%d %b %Y')
        } else {
            this_data$date = as.Date(as.numeric(this_data$Contents), origin = '1899-12-30')
        }
        
        # before binding new data to old data, remove any old dates which the
        # dates within new data spans (new data supersedes old data)
        s1p_df <- subset(s1p_df, !(date %in% this_data$date))
        
        # bind new data to old
        s1p_df    = rbind(s1p_df, this_data)
    }
}

# # plots to check above data
# plot(s1p_df$date, s1p_df$...2)
# plot(s1p_df$date, s1p_df$...5)
# plot(s1p_df$date, s1p_df$...8)
# plot(s1p_df$date, s1p_df$...11)
# plot(s1p_df$date, s1p_df$...14)
# plot(s1p_df$date, s1p_df$...17)
# plot(s1p_df$date, s1p_df$...20)
# plot(s1p_df$date, s1p_df$...23)
# plot(s1p_df$date, s1p_df$...26)
# plot(s1p_df$date, s1p_df$...29)

# next, combine s1p_df with the latest published data in sheets 1b and 1f

# sheet 1b corresponds to the latest published data for England
# sheet 1f corresponds to the latest published data for English regions

# get date ranges from sheets 1b and 1f (check these are correct by hand)
s1b_daterange = s1b$Contents[3]
s1b_dates     = strsplit(s1b_daterange, split = " to ")[[1]]
s1b_start     = as.Date(s1b_dates[1], format = '%d %b %Y')
s1b_end       = as.Date(s1b_dates[2], format = '%d %b %Y')
print(s1b_daterange)
print(s1b_start)
print(s1b_end)
s1f_daterange = s1f$Contents[3]
s1f_dates     = strsplit(s1f_daterange, split = " to ")[[1]]
s1f_start     = as.Date(s1f_dates[1], format = '%d %b %Y')
s1f_end       = as.Date(s1f_dates[2], format = '%d %b %Y')
print(s1f_daterange)
print(s1f_start)
print(s1f_end)

# get actual data from sheets 1b and 1f
s1bd = read_excel(file, sheet = sheet_name_1b, range = cell_rows(5:dim(s1b)[1]))
s1bd_dates = seq(s1b_start, s1b_end, by = 1)
s1bd = s1bd[(1:length(s1bd_dates)),(1:4)]
s1bd$Date = s1bd_dates

s1fd = read_excel(file, sheet = sheet_name_1f, range = cell_rows(5:dim(s1f)[1]))
s1fd_dates = seq(s1f_start, s1f_end, by = 1)
s1fd = s1fd[(1:length(s1fd_dates)+1),c((1:4),(11:13),(20:22),(29:31),(38:40),(47:49),(56:58),(65:67),(74:76))]
s1fd$Date = s1fd_dates

# column bind data from sheets 1b and 1f to match with format used in s1p_df
if (sum(s1bd$Date != s1fd$Date) != 0){
    stop('Differences in dates recorded in sheets 1b and 1f')
} else {
    s1fd$Date = NULL
    finalts = cbind(s1bd, s1fd)
}

# before binding new data in finalts to old data in s1p_df, remove any old dates 
# in s1p_df which already exist in finalts (new data supersedes old data)
s1p_df = subset(s1p_df, !(date %in% finalts$Date))

# adjust columns in finalts to match those already in s1p_df
finalts$publication_date = rep(publicationdate, dim(finalts)[1])
finalts$date = finalts$Date
colnames(finalts) = colnames(s1p_df)

# bind new data to old
s1p_df    = rbind(s1p_df, finalts)
# s1p_df = finalts
if(0){ # plots to check above data
  plot(s1p_df$date, s1p_df$...2)
  plot(s1p_df$date, s1p_df$...5)
  plot(s1p_df$date, s1p_df$...8)
  plot(s1p_df$date, s1p_df$...11)
  plot(s1p_df$date, s1p_df$...14)
  plot(s1p_df$date, s1p_df$...17)
  plot(s1p_df$date, s1p_df$...20)
  plot(s1p_df$date, s1p_df$...23)
  plot(s1p_df$date, s1p_df$...26)
  plot(s1p_df$date, s1p_df$...29)
}

rm(list=setdiff(ls(), "s1p_df"))

# get PHE regions and codes (used by ONS and REACT)
regions = read.csv("./fitting_data/PHE_England_regions.csv")

# initialise data frames to store resized central, lower, upper bound estimates
rszd_swab_prev = data.frame()
rszd_swab_prev_lb = data.frame()
rszd_swab_prev_ub = data.frame()

# for each date, convert central, lower and upper estimates
for (i in 1:length(s1p_df$date)){
    
    # print progress of loop
    print(paste0('Processing date ', i, ' of ', length(s1p_df$date)))
    
    # get list of central estimates to convert (skipping values for England)
    estimates = as.numeric(c(s1p_df[i,5],   # North East
                             s1p_df[i,8],   # North West 
                             s1p_df[i,11],  # Yorkshire and The Humber
                             s1p_df[i,14],  # East Midlands
                             s1p_df[i,17],  # West Midlands
                             s1p_df[i,20],  # East of England
                             s1p_df[i,23],  # London
                             s1p_df[i,26],  # South East
                             s1p_df[i,29])) # South West
    
    # convert estimates to fit nhs regions
    resized_estimates = ogwrangle(estimates, regions$RGN19CD, "e.reg", 
                                   "e.nhser20", "pop2018", "proportion")

    # store converted central estimates
    rszd_swab_prev = rbind(rszd_swab_prev, resized_estimates$estimates)
    
    # get list of lower bound estimates (skipping values for England)
    lb_estimates = as.numeric(c(s1p_df[i,5+1],   # North East
                                s1p_df[i,8+1],   # North West 
                                s1p_df[i,11+1],  # Yorkshire and The Humber
                                s1p_df[i,14+1],  # East Midlands
                                s1p_df[i,17+1],  # West Midlands
                                s1p_df[i,20+1],  # East of England
                                s1p_df[i,23+1],  # London
                                s1p_df[i,26+1],  # South East
                                s1p_df[i,29+1])) # South West
    
    # convert lower bound estimates to fit nhs regions
    resized_lb_estimates = ogwrangle(lb_estimates, regions$RGN19CD, "e.reg", 
                                      "e.nhser20", "pop2018", "proportion")
    
    # store converted lower bound estimates
    rszd_swab_prev_lb = rbind(rszd_swab_prev_lb, 
                               resized_lb_estimates$lb_estimates)
    
    # get list of upper bound estimates (skipping values for England)
    ub_estimates = as.numeric(c(s1p_df[i,5+2],   # North East
                                s1p_df[i,8+2],   # North West 
                                s1p_df[i,11+2],  # Yorkshire and The Humber
                                s1p_df[i,14+2],  # East Midlands
                                s1p_df[i,17+2],  # West Midlands
                                s1p_df[i,20+2],  # East of England
                                s1p_df[i,23+2],  # London
                                s1p_df[i,26+2],  # South East
                                s1p_df[i,29+2])) # South West
    
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

# create data frame to store resized virus prevalence estimates
ONS_CIS_resized <- data.frame(Date = s1p_df$date, 
                              North.East.and.Yorkshire.central = rszd_swab_prev$`North East and Yorkshire`, 
                              North.East.and.Yorkshire.lb = rszd_swab_prev_lb$`North East and Yorkshire`, 
                              North.East.and.Yorkshire.ub = rszd_swab_prev_ub$`North East and Yorkshire`, 
                              North.West.central = rszd_swab_prev$`North West`, 
                              North.West.lb = rszd_swab_prev_lb$`North West`, 
                              North.West.ub = rszd_swab_prev_ub$`North West`, 
                              Midlands.central = rszd_swab_prev$Midlands, 
                              Midlands.lb = rszd_swab_prev_lb$Midlands, 
                              Midlands.ub = rszd_swab_prev_ub$Midlands, 
                              East.of.England.central = rszd_swab_prev$`East of England`, 
                              East.of.England.lb = rszd_swab_prev_lb$`East of England`, 
                              East.of.England.ub = rszd_swab_prev_ub$`East of England`, 
                              London.central = rszd_swab_prev$London, 
                              London.lb = rszd_swab_prev_lb$London, 
                              London.ub = rszd_swab_prev_ub$London, 
                              South.East.central = rszd_swab_prev$`South East`, 
                              South.East.lb = rszd_swab_prev_lb$`South East`, 
                              South.East.ub = rszd_swab_prev_ub$`South East`, 
                              South.West.central = rszd_swab_prev$`South West`, 
                              South.West.lb = rszd_swab_prev_lb$`South West`, 
                              South.West.ub = rszd_swab_prev_ub$`South West`)

# create list of NHS England regional codes
nhs_regions <- c("E40000009", "E40000010", "E40000008", "E40000007", 
                 "E40000003", "E40000005", "E40000006")
names(nhs_regions) <- ogname(nhs_regions)

# initialise final data frame to store virus prevalence data
virus_prev <- data.frame(NHS.region = NULL, Start.date = NULL, End.date = NULL, 
                         Central.estimate = NULL, Lower.bound = NULL, 
                         Upper.bound = NULL, Test = NULL, Median.mean = NULL, 
                         Min.age = NULL, Max.age = NULL, Data.source = NULL, 
                         N.tests = NULL)

# count number of entries of daily modelled virus prevalence in ONS_CIS_resized
num_entries_ONS <- length(ONS_CIS_resized$Date)

for (i in 1:length(nhs_regions)){
    
    # get appropriate string names for each region's estimates in ONS_CIS_resized
    this_region <- str_replace_all(names(nhs_regions)[i], " ", ".")
    central_estim_name <- paste0(this_region, ".central")
    lb_estim_name <- paste0(this_region, ".lb")
    ub_estim_name <- paste0(this_region, ".ub")
    
    # populate entries for this section of the dataframe
    this_region <- rep(names(nhs_regions)[i], num_entries_ONS)
    start_date <- ONS_CIS_resized$Date
    end_date <- ONS_CIS_resized$Date
    central_estimate <- ONS_CIS_resized[, central_estim_name]
    lower_bound <- ONS_CIS_resized[, lb_estim_name]
    upper_bound <- ONS_CIS_resized[, ub_estim_name]
    test <- rep("", num_entries_ONS)
    median_mean <- rep("", num_entries_ONS)
    min_age <- rep("", num_entries_ONS)
    max_age <- rep("", num_entries_ONS)
    data_source <- rep("ONS-CIS", num_entries_ONS)
    N_tests <- rep("", num_entries_ONS)
    
    # make temporary dataframe to bind to master dataframe
    temp_df <- data.frame(NHS.region       = this_region, 
                          Start.date       = start_date, 
                          End.date         = end_date, 
                          Central.estimate = as.numeric(central_estimate), 
                          Lower.bound      = as.numeric(lower_bound), 
                          Upper.bound      = as.numeric(upper_bound), 
                          Test             = test, 
                          Median.mean      = median_mean, 
                          Min.age          = min_age, 
                          Max.age          = max_age, 
                          Data.source      = data_source, 
                          N.tests          = N_tests)
    
    # bind temporary dataframe to master dataframe
    virus_prev <- rbind(virus_prev, temp_df)
}

# optional: add time series for whole of England (change to 'if (0)' to skip)
if (1){
    
    # populate entries for this section of the dataframe
    this_region <- rep('England', num_entries_ONS)
    start_date <- s1p_df$date
    end_date <- s1p_df$date
    central_estimate <- s1p_df$...2
    lower_bound <- s1p_df$...3
    upper_bound <- s1p_df$...4
    test <- rep("", num_entries_ONS)
    median_mean <- rep("", num_entries_ONS)
    min_age <- rep("", num_entries_ONS)
    max_age <- rep("", num_entries_ONS)
    data_source <- rep("ONS-CIS", num_entries_ONS)
    N_tests <- rep("", num_entries_ONS)
    
    # make temporary dataframe to bind to master dataframe
    temp_df <- data.frame(NHS.region       = this_region, 
                          Start.date       = start_date, 
                          End.date         = end_date, 
                          Central.estimate = as.numeric(central_estimate), 
                          Lower.bound      = as.numeric(lower_bound), 
                          Upper.bound      = as.numeric(upper_bound), 
                          Test             = test, 
                          Median.mean      = median_mean, 
                          Min.age          = min_age, 
                          Max.age          = max_age, 
                          Data.source      = data_source, 
                          N.tests          = N_tests)
    
    # bind temporary dataframe to master dataframe
    virus_prev <- rbind(virus_prev, temp_df)
}

# remove any entries with NA dates
virus_prev = virus_prev[!(is.na(virus_prev$Start.date)),]

# # bind new virus prevalence data with existing virus prevalence data
# oldv = read.csv('./fitting_data/virusprev_nhs_regions_20211217122430.csv')
# 
# # remove dates in oldv that exist in virus_prev
# oldv = subset(oldv, !(as.Date(Start.date) %in% virus_prev$Start.date))
# 
# # make sure old data date format matches new
# oldv$Start.date = as.Date(oldv$Start.date)
# oldv$End.date = as.Date(oldv$End.date)
# 
# # bind old and new data together
# virus_prev = rbind(oldv, virus_prev)
# 
# # sort data in order of region and then data (for human readability)
# virus_prev <- virus_prev[order(virus_prev$NHS.region, virus_prev$Start.date),]

# save all virus prevalence estimates
datetime <- str_replace_all(Sys.time(), "[- :BSTGMT]", "")
write.csv(virus_prev, file=paste0("./fitting_data/virusprev_nhs_regions_", 
                                  datetime, ".csv"), row.names = F)

# optional: plots to check data (change to 'if (1)' to execute)
if (0){
    # plot all regions and England as a whole together
    ggplot(data = virus_prev) +
    geom_ribbon(aes(x = as.Date(Start.date), ymin = as.numeric(Lower.bound), 
                    ymax = as.numeric(Upper.bound), fill = NHS.region, 
                    group = NHS.region), alpha = 0.3) +
    labs(x = "Date", y = "Virus prevalence (%)", colour = "NHS England region", fill = "NHS England region") +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    geom_line(aes(x = as.Date(Start.date), y = as.numeric(Central.estimate), 
                  group = NHS.region, colour = NHS.region), size = 0.8) +
    scale_y_continuous(breaks = seq(0, 12, by = 0.5)) +
    scale_color_colorblind() +
    scale_fill_colorblind()
    
    # plot NHS England regions only (remove England)
    ggplot(data = virus_prev[!(virus_prev$NHS.region == 'England'),]) +
        geom_ribbon(aes(x = as.Date(Start.date), ymin = as.numeric(Lower.bound), 
                        ymax = as.numeric(Upper.bound), fill = NHS.region, 
                        group = NHS.region), alpha = 0.3) +
        labs(x = "Date", y = "Virus prevalence (%)", colour = "NHS England region", fill = "NHS England region") +
        scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
        geom_line(aes(x = as.Date(Start.date), y = as.numeric(Central.estimate), 
                      group = NHS.region, colour = NHS.region), size = 0.8) +
        scale_y_continuous(breaks = seq(0, 12, by = 0.5)) +
        scale_color_colorblind() +
        scale_fill_colorblind()
    
    # plot NHS England regions only and remove London (remove England)
    ggplot(data = virus_prev[!(virus_prev$NHS.region %in% c('England', 'London')),]) +
        geom_ribbon(aes(x = as.Date(Start.date), ymin = as.numeric(Lower.bound), 
                        ymax = as.numeric(Upper.bound), fill = NHS.region, 
                        group = NHS.region), alpha = 0.3) +
        labs(x = "Date", y = "Virus prevalence (%)", colour = "NHS England region", fill = "NHS England region") +
        scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
        geom_line(aes(x = as.Date(Start.date), y = as.numeric(Central.estimate), 
                      group = NHS.region, colour = NHS.region), size = 0.8) +
        scale_y_continuous(breaks = seq(0, 12, by = 0.5)) +
        scale_color_colorblind() +
        scale_fill_colorblind()
    
    # plot each region + England as a whole on separate panels
    p = list()
    for (i in 1:length(unique(virus_prev$NHS.region))){
        this_region = unique(virus_prev$NHS.region)[i]
        p[[i]] = ggplot(data = virus_prev[(virus_prev$NHS.region == this_region),]) +
            geom_ribbon(aes(x = as.Date(Start.date), ymin = as.numeric(Lower.bound), 
                            ymax = as.numeric(Upper.bound), alpha = 0.3)) +
            labs(x = "Date", y = "Virus prevalence (%)", title = this_region) +
            scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
            geom_line(aes(x = as.Date(Start.date), 
                          y = as.numeric(Central.estimate)), size = 0.8) +
            scale_y_continuous(breaks = seq(0, 12, by = 0.5)) +
            theme(legend.position = 'none')
    }
    c1 = plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], nrow = 4)
    c2 = plot_grid(p[[5]], p[[6]], p[[7]], p[[8]], nrow = 4)
    plot_grid(c1, c2, nrow = 1)
}
