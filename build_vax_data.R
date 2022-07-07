# build_vax_data.R processes the COVID-19 immunisations linelist from 
# Public Health England and generates schedules for covidm using assumptions on 
# future vaccine rollout (i.e. vaccine supply, uptake and distribution of 
# vaccine products to different age groups)

############################## CHANGES START HERE ##############################

# specify folder and file names here
path_to_data <- "~/Documents/uk_covid_data_sensitive/"
foldername   <- "vaccinations"
filedate     <- "20220504" # change this to match the date the data was built
ref_date     <- '2022-05-04' # this should be in date format and match filedate
sebsfile     <- "vaccination-2022-05-04.rds"

path_to_newcovidvax  <- "."
uk_covid_data_path   <- paste0(path_to_newcovidvax, "/fitting_data/")
cm_path              <- paste0(path_to_newcovidvax, "/covidm_for_fitting/")

# required packages
library(data.table)
library(stringr)
library(lubridate)
library(dplyr)
library(ggplot2)

# turn off scientific notation
options(scipen=999)

date_fitting <- as.character(today())

############################## CHANGES ~END~ HERE ##############################

# load summary vaccination data built from PHE immunisations linelist by Seb
sv <- readRDS(paste0(path_to_data, foldername, '/', sebsfile))

# checks on the summary vaccination data (optional)
checks <- function(vaxdata){
    print(paste0('Total vaccines administered: ', formatC(sum(vaxdata$vaccinated), big.mark=",")))
    print(paste0('Earliest date: ', min(vaxdata$vaccination_date)))
    print(paste0('Latest date: ', max(vaxdata$vaccination_date)))
    print(paste0('Total dates: ', length(unique(vaxdata$vaccination_date))))
    print(paste0(length(vaxdata$vaccination_date[vaxdata$vaccination_date < '2020-12-08']), ' entries with a date before 8th December 2020 (start of rollout)'))
    print(paste0('Comprising ', 100*length(vaxdata$vaccination_date[vaxdata$vaccination_date < '2020-12-08'])/dim(vaxdata)[1], '% of entries'))
    print(paste0(sum(vaxdata$vaccinated[vaxdata$vaccination_date < '2020-12-08']), ' vaccinated on dates before 8th December 2020 (start of rollout)'))
    print(paste0('Comprising ', 100*sum(vaxdata$vaccinated[vaxdata$vaccination_date < '2020-12-08'])/sum(vaxdata$vaccinated), '% of all vaccinations'))
    print(paste0(length(unique(vaxdata$age_group)), ' age groups: ', paste0(unique(vaxdata$age_group), collapse = ", ")))
    print(paste0(length(unique(na.omit(vaxdata$age_group))), ' non-NA age groups (in order): ', paste0(sort(unique(vaxdata$age_group), by = levels(vaxdata$age_group)), collapse = ", ")))
    print(paste0(sum(vaxdata$vaccinated[is.na(vaxdata$age_group) == TRUE]), ' vaccinations wih NA age group'))
    print(paste0('Comprising ', 100*sum(vaxdata$vaccinated[is.na(vaxdata$age_group) == TRUE])/sum(vaxdata$vaccinated), '% of all vaccinations'))
    print(paste0('Dose numbers from ', min(vaxdata$dose_number), ' to ', max(vaxdata$dose_number)))
    print(paste0('Vaccine products: ', paste0(unique(vaxdata$product), collapse = ", ")))
    print(paste0(sum(vaxdata$vaccinated[is.na(vaxdata$product) == TRUE]), ' vaccinations wih NA product'))
    print(paste0('Comprising ', 100*sum(vaxdata$vaccinated[is.na(vaxdata$product) == TRUE])/sum(vaxdata$vaccinated), '% of all vaccinations'))
    print(paste0(sum(vaxdata$vaccinated[vaxdata$product == 'Other'], na.rm = TRUE), ' vaccinations wih `Other` product'))
    print(paste0('Comprising ', 100*sum(vaxdata$vaccinated[vaxdata$product == 'Other'], na.rm = TRUE)/sum(vaxdata$vaccinated), '% of all vaccinations'))
    print(paste0(length(unique(vaxdata$region_of_residence)), ' regions: ', paste0(unique(vaxdata$region_of_residence), collapse = ", ")))
    print(paste0(sum(vaxdata$vaccinated[is.na(vaxdata$region_of_residence) == TRUE]), ' vaccinations wih NA regions'))
    print(paste0('Comprising ', 100*sum(vaxdata$vaccinated[is.na(vaxdata$region_of_residence) == TRUE])/sum(vaxdata$vaccinated), '% of all vaccinations'))
    p1 = hist(vaxdata$dose_number)
    p1
    p2 = qplot(vaxdata$region_of_residence) + theme(axis.text.x = element_text(angle = 90))
    p2
    p3 = qplot(vaxdata$product) + theme(axis.text.x = element_text(angle = 90))
    p3
    plots = list(p1 = p1, p2 = p2, p3 = p3)
    return(plots)
}

oldsv = sv
# keep first doses only
sv = sv[sv$dose_number == 1,]

plots  = checks(oldsv)
plots$p1
plots$p2
plots$p3

plots2 = checks(sv)
plots2$p1
plots2$p2
plots2$p3

# recode erroneous dates.... (assume all vaccinations were really delivered, but dates were recorded incorrectly)

# firstly recode dates occurring before 2020 as 2020 dates
year(sv$vaccination_date[year(sv$vaccination_date) < 2020]) = 2020
# secondly recode dates occurring after this year as this year's dates
year(sv$vaccination_date[year(sv$vaccination_date) > year(today())]) = year(today())
# recode 2020 dates occurring before 8th December 2020 as 2021 dates
year(sv$vaccination_date[sv$vaccination_date < '2020-12-08']) = 2021
# recode dates occurring after reference_date as dates from the previous year
year(sv$vaccination_date[sv$vaccination_date > ref_date]) = year(sv$vaccination_date[sv$vaccination_date > ref_date])-1
# finally, remove any dates which are still before 8th December 2020 (e.g. this could happen when reference_date is < 8th December 2021)
sv = sv[!(sv$vaccination_date < '2020-12-08'),]

# check dates again
print(paste0('Earliest date: ', min(sv$vaccination_date)))
print(paste0('Latest date: ', max(sv$vaccination_date)))

# set up covidm model to get covidm age groups
datapath = function(x) paste0(uk_covid_data_path, x)
cm_force_rebuild = F
cm_build_verbose = T
cm_version = 3
source(paste0(cm_path, "/R/covidm.R"))
popUK = readRDS(datapath("popNHS.rds"))
matricesUK = readRDS(datapath("matricesNHS.rds"))
cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)
nhs_regions = popUK[, unique(name)]
# NUMBER OF REGIONS TO FIT
N_REG = 12;
# Build parameters for NHS regions ###
params = cm_parameters_SEI3R(nhs_regions[1:N_REG], deterministic = T, 
                             date_start = "2020-01-01", 
                             date_end = date_fitting,
                             dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
                             dIp = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
                             dIs = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
                             dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 0.25)$p)
params = cm_split_matrices_ex_in(params, 15)
agegroups <- params$pop[[1]]$group_names
min_ages <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)
oldest_age <- 115
max_ages <- c(4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,oldest_age)
ages <- data.frame(groups = agegroups, min_ages = min_ages, max_ages = max_ages)

# build existing vaccination data (sv) into format suitable for covidm

# get sequence of vaccination dates, from first to last recorded in vax
dates <- seq(min(sv$vaccination_date), max(sv$vaccination_date), by = 1)
num_dates <- length(dates)

# get list of NHS England (and Unknown!) regions of residence in vax
regions <- unique(sv$region_of_residence)

# get list of vaccines in vax -> keep Pfizer, AstraZeneca and Moderna for model purposes
vaccines <- unique(sv$product)
print(vaccines)
`%!in%` <- Negate(`%in%`)
nonvaccs <- vaccines[vaccines %!in% c('PF', 'AZ', 'MD')]
vaccines <- vaccines[vaccines %in% c('PF', 'AZ', 'MD')]

# redistribute doses with vaccine type in 'nonvaccs' (see above) according to 
# the measured split between PF, AZ and MD vaccine products
unknown_vax <- sv[sv$product %in% nonvaccs,]
num_to_correct <- dim(unknown_vax)[1]

# calculate existing split of AZ / PF / MD
total = sum(sv$vaccinated[sv$product %in% vaccines])
aztot = sum(sv$vaccinated[sv$product %in% c('AZ')])
pftot = sum(sv$vaccinated[sv$product %in% c('PF')])
mdtot = sum(sv$vaccinated[sv$product %in% c('MD')])
vax_split = c(pftot/total, aztot/total, mdtot/total)
if (sum(vax_split)!=1){
    stop("Sum of split vector should equal 1")
}

to_distribute = round(vax_split*num_to_correct)
pfsamp = sample(num_to_correct, size = to_distribute[1], replace = FALSE)
remaining = seq(1,num_to_correct, by = 1)
remaining = remaining[!remaining %in% pfsamp]
azsamp = sample(remaining, size = to_distribute[2], replace = FALSE)

product_type_vector = rep('MD', num_to_correct)
product_type_vector[pfsamp] <- "PF"
product_type_vector[azsamp] <- "AZ"

unknown_vax$product <- product_type_vector

# remove unknown vax from main data frame
sv <- sv[!(sv$product %in% nonvaccs),]
# bind unknown vax entries with imputed vaccine product type to main data frame
sv <- rbind(sv, unknown_vax)

# get list of possible doses
doses <- unique(sv$dose_number)

# recode all entries with "" (i.e. unknown) dose number as first doses
if ("" %in% doses){
    num_unknown_doses <- sum(sv$dose_number == "")
    # get unknown doses in separate data frame
    unknown_doses <- sv[sv$dose_number == "",]
    unknown_doses$dose_number <- rep(1, dim(unknown_doses)[1])
    # remove unknown doses from main data frame
    sv <- sv[!(sv$dose_number == ""),]
    # bind updated entries back to main data frame
    sv <- rbind(sv, unknown_doses)
}
# remove "" entry from doses
doses <- doses[!(doses == "")]

# do the same for any NA doses
if (sum(is.na(doses)) > 0){
    num_unknown_doses <- sum(is.na(sv$dose_number) == TRUE)
    # get unknown doses in separate data frame
    unknown_doses <- sv[is.na(sv$dose_number) == TRUE,]
    unknown_doses$dose_number <- rep(1, dim(unknown_doses)[1])
    # remove unknown doses from main data frame
    sv <- sv[!(is.na(sv$dose_number) == TRUE),]
    # bind updated entries back to main data frame
    sv <- rbind(sv, unknown_doses)
}
# remove "" entry from doses
doses <- doses[!is.na(doses) == TRUE]

# initialise dataframe to store vaccine data for covidm
vax_df <- data.frame(region = NULL, age.group = NULL, ages = NULL, 
                     vaccine = NULL, dose = NULL, date = NULL, number = NULL)

# loop through regions
for (i in 1:length(regions)){
    
    this_region <- regions[i]
    print(paste0('Commencing region ', i, ' of ', length(regions), ': ', this_region))
    
    if (is.na(this_region) == TRUE){
        this_region_vax <- sv[is.na(sv$region_of_residence) == TRUE,]
    } else {
        this_region_vax <- sv[sv$region_of_residence == this_region,]
    }
    
    # loop through age groups + 1 (one extra group for NA age entries)
    for (j in 1:(length(agegroups)+1)){
        
        print(paste0('Commencing age group ', j, ' of ', length(agegroups)+1))
        
        if (j == length(agegroups)+1){
            
            this_age_group <- NA
            this_age_group_vax <- this_region_vax[is.na(this_region_vax$age_group) == T,]
            
            print('Age group is NA ages')
            
        } else {
            
            this_age_group <- agegroups[j]
            this_age_group_vax <- this_region_vax[this_region_vax$age_group == this_age_group,]
            
            print(paste0('Age group is ', this_age_group))
            
        }
        
        # loop through vaccines
        for (k in 1:length(vaccines)){
            
            this_vaccine <- vaccines[k]
            print(paste0('Commencing vaccine product ', k, ' of ', length(vaccines), ': ', this_vaccine))
            this_vaccine_vax <- this_age_group_vax[this_age_group_vax$product == this_vaccine,]
            
            # loop through doses
            for (l in 1:length(doses)){
                
                print(paste0('Commencing dose number ', l, ' of ', length(doses)))
                this_dose <- doses[l]
                this_dose_vax <- this_vaccine_vax[this_vaccine_vax$dose_number == this_dose,]
                
                # initialise vector to store number of vaccines delivered on each date
                deliveries <- rep(0, length(dates))
                
                # loop through dates
                for (m in 1:length(dates)){
                    
                    # calculate number of doses delivered for this region, age group, 
                    # vaccine, dose and date
                    this_date <- dates[m]
                    this_date_vax <- this_dose_vax[this_dose_vax$vaccination_date == this_date,]
                    if (dim(this_date_vax)[1] > 1){
                        deliveries[m] <- sum(this_date_vax$vaccinated)
                    } else if (dim(this_date_vax)[1] == 0){
                        deliveries[m] <- 0
                    } else {
                        deliveries[m] <- this_date_vax$vaccinated
                    }
                    rm(this_date_vax)
                }
                
                # populate dataframe for this region, age group, vaccine, dose and date
                this_df <- data.frame(region = rep(this_region, num_dates),
                                      age.group = rep(j, num_dates),
                                      ages = rep(this_age_group, num_dates), 
                                      vaccine = rep(this_vaccine, num_dates), 
                                      dose = rep(this_dose, num_dates), 
                                      date = dates,
                                      number = deliveries)
                
                rm(this_dose, this_dose_vax, deliveries)
                
                # bind dataframe to master dataframe
                vax_df <- rbind(vax_df, this_df)
                rm(this_df)
            }
            
            rm(this_vaccine, this_vaccine_vax)
        }
        
        rm(this_age_group, this_age_group_vax)
        rm(lower_age, upper_age)
    }
    
    rm(this_region, this_region_vax)
}

# check that total doses in vax_df matches total doses in vax
if (sum(sv$vaccinated) != sum(vax_df$number)){
    warning("Total doses in sv does not match those in processed dataframe")
}

# save processed linelist output as .csv
datetime <- str_replace_all(Sys.time(), "[ :GMTBST-]", "")
write.csv(vax_df, file = paste0(path_to_data, foldername, "/", filedate, "-", 
                                datetime, "-processedvaccinedata.csv"), 
          row.names = FALSE)

if (0){
    vax_df = read.csv('~/Documents/uk_covid_data_sensitive/vaccinations/20220504-20220505174155-processedvaccinedata.csv')
}

# redistribute vaccination entries with 'NA' region and 'NA' age (or both)
# into appropriate NHS England regions and covidm age groups

if (0){ # testing code on small subset of data
    # sample at random 5000 rows from data frame
    idxs = round(runif(5000, 0, dim(vax_df)[1]))
    vax_df = vax_df[idxs,]
}

# rename regions from codes to names
vax_df = as.data.table(vax_df)
vax_df[region == 'E40000010', new_region := 'North West']
vax_df[region == 'E40000008', new_region := 'Midlands']
vax_df[region == 'E40000009', new_region := 'North East and Yorkshire']
vax_df[region == 'E40000007', new_region := 'East of England']
vax_df[region == 'E40000003', new_region := 'London']
vax_df[region == 'E40000006', new_region := 'South West']
vax_df[region == 'E40000005', new_region := 'South East']

vax_df[region == 'E40000010', region := 'North West']
vax_df[region == 'E40000008', region := 'Midlands']
vax_df[region == 'E40000009', region := 'North East and Yorkshire']
vax_df[region == 'E40000007', region := 'East of England']
vax_df[region == 'E40000003', region := 'London']
vax_df[region == 'E40000006', region := 'South West']
vax_df[region == 'E40000005', region := 'South East']

vax_df[is.na(region) == TRUE, new_region := NA]

# next: calculate distribution of vaccines into age groups
age_dist = data.table(age_group = NULL, total = NULL)
for (i in 1:(length(agegroups))){
    group = agegroups[i]
    total = sum(vax_df$number[vax_df$age.group == i])
    new_entry = data.table(age_group = group, 
                           total = total)
    age_dist = rbind(age_dist, new_entry)
}
age_dist$prop = age_dist$total / sum(age_dist$total)
age_dist$age_group = factor(age_dist$age_group, levels = levels(sv$age_group))
age_dist$cumprobs <- cumsum(age_dist$prop)

qplot(as.factor(age_dist$age_group), age_dist$prop) + theme(axis.text.x = element_text(angle = 90))

# remove any entries with NA number
vax_df = vax_df[!(is.na(vax_df$number) == TRUE),]

counter = 0

if (1){
    # distribute vaccines with known region and unknown age
    krua_vaxp <- vax_df[is.na(vax_df$ages) == TRUE,]
    known_regions <- regions[is.na(regions) == FALSE]
    for (i in 1:length(known_regions)){
        this_region = known_regions[i]
        krua_vaxp_tr <- krua_vaxp[krua_vaxp$new_region == this_region,]
        # get non-zero vaccine entries
        krua_vaxp_tr_n0 <- krua_vaxp_tr[is.na(krua_vaxp_tr$number) == FALSE,]
        krua_vaxp_tr_n0 <- krua_vaxp_tr_n0[krua_vaxp_tr_n0$number > 0,]
        if (dim(krua_vaxp_tr_n0)[1] > 0){
            for (o in 1:dim(krua_vaxp_tr_n0)[1]){
                
                num_to_distribute = krua_vaxp_tr_n0$number[o]
                counter = counter + num_to_distribute
                print(paste0('Known region unknown age counter: ', counter))
                
                # loop through number to distribute 
                rand <- runif(num_to_distribute,0,1)
                idx <- rep(1,num_to_distribute)
                # recoded_ages <- rep(NA,num_to_distribute)
                for (m in 1:num_to_distribute){
                    for (l in 1:(length(age_dist$cumprobs)-1)){
                        if (rand[m] > age_dist$cumprobs[l]){
                            idx[m] <- idx[m] + 1
                        }
                    }
                    
                    # now we need to add krua_vaxp_tr_n0$number[m] to the correct age group in the master df
                    vax_df$number[vax_df$new_region == this_region && 
                                      vax_df$age.group == idx[m] &&
                                      vax_df$vaccine == krua_vaxp_tr_n0$vaccine[m] &&
                                      vax_df$dose == krua_vaxp_tr_n0$dose[m] && 
                                      vax_df$date == krua_vaxp_tr_n0$date[m]] =
                        vax_df$number[vax_df$new_region == this_region && 
                                          vax_df$age.group == idx[m] &&
                                          vax_df$vaccine == krua_vaxp_tr_n0$vaccine[m] &&
                                          vax_df$dose == krua_vaxp_tr_n0$dose[m] && 
                                          vax_df$date == krua_vaxp_tr_n0$date[m]] + 1
                }
                
            }
        }
    }
    
    # distribute vaccines with unknown region and known age
    urka_vaxp <- vax_df[is.na(vax_df$region) == TRUE,]
    urka_vaxp <- urka_vaxp[is.na(urka_vaxp$ages) == FALSE,]
    urka_vaxp_n0 <- urka_vaxp[!(is.na(urka_vaxp$number) == TRUE),]
    urka_vaxp_n0 <- urka_vaxp_n0[urka_vaxp_n0$number > 0,]
    # loop through age groups recorded in urka_vaxp_n0
    age_groups_to_loop <- unique(urka_vaxp_n0$ages)
    known_regions = unique(vax_df$new_region)
    known_regions = known_regions[!(is.na(known_regions))]
    if (length(age_groups_to_loop) == 0){
        break
    } else if (length(age_groups_to_loop) >= 1){
        for (i in 1:length(age_groups_to_loop)){
            this_age_group <- age_groups_to_loop[i]
            # calculate population sizes for this age group across NHS England regions
            pop_vec <- NULL
            if (this_age_group == "75+"){
                print(paste0("Unknown region, known age group: ", this_age_group))
                all_pops <- popUK[popUK$age %in% c("75-79","80-84","85-89","90+"),]
                for (j in 1:length(known_regions)){
                    these_pops <- all_pops[all_pops$name == known_regions[j],]
                    these_pops$all <- these_pops$f + these_pops$m
                    pop_vec[j] <- sum(these_pops$all)
                }
            } else {
                print(paste0('Unknown region, known age group: ', this_age_group))
                all_pops <- popUK[popUK$age == this_age_group,]
                for (j in 1:length(known_regions)){
                    these_pops <- all_pops[all_pops$name == known_regions[j],]
                    pop_vec[j] <- these_pops$f + these_pops$m
                }
            }
            # print(pop_vec)
            # normalise population size vector and calculate cumulative sum for allocation
            cs_norm_pop <- cumsum(pop_vec / sum(pop_vec))
            # get vaccines to distribute for this age group
            to_distribute <- urka_vaxp_n0[urka_vaxp_n0$ages == this_age_group,]
            # loop through number of vaccine dates with vaccines to distribute
            for (j in 1:length(to_distribute$date)){
                this_vaccine <- to_distribute$vaccine[j]
                this_dose <- to_distribute$dose[j]
                this_date <- to_distribute$date[j]
                num_to_distribute <- to_distribute$number[j]
                counter = counter + num_to_distribute
                print(paste0('Unknown region, known age counter: ', counter))
                # print(counter)
                # loop through each vaccine of this category that needs to be distributed
                for (k in 1:num_to_distribute){
                    # each dose needs to be randomly allocated to a region
                    rand <- runif(1,0,1)
                    idx <- 1
                    for (l in 1:(length(cs_norm_pop)-1)){
                        if (rand > cs_norm_pop[l]){
                            idx <- idx + 1
                        }
                    }
                    region_to_go <- known_regions[idx]
                    # allocate dose to correct region, age group, vaccine, dose, date
                    row_idx <- which(vax_df$new_region == region_to_go & 
                                         vax_df$ages == this_age_group & 
                                         vax_df$vaccine == this_vaccine & 
                                         vax_df$dose == this_dose & 
                                         vax_df$date == this_date)
                    if (length(row_idx) > 1){
                        stop("More than one location for this vaccine to go: code fix required")
                    }
                    vax_df$number[row_idx] <- vax_df$number[row_idx] + 1
                }
            }
        }
    }
    
    
    # distribute vaccines with unknown region and unknown age
    agegroups <- ages$groups
    urua_vaxp <- vax_df[is.na(vax_df$region) == TRUE,]
    urua_vaxp <- urua_vaxp[is.na(urua_vaxp$ages) == TRUE,]
    urua_vaxp_n0 <- urua_vaxp[urua_vaxp$number > 0,]
    # pre calculate probabilities for regions by age group
    cs_norm_pop_allages <- list()
    cumprobs = age_dist$cumprobs
    for (i in 1:length(agegroups)){
        this_age_group <- agegroups[i]
        # calculate population sizes for this age group across NHS England regions
        pop_vec <- NULL
        if (this_age_group == "75+"){
            all_pops <- popUK[popUK$age %in% c("75-79","80-84","85-89","90+"),]
            for (j in 1:length(known_regions)){
                these_pops <- all_pops[all_pops$name == known_regions[j],]
                these_pops$all <- these_pops$f + these_pops$m
                pop_vec[j] <- sum(these_pops$all)
            }
        } else {
            all_pops <- popUK[popUK$age == this_age_group,]
            for (j in 1:length(known_regions)){
                these_pops <- all_pops[all_pops$name == known_regions[j],]
                pop_vec[j] <- these_pops$f + these_pops$m
            }
        }
        # normalise population size vector and calculate cumulative sum for allocation
        cs_norm_pop <- cumsum(pop_vec / sum(pop_vec))
        # store vector in list at appropriate index
        cs_norm_pop_allages[[i]] <- cs_norm_pop
    }
    # loop through number of entries recorded in urua_vaxp_n0
    for (i in 1:dim(urua_vaxp_n0)[1]){
        this_vaccine <- urua_vaxp_n0$vaccine[i]
        this_dose <- urua_vaxp_n0$dose[i]
        this_date <- urua_vaxp_n0$date[i]
        num_to_allocate <- urua_vaxp_n0$number[i]
        counter = counter + num_to_allocate
        # loop through each dose that needs allocating
        print(paste0('Unknown region, unknown age counter: ', counter))
        for (j in 1:num_to_allocate){
            # randomly select an age group for this dose to be allocated to
            rand <- runif(1,0,1)
            idx <- 1
            for (k in 1:(length(cumprobs))){
                if (rand > cumprobs[k]){
                    idx <- idx + 1
                }
            }
            age_group_to_go <- agegroups[idx]
            cs_norm_pop <- cs_norm_pop_allages[[idx]]
            # select a region for this dose to be allocated to
            rand <- runif(1,0,1)
            idx <- 1
            for (l in 1:(length(cs_norm_pop)-1)){
                if (rand > cs_norm_pop[l]){
                    idx <- idx + 1
                }
            }
            region_to_go <- known_regions[idx]
            # allocate dose to correct region, age group, vaccine, dose, date
            row_idx <- which(vax_df$new_region == region_to_go & 
                                 vax_df$ages == age_group_to_go & 
                                 vax_df$vaccine == this_vaccine & 
                                 vax_df$dose == this_dose & 
                                 vax_df$date == this_date)
            if (length(row_idx) > 1){
                stop("More than one location for this vaccine to go: code fix required")
            }
            vax_df$number[row_idx] <- vax_df$number[row_idx] + 1
        }
    }
    
    # finally, remove entries within vaxp that have region = "Unknown" and age = NA
    final_vaxp <- vax_df[!(is.na(vax_df$region) == TRUE),]
    final_vaxp <- final_vaxp[!((is.na(final_vaxp$ages))==TRUE),]
    # vax_df = final_vaxp
    
    # save processed dataframe as .csv
    datetime <- str_replace_all(Sys.time(), "[ :GMTBST-]", "")
    write.csv(final_vaxp, file = paste0(path_to_data,
                                        foldername, "/", filedate, "-", datetime, 
                                        "-processedvaccinedata-reallocated.csv"), 
              row.names = FALSE)
}

if(0){
    final_vaxp = read.csv('~/Documents/uk_covid_data_sensitive/vaccinations/20220504-20220505195127-processedvaccinedata-reallocated.csv')
}

vax_df_old_copy = vax_df
vax_df = final_vaxp

# calculate actual uptake across England by age
source('./vax_funcs.R')
actuals = calculate_uptake(vax_df, popUK)
first_dose_coverage = actuals[actuals$dose == 1,]
print(first_dose_coverage)

# results below calculated on 5th May 2022 (see coverage column)

# age.group dose       date cum_doses   age     pop     coverage
# 1:         1    1 2022-05-02       168   0-4 3239447 0.0000518607
# 2:         2    1 2022-05-02    167035   5-9 3539458 0.0471922537
# 3:         3    1 2022-05-02   1481081 10-14 3435579 0.4311008421
# 4:         4    1 2022-05-02   2604922 15-19 3115871 0.8360172806
# 5:         5    1 2022-05-02   2928173 20-24 3472522 0.8432410219
# 6:         6    1 2022-05-02   3185664 25-29 3771493 0.8446692066
# 7:         7    1 2022-05-02   3453348 30-34 3824652 0.9029182263
# 8:         8    1 2022-05-02   3420699 35-39 3738209 0.9150636040
# 9:         9    1 2022-05-02   3310558 40-44 3476303 0.9523214749
# 10:        10    1 2022-05-02   3399008 45-49 3638639 0.9341426836
# 11:        11    1 2022-05-02   3760555 50-54 3875351 0.9703779090
# 12:        12    1 2022-05-02   3688624 55-59 3761782 0.9805523021
# 13:        13    1 2022-05-02   3183170 60-64 3196813 0.9957323121
# 14:        14    1 2022-05-02   2718777 65-69 2784300 0.9764669755
# 15:        15    1 2022-05-02   2757870 70-74 2814128 0.9800087274
# 16:        16    1 2022-05-02   4734340   75+ 4865591 0.9730246541

# assume existing first dose uptake/coverage limit for individuals aged 15+
uptake         = rep(0, 16)
uptake[4:16]   = pmin(first_dose_coverage$coverage[4:16],1) # setting existing uptake for ages 15+ (>80% for all age groups)
# assume 80% uptake limit in individuals aged 5-14
uptake[2:3]    = 0.8
uptake

# calculate vaccine supply (first doses) by week over time
vax_supply = first_dose_supply(vax_df)
vax_supply$rollmean = zoo::rollmean(vax_supply$V1, 7, fill = NA)
# plot rolling 7-day average number of first doses
plot(as.Date(vax_supply$date), vax_supply$rollmean, type = 'l')

# calculate mean number of vaccines delivered (PER WEEK) in 2022
supp22 = mean(vax_supply$rollmean[vax_supply$date > '2021-12-31'], na.rm = TRUE)*7
print(supp22)

# calculate done on 5th May
# [1] 69465.87

# calculate vector containing future dose supply: each element corresponds to 
# the total number of vaccine doses for *one week* in England

# assume 150000 first doses available per week for 35 weeks (takes us from 2nd May 2022 to 2nd January 2023)
num_weeks = 35
doses = c(rep(150000, num_weeks))
dates = seq(as.Date(max(vax_df$date))+1, as.Date(max(vax_df$date))+1+((num_weeks*7)-1), by = 1)
ddoses <- NULL
for (i in 1:(length(doses))){
    ddoses[(1+(i-1)*7):(i*7)] <- doses[i]/7
}
last_delivery <- as.Date(max(vax_df$date))
idx <- which (dates == last_delivery)
if (length(idx) > 0){
    dates <- dates[-(1:idx)]
    ddoses <- ddoses[-(1:idx)]
}
final_doses <- NULL
for (i in 1:ceiling(length(ddoses)/7)){
    if (i == ceiling(length(ddoses)/7)){
        final_doses[i] = sum(ddoses[(1+(i-1)*7):length(ddoses)])
    } else {
        final_doses[i] = sum(ddoses[(1+(i-1)*7):(i*7)])
    }
}

# calculate data frame containing the proportion of vaccine product types to be
# administered by age group in covidm (currently AstraZeneca, Pfizer, Moderna)

# calculate proportion of product type delivered so far (for all doses and for first doses only)
props          <- products_delivered(vax_df, age_group_limit = 11) # age_group_limit = 11 corresponds to ages 50+ only 

prop_AZ        <- rep(0,16)    # no AZ for individuals <40 years old
prop_AZ[9:16]  <- 0.6          # 60% AZ for individuals aged 40+
prop_PZ        <- rep(0.75,16) # 75% Pfizer for individuals <40 years old
prop_PZ[9:16]  <- 0.3          # 30% Pfizer for individuals aged 40+
prop_MD        <- rep(0.25,16) # 25% Moderna for individuals <40 years old
prop_MD[9:16]  <- 0.1          # 10% Moderna for individuals aged 40+
products <- data.frame(ages = agegroups, AZ = prop_AZ, PZ = prop_PZ, MD = prop_MD)

# define second_doses as TRUE if doses vector contains both first and second 
# doses and as FALSE if doses vector contains first doses only
second_doses <- FALSE

# generate default schedule

# check input values are suitable
if (length(uptake) != 16){
    stop('Uptake vector should contain 16 values for 16 age groups')
} else if (max(uptake) > 1){
    stop('Values in uptake vector must not exceed 1')
} else if (min(uptake) < 0){
    stop('Values in uptake vector must be non-negative')
} else if (min(doses) < 0){
    stop('Values in doses vector must be non-negative')
} else if (dim(products)[1] != 16){
    stop('Products data frame should contain 16 rows for all age groups')
}

# calculate start and end date for future vaccine schedule
start_date    <- as.Date(max(vax_df$date))+1
end_date      <- start_date + 7*(length(doses))

# make sure dates in vax_df are Date format
vax_df$date = as.Date(vax_df$date)

# remove 'new_region column header from vax_df before next section
vax_df$new_region = NULL

if (second_doses == TRUE){
    print('Second doses = TRUE')
    stop('Need code for second_doses == TRUE -> see vax_funcs.R')
} else if (second_doses == FALSE){
    print('Second doses = FALSE')
    
    # calculate list of dates and number of doses *per day*
    vpdates <- dates
    vpdoses <- ddoses
    
    # calculate normalised population size vector for NHS England regions
    regions <- unique(vax_df$region)
    pop_vec <- NULL
    for (i in 1:length(regions)){
        this_region <- regions[i]
        all_pops <- popUK[popUK$name == this_region,]
        pop_vec[i] <- sum(all_pops$f) + sum(all_pops$m)
    }
    pop_vec <- pop_vec / sum(pop_vec)
    
    # calculate total doses administered by age group, add to vax_df
    # calculate number of individuals left to dose in each age group, 
    # vaccine group and dose group
    age_groups <- unique(vax_df$ages)
    final_df = data.frame(region = NULL, age.group = NULL, ages = NULL, 
                          vaccine = NULL, dose = NULL, date = NULL, 
                          number = NULL, cum_doses = NULL, pop_size = NULL, 
                          max_uptake = NULL, left_to_dose = NULL, 
                          left_to_dose_ut = NULL)
    vaccines <- unique(vax_df$vaccine)
    doses <- unique(vax_df$dose)
    for (i in 1:length(regions)){
        this_region <- regions[i]
        this_region_data <- vax_df[vax_df$region == this_region,]
        this_region_pop <- popUK[popUK$name == this_region,]
        for (j in 1:length(age_groups)){
            this_age_group <- age_groups[j]
            this_data <- this_region_data[this_region_data$ages == this_age_group,]
            
            # note that covidm population sizes are listed in 1000's
            if (this_age_group == "75+"){
                all_pops <- this_region_pop[this_region_pop$age %in% c("75-79","80-84","85-89","90+"),]
                pop_size <- (sum(all_pops$f) + sum(all_pops$m)) * 1000
            } else {
                this_region_popage <- this_region_pop[this_region_pop$age == this_age_group,]
                pop_size <- (this_region_popage$f + this_region_popage$m) * 1000
            }
            
            inner_df = data.frame(region = NULL, age.group = NULL, 
                                  ages = NULL, vaccine = NULL, dose = NULL, 
                                  date = NULL, number = NULL, 
                                  cum_doses = NULL, pop_size = NULL, 
                                  max_uptake = NULL, left_to_dose = NULL, 
                                  left_to_dose_ut = NULL)
            
            # loop through vaccine type and dose number to calculate cumulative doses
            for (k in 1:length(vaccines)){
                this_vax <- vaccines[k]
                this_data_v <- this_data[this_data$vaccine == this_vax,]
                for (l in 1:length(doses)){
                    this_dose <- doses[l]
                    this_data_vd <- this_data_v[this_data_v$dose == this_dose,]
                    this_data_vd$cum_doses <- cumsum(this_data_vd$number)
                    this_data_vd$pop_size <- rep(round(pop_size), length(this_data_vd$cum_doses))
                    this_data_vd$max_uptake <- round(this_data_vd$pop_size*uptake[j])
                    inner_df <- rbind(inner_df, this_data_vd)
                }
            }
            # calculate the number of people left to dose in each age group ACROSS
            # vaccine products (i.e. AZ+PZ+Moderna), for each dose (first and second)
            inner2_df = data.frame(region = NULL, age.group = NULL, ages = NULL, 
                                   vaccine = NULL, dose = NULL, date = NULL, 
                                   number = NULL, cum_doses = NULL, pop_size = NULL, 
                                   max_uptake = NULL, left_to_dose = NULL, 
                                   left_to_dose_ut = NULL)
            for (k in 1:length(doses)){
                this_dose <- doses[k]
                this_dose_data <- inner_df[inner_df$dose == this_dose,]
                dates_to_loop <- unique(this_dose_data$date)
                left_to_dose <- NULL
                left_to_dose_ut <- NULL
                for (l in 1:length(dates_to_loop)){
                    this_date <- dates_to_loop[l]
                    this_date_data <- this_dose_data[this_dose_data$date == this_date,]
                    left_to_dose[l] <- this_date_data$pop_size[1] - sum(this_date_data$cum_doses)
                    left_to_dose_ut[l] <- this_date_data$max_uptake[1] - sum(this_date_data$cum_doses) 
                }
                final_this_dose_date <- cbind(this_dose_data, 
                                              left_to_dose = rep(left_to_dose, length(vaccines)),
                                              left_to_dose_ut = rep(left_to_dose_ut, length(vaccines)))
                inner2_df <- rbind(inner2_df, final_this_dose_date)
            }
            final_df <- rbind(final_df, inner2_df)
        }
    }
    
    # remove dose > 1 information in data frame (not required)
    final_df <- final_df[final_df$dose == 1,]
    final_df$date = as.Date(final_df$date)
    
    # calculate age distribution of vaccines delivered
    
    # get list of age groups in covidm (for England)
    # oldest_age <- max(vax$age, na.rm = TRUE)
    # min_ages <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)
    # max_ages <- c(4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,oldest_age)
    # ages <- data.frame(groups = agegroups, min_ages = min_ages, max_ages = max_ages)
    # breaks = c(0,4.5,9.5,14.5,19.5,24.5,29.5,34.5,39.5,44.5,49.5,54.5,59.5,64.5,69.5,74.5,oldest_age)
    # hist <- hist(vax_df$age, freq = FALSE, breaks)
    # nbins <- length(breaks) - 1
    # binwidths <- diff(breaks)
    agegroup_probs <- age_dist$prop
    
    # loop through future dates to add
    leftover_doses <- 0
    for (i in 1:length(vpdates)){
        
        date_today  <- vpdates[i]
        date_yesterday <- date_today - 1
        doses_today <- vpdoses[i]
        doses_per_region <- doses_today * pop_vec
        
        print(date_today)
        
        # now distribute doses per region to each region
        for (j in 1:length(regions)){
            
            this_region <- regions[j]
            print(this_region)
            these_doses <- doses_per_region[j]
            
            # distribute these_doses into age groups for this_region
            
            # IMPORTANT: we assume that the age distribution of vaccines to be delivered
            # follows the existing age distribution in the vaccines delivered already,
            # starting with the oldest age group to the youngest. Any leftover doses 
            # get carried over to doses for the next age group down, the 
            # next region along, or the following date along, or are otherwise
            # recorded as a leftover dose
            doses_per_group <- round(these_doses * agegroup_probs)
            AZ_doses <- round(doses_per_group*products$AZ)
            PZ_doses <- round(doses_per_group*products$PZ)
            MD_doses <- round(doses_per_group*products$MD)
            
            # get data for this region only
            dfr <- final_df[final_df$region == this_region,]
            
            # check if these doses *can* be allocated (cannot exceed either 100% of
            # population sizes or vaccine uptake limits e.g. 95% of population size)
            # looping from oldest age group to youngest
            for (k in length(doses_per_group):1){
                
                # get data for this age group only
                this_age_group <- age_groups[k]
                print(this_age_group)
                dfra <- dfr[dfr$ages == this_age_group,]
                
                # get data from last date of vaccines already delivered
                first_dose_data <- dfra[dfra$dose == 1,]
                fdAZ_data <- first_dose_data[first_dose_data$vaccine == "AZ",]
                fdAZ_data_lastday <- fdAZ_data[fdAZ_data$date == date_yesterday,]
                fdPZ_data <- first_dose_data[first_dose_data$vaccine == "PF",]
                fdPZ_data_lastday <- fdPZ_data[fdPZ_data$date == date_yesterday,]
                fdMD_data <- first_dose_data[first_dose_data$vaccine == "MD",]
                fdMD_data_lastday <- fdMD_data[fdMD_data$date == date_yesterday,]
                
                AZ_doses_to_deliver <- max(AZ_doses[k],0)
                PZ_doses_to_deliver <- max(PZ_doses[k],0)
                MD_doses_to_deliver <- max(MD_doses[k],0)
                
                # check if product types need redistributing here (e.g. AZ vaccines from older
                # groups can get carried over down to younger groups, but younger groups do not
                # receive AZ vaccines in the real world. We assume vaccine supply is fixed, but
                # we can switch across product types without constraints on supply)
                product_split_this_age_group = products[k,c(2:4)]
                all_doses_to_deliver = sum(AZ_doses_to_deliver,
                                           PZ_doses_to_deliver,
                                           MD_doses_to_deliver)
                AZ_doses_to_deliver = round(all_doses_to_deliver*product_split_this_age_group$AZ)
                PZ_doses_to_deliver = round(all_doses_to_deliver*product_split_this_age_group$PZ)
                MD_doses_to_deliver = round(all_doses_to_deliver*product_split_this_age_group$MD)
                
                # calculate first and second doses to (attempt to) deliver
                fdtd_AZ <- AZ_doses_to_deliver
                fdtd_PZ <- PZ_doses_to_deliver
                fdtd_MD <- MD_doses_to_deliver
                
                # limit doses to the uptake threshold of each age group
                
                # deliver as many as possible first doses for AZ, Pfizer, and Moderna,
                # limiting the number of doses to uptake threshold
                if (fdAZ_data_lastday$left_to_dose_ut > 0){
                    if (fdtd_AZ + fdtd_PZ + fdtd_MD > fdAZ_data_lastday$left_to_dose_ut){
                        # split left_to_dose_ut into proportions of AZ, Pfizer and Moderna
                        final_fdtd_AZ <- round(products$AZ[k] *  fdAZ_data_lastday$left_to_dose_ut)
                        final_fdtd_PZ <- round(products$PZ[k] *  fdAZ_data_lastday$left_to_dose_ut)
                        final_fdtd_MD <- round(products$MD[k] *  fdAZ_data_lastday$left_to_dose_ut)
                        # shift all remaining second doses to...
                        if (k > 1){
                            # next age group down
                            AZ_doses[k-1] <- AZ_doses[k-1] + (fdtd_AZ - final_fdtd_AZ)
                            PZ_doses[k-1] <- PZ_doses[k-1] + (fdtd_PZ - final_fdtd_PZ)
                            MD_doses[k-1] <- MD_doses[k-1] + (fdtd_MD - final_fdtd_MD)
                        } else if (j < length(regions)){
                            # next region along
                            doses_per_region[j+1] <- doses_per_region[j+1] + 
                                (fdtd_AZ - final_fdtd_AZ) +
                                (fdtd_PZ - final_fdtd_PZ) +
                                (fdtd_MD - final_fdtd_MD)
                        } else if (i < length(vpdates)){
                            # next date along
                            vpdoses[i+1] <- vpdoses[i+1] + 
                                (fdtd_AZ - final_fdtd_AZ) +
                                (fdtd_PZ - final_fdtd_PZ) +
                                (fdtd_MD - final_fdtd_MD)
                        } else {
                            # or record leftover doses
                            leftover_doses <- leftover_doses + 
                                (fdtd_AZ - final_fdtd_AZ) +
                                (fdtd_PZ - final_fdtd_PZ) +
                                (fdtd_MD - final_fdtd_MD)
                        }
                        # deliver final first doses for AZ, Pfizer and Moderna
                        new_entries <- data.frame(region = rep(this_region, 3), 
                                                  age.group = rep(k, 3), 
                                                  ages = rep(this_age_group, 3), 
                                                  vaccine = c("AZ", "PF", "MD"),
                                                  dose = rep(1, 3), 
                                                  date = rep(as.Date(date_today), 3),
                                                  number = c(final_fdtd_AZ, final_fdtd_PZ, final_fdtd_MD),
                                                  cum_doses = c(fdAZ_data_lastday$cum_doses + final_fdtd_AZ, 
                                                                fdPZ_data_lastday$cum_doses + final_fdtd_PZ,
                                                                fdMD_data_lastday$cum_doses + final_fdtd_MD),
                                                  pop_size = rep(fdAZ_data_lastday$pop_size, 3),
                                                  max_uptake = rep(fdAZ_data_lastday$max_uptake, 3),
                                                  left_to_dose = rep(fdAZ_data_lastday$left_to_dose - fdAZ_data_lastday$left_to_dose_ut, 3),
                                                  left_to_dose_ut = rep(0, 3))
                        # bind new entries to final data frame
                        final_df <- rbind(final_df, new_entries)
                    } else { # fdtd_AZ + fdtd_PZ + fdtd_MD <= fdAZ_data_lastday$left_to_dose_ut
                        # deliver first doses for AZ, Pfizer and Moderna
                        new_entries <- data.frame(region = rep(this_region, 3), 
                                                  age.group = rep(k, 3), 
                                                  ages = rep(this_age_group, 3), 
                                                  vaccine = c("AZ", "PF", "MD"),
                                                  dose = rep(1, 3), 
                                                  date = rep(as.Date(date_today), 3),
                                                  number = c(fdtd_AZ, fdtd_PZ, fdtd_MD),
                                                  cum_doses = c(fdAZ_data_lastday$cum_doses + fdtd_AZ, 
                                                                fdPZ_data_lastday$cum_doses + fdtd_PZ,
                                                                fdMD_data_lastday$cum_doses + fdtd_MD),
                                                  pop_size = rep(fdAZ_data_lastday$pop_size, 3),
                                                  max_uptake = rep(fdAZ_data_lastday$max_uptake, 3),
                                                  left_to_dose = rep(fdAZ_data_lastday$left_to_dose - (fdtd_AZ+fdtd_PZ+fdtd_MD), 3),
                                                  left_to_dose_ut = rep(fdAZ_data_lastday$left_to_dose_ut - (fdtd_AZ+fdtd_PZ+fdtd_MD), 3))
                        # bind new entries to final data frame
                        final_df <- rbind(final_df, new_entries)
                    }
                } else { 
                    # fdAZ_data_lastday$left_to_dose_ut <= 0 
                    # record zero first doses getting delivered for AZ, Pfizer and Moderna
                    new_entries <- data.frame(region = rep(this_region, 3), 
                                              age.group = rep(k, 3), 
                                              ages = rep(this_age_group, 3), 
                                              vaccine = c("AZ", "PF", "MD"),
                                              dose = rep(1, 3), 
                                              date = rep(as.Date(date_today), 3),
                                              number = rep(0, 3),
                                              cum_doses = c(fdAZ_data_lastday$cum_doses, 
                                                            fdPZ_data_lastday$cum_doses,
                                                            fdMD_data_lastday$cum_doses),
                                              pop_size = rep(fdAZ_data_lastday$pop_size, 3),
                                              max_uptake = rep(fdAZ_data_lastday$max_uptake, 3),
                                              left_to_dose = rep(fdAZ_data_lastday$left_to_dose, 3),
                                              left_to_dose_ut = rep(fdAZ_data_lastday$left_to_dose_ut, 3))
                    # bind new entries to final data frame
                    final_df <- rbind(final_df, new_entries)
                    
                    
                    # shift all remaining second doses to...
                    if (k > 1){
                        # next age group down
                        AZ_doses[k-1] <- AZ_doses[k-1] + fdtd_AZ
                        PZ_doses[k-1] <- PZ_doses[k-1] + fdtd_PZ
                        MD_doses[k-1] <- MD_doses[k-1] + fdtd_MD
                    } else if (j < length(regions)){
                        # next region along
                        doses_per_region[j+1] <- doses_per_region[j+1] + fdtd_AZ + fdtd_PZ + fdtd_MD
                    } else if (i < length(vpdates)){
                        # next date along
                        vpdoses[i+1] <- vpdoses[i+1] + fdtd_AZ + fdtd_PZ + fdtd_MD
                    } else {
                        # or record leftover doses
                        leftover_doses <- leftover_doses + fdtd_AZ + fdtd_PZ + fdtd_MD
                    }
                }
            }
        }
    }
    
    # print leftover doses
    print(leftover_doses)
    
    # sort final_df data frame rows by region, age group, vaccine, date
    vax_schedule <- final_df[order(final_df$region, final_df$age.group, 
                                   final_df$vaccine, final_df$date),]
    
}

# save schedule containing data on delivered and projected deliveries in future
datetime <- str_replace_all(Sys.time(), "[ :GMTBST-]", "")
write.csv(vax_schedule, file = paste0(path_to_data,
                                      foldername, "/", filedate, "-", datetime, 
                                      "-delivered_and_projected_vaccines.csv"), 
          row.names = FALSE)

# convert schedule

# add lagged dates
first_doses <- vax_schedule[vax_schedule$dose == 1,]
second_doses <- vax_schedule[vax_schedule$dose == 2,]
first_doses$lag_date <- as.Date(as.Date(first_doses$date) + 28)
second_doses$lag_date <- as.Date(as.Date(second_doses$date) + 14)
vax_schedule2 <- rbind(first_doses, second_doses)

# NHS England popset
popset = c(1, 3, 4, 5, 6, 9, 10)
covidm_regions = NULL
for (i in popset){
    covidm_regions[i] <- params$pop[[i]]$name  
}
covidm_regions <- covidm_regions[!is.na(covidm_regions)]


lag_dates <- as.Date(sort(unique(vax_schedule2$lag_date)))
va1 <- list()
va2 <- list()
vb1 <- list()
vb2 <- list()
ret <- list()

for (i in 1:length(covidm_regions)){
    this_region <- covidm_regions[i]
    print(this_region)
    this_schedu <- vax_schedule2[vax_schedule2$region == this_region,]
    for (j in 1:length(lag_dates)){
        
        this_date <- lag_dates[j]
        print(this_date)
        this_schedu_date <- this_schedu[this_schedu$lag_date == this_date,]
        
        AZ <- this_schedu_date[this_schedu_date$vaccine == "AZ",]
        PZ <- this_schedu_date[this_schedu_date$vaccine == "PF",]
        MD <- this_schedu_date[this_schedu_date$vaccine == "MD",]
        AZ1 <- AZ[AZ$dose == 1,]
        AZ1 <- AZ1[order(AZ1$age.group),]
        AZ2 <- AZ[AZ$dose == 2,]
        AZ2 <- AZ2[order(AZ2$age.group),]
        PZ1 <- PZ[PZ$dose == 1,]
        PZ1 <- PZ1[order(PZ1$age.group),]
        PZ2 <- PZ[PZ$dose == 2,]
        PZ2 <- PZ2[order(PZ2$age.group),]
        MD1 <- MD[MD$dose == 1,]
        MD1 <- MD1[order(MD1$age.group),]
        MD2 <- MD[MD$dose == 2,]
        MD2 <- MD2[order(MD2$age.group),]
        
        if (dim(AZ1)[1] == 0){
            va1[[j]] <- rep(0, 16)
            print('AZ dose 1:')
            print(paste0('All ', 0))
        } else {
            va1[[j]] <- AZ1$number
            print('AZ dose 1:')
            print(AZ1$number)
        }
        if (dim(AZ2)[1] == 0){
            va2[[j]] <- rep(0,16)
            print('AZ dose 2:')
            print(paste0('All ', 0))
        } else {
            va2[[j]] <- AZ2$number
            print('AZ dose 2:')
            print(AZ2$number)
        }
        if (dim(PZ1)[1] == 0 && dim(MD1)[1] == 0){
            vb1[[j]] <- rep(0, 16)
            print('Pfizer & Moderna dose 1:')
            print(paste0('All ', 0))
        } else {
            vb1[[j]] <- PZ1$number + MD1$number
            
            print('Pfizer & Moderna dose 1:')
            print(PZ1$number + MD1$number)
        }
        if (dim(PZ2)[1] == 0 && dim(MD2)[1] == 0){
            vb2[[j]] <- rep(0,16)
            print('Pfizer & Moderna dose 2:')
            print(paste0('All ', 0))
        } else {
            vb2[[j]] <- PZ2$number + MD2$number
            print('Pfizer & Moderna dose 2:')
            print(PZ2$number + MD2$number)
        }
        
    }
    idx <- popset[i]
    ret[[idx]] = list(
        vt = lag_dates,
        va1 = va1,
        va2 = va2,
        vb1 = vb1,
        vb2 = vb2
    )
}

# save final schedule as .rds file
datetime <- str_replace_all(Sys.time(), "[ :GMTBST-]", "")
saveRDS(ret, paste0(path_to_data, foldername, "/vax-covidm", datetime, ".rds"))
saveRDS(ret, paste0(uk_covid_data_path, "vax-covidm", datetime, ".rds"))

## create alternative vaccine schedule for paper (50% uptake limit for kids)

name = '5plus_50percent'

# assume existing first dose uptake/coverage limit for individuals aged 15+
uptake         = rep(0, 16)
uptake[4:16]   = pmin(first_dose_coverage$coverage[4:16],1) # setting existing uptake for ages 15+ (>80% for all age groups)
# assume 50% uptake limit in individuals aged 5-14
uptake[2:3]    = 0.5
uptake

# calculate vaccine supply (first doses) by week over time
vax_supply = first_dose_supply(vax_df)
vax_supply$rollmean = zoo::rollmean(vax_supply$V1, 7, fill = NA)
# plot rolling 7-day average number of first doses
plot(as.Date(vax_supply$date), vax_supply$rollmean, type = 'l')

# calculate mean number of vaccines delivered (PER WEEK) in 2022
supp22 = mean(vax_supply$rollmean[vax_supply$date > '2021-12-31'], na.rm = TRUE)*7
print(supp22)

# calculate done on 29th March
# [1] 72799.49

# calculate vector containing future dose supply: each element corresponds to 
# the total number of vaccine doses for *one week* in England

# assume 150000 first doses available per week for 35 weeks
doses = c(rep(150000, num_weeks))
dates = seq(as.Date(max(vax_df$date))+1, as.Date(max(vax_df$date))+1+((num_weeks*7)-1), by = 1)
ddoses <- NULL
for (i in 1:(length(doses))){
    ddoses[(1+(i-1)*7):(i*7)] <- doses[i]/7
}
last_delivery <- as.Date(max(vax_df$date))
idx <- which (dates == last_delivery)
if (length(idx) > 0){
    dates <- dates[-(1:idx)]
    ddoses <- ddoses[-(1:idx)]
}
final_doses <- NULL
for (i in 1:ceiling(length(ddoses)/7)){
    if (i == ceiling(length(ddoses)/7)){
        final_doses[i] = sum(ddoses[(1+(i-1)*7):length(ddoses)])
    } else {
        final_doses[i] = sum(ddoses[(1+(i-1)*7):(i*7)])
    }
}

# calculate data frame containing the proportion of vaccine product types to be
# administered by age group in covidm (currently AstraZeneca, Pfizer, Moderna)

# calculate proportion of product type delivered so far (for all doses and for first doses only)
props          <- products_delivered(vax_df, age_group_limit = 11) # age_group_limit = 11 corresponds to ages 50+ only 

prop_AZ        <- rep(0,16)    # no AZ for individuals <40 years old
prop_AZ[9:16]  <- 0.6          # 60% AZ for individuals aged 40+
prop_PZ        <- rep(0.75,16) # 75% Pfizer for individuals <40 years old
prop_PZ[9:16]  <- 0.3          # 30% Pfizer for individuals aged 40+
prop_MD        <- rep(0.25,16) # 25% Moderna for individuals <40 years old
prop_MD[9:16]  <- 0.1          # 10% Moderna for individuals aged 40+
products <- data.frame(ages = agegroups, AZ = prop_AZ, PZ = prop_PZ, MD = prop_MD)

# define second_doses as TRUE if doses vector contains both first and second 
# doses and as FALSE if doses vector contains first doses only
second_doses <- FALSE

# generate default schedule

# check input values are suitable
if (length(uptake) != 16){
    stop('Uptake vector should contain 16 values for 16 age groups')
} else if (max(uptake) > 1){
    stop('Values in uptake vector must not exceed 1')
} else if (min(uptake) < 0){
    stop('Values in uptake vector must be non-negative')
} else if (min(doses) < 0){
    stop('Values in doses vector must be non-negative')
} else if (dim(products)[1] != 16){
    stop('Products data frame should contain 16 rows for all age groups')
}

# calculate start and end date for future vaccine schedule
start_date    <- as.Date(max(vax_df$date))+1
end_date      <- start_date + 7*(length(doses))

# make sure dates in vax_df are Date format
vax_df$date = as.Date(vax_df$date)

# remove 'new_region column header from vax_df before next section
vax_df$new_region = NULL

if (second_doses == TRUE){
    print('Second doses = TRUE')
    stop('Need code for second_doses == TRUE -> see vax_funcs.R')
} else if (second_doses == FALSE){
    print('Second doses = FALSE')
    
    # calculate list of dates and number of doses *per day*
    vpdates <- dates
    vpdoses <- ddoses
    
    # calculate normalised population size vector for NHS England regions
    regions <- unique(vax_df$region)
    pop_vec <- NULL
    for (i in 1:length(regions)){
        this_region <- regions[i]
        all_pops <- popUK[popUK$name == this_region,]
        pop_vec[i] <- sum(all_pops$f) + sum(all_pops$m)
    }
    pop_vec <- pop_vec / sum(pop_vec)
    
    # calculate total doses administered by age group, add to vax_df
    # calculate number of individuals left to dose in each age group, 
    # vaccine group and dose group
    age_groups <- unique(vax_df$ages)
    final_df = data.frame(region = NULL, age.group = NULL, ages = NULL, 
                          vaccine = NULL, dose = NULL, date = NULL, 
                          number = NULL, cum_doses = NULL, pop_size = NULL, 
                          max_uptake = NULL, left_to_dose = NULL, 
                          left_to_dose_ut = NULL)
    vaccines <- unique(vax_df$vaccine)
    doses <- unique(vax_df$dose)
    for (i in 1:length(regions)){
        this_region <- regions[i]
        this_region_data <- vax_df[vax_df$region == this_region,]
        this_region_pop <- popUK[popUK$name == this_region,]
        for (j in 1:length(age_groups)){
            this_age_group <- age_groups[j]
            this_data <- this_region_data[this_region_data$ages == this_age_group,]
            
            # note that covidm population sizes are listed in 1000's
            if (this_age_group == "75+"){
                all_pops <- this_region_pop[this_region_pop$age %in% c("75-79","80-84","85-89","90+"),]
                pop_size <- (sum(all_pops$f) + sum(all_pops$m)) * 1000
            } else {
                this_region_popage <- this_region_pop[this_region_pop$age == this_age_group,]
                pop_size <- (this_region_popage$f + this_region_popage$m) * 1000
            }
            
            inner_df = data.frame(region = NULL, age.group = NULL, 
                                  ages = NULL, vaccine = NULL, dose = NULL, 
                                  date = NULL, number = NULL, 
                                  cum_doses = NULL, pop_size = NULL, 
                                  max_uptake = NULL, left_to_dose = NULL, 
                                  left_to_dose_ut = NULL)
            
            # loop through vaccine type and dose number to calculate cumulative doses
            for (k in 1:length(vaccines)){
                this_vax <- vaccines[k]
                this_data_v <- this_data[this_data$vaccine == this_vax,]
                for (l in 1:length(doses)){
                    this_dose <- doses[l]
                    this_data_vd <- this_data_v[this_data_v$dose == this_dose,]
                    this_data_vd$cum_doses <- cumsum(this_data_vd$number)
                    this_data_vd$pop_size <- rep(round(pop_size), length(this_data_vd$cum_doses))
                    this_data_vd$max_uptake <- round(this_data_vd$pop_size*uptake[j])
                    inner_df <- rbind(inner_df, this_data_vd)
                }
            }
            # calculate the number of people left to dose in each age group ACROSS
            # vaccine products (i.e. AZ+PZ+Moderna), for each dose (first and second)
            inner2_df = data.frame(region = NULL, age.group = NULL, ages = NULL, 
                                   vaccine = NULL, dose = NULL, date = NULL, 
                                   number = NULL, cum_doses = NULL, pop_size = NULL, 
                                   max_uptake = NULL, left_to_dose = NULL, 
                                   left_to_dose_ut = NULL)
            for (k in 1:length(doses)){
                this_dose <- doses[k]
                this_dose_data <- inner_df[inner_df$dose == this_dose,]
                dates_to_loop <- unique(this_dose_data$date)
                left_to_dose <- NULL
                left_to_dose_ut <- NULL
                for (l in 1:length(dates_to_loop)){
                    this_date <- dates_to_loop[l]
                    this_date_data <- this_dose_data[this_dose_data$date == this_date,]
                    left_to_dose[l] <- this_date_data$pop_size[1] - sum(this_date_data$cum_doses)
                    left_to_dose_ut[l] <- this_date_data$max_uptake[1] - sum(this_date_data$cum_doses) 
                }
                final_this_dose_date <- cbind(this_dose_data, 
                                              left_to_dose = rep(left_to_dose, length(vaccines)),
                                              left_to_dose_ut = rep(left_to_dose_ut, length(vaccines)))
                inner2_df <- rbind(inner2_df, final_this_dose_date)
            }
            final_df <- rbind(final_df, inner2_df)
        }
    }
    
    # remove dose > 1 information in data frame (not required)
    final_df <- final_df[final_df$dose == 1,]
    final_df$date = as.Date(final_df$date)
    
    # calculate age distribution of vaccines delivered
    
    # get list of age groups in covidm (for England)
    # oldest_age <- max(vax$age, na.rm = TRUE)
    # min_ages <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)
    # max_ages <- c(4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,oldest_age)
    # ages <- data.frame(groups = agegroups, min_ages = min_ages, max_ages = max_ages)
    # breaks = c(0,4.5,9.5,14.5,19.5,24.5,29.5,34.5,39.5,44.5,49.5,54.5,59.5,64.5,69.5,74.5,oldest_age)
    # hist <- hist(vax_df$age, freq = FALSE, breaks)
    # nbins <- length(breaks) - 1
    # binwidths <- diff(breaks)
    agegroup_probs <- age_dist$prop
    
    # loop through future dates to add
    leftover_doses <- 0
    for (i in 1:length(vpdates)){
        
        date_today  <- vpdates[i]
        date_yesterday <- date_today - 1
        doses_today <- vpdoses[i]
        doses_per_region <- doses_today * pop_vec
        
        print(date_today)
        
        # now distribute doses per region to each region
        for (j in 1:length(regions)){
            
            this_region <- regions[j]
            print(this_region)
            these_doses <- doses_per_region[j]
            
            # distribute these_doses into age groups for this_region
            
            # IMPORTANT: we assume that the age distribution of vaccines to be delivered
            # follows the existing age distribution in the vaccines delivered already,
            # starting with the oldest age group to the youngest. Any leftover doses 
            # get carried over to doses for the next age group down, the 
            # next region along, or the following date along, or are otherwise
            # recorded as a leftover dose
            doses_per_group <- round(these_doses * agegroup_probs)
            AZ_doses <- round(doses_per_group*products$AZ)
            PZ_doses <- round(doses_per_group*products$PZ)
            MD_doses <- round(doses_per_group*products$MD)
            
            # get data for this region only
            dfr <- final_df[final_df$region == this_region,]
            
            # check if these doses *can* be allocated (cannot exceed either 100% of
            # population sizes or vaccine uptake limits e.g. 95% of population size)
            # looping from oldest age group to youngest
            for (k in length(doses_per_group):1){
                
                # get data for this age group only
                this_age_group <- age_groups[k]
                print(this_age_group)
                dfra <- dfr[dfr$ages == this_age_group,]
                
                # get data from last date of vaccines already delivered
                first_dose_data <- dfra[dfra$dose == 1,]
                fdAZ_data <- first_dose_data[first_dose_data$vaccine == "AZ",]
                fdAZ_data_lastday <- fdAZ_data[fdAZ_data$date == date_yesterday,]
                fdPZ_data <- first_dose_data[first_dose_data$vaccine == "PF",]
                fdPZ_data_lastday <- fdPZ_data[fdPZ_data$date == date_yesterday,]
                fdMD_data <- first_dose_data[first_dose_data$vaccine == "MD",]
                fdMD_data_lastday <- fdMD_data[fdMD_data$date == date_yesterday,]
                
                AZ_doses_to_deliver <- max(AZ_doses[k],0)
                PZ_doses_to_deliver <- max(PZ_doses[k],0)
                MD_doses_to_deliver <- max(MD_doses[k],0)
                
                # check if product types need redistributing here (e.g. AZ vaccines from older
                # groups can get carried over down to younger groups, but younger groups do not
                # receive AZ vaccines in the real world. We assume vaccine supply is fixed, but
                # we can switch across product types without constraints on supply)
                product_split_this_age_group = products[k,c(2:4)]
                all_doses_to_deliver = sum(AZ_doses_to_deliver,
                                           PZ_doses_to_deliver,
                                           MD_doses_to_deliver)
                AZ_doses_to_deliver = round(all_doses_to_deliver*product_split_this_age_group$AZ)
                PZ_doses_to_deliver = round(all_doses_to_deliver*product_split_this_age_group$PZ)
                MD_doses_to_deliver = round(all_doses_to_deliver*product_split_this_age_group$MD)
                
                # calculate first and second doses to (attempt to) deliver
                fdtd_AZ <- AZ_doses_to_deliver
                fdtd_PZ <- PZ_doses_to_deliver
                fdtd_MD <- MD_doses_to_deliver
                
                # limit doses to the uptake threshold of each age group
                
                # deliver as many as possible first doses for AZ, Pfizer, and Moderna,
                # limiting the number of doses to uptake threshold
                if (fdAZ_data_lastday$left_to_dose_ut > 0){
                    if (fdtd_AZ + fdtd_PZ + fdtd_MD > fdAZ_data_lastday$left_to_dose_ut){
                        # split left_to_dose_ut into proportions of AZ, Pfizer and Moderna
                        final_fdtd_AZ <- round(products$AZ[k] *  fdAZ_data_lastday$left_to_dose_ut)
                        final_fdtd_PZ <- round(products$PZ[k] *  fdAZ_data_lastday$left_to_dose_ut)
                        final_fdtd_MD <- round(products$MD[k] *  fdAZ_data_lastday$left_to_dose_ut)
                        # shift all remaining second doses to...
                        if (k > 1){
                            # next age group down
                            AZ_doses[k-1] <- AZ_doses[k-1] + (fdtd_AZ - final_fdtd_AZ)
                            PZ_doses[k-1] <- PZ_doses[k-1] + (fdtd_PZ - final_fdtd_PZ)
                            MD_doses[k-1] <- MD_doses[k-1] + (fdtd_MD - final_fdtd_MD)
                        } else if (j < length(regions)){
                            # next region along
                            doses_per_region[j+1] <- doses_per_region[j+1] + 
                                (fdtd_AZ - final_fdtd_AZ) +
                                (fdtd_PZ - final_fdtd_PZ) +
                                (fdtd_MD - final_fdtd_MD)
                        } else if (i < length(vpdates)){
                            # next date along
                            vpdoses[i+1] <- vpdoses[i+1] + 
                                (fdtd_AZ - final_fdtd_AZ) +
                                (fdtd_PZ - final_fdtd_PZ) +
                                (fdtd_MD - final_fdtd_MD)
                        } else {
                            # or record leftover doses
                            leftover_doses <- leftover_doses + 
                                (fdtd_AZ - final_fdtd_AZ) +
                                (fdtd_PZ - final_fdtd_PZ) +
                                (fdtd_MD - final_fdtd_MD)
                        }
                        # deliver final first doses for AZ, Pfizer and Moderna
                        new_entries <- data.frame(region = rep(this_region, 3), 
                                                  age.group = rep(k, 3), 
                                                  ages = rep(this_age_group, 3), 
                                                  vaccine = c("AZ", "PF", "MD"),
                                                  dose = rep(1, 3), 
                                                  date = rep(as.Date(date_today), 3),
                                                  number = c(final_fdtd_AZ, final_fdtd_PZ, final_fdtd_MD),
                                                  cum_doses = c(fdAZ_data_lastday$cum_doses + final_fdtd_AZ, 
                                                                fdPZ_data_lastday$cum_doses + final_fdtd_PZ,
                                                                fdMD_data_lastday$cum_doses + final_fdtd_MD),
                                                  pop_size = rep(fdAZ_data_lastday$pop_size, 3),
                                                  max_uptake = rep(fdAZ_data_lastday$max_uptake, 3),
                                                  left_to_dose = rep(fdAZ_data_lastday$left_to_dose - fdAZ_data_lastday$left_to_dose_ut, 3),
                                                  left_to_dose_ut = rep(0, 3))
                        # bind new entries to final data frame
                        final_df <- rbind(final_df, new_entries)
                    } else { # fdtd_AZ + fdtd_PZ + fdtd_MD <= fdAZ_data_lastday$left_to_dose_ut
                        # deliver first doses for AZ, Pfizer and Moderna
                        new_entries <- data.frame(region = rep(this_region, 3), 
                                                  age.group = rep(k, 3), 
                                                  ages = rep(this_age_group, 3), 
                                                  vaccine = c("AZ", "PF", "MD"),
                                                  dose = rep(1, 3), 
                                                  date = rep(as.Date(date_today), 3),
                                                  number = c(fdtd_AZ, fdtd_PZ, fdtd_MD),
                                                  cum_doses = c(fdAZ_data_lastday$cum_doses + fdtd_AZ, 
                                                                fdPZ_data_lastday$cum_doses + fdtd_PZ,
                                                                fdMD_data_lastday$cum_doses + fdtd_MD),
                                                  pop_size = rep(fdAZ_data_lastday$pop_size, 3),
                                                  max_uptake = rep(fdAZ_data_lastday$max_uptake, 3),
                                                  left_to_dose = rep(fdAZ_data_lastday$left_to_dose - (fdtd_AZ+fdtd_PZ+fdtd_MD), 3),
                                                  left_to_dose_ut = rep(fdAZ_data_lastday$left_to_dose_ut - (fdtd_AZ+fdtd_PZ+fdtd_MD), 3))
                        # bind new entries to final data frame
                        final_df <- rbind(final_df, new_entries)
                    }
                } else { 
                    # fdAZ_data_lastday$left_to_dose_ut <= 0 
                    # record zero first doses getting delivered for AZ, Pfizer and Moderna
                    new_entries <- data.frame(region = rep(this_region, 3), 
                                              age.group = rep(k, 3), 
                                              ages = rep(this_age_group, 3), 
                                              vaccine = c("AZ", "PF", "MD"),
                                              dose = rep(1, 3), 
                                              date = rep(as.Date(date_today), 3),
                                              number = rep(0, 3),
                                              cum_doses = c(fdAZ_data_lastday$cum_doses, 
                                                            fdPZ_data_lastday$cum_doses,
                                                            fdMD_data_lastday$cum_doses),
                                              pop_size = rep(fdAZ_data_lastday$pop_size, 3),
                                              max_uptake = rep(fdAZ_data_lastday$max_uptake, 3),
                                              left_to_dose = rep(fdAZ_data_lastday$left_to_dose, 3),
                                              left_to_dose_ut = rep(fdAZ_data_lastday$left_to_dose_ut, 3))
                    # bind new entries to final data frame
                    final_df <- rbind(final_df, new_entries)
                    
                    
                    # shift all remaining second doses to...
                    if (k > 1){
                        # next age group down
                        AZ_doses[k-1] <- AZ_doses[k-1] + fdtd_AZ
                        PZ_doses[k-1] <- PZ_doses[k-1] + fdtd_PZ
                        MD_doses[k-1] <- MD_doses[k-1] + fdtd_MD
                    } else if (j < length(regions)){
                        # next region along
                        doses_per_region[j+1] <- doses_per_region[j+1] + fdtd_AZ + fdtd_PZ + fdtd_MD
                    } else if (i < length(vpdates)){
                        # next date along
                        vpdoses[i+1] <- vpdoses[i+1] + fdtd_AZ + fdtd_PZ + fdtd_MD
                    } else {
                        # or record leftover doses
                        leftover_doses <- leftover_doses + fdtd_AZ + fdtd_PZ + fdtd_MD
                    }
                }
            }
        }
    }
    
    # print leftover doses
    print(leftover_doses)
    
    # sort final_df data frame rows by region, age group, vaccine, date
    vax_schedule <- final_df[order(final_df$region, final_df$age.group, 
                                   final_df$vaccine, final_df$date),]
    
}

# save schedule containing data on delivered and projected deliveries in future
datetime <- str_replace_all(Sys.time(), "[ :GMTBST-]", "")
write.csv(vax_schedule, file = paste0(path_to_data, foldername, "/", filedate, "-", 
                                  datetime, "-delivered_and_projected_vaccines_", 
                                  name, ".csv"), row.names = FALSE)

# convert schedule

# add lagged dates
first_doses <- vax_schedule[vax_schedule$dose == 1,]
second_doses <- vax_schedule[vax_schedule$dose == 2,]
first_doses$lag_date <- as.Date(as.Date(first_doses$date) + 28)
second_doses$lag_date <- as.Date(as.Date(second_doses$date) + 14)
vax_schedule2 <- rbind(first_doses, second_doses)

# NHS England popset
popset = c(1, 3, 4, 5, 6, 9, 10)
covidm_regions = NULL
for (i in popset){
    covidm_regions[i] <- params$pop[[i]]$name  
}
covidm_regions <- covidm_regions[!is.na(covidm_regions)]


lag_dates <- as.Date(sort(unique(vax_schedule2$lag_date)))
va1 <- list()
va2 <- list()
vb1 <- list()
vb2 <- list()
ret <- list()

for (i in 1:length(covidm_regions)){
    this_region <- covidm_regions[i]
    print(this_region)
    this_schedu <- vax_schedule2[vax_schedule2$region == this_region,]
    for (j in 1:length(lag_dates)){
        
        this_date <- lag_dates[j]
        print(this_date)
        this_schedu_date <- this_schedu[this_schedu$lag_date == this_date,]
        
        AZ <- this_schedu_date[this_schedu_date$vaccine == "AZ",]
        PZ <- this_schedu_date[this_schedu_date$vaccine == "PF",]
        MD <- this_schedu_date[this_schedu_date$vaccine == "MD",]
        AZ1 <- AZ[AZ$dose == 1,]
        AZ1 <- AZ1[order(AZ1$age.group),]
        AZ2 <- AZ[AZ$dose == 2,]
        AZ2 <- AZ2[order(AZ2$age.group),]
        PZ1 <- PZ[PZ$dose == 1,]
        PZ1 <- PZ1[order(PZ1$age.group),]
        PZ2 <- PZ[PZ$dose == 2,]
        PZ2 <- PZ2[order(PZ2$age.group),]
        MD1 <- MD[MD$dose == 1,]
        MD1 <- MD1[order(MD1$age.group),]
        MD2 <- MD[MD$dose == 2,]
        MD2 <- MD2[order(MD2$age.group),]
        
        if (dim(AZ1)[1] == 0){
            va1[[j]] <- rep(0, 16)
            print('AZ dose 1:')
            print(paste0('All ', 0))
        } else {
            va1[[j]] <- AZ1$number
            print('AZ dose 1:')
            print(AZ1$number)
        }
        if (dim(AZ2)[1] == 0){
            va2[[j]] <- rep(0,16)
            print('AZ dose 2:')
            print(paste0('All ', 0))
        } else {
            va2[[j]] <- AZ2$number
            print('AZ dose 2:')
            print(AZ2$number)
        }
        if (dim(PZ1)[1] == 0 && dim(MD1)[1] == 0){
            vb1[[j]] <- rep(0, 16)
            print('Pfizer & Moderna dose 1:')
            print(paste0('All ', 0))
        } else {
            vb1[[j]] <- PZ1$number + MD1$number
            
            print('Pfizer & Moderna dose 1:')
            print(PZ1$number + MD1$number)
        }
        if (dim(PZ2)[1] == 0 && dim(MD2)[1] == 0){
            vb2[[j]] <- rep(0,16)
            print('Pfizer & Moderna dose 2:')
            print(paste0('All ', 0))
        } else {
            vb2[[j]] <- PZ2$number + MD2$number
            print('Pfizer & Moderna dose 2:')
            print(PZ2$number + MD2$number)
        }
        
    }
    idx <- popset[i]
    ret[[idx]] = list(
        vt = lag_dates,
        va1 = va1,
        va2 = va2,
        vb1 = vb1,
        vb2 = vb2
    )
}

# save final schedule as .rds file
datetime <- str_replace_all(Sys.time(), "[ :GMTBST-]", "")
saveRDS(ret, paste0(path_to_data, foldername, "/vax-covidm", datetime, name, ".rds"))
saveRDS(ret, paste0(uk_covid_data_path, "vax-covidm", datetime, name, ".rds"))

# 
# # create alternative vaccine schedules for paper using same datetime
# #####
# #####
# name = '12plus_50percent'
# uptake[3]      = 0.5*(3/5) # ages 12-14 out of 10-14 age group at 50% uptake
# uptake[4]      = 0.5       # 50% uptake limit for ages 15-19
# schedule <- generate_schedule(ll, ll4, uptake, final_doses, products, second_doses, popUK)
# write.csv(schedule, file = paste0(path_to_data,
#                                   foldername, "/", filedate, "-", datetime, 
#                                   "-delivered_and_projected_vaccines_", name, ".csv"), 
#           row.names = FALSE)
# final <- convert_schedule(schedule, params)
# saveRDS(final, paste0(path_to_data, foldername, "/vax-covidm", datetime, name, ".rds"))
# saveRDS(final, paste0(uk_covid_data_path, "vax-covidm", datetime, name, ".rds"))
# #####
# #####
# name = '5plus_50percent'
# uptake[2:4]      = 0.5       # 50% uptake limit for ages 5-19
# schedule <- generate_schedule(ll, ll4, uptake, final_doses, products, second_doses, popUK)
# write.csv(schedule, file = paste0(path_to_data,
#                                   foldername, "/", filedate, "-", datetime, 
#                                   "-delivered_and_projected_vaccines_", name, ".csv"), 
#           row.names = FALSE)
# final <- convert_schedule(schedule, params)
# saveRDS(final, paste0(path_to_data, foldername, "/vax-covidm", datetime, name, ".rds"))
# saveRDS(final, paste0(uk_covid_data_path, "vax-covidm", datetime, name, ".rds"))
# #####
# #####
# name = '5plus_80percent'
# uptake[2:4]      = 0.8       # 80% uptake limit for ages 5-19
# schedule <- generate_schedule(ll, ll4, uptake, final_doses, products, second_doses, popUK)
# write.csv(schedule, file = paste0(path_to_data,
#                                   foldername, "/", filedate, "-", datetime, 
#                                   "-delivered_and_projected_vaccines_", name, ".csv"), 
#           row.names = FALSE)
# final <- convert_schedule(schedule, params)
# saveRDS(final, paste0(path_to_data, foldername, "/vax-covidm", datetime, name, ".rds"))
# saveRDS(final, paste0(uk_covid_data_path, "vax-covidm", datetime, name, ".rds"))
