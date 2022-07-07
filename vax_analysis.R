# script to calculate quantities of interest from vaccines already delivered

# ---------------------- changeable items here --------------------------------
path_to_data <- "~/Documents/uk_covid_data_sensitive/phe/"
foldername   <- "2021-10-04"
filedate     <- "20211003" # change this to match the PHE reference date
linelist     <- "20211003 immunisations SPIM.csv"
pre_realloc  <- "20210826-20210827234220-processedvaccinedata.csv"
vax_file     <- "20210826-20210828003412-processedvaccinedata-reallocated.csv"
vax_project  <- "20210701-20210703005516-delivered_and_projected_vaccines.csv"
path_to_newcovidvax  <- "."
uk_covid_data_path   <- paste0(path_to_newcovidvax, "/fitting_data/")
cm_path              <- paste0(path_to_newcovidvax, "/covidm_for_fitting/")
setwd("~/OneDrive_LSHTM/GitHub/newcovid3/")
# ---------------------- end changeable items  --------------------------------

# load processed vaccine files
vax_delivered <- read.csv(paste0(path_to_data, foldername, "/", vax_file), row.names = NULL)
vax_projected <- read.csv(paste0(path_to_data, foldername, "/", vax_project), row.names = NULL)

# 1. calculate proportion of vaccine type to individuals aged 50+ in England
source(paste0(path_to_newcovidvax, '/vax_funcs.R'))
props <- products_delivered(vax_delivered, age_group_limit = 11)

# 2. calculate date by which all adults have received their second vaccine dose
youngest_adults <- vax_projected[vax_projected$ages == "15-19",]
adults2024 <- vax_projected[vax_projected$ages == "20-24",]
adults2529 <- vax_projected[vax_projected$ages == "25-29",]
for (i in unique(youngest_adults$region)){
    print(i)
    this_region <- youngest_adults[youngest_adults$region == i,]
    second_doses <- this_region[this_region$dose == "Second",]
    left_to_dose_0 <- second_doses[second_doses$left_to_dose_ut == 0,]
    AZ <- left_to_dose_0[left_to_dose_0$vaccine == "AZ",]
    AZ_all_dose_date <- min(AZ$date)
    PF <- left_to_dose_0[left_to_dose_0$vaccine == "PF",]
    PF_all_dose_date <- min(PF$date)
    MD <- left_to_dose_0[left_to_dose_0$vaccine == "MD",]
    MD_all_dose_date <- min(MD$date)
    all_dose_date <- max(AZ_all_dose_date, PF_all_dose_date, MD_all_dose_date)
    print(all_dose_date)
    
    this_region <- adults2024[adults2024$region == i,]
    second_doses <- this_region[this_region$dose == "Second",]
    left_to_dose_0 <- second_doses[second_doses$left_to_dose_ut == 0,]
    AZ <- left_to_dose_0[left_to_dose_0$vaccine == "AZ",]
    AZ_all_dose_date <- min(AZ$date)
    PF <- left_to_dose_0[left_to_dose_0$vaccine == "PF",]
    PF_all_dose_date <- min(PF$date)
    MD <- left_to_dose_0[left_to_dose_0$vaccine == "MD",]
    MD_all_dose_date <- min(MD$date)
    all_dose_date <- max(AZ_all_dose_date, PF_all_dose_date, MD_all_dose_date)
    print(all_dose_date)
    
    this_region <- adults2529[adults2529$region == i,]
    second_doses <- this_region[this_region$dose == "Second",]
    left_to_dose_0 <- second_doses[second_doses$left_to_dose_ut == 0,]
    AZ <- left_to_dose_0[left_to_dose_0$vaccine == "AZ",]
    AZ_all_dose_date <- min(AZ$date)
    PF <- left_to_dose_0[left_to_dose_0$vaccine == "PF",]
    PF_all_dose_date <- min(PF$date)
    MD <- left_to_dose_0[left_to_dose_0$vaccine == "MD",]
    MD_all_dose_date <- min(MD$date)
    all_dose_date <- max(AZ_all_dose_date, PF_all_dose_date, MD_all_dose_date)
    print(all_dose_date)
    
}

# 3. calculate date by which all adults have received their first vaccine dose 
vax_projected <- read.csv(paste0(path_to_data, foldername, "/", vax_project), row.names = NULL)
youngest_adults <- vax_projected[vax_projected$ages == "15-19",]
adults2024 <- vax_projected[vax_projected$ages == "20-24",]
adults2529 <- vax_projected[vax_projected$ages == "25-29",]
for (i in unique(youngest_adults$region)){
    print(i)
    this_region <- youngest_adults[youngest_adults$region == i,]
    second_doses <- this_region[this_region$dose == "First",]
    left_to_dose_0 <- second_doses[second_doses$left_to_dose_ut == 0,]
    AZ <- left_to_dose_0[left_to_dose_0$vaccine == "AZ",]
    AZ_all_dose_date <- min(AZ$date)
    PF <- left_to_dose_0[left_to_dose_0$vaccine == "PF",]
    PF_all_dose_date <- min(PF$date)
    MD <- left_to_dose_0[left_to_dose_0$vaccine == "MD",]
    MD_all_dose_date <- min(MD$date)
    all_dose_date <- max(AZ_all_dose_date, PF_all_dose_date, MD_all_dose_date)
    print(all_dose_date)
    
    this_region <- adults2024[adults2024$region == i,]
    second_doses <- this_region[this_region$dose == "Second",]
    left_to_dose_0 <- second_doses[second_doses$left_to_dose_ut == 0,]
    AZ <- left_to_dose_0[left_to_dose_0$vaccine == "AZ",]
    AZ_all_dose_date <- min(AZ$date)
    PF <- left_to_dose_0[left_to_dose_0$vaccine == "PF",]
    PF_all_dose_date <- min(PF$date)
    MD <- left_to_dose_0[left_to_dose_0$vaccine == "MD",]
    MD_all_dose_date <- min(MD$date)
    all_dose_date <- max(AZ_all_dose_date, PF_all_dose_date, MD_all_dose_date)
    print(all_dose_date)
    
    this_region <- adults2529[adults2529$region == i,]
    second_doses <- this_region[this_region$dose == "Second",]
    left_to_dose_0 <- second_doses[second_doses$left_to_dose_ut == 0,]
    AZ <- left_to_dose_0[left_to_dose_0$vaccine == "AZ",]
    AZ_all_dose_date <- min(AZ$date)
    PF <- left_to_dose_0[left_to_dose_0$vaccine == "PF",]
    PF_all_dose_date <- min(PF$date)
    MD <- left_to_dose_0[left_to_dose_0$vaccine == "MD",]
    MD_all_dose_date <- min(MD$date)
    all_dose_date <- max(AZ_all_dose_date, PF_all_dose_date, MD_all_dose_date)
    print(all_dose_date)
    
}

# calculate overall vaccine coverage

# set up covidm
datapath = function(x) paste0(uk_covid_data_path, x)
cm_force_rebuild = F
cm_build_verbose = T
cm_version = 2
source(paste0(cm_path, "/R/covidm.R"))
popUK = readRDS(datapath("popNHS.rds"))
matricesUK = readRDS(datapath("matricesNHS.rds"))
cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)
nhs_regions = popUK[, unique(name)]


setDT(vax_projected)
# add column for vaccine coverage by vaccine product (in each age group, region, dose)
vax_projected[, coverage_thisvax := cum_doses / pop_size]
# add column for cumulative doses ACROSS VACCINE TYPE by region, age group
vax_projected[, cum_doses_allvax := max_uptake - left_to_dose_ut]
# add column for vaccine coverage across vaccine products (in each age group, region, dose)
vax_projected[, coverage_allvax := cum_doses_allvax / pop_size]

# get vaccine coverage by region but ACROSS all age groups
agegroups <- unique(vax_projected$age.group)
regions   <- unique(vax_projected$region)
dates     <- unique(vax_projected$date)
vaccines  <- unique(vax_projected$vaccine)
doses     <- unique(vax_projected$dose)
regional_coverage <- data.frame(region = NULL, vaccine = NULL, dose = NULL, 
                                date = NULL, number = NULL, pop_size = NULL,
                                coverage_thisvax = NULL)
for (i in regions){
    print(i)
    regional <- vax_projected[vax_projected$region == i,]
    for (j in dates){
        regional_date <- regional[regional$date == j,]
        for (k in vaccines){
            regional_date_vax <- regional_date[regional_date$vaccine == k,]
            for (l in doses){
                regional_date_vax_dose <- regional_date_vax[regional_date_vax$dose == l,]
                numdose <- sum(regional_date_vax_dose$cum_doses)
                popsize <- sum(regional_date_vax_dose$pop_size)
                cov     <- numdose / popsize
                new_entry <- data.frame(region = i, vaccine = k, dose = l, date = j,
                                        number = numdose, pop_size = popsize, 
                                        coverage_thisvax = cov)
                regional_coverage <- rbind(regional_coverage, new_entry)
            }
        }
    }
}

regional_coverage_all <- data.frame(region = NULL, vaccine = NULL, dose = NULL, 
                                    date = NULL, number = NULL, pop_size = NULL,
                                    coverage_thisvax = NULL, coverage_allvax = NULL)
# now, also calculate regional coverage across vaccines
for (i in regions){
    regional <- regional_coverage[regional_coverage$region == i,]
    for (j in dates){
        regional_date <- regional[regional$date == j,]
        for (l in doses){
            regional_date_dose <- regional_date[regional_date$dose == l,]
            regional_date_dose$coverage_allvax = sum(regional_date_dose$coverage_thisvax)
            regional_coverage_all <- rbind(regional_coverage_all, regional_date_dose)
        }
    }
}

# get vaccine coverage by age group but ACROSS all regions
agegroup_coverage <- data.frame(age.group = NULL, vaccine = NULL, dose = NULL, 
                                date = NULL, number = NULL, pop_size = NULL,
                                coverage_thisvax = NULL)
for (i in agegroups){
    print(i)
    by_age <- vax_projected[vax_projected$age.group == i,]
    for (j in dates){
        by_age_date <- by_age[by_age$date == j,]
        for (k in vaccines){
            by_age_date_vax <- by_age_date[by_age_date$vaccine == k,]
            for (l in doses){
                by_age_date_vax_dose <- by_age_date_vax[by_age_date_vax$dose == l,]
                numdose <- sum(by_age_date_vax_dose$cum_doses)
                popsize <- sum(by_age_date_vax_dose$pop_size)
                cov     <- numdose / popsize
                new_entry <- data.frame(age.group = i, vaccine = k, dose = l, date = j,
                                        number = numdose, pop_size = popsize, 
                                        coverage_thisvax = cov)
                agegroup_coverage <- rbind(agegroup_coverage, new_entry)
            }
        }
    }
}
agegroup_coverage_all <- data.frame(age.group = NULL, vaccine = NULL, dose = NULL, 
                                    date = NULL, number = NULL, pop_size = NULL,
                                    coverage_thisvax = NULL, coverage_allvax = NULL)
# now, also calculate age group coverage across vaccines
for (i in agegroups){
    age_group <- agegroups[i]
    print(age_group)
    by_age <- agegroup_coverage[agegroup_coverage$age.group == age_group,]
    for (j in dates){
        by_age_date <- by_age[by_age$date == j,]
        for (l in doses){
            by_age_date_dose <- by_age_date[by_age_date$dose == l,]
            by_age_date_dose$coverage_allvax = sum(by_age_date_dose$coverage_thisvax)
            agegroup_coverage_all <- rbind(agegroup_coverage_all, by_age_date_dose)
        }
    }
}

# remove MD, PF entries (duplicates) and coverage_thisvax
agegroup_coverage_all <- agegroup_coverage_all[!(agegroup_coverage_all$vaccine == 'MD'),]
agegroup_coverage_all <- agegroup_coverage_all[!(agegroup_coverage_all$vaccine == 'PF'),]
agegroup_coverage_all$coverage_thisvax <- NULL
agegroup_coverage_all$vaccine <- NULL
first_doses <- agegroup_coverage_all[agegroup_coverage_all$dose == "First",]
second_doses <- agegroup_coverage_all[agegroup_coverage_all$dose == "Second",]

first_doses$date <- as.Date(first_doses$date)

first_doses <- first_doses[first_doses$age.group > 3,]
first_doses <- first_doses[!(is.na(first_doses$age.group == T)),]

ggplot(first_doses, aes(date, coverage_allvax, colour=as.factor(age.group))) + 
    geom_line(lwd = 0.8) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") + 
    labs(x = "Date", y = "Coverage", colour = "Age group") +
    scale_color_discrete(labels = c("15-19", "20-24", "25-29", "30-34", "35-39",
                                    "40-44", "45-49", "50-54", "55-59", "60-64",
                                    "65-69", "70-74", "75+"))

# try the same as above but plotting first and second dose coverage
second_doses <- second_doses[second_doses$age.group > 3,]
second_doses <- second_doses[!(is.na(second_doses$age.group == T)),]
first_doses$dose <- rep("First", dim(first_doses)[1])
second_doses$dose <- rep("Second", dim(second_doses)[1])
all_doses <- rbind(first_doses, second_doses)
# remove doses after a certain point
all_doses <- all_doses[all_doses$date < "2021-07-02",]

colours2 <- c("#FABAFB","#E7ABE8","#D59DD5","#C28EC2","#B081AF","#9E739D","#8D668C","#7C597A","#6B4C6A","#5B4059","#4B3449","#3C283A","#2D1D2C")

p1 <- ggplot(all_doses, aes(date, coverage_allvax, colour=as.factor(age.group))) + 
    geom_line(lwd = 0.8, aes(linetype = dose)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") + 
    scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1)) +
    labs(x = "", y = "Vaccination coverage", colour = "Age group", linetype = "Dose") +
    #scale_color_discrete(labels = c("15-19", "20-24", "25-29", "30-34", "35-39","40-44", "45-49", "50-54", "55-59", "60-64","65-69", "70-74", "75+")) +
    theme_cowplot(font_size = 9) + theme(strip.background = element_blank()) +
    background_grid() + 
    scale_color_manual(values = colours2, labels = c("15-19", "20-24", "25-29", "30-34", "35-39",
                                                     "40-44", "45-49", "50-54", "55-59", "60-64",
                                                     "65-69", "70-74", "75+"))

# next try plotting vaccine coverage by region
regional_coverage_all <- regional_coverage_all[!(regional_coverage_all$vaccine == 'MD'),]
regional_coverage_all <- regional_coverage_all[!(regional_coverage_all$vaccine == 'PF'),]
regional_coverage_all$coverage_thisvax <- NULL
regional_coverage_all$vaccine <- NULL

toplot <- regional_coverage_all[regional_coverage_all$date < "2021-07-02",]
toplot$date <- as.Date(toplot$date)
p2 <- ggplot(toplot, aes(date, coverage_allvax, colour=as.factor(region))) + 
    geom_line(lwd = 0.8, aes(linetype = dose)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") + 
    scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1)) +
    labs(x = "", y = "Vaccination coverage", colour = "NHS England region", linetype = "Dose") +
    theme_cowplot(font_size = 9) + theme(strip.background = element_blank()) +
    background_grid() +
    scale_color_colorblind()

# make a side by side plot showing coverage by age group and region
both <- cowplot::plot_grid(p1, p2, ncol=1, align = "v", labels = c("A","B"))
ggsave("./output/figures/vax_coverage_overview_20210705_ax_CBfriendly.png", both, width = 15, height = 15, units = "cm")
ggsave("./output/figures/vax_coverage_overview_20210705_ax_CBfriendly.pdf", both, width = 15, height = 15, units = "cm", useDingbats = FALSE)

# new plot 10th July 2021 - vaccine coverage without imputing missing values and with colourblind friendly colour palette

vax_processed <- read.csv(paste0(path_to_data, foldername, "/", pre_realloc), row.names = NULL)
setDT(vax_processed)
# remove region == "" and ages == NA
vax_processed <- vax_processed[!(vax_processed$region == ""),]
vax_processed <- vax_processed[!(is.na(vax_processed$ages) == TRUE),]

# need to calculate additional elements such as cum_doses and max_uptake

# assumptions on uptake for roadmap modelling step 4 encore (26th June 2021)
vac_uptake        <- rep(0, 16)
vac_uptake[4]     <- 0.8*(2/5) # 80% uptake for ages 18-19 in 15-19 age group
vac_uptake[5:8]   <- 0.8  # 80% uptake for ages 20-39 inclusive
# we need to calculate uptake limits for individuals aged 40 and above using
# the numbers of first doses delivered already divided by the total population
# loop through regions
num_delivered <- 0
num_delivered_d2 <- 0
agegroups <- unique(vax_processed$ages)
for (i in unique(vax_processed$region)){
    
    print(i)
    this_region <- vax_processed[vax_processed$region == i,]
    age_groups_40_above <- agegroups[9:16]
    # loop through age groups over 40
    for (j in age_groups_40_above){
        
        print(j)
        this_region_age <- this_region[this_region$ages == j,]
        this_region_age_d1 <- this_region_age[this_region_age$dose == "First",]
        this_region_age_d2 <- this_region_age[this_region_age$dose == "Second",]
        num_delivered <- num_delivered + sum(this_region_age_d1$number)
        num_delivered_d2 <- num_delivered_d2 + sum(this_region_age_d2$number)
        
    }
}

# set up covidm data path
datapath = function(x) paste0(uk_covid_data_path, x)
popUK = readRDS(datapath("popNHS.rds"))

popUK_ages <- unique(popUK$age)
popUK_ages_above_40 <- popUK_ages[9:19]
total_pop_over_40 <- 0

for (i in unique(vax_processed$region)){
    
    print(i)
    this_region <- popUK[popUK$name == i,]
    
    # loop through popUK age groups over 40
    for (j in popUK_ages_above_40){
        
        print(j)
        this_region_age <- this_region[this_region$age == j,]
        total_pop_over_40 <- total_pop_over_40 + (this_region_age$f*1000) + (this_region_age$m*1000)
        
    }
}

# using population size for 40+ year olds in England, divide num_delivered by 
# that population size for the vac_uptake[9:16] entry
vac_uptake[9:16] <- min(num_delivered / total_pop_over_40,1)
print(vac_uptake)

final_df = data.frame(region = NULL, age.group = NULL, ages = NULL, 
                      vaccine = NULL, dose = NULL, date = NULL, number = NULL,
                      cum_doses = NULL, pop_size = NULL, max_uptake = NULL, 
                      left_to_dose = NULL, left_to_dose_ut = NULL)
vaccines <- unique(vax_processed$vaccine)
doses <- unique(vax_processed$dose)
# remove NA ages from vax_processed data frame and age_groups
vax_processed <- vax_processed[!(is.na(vax_processed$age.group) == TRUE),]
age_groups <- unique(vax_processed$ages)
age_groups <- age_groups[!(is.na(age_groups)==TRUE)]

# rename regions from codes to names
vax_processed[vax_processed$region == "E40000010",]$region <- "North West" 
vax_processed[vax_processed$region == "E40000008",]$region <- "Midlands" 
vax_processed[vax_processed$region == "E40000009",]$region <- "North East and Yorkshire" 
vax_processed[vax_processed$region == "E40000007",]$region <- "East of England" 
vax_processed[vax_processed$region == "E40000003",]$region <- "London" 
vax_processed[vax_processed$region == "E40000006",]$region <- "South West" 
vax_processed[vax_processed$region == "E40000005",]$region <- "South East" 

regions <- unique(vax_processed$region)
for (i in 1:length(regions)){
    this_region <- regions[i]
    print(this_region)
    this_region_data <- vax_processed[vax_processed$region == this_region,]
    this_region_pop <- popUK[popUK$name == this_region,]
    for (j in 1:length(age_groups)){
        this_age_group <- age_groups[j]
        this_data <- this_region_data[this_region_data$ages == this_age_group,]
        print(this_age_group)
        # note that covidm population sizes are listed in 1000's
        if (this_age_group == "75+"){
            all_pops <- this_region_pop[this_region_pop$age %in% c("75-79","80-84","85-89","90+"),]
            pop_size <- (sum(all_pops$f) + sum(all_pops$m)) * 1000
        } else {
            this_region_popage <- this_region_pop[this_region_pop$age == this_age_group,]
            pop_size <- (this_region_popage$f + this_region_popage$m) * 1000
        }
        
        inner_df = data.frame(region = NULL, age.group = NULL, ages = NULL, 
                              vaccine = NULL, dose = NULL, date = NULL, 
                              number = NULL, cum_doses = NULL, pop_size = NULL, 
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
                this_data_vd$max_uptake <- round(this_data_vd$pop_size*vac_uptake[j])
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

vax_processed <- final_df
# add column for vaccine coverage by vaccine product (in each age group, region, dose)
vax_processed[, coverage_thisvax := cum_doses / pop_size]
# add column for cumulative doses ACROSS VACCINE TYPE by region, age group
vax_processed[, cum_doses_allvax := max_uptake - left_to_dose_ut]
# add column for vaccine coverage across vaccine products (in each age group, region, dose)
vax_processed[, coverage_allvax := cum_doses_allvax / pop_size]

# get vaccine coverage by region but ACROSS all age groups
agegroups <- unique(vax_processed$age.group)
regions   <- unique(vax_processed$region)
dates     <- unique(vax_processed$date)
vaccines  <- unique(vax_processed$vaccine)
doses     <- unique(vax_processed$dose)
regional_coverage <- data.frame(region = NULL, vaccine = NULL, dose = NULL, 
                                date = NULL, number = NULL, pop_size = NULL,
                                coverage_thisvax = NULL)
for (i in regions){
    print(i)
    regional <- vax_processed[vax_processed$region == i,]
    for (j in dates){
        regional_date <- regional[regional$date == j,]
        for (k in vaccines){
            regional_date_vax <- regional_date[regional_date$vaccine == k,]
            for (l in doses){
                regional_date_vax_dose <- regional_date_vax[regional_date_vax$dose == l,]
                numdose <- sum(regional_date_vax_dose$cum_doses)
                popsize <- sum(regional_date_vax_dose$pop_size)
                cov     <- numdose / popsize
                new_entry <- data.frame(region = i, vaccine = k, dose = l, date = j,
                                        number = numdose, pop_size = popsize, 
                                        coverage_thisvax = cov)
                regional_coverage <- rbind(regional_coverage, new_entry)
            }
        }
    }
}

regional_coverage_all <- data.frame(region = NULL, vaccine = NULL, dose = NULL, 
                                    date = NULL, number = NULL, pop_size = NULL,
                                    coverage_thisvax = NULL, coverage_allvax = NULL)
# now, also calculate regional coverage across vaccines
for (i in regions){
    regional <- regional_coverage[regional_coverage$region == i,]
    for (j in dates){
        regional_date <- regional[regional$date == j,]
        for (l in doses){
            regional_date_dose <- regional_date[regional_date$dose == l,]
            regional_date_dose$coverage_allvax = sum(regional_date_dose$coverage_thisvax)
            regional_coverage_all <- rbind(regional_coverage_all, regional_date_dose)
        }
    }
}

# get vaccine coverage by age group but ACROSS all regions
agegroup_coverage <- data.frame(age.group = NULL, vaccine = NULL, dose = NULL, 
                                date = NULL, number = NULL, pop_size = NULL,
                                coverage_thisvax = NULL)
for (i in agegroups){
    print(i)
    by_age <- vax_processed[vax_processed$age.group == i,]
    for (j in dates){
        by_age_date <- by_age[by_age$date == j,]
        for (k in vaccines){
            by_age_date_vax <- by_age_date[by_age_date$vaccine == k,]
            for (l in doses){
                by_age_date_vax_dose <- by_age_date_vax[by_age_date_vax$dose == l,]
                numdose <- sum(by_age_date_vax_dose$cum_doses)
                popsize <- sum(by_age_date_vax_dose$pop_size)
                cov     <- numdose / popsize
                new_entry <- data.frame(age.group = i, vaccine = k, dose = l, date = j,
                                        number = numdose, pop_size = popsize, 
                                        coverage_thisvax = cov)
                agegroup_coverage <- rbind(agegroup_coverage, new_entry)
            }
        }
    }
}
agegroup_coverage_all <- data.frame(age.group = NULL, vaccine = NULL, dose = NULL, 
                                    date = NULL, number = NULL, pop_size = NULL,
                                    coverage_thisvax = NULL, coverage_allvax = NULL)
# now, also calculate age group coverage across vaccines
for (i in agegroups){
    age_group <- agegroups[i]
    print(age_group)
    by_age <- agegroup_coverage[agegroup_coverage$age.group == age_group,]
    for (j in dates){
        by_age_date <- by_age[by_age$date == j,]
        for (l in doses){
            by_age_date_dose <- by_age_date[by_age_date$dose == l,]
            by_age_date_dose$coverage_allvax = sum(by_age_date_dose$coverage_thisvax)
            agegroup_coverage_all <- rbind(agegroup_coverage_all, by_age_date_dose)
        }
    }
}

# remove MD, PF entries (duplicates) and coverage_thisvax
agegroup_coverage_all <- agegroup_coverage_all[!(agegroup_coverage_all$vaccine == 'MD'),]
agegroup_coverage_all <- agegroup_coverage_all[!(agegroup_coverage_all$vaccine == 'PF'),]
agegroup_coverage_all$coverage_thisvax <- NULL
agegroup_coverage_all$vaccine <- NULL
first_doses <- agegroup_coverage_all[agegroup_coverage_all$dose == "First",]
second_doses <- agegroup_coverage_all[agegroup_coverage_all$dose == "Second",]

first_doses$date <- as.Date(first_doses$date)

first_doses <- first_doses[first_doses$age.group > 3,]
first_doses <- first_doses[!(is.na(first_doses$age.group == T)),]

ggplot(first_doses, aes(date, coverage_allvax, colour=as.factor(age.group))) + 
    geom_line(lwd = 0.8) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") + 
    labs(x = "Date", y = "Coverage", colour = "Age group") +
    scale_color_discrete(labels = c("15-19", "20-24", "25-29", "30-34", "35-39",
                                    "40-44", "45-49", "50-54", "55-59", "60-64",
                                    "65-69", "70-74", "75+"))

# try the same as above but plotting first and second dose coverage
second_doses <- second_doses[second_doses$age.group > 3,]
second_doses <- second_doses[!(is.na(second_doses$age.group == T)),]
first_doses$dose <- rep("First", dim(first_doses)[1])
second_doses$dose <- rep("Second", dim(second_doses)[1])
all_doses <- rbind(first_doses, second_doses)
# remove doses after a certain point
all_doses <- all_doses[all_doses$date < "2021-07-02",]
library(ggthemes)

colours <- c("#A275A0","#996E97","#90688E","#886186","#7F5B7D","#775575","#6E4E6D","#664865","#5E425D","#563C55","#4E374D","#463145","#3F2B3E")

colours2 <- c("#FABAFB","#E7ABE8","#D59DD5","#C28EC2","#B081AF","#9E739D","#8D668C","#7C597A","#6B4C6A","#5B4059","#4B3449","#3C283A","#2D1D2C")

p1 <- ggplot(all_doses, aes(date, coverage_allvax, colour=as.factor(age.group))) + 
    geom_line(lwd = 0.8, aes(linetype = dose)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") + 
    scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1)) +
    labs(x = "", y = "Vaccination coverage", colour = "Age group", linetype = "Dose") +
    #scale_color_discrete(labels = c("15-19", "20-24", "25-29", "30-34", "35-39","40-44", "45-49", "50-54", "55-59", "60-64","65-69", "70-74", "75+")) +
    theme_cowplot(font_size = 9) + theme(strip.background = element_blank()) +
    background_grid() + 
    scale_color_manual(values = colours2, labels = c("15-19", "20-24", "25-29", "30-34", "35-39",
                                                    "40-44", "45-49", "50-54", "55-59", "60-64",
                                                    "65-69", "70-74", "75+"))
p1
# next try plotting vaccine coverage by region
regional_coverage_all <- regional_coverage_all[!(regional_coverage_all$vaccine == 'MD'),]
regional_coverage_all <- regional_coverage_all[!(regional_coverage_all$vaccine == 'PF'),]
regional_coverage_all$coverage_thisvax <- NULL
regional_coverage_all$vaccine <- NULL

toplot <- regional_coverage_all[regional_coverage_all$date < "2021-07-02",]
toplot$date <- as.Date(toplot$date)
p2 <- ggplot(toplot, aes(date, coverage_allvax, colour=as.factor(region))) + 
    geom_line(lwd = 0.8, aes(linetype = dose)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") + 
    scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1)) +
    labs(x = "", y = "Vaccination coverage", colour = "NHS England region", linetype = "Dose") +
    theme_cowplot(font_size = 9) + theme(strip.background = element_blank()) +
    background_grid() +
    scale_color_colorblind()
p2
# make a side by side plot showing coverage by age group and region
both <- cowplot::plot_grid(p1, p2, ncol=1, align = "v", labels = c("A","B"))
ggsave("./output/figures/vax_coverage_overview_20210705_ax_noimpute_CBfriendly.png", both, width = 15, height = 15, units = "cm")
ggsave("./output/figures/vax_coverage_overview_20210705_ax_noimpute_CBfriendly.pdf", both, width = 15, height = 15, units = "cm", useDingbats = FALSE)


# ---------------------


# 4. calculate delays between first and second doses
library(data.table)
ll <- fread(paste0(path_to_data, foldername, "/", linelist))

library(dplyr)

# count number of times each patient_pseudo_id occurs in ll
person_freq = count(ll, ll$patient_pseudo_id)

# only keep entries where a patient_pseudo_id occurs twice (so two doses)
person_freq2 <- person_freq[person_freq$n == 2,]

# get entries in ll where patient_pseudo_id occurs in person_freq2
final_ll <- ll[ll$patient_pseudo_id %in% person_freq2$`ll$patient_pseudo_id`,]

# calculate date in suitable format
final_ll$date <- as.Date(final_ll$vaccination_date, format = "%d%b%Y")

# make sure final_ll is sorted in order of patient_pseudo_id and date
final_ll <- final_ll[order(final_ll$patient_pseudo_id, final_ll$date),]

# split final_ll into first and second doses (use odd and even row numbers)
first_doses <- final_ll[seq(1,dim(final_ll)[1]-1,2),]
secnd_doses <- final_ll[seq(2,dim(final_ll)[1],2),]
# first_doses <- final_ll[final_ll$dose_number == "First",]
# secnd_doses <- final_ll[final_ll$dose_number == "Second",]

# calculate delay between first and second doses for all individuals
dose_delay <- as.numeric(secnd_doses$date - first_doses$date)

# store dose_delay as column for secnd_doses data.frame
secnd_doses$dose_delay <- dose_delay
mean(secnd_doses$dose_delay)

# look more closely at short delays (<50 days)
short_delays <- secnd_doses[secnd_doses$dose_delay < 50,]
hist(short_delays$dose_delay)

# look more closely at earliest second doses (<= 26 January 2021)
# N.B JCVI issued change in dose delay guidance on 26th January 2021 here
# https://www.gov.uk/government/publications/prioritising-the-first-covid-19-vaccine-dose-jcvi-statement/optimising-the-covid-19-vaccination-programme-for-maximum-short-term-impact
early_doses <- secnd_doses[secnd_doses$date <= "2021-01-26",]
print(mean(early_doses$dose_delay))
hist(early_doses$dose_delay, main = 'COVID-19 vaccine dose 1:2 delay \n up to 26th January 2021')

# look more closely at later second doses (> 26 January 2021)
late_doses <- secnd_doses[secnd_doses$date > "2021-01-26",]
print(mean(late_doses$dose_delay))
hist(late_doses$dose_delay, main = 'COVID-19 vaccine dose 1:2 delay \n after 26th January 2021', breaks = 100, xlim = c(0,150))

# look more closely at delays for younger individuals
young_doses <- secnd_doses[secnd_doses$age < 50,]
print(mean(young_doses$dose_delay))
hist(young_doses$dose_delay, main = 'COVID-19 vaccine dose 1:2 delay \n individuals <50', breaks = 100, xlim = c(0,150))

# and older individuals
older_doses <- secnd_doses[secnd_doses$age >= 50,]
print(mean(older_doses$dose_delay))
hist(older_doses$dose_delay, main = 'COVID-19 vaccine dose 1:2 delay \n individuals 50+', breaks = 100, xlim = c(0,150))
