# function to check contents of PHE immunisations linelist
check_ll <- function(ll, reference_date){
    
    entries <- dim(ll)[1]
    print(paste0('Immunisations linelist has ', entries, ' entries'))
    
    ll$actual_date <- as.Date(ll$vaccination_date, format = "%d%b%Y")
    print(paste0('Dates range between ', min(ll$actual_date), ' and ', max(ll$actual_date)))
    
    # recode dates occurring before 2020 as 2020 dates
    year(ll$actual_date[year(ll$actual_date) < 2020]) = 2020
    # recode dates occurring after this year as this year's dates
    year(ll$actual_date[year(ll$actual_date) > year(today())]) = year(today())
    # recode 2020 dates occurring before 8th December 2020 as 2021 dates
    year(ll$actual_date[ll$actual_date < '2020-12-08']) = 2021
    # recode dates occurring after reference_date as dates from the previous year
    year(ll$actual_date[ll$actual_date > reference_date]) = year(ll$actual_date[ll$actual_date > reference_date])-1
    # finally, remove any dates which are still before 8th December 2020 (e.g. this could happen when reference_date is < 8th December 2021)
    ll = ll[!(ll$actual_date < '2020-12-08'),]
    
    # reference_date input refers to the maximum date in the vaccinations linelist
    # the earliest date in the vaccinations linelist should be 8th December 2020
    
    # 
    # erroneous_dates <- sum(ll$actual_date < "2020-12-08")
    # if (erroneous_dates > 0){
    #     erroneous_proportion <- erroneous_dates / entries
    #     print(paste0('The linelist has ', erroneous_dates, ' entries occurring before vaccination rollout started in England on 8th December 2020'))
    #     print(paste0('These erroneous entries comprise ', erroneous_proportion*100, '% of all entries'))
    #     print('Removing entries occurring after today in 2020 and before 8th December 2020')
    #     today_last_year <- as.Date(as.numeric(today())-365, origin = "1970-01-01")
    #     ll <- ll[!(ll$actual_date > today_last_year & ll$actual_date < "2020-12-08"),]
    #     print('Recoding entries occuring before today in 2020 as 2021 deliveries')
    #     ll$actual_date[ll$actual_date < today_last_year] <- as.Date(as.numeric(ll$actual_date[ll$actual_date < today_last_year]) + 365, origin = "1970-01-01")
    # }
    
    return(ll)
}

# function to process PHE immunisations linelist into format suitable for covidm
process_ll <- function(ll, ages){
    
    ll$vaccination_date <- ll$actual_date

    agegroups <- ages$groups
    
    # get sequence of vaccination dates, from first to last recorded in vax
    dates <- seq(as.Date(min(ll$vaccination_date)), 
                 as.Date(max(ll$vaccination_date)), by = 1)
    num_dates <- length(dates)
    
    # get list of NHS England (and Unknown!) regions of residence in vax
    regions <- unique(ll$region_of_residence)
    
    # get list of vaccines in vax -> remove "" entries for vaccine type
    vaccines <- unique(ll$product_display_type)
    vaccines <- vaccines[1:3]
    print(vaccines)
    # if ("NA" %in% vaccines){
    #   vaccines <- vaccines[-(which(is.na(vaccines)==TRUE))] 
    # }
    
    # redistribute doses with "" vaccine type following an 80/20 AZ/Pfizer split
    unknown_vax <- ll[ll$product_display_type == "",]
    num_to_correct <- dim(unknown_vax)[1]
    # define AZ / PF / MD split to distribute unknown vaccines (must sum to 1)
    vax_split <- c(0.8, 0.2, 0)
    if (sum(vax_split)!=1){
        stop("Sum of split vector should equal 1")
    }
    
    # use the Pfizer proportion and sample required number without replacement
    # DO THIS A DIFFERENT WAY IF MD PROPORTION IS NON-ZERO!
    PF_to_distribute <- round(vax_split[2]*num_to_correct)  
    samples <- sample(num_to_correct, size = PF_to_distribute, replace = FALSE)
    product_type_vector <- rep("AZ", num_to_correct)
    product_type_vector[samples] <- "PF"
    unknown_vax$product_display_type <- product_type_vector
    
    # remove unknown vax from main data frame
    ll <- ll[!(ll$product_display_type == ""),]
    # bind unknown vax entries with imputed vaccine product type to main data frame
    ll <- rbind(ll, unknown_vax)
    
    # get list of possible doses in vax
    doses <- unique(ll$dose_number)
    # remove "" entry from doses
    doses <- doses[!(doses == "")]
    
    # recode all entries with "" (i.e. unknown) dose number as first doses
    num_unknown_doses <- sum(ll$dose_number == "")
    
    if (num_unknown_doses > 0){
        # get unknown doses in separate data frame
        unknown_doses <- ll[ll$dose_number == "",]
        unknown_doses$dose_number <- rep("First", dim(unknown_doses)[1])
        # remove unknown doses from main data frame
        ll <- ll[!(ll$dose_number == ""),]
        # bind updated entries back to main data frame
        ll <- rbind(ll, unknown_doses)
    }
    
    # get list of possible ages in vax (not required below but checking for NAs)
    agelist <- unique(ll$age)
    
    # initialise dataframe to store vaccine data for covidm
    vax_df <- data.frame(region = NULL, age.group = NULL, ages = NULL, 
                         vaccine = NULL, dose = NULL, date = NULL, number = NULL)
    
    # loop through regions
    for (i in 1:length(regions)){
        
        this_region <- regions[i]
        print(paste0('Commencing region ', i, ' of ', length(regions), ': ', this_region))
        this_region_vax <- ll[ll$region_of_residence == this_region,]
        
        # loop through age groups + 1 (one extra group for NA age entries)
        for (j in 1:(length(agegroups)+1)){
            
            print(paste0('Commencing age group ', j, ' of ', length(agegroups)+1))
            
            if (j == length(agegroups)+1){
                
                this_age_group <- NA
                this_age_group_vax <- this_region_vax[is.na(this_region_vax$age) == T,]
                
                print('Age group is NA ages')
                
            } else {
                
                this_age_group <- agegroups[j]
                lower_age <- ages$min_ages[j]
                upper_age <- ages$max_ages[j]
                this_age_group_vax <- this_region_vax[is.na(this_region_vax$age) == F,]
                this_age_group_vax <- this_age_group_vax[this_age_group_vax$age >= lower_age,]
                this_age_group_vax <- this_age_group_vax[this_age_group_vax$age <= upper_age,]
                
                print(paste0('Age group is ', this_age_group))
                
            }
            
            # loop through vaccines
            for (k in 1:length(vaccines)){
                
                this_vaccine <- vaccines[k]
                print(paste0('Commencing vaccine product ', k, ' of ', length(vaccines), ': ', this_vaccine))
                this_vaccine_vax <- this_age_group_vax[this_age_group_vax$product_display_type == this_vaccine,]
                
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
                        deliveries[m] <- dim(this_date_vax)[1]
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
    if (dim(ll)[1] != sum(vax_df$number)){
        warning("Total doses in linelist does not match total doses in processed dataframe")
    } 
    
    return(vax_df)
}

# function to redistribute entries with unknown region or age
redistribute_entries <- function(vaxp, ages, vax, popUK){
    
    # rename regions from codes to names
    NW <- vaxp[vaxp$region == "E40000010",]
    NW$region <- rep("North West", dim(NW)[1])
    MD <- vaxp[vaxp$region == "E40000008",]
    MD$region <- rep("Midlands", dim(MD)[1])
    NEY <- vaxp[vaxp$region == "E40000009",]
    NEY$region <- rep("North East and Yorkshire", dim(NEY)[1])
    EE <- vaxp[vaxp$region == "E40000007",]
    EE$region <- rep("East of England", dim(EE)[1])
    LDN <- vaxp[vaxp$region == "E40000003",]
    LDN$region <- rep("London", dim(LDN)[1])
    SW <- vaxp[vaxp$region == "E40000006",]
    SW$region <- rep("South West", dim(SW)[1])
    SE <- vaxp[vaxp$region == "E40000005",]
    SE$region <- rep("South East", dim(SE)[1])
    unknown <- vaxp[vaxp$region == "",]
    
    rm(vaxp)
    vaxp <- rbind(NW, MD)
    vaxp <- rbind(vaxp, NEY)
    vaxp <- rbind(vaxp, EE)
    vaxp <- rbind(vaxp, LDN)
    vaxp <- rbind(vaxp, SW)
    vaxp <- rbind(vaxp, SE)
    vaxp <- rbind(vaxp, unknown)
    
    # calculate age distribution of vaccines delivered
    oldest_age <- max(vax$age, na.rm = TRUE)
    breaks = c(0,4.5,9.5,14.5,19.5,24.5,29.5,34.5,39.5,44.5,49.5,54.5,59.5,64.5,
               69.5,74.5,oldest_age)
    hist <- hist(vax$age, freq = FALSE, breaks)
    nbins <- length(breaks) - 1
    binwidths <- diff(breaks)
    agegroup_probs <- binwidths * hist$density
    cumprobs <- cumsum(agegroup_probs)
    
    # get list of NHS England (and Unknown!) regions of residence in vaxp
    regions <- unique(vaxp$region)
    
    # distribute vaccines with known region and unknown age
    krua_vaxp <- vaxp[is.na(vaxp$ages) == TRUE,]
    known_regions <- regions[! (regions %in% "")]
    # for each region
    for (i in 1:length(known_regions)){
        krua_vaxp_tr <- krua_vaxp[krua_vaxp$region == known_regions[i],]
        # get non-zero vaccine entries
        krua_vaxp_tr_n0 <- krua_vaxp_tr[krua_vaxp_tr$number > 0,]
        if (dim(krua_vaxp_tr_n0)[1] > 0){
            # need to distribute vaccines
            stop('Add code to distribute vaccines with known region unknown age')
        }
    }
    
    # distribute vaccines with unknown region and known age
    urka_vaxp <- vaxp[vaxp$region == "",]
    urka_vaxp <- urka_vaxp[is.na(urka_vaxp$ages) == FALSE,]
    urka_vaxp_n0 <- urka_vaxp[urka_vaxp$number > 0,]
    # loop through age groups recorded in urka_vaxp_n0
    age_groups_to_loop <- unique(urka_vaxp_n0$ages)
    for (i in 1:length(age_groups_to_loop)){
        this_age_group <- age_groups_to_loop[i]
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
        # get vaccines to distribute for this age group
        to_distribute <- urka_vaxp_n0[urka_vaxp_n0$ages == this_age_group,]
        # loop through number of vaccine dates with vaccines to distribute
        for (j in 1:length(to_distribute$date)){
            this_vaccine <- to_distribute$vaccine[j]
            this_dose <- to_distribute$dose[j]
            this_date <- to_distribute$date[j]
            num_to_distribute <- to_distribute$number[j]
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
                row_idx <- which(vaxp$region == region_to_go & 
                                     vaxp$ages == this_age_group & 
                                     vaxp$vaccine == this_vaccine & 
                                     vaxp$dose == this_dose & vaxp$date == this_date)
                if (length(row_idx) > 1){
                    stop("More than one location for this vaccine to go: code fix required")
                }
                vaxp$number[row_idx] <- vaxp$number[row_idx] + 1
            }
        }
    }
    
    # distribute vaccines with unknown region and unknown age
    agegroups <- ages$groups
    urua_vaxp <- vaxp[vaxp$region == "",]
    urua_vaxp <- urua_vaxp[is.na(urua_vaxp$ages) == TRUE,]
    urua_vaxp_n0 <- urua_vaxp[urua_vaxp$number > 0,]
    # pre calculate probabilities for regions by age group
    cs_norm_pop_allages <- list()
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
        # loop through each dose that needs allocating
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
            row_idx <- which(vaxp$region == region_to_go & 
                                 vaxp$ages == age_group_to_go & 
                                 vaxp$vaccine == this_vaccine & 
                                 vaxp$dose == this_dose & vaxp$date == this_date)
            if (length(row_idx) > 1){
                stop("More than one location for this vaccine to go: code fix required")
            }
            vaxp$number[row_idx] <- vaxp$number[row_idx] + 1
        }
    }
    
    # finally, remove entries within vaxp that have region = "Unknown" and age = NA
    final_vaxp <- vaxp[!(vaxp$region == ""),]
    final_vaxp <- final_vaxp[!((is.na(final_vaxp$ages))==TRUE),]
    
    return(final_vaxp)
}

# function to generate schedule of vaccinations (existing + future vaccinations)
generate_schedule <- function(vax, vax_df, uptake, doses, products, second_doses, popUK){
    
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
    
    # if second_doses == TRUE then doses contains first *and* second doses
    if (second_doses == TRUE){
        
        # calculate list of dates and number of doses *per day*
        vpdates <- seq(start_date, end_date, by = "day")
        vpdoses <- rep(doses[1]/7, length(vpdates))
        
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
        
        # calculate age distribution of vaccines delivered
        
        # get list of age groups in covidm (for England)
        oldest_age <- max(vax$age, na.rm = TRUE)
        min_ages <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)
        max_ages <- c(4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,oldest_age)
        ages <- data.frame(groups = agegroups, min_ages = min_ages, max_ages = max_ages)
        breaks = c(0,4.5,9.5,14.5,19.5,24.5,29.5,34.5,39.5,44.5,49.5,54.5,59.5,64.5,
                   69.5,74.5,oldest_age)
        hist <- hist(vax$age, freq = FALSE, breaks)
        nbins <- length(breaks) - 1
        binwidths <- diff(breaks)
        agegroup_probs <- binwidths * hist$density
        
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
                these_doses <- doses_per_region[j]
                
                # distribute these_doses into age groups for this_region
                
                # IMPORTANT: we assume that the age distribution of vaccines to be delivered
                # follows the existing age distribution in the vaccines delivered already,
                # starting with the oldest age group to the youngest. Any leftover doses 
                # get carried over to equivalent second doses, doses for the next age group 
                # down, the next region along, or the following date along, or are otherwise
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
                    first_dose_data <- dfra[dfra$dose == "First",]
                    fdAZ_data <- first_dose_data[first_dose_data$vaccine == "AZ",]
                    fdAZ_data_lastday <- fdAZ_data[fdAZ_data$date == date_yesterday,]
                    fdPZ_data <- first_dose_data[first_dose_data$vaccine == "PF",]
                    fdPZ_data_lastday <- fdPZ_data[fdPZ_data$date == date_yesterday,]
                    fdMD_data <- first_dose_data[first_dose_data$vaccine == "MD",]
                    fdMD_data_lastday <- fdMD_data[fdMD_data$date == date_yesterday,]
                    
                    second_dose_data <- dfra[dfra$dose == "Second",]
                    sdAZ_data <- second_dose_data[second_dose_data$vaccine == "AZ",]
                    sdAZ_data_lastday <- sdAZ_data[sdAZ_data$date == date_yesterday,]
                    sdPZ_data <- second_dose_data[second_dose_data$vaccine == "PF",]
                    sdPZ_data_lastday <- sdPZ_data[sdPZ_data$date == date_yesterday,]
                    sdMD_data <- second_dose_data[second_dose_data$vaccine == "MD",]
                    sdMD_data_lastday <- sdMD_data[sdMD_data$date == date_yesterday,]
                    
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
                    
                    # check how many second doses to deliver 
                    # this has changed slightly to assume a dosing gap of 8 weeks for individuals
                    # aged 50 and over, and 11 weeks for individuals less than 50 years old
                    # then, appropriately adjust first doses to deliver if necessary
                    if (k >= 11){
                        twelveweeks_prior <- date_today - 7*8
                    } else {
                        twelveweeks_prior <- date_today - 7*11
                    }
                    
                    if (twelveweeks_prior %in% first_dose_data$date){
                        first_doses_then <- first_dose_data[first_dose_data$date == twelveweeks_prior,]
                        
                        first_doses_then_AZ <- first_doses_then[first_doses_then$vaccine == "AZ",]$cum_doses
                        last_day_second_doses_AZ <- sdAZ_data_lastday$cum_doses
                        sdtd_AZ <- first_doses_then_AZ - last_day_second_doses_AZ
                        if (sdtd_AZ < 0){
                            sdtd_AZ <- 0
                        }
                        fdtd_AZ <- fdtd_AZ - sdtd_AZ
                        if (fdtd_AZ < 0){
                            fdtd_AZ <- 0
                        }
                        
                        first_doses_then_PZ <- first_doses_then[first_doses_then$vaccine == "PF",]$cum_doses
                        last_day_second_doses_PZ <- sdPZ_data_lastday$cum_doses
                        sdtd_PZ <- first_doses_then_PZ - last_day_second_doses_PZ
                        if (sdtd_PZ < 0){
                            sdtd_PZ <- 0
                        }
                        fdtd_PZ <- fdtd_PZ - sdtd_PZ
                        if (fdtd_PZ < 0){
                            fdtd_PZ <- 0
                        }
                        
                        first_doses_then_MD <- first_doses_then[first_doses_then$vaccine == "MD",]$cum_doses
                        last_day_second_doses_MD <- sdMD_data_lastday$cum_doses
                        sdtd_MD <- first_doses_then_MD - last_day_second_doses_MD
                        if (sdtd_MD < 0){
                            sdtd_MD <- 0
                        }
                        fdtd_MD <- fdtd_MD - sdtd_MD
                        if (fdtd_MD < 0){
                            fdtd_MD <- 0
                        }
                        
                    } else {
                        sdtd_AZ <- 0
                        sdtd_PZ <- 0
                        sdtd_MD <- 0
                    }
                    
                    
                    # limit doses to the uptake threshold of each age group
                    
                    # deliver as many as possible first doses for AZ, Pfizer, and Moderna,
                    # limiting the number of doses to uptake threshold
                    if (fdAZ_data_lastday$left_to_dose_ut > 0){
                        if (fdtd_AZ + fdtd_PZ + fdtd_MD > fdAZ_data_lastday$left_to_dose_ut){
                            # split left_to_dose_ut into proportions of AZ, Pfizer and Moderna
                            final_fdtd_AZ <- round(products$AZ[k] *  fdAZ_data_lastday$left_to_dose_ut)
                            final_fdtd_PZ <- round(products$PZ[k] *  fdAZ_data_lastday$left_to_dose_ut)
                            final_fdtd_MD <- round(products$MD[k] *  fdAZ_data_lastday$left_to_dose_ut)
                            # shift remaining first doses to second doses
                            sdtd_AZ <- sdtd_AZ + (fdtd_AZ - final_fdtd_AZ)
                            sdtd_PZ <- sdtd_PZ + (fdtd_PZ - final_fdtd_PZ)
                            sdtd_MD <- sdtd_MD + (fdtd_MD - final_fdtd_MD)
                            if (sdtd_AZ < 0){
                                sdtd_AZ <- 0
                            }
                            if (sdtd_PZ < 0){
                                sdtd_PZ <- 0
                            }
                            if (sdtd_MD < 0){
                                sdtd_MD <- 0
                            }
                            # deliver final first doses for AZ, Pfizer and Moderna
                            new_entries <- data.frame(region = rep(this_region, 3), 
                                                      age.group = rep(k, 3), 
                                                      ages = rep(this_age_group, 3), 
                                                      vaccine = c("AZ", "PF", "MD"),
                                                      dose = rep("First", 3), 
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
                                                      dose = rep("First", 3), 
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
                                                  dose = rep("First", 3), 
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
                        # shift all first doses to equivalent second doses
                        sdtd_AZ <- sdtd_AZ + fdtd_AZ
                        sdtd_PZ <- sdtd_PZ + fdtd_PZ
                        sdtd_MD <- sdtd_MD + fdtd_MD
                    }
                    
                    # deliver as many as possible second doses for AZ, Pfizer and Moderna,
                    # limiting the number of doses to relevant uptake threshold OR the value
                    # of the relevant left_to_dose value for first doses 12 weeks ago
                    left_to_dose_then <- min(first_dose_data$left_to_dose_ut)
                    if (sdAZ_data_lastday$left_to_dose_ut > left_to_dose_then){
                        diff <- sdAZ_data_lastday$left_to_dose_ut - left_to_dose_then
                        if (diff < 0){
                            diff <- 0
                        }
                        if (sdtd_AZ + sdtd_PZ + sdtd_MD > diff){
                            # split diff into specified proportions of AZ, Pfizer and Moderna
                            final_sdtd_AZ <- round(products$AZ[k] * diff)
                            if (final_sdtd_AZ < 0){
                                final_sdtd_AZ <- 0
                            }
                            final_sdtd_PZ <- round(products$PZ[k] * diff)
                            if (final_sdtd_PZ < 0){
                                final_sdtd_PZ <- 0
                            }
                            final_sdtd_MD <- round(products$MD[k] * diff)
                            if (final_sdtd_MD < 0){
                                final_sdtd_MD <- 0
                            }
                            # shift remaining second doses to...
                            if (k > 1){
                                # next age group down
                                AZ_doses[k-1] <- AZ_doses[k-1] + (sdtd_AZ - final_sdtd_AZ)
                                PZ_doses[k-1] <- PZ_doses[k-1] + (sdtd_PZ - final_sdtd_PZ)
                                MD_doses[k-1] <- MD_doses[k-1] + (sdtd_MD - final_sdtd_MD)
                            } else if (j < length(regions)){
                                # next region along
                                doses_per_region[j+1] <- doses_per_region[j+1] + (sdtd_AZ + sdtd_PZ + sdtd_MD - diff)
                            } else if (i < length(vpdates)){
                                # next date along
                                vpdoses[i+1] <- vpdoses[i+1] + (sdtd_AZ + sdtd_PZ + sdtd_MD - diff)
                            } else {
                                # or record leftover doses
                                leftover_doses <- leftover_doses + (sdtd_AZ + sdtd_PZ + sdtd_MD - diff)
                            }
                            # deliver final second doses for AZ, Pfizer and Moderna
                            new_entries <- data.frame(region = rep(this_region, 3), 
                                                      age.group = rep(k, 3), 
                                                      ages = rep(this_age_group, 3), 
                                                      vaccine = c("AZ", "PF", "MD"),
                                                      dose = rep("Second", 3), 
                                                      date = rep(as.Date(date_today), 3),
                                                      number = c(final_sdtd_AZ, final_sdtd_PZ, final_sdtd_MD),
                                                      cum_doses = c(sdAZ_data_lastday$cum_doses + final_sdtd_AZ, 
                                                                    sdPZ_data_lastday$cum_doses + final_sdtd_PZ,
                                                                    sdMD_data_lastday$cum_doses + final_sdtd_MD),
                                                      pop_size = rep(sdAZ_data_lastday$pop_size, 3),
                                                      max_uptake = rep(sdAZ_data_lastday$max_uptake, 3),
                                                      left_to_dose = rep(sdAZ_data_lastday$left_to_dose - diff, 3),
                                                      left_to_dose_ut = rep(sdAZ_data_lastday$left_to_dose_ut - diff, 3))
                            # bind new entries to final data frame
                            final_df <- rbind(final_df, new_entries)
                        } else { # sdtd_AZ + sdtd_PZ + sdtd_MD <= diff
                            # deliver second doses for AZ, Pfizer and Moderna
                            new_entries <- data.frame(region = rep(this_region, 3), 
                                                      age.group = rep(k, 3), 
                                                      ages = rep(this_age_group, 3), 
                                                      vaccine = c("AZ", "PF", "MD"),
                                                      dose = rep("Second", 3), 
                                                      date = rep(as.Date(date_today), 3),
                                                      number = c(sdtd_AZ, sdtd_PZ, sdtd_MD),
                                                      cum_doses = c(sdAZ_data_lastday$cum_doses + sdtd_AZ, 
                                                                    sdPZ_data_lastday$cum_doses + sdtd_PZ,
                                                                    sdMD_data_lastday$cum_doses + sdtd_MD),
                                                      pop_size = rep(sdAZ_data_lastday$pop_size, 3),
                                                      max_uptake = rep(sdAZ_data_lastday$max_uptake, 3),
                                                      left_to_dose = rep(sdAZ_data_lastday$left_to_dose - (sdtd_AZ+sdtd_PZ+sdtd_MD), 3),
                                                      left_to_dose_ut = rep(sdAZ_data_lastday$left_to_dose_ut - (sdtd_AZ+sdtd_PZ+sdtd_MD), 3))
                            # bind new entries to final data frame
                            final_df <- rbind(final_df, new_entries)
                        }
                    } else {
                        # sdAZ_data_lastday$left_to_dose_ut <= left_to_dose_then
                        # shift all remaining second doses to...
                        if (k > 1){
                            # next age group down
                            AZ_doses[k-1] <- AZ_doses[k-1] + sdtd_AZ
                            PZ_doses[k-1] <- PZ_doses[k-1] + sdtd_PZ
                            MD_doses[k-1] <- MD_doses[k-1] + sdtd_MD
                        } else if (j < length(regions)){
                            # next region along
                            doses_per_region[j+1] <- doses_per_region[j+1] + sdtd_AZ + sdtd_PZ + sdtd_MD
                        } else if (i < length(vpdates)){
                            # next date along
                            vpdoses[i+1] <- vpdoses[i+1] + sdtd_AZ + sdtd_PZ + sdtd_MD
                        } else {
                            # or record leftover doses
                            leftover_doses <- leftover_doses + sdtd_AZ + sdtd_PZ + sdtd_MD
                        }
                        # deliver zero second doses for AZ and Pfizer
                        new_entries <- data.frame(region = rep(this_region, 3), 
                                                  age.group = rep(k, 3), 
                                                  ages = rep(this_age_group, 3), 
                                                  vaccine = c("AZ", "PF", "MD"),
                                                  dose = rep("Second", 3), 
                                                  date = rep(as.Date(date_today), 3),
                                                  number = rep(0, 3),
                                                  cum_doses = c(sdAZ_data_lastday$cum_doses, 
                                                                sdPZ_data_lastday$cum_doses,
                                                                sdMD_data_lastday$cum_doses),
                                                  pop_size = rep(sdAZ_data_lastday$pop_size, 3),
                                                  max_uptake = rep(sdAZ_data_lastday$max_uptake, 3),
                                                  left_to_dose = rep(sdAZ_data_lastday$left_to_dose, 3),
                                                  left_to_dose_ut = rep(sdAZ_data_lastday$left_to_dose_ut, 3))
                        # bind new entries to final data frame
                        final_df <- rbind(final_df, new_entries)
                    }
                    
                    
                }
            }
        }
        
        # sort final_df data frame rows by region, age group, vaccine, dose, date
        vax_schedule <- final_df[order(final_df$region, final_df$age.group, 
                                       final_df$vaccine, final_df$dose, 
                                       final_df$date),]
        
    } else if (second_doses == FALSE){
        
        # second_doses == FALSE so doses contains first doses only
        
        # calculate list of dates and number of doses *per day*
        vpdates <- seq(start_date, end_date, by = "day")
        vpdoses <- rep(doses[1]/7, length(vpdates))
        
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
        oldest_age <- max(vax$age, na.rm = TRUE)
        min_ages <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75)
        max_ages <- c(4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,oldest_age)
        ages <- data.frame(groups = agegroups, min_ages = min_ages, max_ages = max_ages)
        breaks = c(0,4.5,9.5,14.5,19.5,24.5,29.5,34.5,39.5,44.5,49.5,54.5,59.5,64.5,
                   69.5,74.5,oldest_age)
        hist <- hist(vax$age, freq = FALSE, breaks)
        nbins <- length(breaks) - 1
        binwidths <- diff(breaks)
        agegroup_probs <- binwidths * hist$density
        
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
        
    } else {
        stop('second_doses should be TRUE or FALSE')
    }
    
    return(vax_schedule)
}

# function to convert existing data frame based schedule into list for covidm
convert_schedule <- function(vax_schedule, params){
    
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
    
    return(ret)
}

# function to calculate the proportion of each vaccine type delivered to date
products_delivered <- function(vax_delivered, age_group_limit = 11){
    
    if (length(unique(vax_delivered$vaccine)) != 3){
        stop('Input data does not have 3 vaccine products')
    }
    
    # get subset of data only
    vax_subset <- vax_delivered[(vax_delivered$age.group >= age_group_limit),]
    
    AZ <- vax_subset[(vax_subset$vaccine == "AZ"),]
    if (dim(AZ)[1] == 0){
        stop('No entries in AZ data frame - check vaccine name in input data')
    }
    PZ <- vax_subset[(vax_subset$vaccine == "PF"),]
    if (dim(PZ)[1] == 0){
        stop('No entries in PZ data frame - check vaccine name in input data')
    }
    MD <- vax_subset[(vax_subset$vaccine == "MD"),]
    if (dim(MD)[1] == 0){
        stop('No entries in MD data frame - check vaccine name in input data')
    }
    
    total_AZ <- sum(AZ$number)
    total_PZ <- sum(PZ$number)
    total_MD <- sum(MD$number)
    total_all <- total_AZ + total_PZ + total_MD
    
    prop_AZ <- total_AZ / total_all
    prop_PZ <- total_PZ / total_all
    prop_MD <- total_MD / total_all
    
    AZd1 <- AZ[(AZ$dose == 1),]
    if (dim(AZd1)[1] == 0){
        stop('No entries in AZd1 data frame - check dose name in input data')
    }
    PZd1 <- PZ[(PZ$dose == 1),]
    if (dim(PZd1)[1] == 0){
        stop('No entries in PZd1 data frame - check dose name in input data')
    }
    MDd1 <- MD[(MD$dose == 1),]
    if (dim(MDd1)[1] == 0){
        stop('No entries in MDd1 data frame - check dose name in input data')
    }
    
    total_AZd1 <- sum(AZd1$number)
    total_PZd1 <- sum(PZd1$number)
    total_MDd1 <- sum(MDd1$number)
    total_alld1 <- total_AZd1 + total_PZd1 + total_MDd1
    
    prop_AZd1 <- total_AZd1 / total_alld1
    prop_PZd1 <- total_PZd1 / total_alld1
    prop_MDd1 <- total_MDd1 / total_alld1
    
    to_return <- data.frame(doses = c('All', 'First only'), 
                            AZ = c(prop_AZ, prop_AZd1),
                            PZ = c(prop_PZ, prop_PZd1),
                            MD = c(prop_MD, prop_MDd1))
    
    return(to_return)
}

# function to calculate England-wide uptake based on existed vaccines delivered
calculate_uptake <- function(vax_delivered, popUK){
    
    setDT(vax_delivered)
    england = vax_delivered[, sum(number), by = .(age.group, dose, date)]
    england[, cum_doses := cumsum(V1), by = .(age.group, dose)]
    england = england[england$date == max(england$date)]
    
    popeng = popUK[popUK$name == 'England',]
    popeng$total = (popeng$f + popeng$m)*1000
    popeng$covidm_age = popeng$age
    popeng[age %in% c('75-79', '80-84', '85-89', '90+'), covidm_age := '75+']
    popeng[, covidm_pop := sum(total), by = covidm_age]
    popeng = data.table(age = popeng$covidm_age[1:16], pop = popeng$covidm_pop[1:16])
    popeng$age.group = (1:16)
    
    merged = merge(england, popeng, # dataset names
          by = 'age.group', # by common variable
          all = FALSE
    )
    merged[, coverage := cum_doses / pop]
    merged$V1 = NULL
   
    return(merged)
}

# function to calculate supply of first vaccine doses over time
first_dose_supply <- function(vax_delivered){
    
    # sum data over England and over age groups and get first doses only
    setDT(vax_delivered)
    england = vax_delivered[, sum(number), by = .(dose, date)]
    england = england[england$dose == 1,]

    return(england)    
}
