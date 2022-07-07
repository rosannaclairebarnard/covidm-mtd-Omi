library(lubridate)
library(data.table)
library(stringr)
library(ggplot2)

###############################################################################
# writing a function to plot vaccine schedule output from an .rds file

process_vax_schedule <- function(path_to_schedule){
    
    vax_schedule = readRDS(path_to_schedule)
    num_to_loop = length(vax_schedule[[1]]$va1)
    region_idxs = c(1,3,4,5,6,9,10)
    master_df = data.table(date = NULL, 
                           region = NULL, 
                           age.group = NULL, 
                           vaccine = NULL,
                           dose = NULL,
                           num = NULL)
    for (i in 1:num_to_loop){
        day = vax_schedule[[1]]$vt[i]
        for (j in 1:16){
            age_sum = rep(0,4)
            for (k in region_idxs){
                num_va1 = vax_schedule[[k]]$va1[[i]][j]
                num_va2 = vax_schedule[[k]]$va2[[i]][j]
                num_vb1 = vax_schedule[[k]]$vb1[[i]][j]
                num_vb2 = vax_schedule[[k]]$vb2[[i]][j]
                df = data.table(date = rep(day,4), 
                                region = rep(k,4), 
                                age.group = rep(j,4),
                                vaccine = c('AZ', 'AZ', 'Pfizer', 'Pfizer'),
                                dose = c('First', 'Second', 'First', 'Second'),
                                num = c(num_va1, num_va2, num_vb1, num_vb2))
                master_df = rbind(master_df, df)
                age_sum[1] = age_sum[1] + num_va1
                age_sum[2] = age_sum[2] + num_va2
                age_sum[3] = age_sum[3] + num_vb1
                age_sum[4] = age_sum[4] + num_vb2
            }
        }
    }
    return(master_df)   
}

###############################################################################

# process and save comparison of vaccination schedules (5th May 2022)

scheda = process_vax_schedule('./fitting_data/vax-covidm20220505205235.rds')
scheda$title = '5plus_80%'

schedb = process_vax_schedule('./fitting_data/vax-covidm202205052222505plus_50percent.rds')
schedb$title = '5plus_50%'

two_scheds = rbind(scheda, schedb)

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
write.csv(two_scheds, paste0('./output/paper/apr22/vaccine_schedules_', datetime, '.csv'), row.names = FALSE)

eng_sched = two_scheds[, .(num = sum(num)), by = .(date, age.group, dose, title)]

eng_sched$age[eng_sched$age.group == 1]  = '0-4'
eng_sched$age[eng_sched$age.group == 2]  = '5-9' 
eng_sched$age[eng_sched$age.group == 3]  = '10-14'
eng_sched$age[eng_sched$age.group == 4]  = '15-19'
eng_sched$age[eng_sched$age.group == 5]  = '20-24'
eng_sched$age[eng_sched$age.group == 6]  = '25-29'
eng_sched$age[eng_sched$age.group == 7]  = '30-34'
eng_sched$age[eng_sched$age.group == 8]  = '35-39'
eng_sched$age[eng_sched$age.group == 9]  = '40-44'
eng_sched$age[eng_sched$age.group == 10] = '45-49'
eng_sched$age[eng_sched$age.group == 11] = '50-54'
eng_sched$age[eng_sched$age.group == 12] = '55-59'
eng_sched$age[eng_sched$age.group == 13] = '60-64'
eng_sched$age[eng_sched$age.group == 14] = '65-69'
eng_sched$age[eng_sched$age.group == 15] = '70-74'
eng_sched$age[eng_sched$age.group == 16] = '75+' 
eng_sched$age = factor(eng_sched$age, levels = c('0-4'   ,
                                                 '5-9' ,
                                                 '10-14',
                                                 '15-19',
                                                 '20-24',
                                                 '25-29',
                                                 '30-34',
                                                 '35-39',
                                                 '40-44',
                                                 '45-49',
                                                 '50-54',
                                                 '55-59',
                                                 '60-64',
                                                 '65-69',
                                                 '70-74',
                                                 '75+' ))


path_to_newcovidvax  <- "."
uk_covid_data_path   <- paste0(path_to_newcovidvax, "/fitting_data/")
datapath = function(x) paste0(uk_covid_data_path, x)
popUK = readRDS(datapath("popNHS.rds"))
engpop = popUK[popUK$name == 'England',]
engpop$total = (engpop$f + engpop$m)*1000
over75s = sum(engpop$total[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
over75sf = sum(engpop$f[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
over75sm = sum(engpop$m[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
new_row = data.table(name = 'England',
                     age = '75+',
                     f = over75sf,
                     m = over75sm,
                     location_type = 5,
                     country_code = 2222,
                     total = over75s)
engpop = rbind(engpop, new_row)
engpop = engpop[!(engpop$age %in% c('75-79', '80-84', '85-89', '90+')),]
engpop$age.group = NULL
engpop$age.group[engpop$age == '0-4'] = 1
engpop$age.group[engpop$age == '5-9'] = 2
engpop$age.group[engpop$age == '10-14'] = 3
engpop$age.group[engpop$age == '15-19'] = 4
engpop$age.group[engpop$age == '20-24'] = 5
engpop$age.group[engpop$age == '25-29'] = 6
engpop$age.group[engpop$age == '30-34'] = 7
engpop$age.group[engpop$age == '35-39'] = 8
engpop$age.group[engpop$age == '40-44'] = 9
engpop$age.group[engpop$age == '45-49'] = 10
engpop$age.group[engpop$age == '50-54'] = 11
engpop$age.group[engpop$age == '55-59'] = 12
engpop$age.group[engpop$age == '60-64'] = 13
engpop$age.group[engpop$age == '65-69'] = 14
engpop$age.group[engpop$age == '70-74'] = 15
engpop$age.group[engpop$age == '75+'] = 16

all = merge(eng_sched[, .(date, age.group, dose, title, num, age)], engpop[, .(total, age.group)], by = "age.group")

all[,cumulative := cumsum(num), by = .(age.group, dose, title)]
all$coverage = (all$cumulative / all$total)*100


all$title[all$title == '5plus_80%'] = '5 plus, 80% uptake'
all$title[all$title == '5plus_50%'] = '5 plus, 50% uptake'
all$title = factor(all$title, levels = c( '5 plus, 80% uptake',
                                          '5 plus, 50% uptake'))

ggplot(all[all$dose == 'First' & all$date > '2022-01-01' & all$date <= '2022-12-31' & all$age.group %in% c(2,3,4),]) + 
    geom_line(aes(date, coverage, colour=factor(age), group = age)) +
    facet_wrap(~title, ncol = 1) +
    scale_x_date(date_breaks = '2 months', date_labels = "%m/%y", expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = 'Vaccination coverage (%)', colour = 'Age group') +
    cowplot::theme_cowplot(font_size = 10) +
    theme(strip.background = element_blank(), legend.position = c(0.8, 0.9),
          legend.background = element_rect(fill = "white", # Background
                                           colour = 1),
          legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) +
    geom_vline(aes(xintercept = as.Date('2020-01-01')+765), size = .5, linetype = "dotted")

ggplot(all[all$dose == 'First' & all$date <= '2022-05-06',]) + 
    geom_line(aes(date, coverage, colour=factor(age), group = age)) +
    facet_wrap(~title, ncol = 1) +
    scale_x_date(date_breaks = '2 months', date_labels = "%m/%y", expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = 'Vaccination coverage (%)', colour = 'Age group') +
    cowplot::theme_cowplot(font_size = 10) +
    theme(strip.background = element_blank(), legend.position = c(0.8, 0.9),
          legend.background = element_rect(fill = "white", # Background
                                           colour = 1),
          legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) +
    geom_vline(aes(xintercept = as.Date('2020-01-01')+765), size = .5, linetype = "dotted")

###############################################################################

# process and save comparison of new vs. old vaccinations schedules (March 2022)
four_old_scheds = read.csv('./output/paper/jan22/vaccine_schedules_20220221131941.csv')
unique(four_old_scheds$title)
old_sched = four_old_scheds[four_old_scheds$title == '12plus_80%',]

# process new schedule
new_sched = process_vax_schedule('./fitting_data/vax-covidm20220330102111.rds')
new_sched$title = '5plus_80%'

old_sched$date = as.Date(old_sched$date)
two_scheds = rbind(old_sched, new_sched)

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
write.csv(two_scheds, paste0('./output/paper/marchapril22/vaccine_schedules_', datetime, '.csv'), row.names = FALSE)

old = ggplot(old_sched[old_sched$vaccine == 'Pfizer' & old_sched$date > '2022-04-01',]) + 
    geom_line(aes(as.Date(date), num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) +
    facet_wrap(~region, ncol = 1, scales = "free") 
old

new = ggplot(new_sched[new_sched$vaccine == 'Pfizer' & new_sched$date > '2022-04-01',]) + 
    geom_line(aes(as.Date(date), num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) +
    facet_wrap(~region, ncol = 1, scales = "free")
new

cowplot::plot_grid(old + scale_x_date(limits = as.Date(c('2022-03-21', '2022-09-30'))), 
                   new + scale_x_date(limits = as.Date(c('2022-03-21', '2022-09-30'))), ncol = 2)

old2 = ggplot(old_sched[old_sched$vaccine == 'Pfizer' & old_sched$date < '2022-03-21',]) + 
    geom_line(aes(as.Date(date), num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) +
    facet_wrap(~region, ncol = 1, scales = "free") 
old2

new2 = ggplot(new_sched[new_sched$vaccine == 'Pfizer' & new_sched$date < '2022-03-21',]) + 
    geom_line(aes(as.Date(date), num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) +
    facet_wrap(~region, ncol = 1, scales = "free")
new2

cowplot::plot_grid(old2 + scale_x_date(limits = as.Date(c('2021-01-01', '2022-03-21'))), 
                   new2 + scale_x_date(limits = as.Date(c('2021-01-01', '2022-03-21'))), ncol = 2)

ggplot(new_sched[new_sched$vaccine == 'AZ' & new_sched$date > '2022-04-01',]) + 
    geom_line(aes(as.Date(date), num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) +
    facet_wrap(~region, ncol = 1, scales = "free")

setDT(two_scheds)
eng_sched = two_scheds[, .(num = sum(num)), by = .(date, age.group, dose, title)]

eng_sched$age[eng_sched$age.group == 1]  = '0-4'
eng_sched$age[eng_sched$age.group == 2]  = '5-9' 
eng_sched$age[eng_sched$age.group == 3]  = '10-14'
eng_sched$age[eng_sched$age.group == 4]  = '15-19'
eng_sched$age[eng_sched$age.group == 5]  = '20-24'
eng_sched$age[eng_sched$age.group == 6]  = '25-29'
eng_sched$age[eng_sched$age.group == 7]  = '30-34'
eng_sched$age[eng_sched$age.group == 8]  = '35-39'
eng_sched$age[eng_sched$age.group == 9]  = '40-44'
eng_sched$age[eng_sched$age.group == 10] = '45-49'
eng_sched$age[eng_sched$age.group == 11] = '50-54'
eng_sched$age[eng_sched$age.group == 12] = '55-59'
eng_sched$age[eng_sched$age.group == 13] = '60-64'
eng_sched$age[eng_sched$age.group == 14] = '65-69'
eng_sched$age[eng_sched$age.group == 15] = '70-74'
eng_sched$age[eng_sched$age.group == 16] = '75+' 
eng_sched$age = factor(eng_sched$age, levels = c('0-4'   ,
                                                 '5-9' ,
                                                 '10-14',
                                                 '15-19',
                                                 '20-24',
                                                 '25-29',
                                                 '30-34',
                                                 '35-39',
                                                 '40-44',
                                                 '45-49',
                                                 '50-54',
                                                 '55-59',
                                                 '60-64',
                                                 '65-69',
                                                 '70-74',
                                                 '75+' ))

path_to_newcovidvax  <- "."
uk_covid_data_path   <- paste0(path_to_newcovidvax, "/fitting_data/")
datapath = function(x) paste0(uk_covid_data_path, x)
popUK = readRDS(datapath("popNHS.rds"))
engpop = popUK[popUK$name == 'England',]
engpop$total = (engpop$f + engpop$m)*1000
over75s = sum(engpop$total[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
over75sf = sum(engpop$f[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
over75sm = sum(engpop$m[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
new_row = data.table(name = 'England',
                     age = '75+',
                     f = over75sf,
                     m = over75sm,
                     location_type = 5,
                     country_code = 2222,
                     total = over75s)
engpop = rbind(engpop, new_row)
engpop = engpop[!(engpop$age %in% c('75-79', '80-84', '85-89', '90+')),]
engpop$age.group = NULL
engpop$age.group[engpop$age == '0-4'] = 1
engpop$age.group[engpop$age == '5-9'] = 2
engpop$age.group[engpop$age == '10-14'] = 3
engpop$age.group[engpop$age == '15-19'] = 4
engpop$age.group[engpop$age == '20-24'] = 5
engpop$age.group[engpop$age == '25-29'] = 6
engpop$age.group[engpop$age == '30-34'] = 7
engpop$age.group[engpop$age == '35-39'] = 8
engpop$age.group[engpop$age == '40-44'] = 9
engpop$age.group[engpop$age == '45-49'] = 10
engpop$age.group[engpop$age == '50-54'] = 11
engpop$age.group[engpop$age == '55-59'] = 12
engpop$age.group[engpop$age == '60-64'] = 13
engpop$age.group[engpop$age == '65-69'] = 14
engpop$age.group[engpop$age == '70-74'] = 15
engpop$age.group[engpop$age == '75+'] = 16

all = merge(eng_sched[, .(date, age.group, dose, title, num, age)], engpop[, .(total, age.group)], by = "age.group")

all[,cumulative := cumsum(num), by = .(age.group, dose, title)]
all$coverage = (all$cumulative / all$total)*100


all$title[all$title == '12plus_80%'] = '12 plus, 80% uptake'
all$title[all$title == '5plus_80%'] = '5 plus, 80% uptake'
all$title = factor(all$title, levels = c( '12 plus, 80% uptake',
                                          '5 plus, 80% uptake'))

ggplot(all[all$dose == 'First' & all$date > '2022-03-20' & all$date < '2022-12-01' & all$age.group %in% c(2,3,4),]) + 
    geom_line(aes(as.Date(date), coverage, colour=factor(age), group = age)) +
    facet_wrap(~title, ncol = 1) +
    scale_x_date(date_breaks = '2 months', date_labels = "%m/%y", expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = 'Vaccination coverage (%)', colour = 'Age group') +
    cowplot::theme_cowplot(font_size = 10) +
    theme(strip.background = element_blank(), legend.position = c(0.8, 0.9),
          legend.background = element_rect(fill = "white", # Background
                                           colour = 1),
          legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) +
    geom_vline(aes(xintercept = as.Date('2020-01-01')+765), size = .5, linetype = "dotted")



# process and save summary of existing vaccination schedules (February 2022)

sched1 = process_vax_schedule('./fitting_data/vax-covidm20220117203900.rds')
sched2 = process_vax_schedule('./fitting_data/vax-covidm2022011720390012plus_50percent.rds')
sched3 = process_vax_schedule('./fitting_data/vax-covidm202201172039005plus_50percent.rds')
sched4 = process_vax_schedule('./fitting_data/vax-covidm202201172039005plus_80percent.rds')

sched1$title = '12plus_80%'
sched2$title = '12plus_50%'
sched3$title = '5plus_50%'
sched4$title = '5plus_80%'

four_scheds = rbind(sched1, sched2, sched3, sched4)
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
write.csv(four_scheds, paste0('./output/paper/jan22/vaccine_schedules_', datetime, '.csv'), row.names = FALSE)

# calculate sums over all NHS England regions

eng_sched = four_scheds[, .(num = sum(num)), by = .(date, age.group, dose, title)]

eng_sched$age[eng_sched$age.group == 1]  = '0-4'
eng_sched$age[eng_sched$age.group == 2]  = '5-9' 
eng_sched$age[eng_sched$age.group == 3]  = '10-14'
eng_sched$age[eng_sched$age.group == 4]  = '15-19'
eng_sched$age[eng_sched$age.group == 5]  = '20-24'
eng_sched$age[eng_sched$age.group == 6]  = '25-29'
eng_sched$age[eng_sched$age.group == 7]  = '30-34'
eng_sched$age[eng_sched$age.group == 8]  = '35-39'
eng_sched$age[eng_sched$age.group == 9]  = '40-44'
eng_sched$age[eng_sched$age.group == 10] = '45-49'
eng_sched$age[eng_sched$age.group == 11] = '50-54'
eng_sched$age[eng_sched$age.group == 12] = '55-59'
eng_sched$age[eng_sched$age.group == 13] = '60-64'
eng_sched$age[eng_sched$age.group == 14] = '65-69'
eng_sched$age[eng_sched$age.group == 15] = '70-74'
eng_sched$age[eng_sched$age.group == 16] = '75+' 
eng_sched$age = factor(eng_sched$age, levels = c('0-4'   ,
                                                 '5-9' ,
                                                 '10-14',
                                                 '15-19',
                                                 '20-24',
                                                 '25-29',
                                                 '30-34',
                                                 '35-39',
                                                 '40-44',
                                                 '45-49',
                                                 '50-54',
                                                 '55-59',
                                                 '60-64',
                                                 '65-69',
                                                 '70-74',
                                                 '75+' ))


path_to_newcovidvax  <- "."
uk_covid_data_path   <- paste0(path_to_newcovidvax, "/fitting_data/")
datapath = function(x) paste0(uk_covid_data_path, x)
popUK = readRDS(datapath("popNHS.rds"))
engpop = popUK[popUK$name == 'England',]
engpop$total = (engpop$f + engpop$m)*1000
over75s = sum(engpop$total[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
over75sf = sum(engpop$f[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
over75sm = sum(engpop$m[engpop$age %in% c('75-79', '80-84', '85-89', '90+')])
new_row = data.table(name = 'England',
                     age = '75+',
                     f = over75sf,
                     m = over75sm,
                     location_type = 5,
                     country_code = 2222,
                     total = over75s)
engpop = rbind(engpop, new_row)
engpop = engpop[!(engpop$age %in% c('75-79', '80-84', '85-89', '90+')),]
engpop$age.group = NULL
engpop$age.group[engpop$age == '0-4'] = 1
engpop$age.group[engpop$age == '5-9'] = 2
engpop$age.group[engpop$age == '10-14'] = 3
engpop$age.group[engpop$age == '15-19'] = 4
engpop$age.group[engpop$age == '20-24'] = 5
engpop$age.group[engpop$age == '25-29'] = 6
engpop$age.group[engpop$age == '30-34'] = 7
engpop$age.group[engpop$age == '35-39'] = 8
engpop$age.group[engpop$age == '40-44'] = 9
engpop$age.group[engpop$age == '45-49'] = 10
engpop$age.group[engpop$age == '50-54'] = 11
engpop$age.group[engpop$age == '55-59'] = 12
engpop$age.group[engpop$age == '60-64'] = 13
engpop$age.group[engpop$age == '65-69'] = 14
engpop$age.group[engpop$age == '70-74'] = 15
engpop$age.group[engpop$age == '75+'] = 16

all = merge(eng_sched[, .(date, age.group, dose, title, num, age)], engpop[, .(total, age.group)], by = "age.group")

all[,cumulative := cumsum(num), by = .(age.group, dose, title)]
all$coverage = (all$cumulative / all$total)*100


all$title[all$title == '12plus_80%'] = '12 plus, 80% uptake'
all$title[all$title == '12plus_50%'] = '12 plus, 50% uptake'
all$title[all$title == '5plus_80%'] = '5 plus, 80% uptake'
all$title[all$title == '5plus_50%'] = '5 plus, 50% uptake'
all$title = factor(all$title, levels = c( '12 plus, 80% uptake',
                                          '12 plus, 50% uptake',
                                          '5 plus, 80% uptake',
                                          '5 plus, 50% uptake'))

ggplot(all[all$dose == 'First' & all$date > '2022-01-01' & all$date < '2022-09-30' & all$age.group %in% c(2,3,4),]) + 
    geom_line(aes(date, coverage, colour=factor(age), group = age)) +
    facet_wrap(~title, ncol = 1) +
    scale_x_date(date_breaks = '2 months', date_labels = "%m/%y", expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = 'Vaccination coverage (%)', colour = 'Age group') +
    cowplot::theme_cowplot(font_size = 10) +
    theme(strip.background = element_blank(), legend.position = c(0.8, 0.9),
          legend.background = element_rect(fill = "white", # Background
                                           colour = 1),
          legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) +
    geom_vline(aes(xintercept = as.Date('2020-01-01')+765), size = .5, linetype = "dotted")
    
    
    
    
    


# calculate cumulative doses by different variables

# master_df_new[,cum_doses := cumsum(num), by =.(age.group, region, vaccine, dose)]
# master_df_new[,overall_doses := cumsum(num), by =.(region, dose)]
# 
# master_df_old[,cum_doses := cumsum(num), by =.(age.group, region, vaccine, dose)]
# master_df_old[,overall_doses := cumsum(num), by =.(region, dose)]
# 
# overall_df = rbind(master_df_old, master_df_new)

###############################################################################
# plot output showing different vaccine schedules

ggplot(four_scheds[four_scheds$region == 1 & four_scheds$vaccine == 'Pfizer' & four_scheds$date > '2022-01-01',]) + 
    geom_line(aes(date, num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) +
    facet_wrap(~title, scales = "free") +
    scale_x_date(limits = as.Date(c('2022-01-01','2022-09-30')))

ggplot(eng_sched[eng_sched$dose == 'First' & eng_sched$date > '2022-01-01',]) + 
    geom_line(aes(date, num, colour=factor(age), group = age)) +
    facet_wrap(~title, scales = "free") +
    scale_x_date(limits = as.Date(c('2022-01-01','2022-09-30')))

###############################################################################

plot1 <- function(master_df){
    ggplot(master_df[master_df$vaccine == 'Pfizer',]) + 
        geom_line(aes(date, num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) + 
        facet_wrap(~region, scales = "free")
}
plot1(sched1)
sched_list = list(sched1, sched2, sched3, sched4)
compare_schedules <- function(schedule_list, cut_off_date){
    master_df = data.table(date = NULL, 
                           region = NULL,
                           age.group = NULL,
                           vaccine = NULL,
                           dose = NULL,
                           num = NULL,
                           schedule = NULL)
    num_to_loop = length(schedule_list)
    for (i in 1:num_to_loop){
        this_sched = schedule_list[[i]]
        this_sched$schedule = rep(i, dim(this_sched)[1])
        master_df = rbind(master_df, this_sched)
    }
    # now plot with facet_wrap according to schedule number and region
    ggplot(master_df[master_df$vaccine == 'Pfizer' & master_df$date > cut_off_date & master_df$age.group < 5,]) + 
        geom_line(aes(date, num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) + 
        facet_wrap(schedule~region, scales = "free", ncol = 7)
    
}
compare_schedules(sched_list,as.Date('2022-01-01'))
###############################################################################

vax_schedule = readRDS('./fitting_data/vax-covidm20210925084913.rds')


num_to_loop = length(vax_schedule[[1]]$va1)

region_idxs = c(1,3,4,5,6,9,10)

master_df = data.table(date = NULL, 
                       region = NULL, 
                       age.group = NULL, 
                       vaccine = NULL,
                       dose = NULL,
                       num = NULL)

for (i in 1:num_to_loop){
    day = vax_schedule[[1]]$vt[i]
    for (j in 1:16){
        age_sum = rep(0,4)
        for (k in region_idxs){
            num_va1 = vax_schedule[[k]]$va1[[i]][j]
            num_va2 = vax_schedule[[k]]$va2[[i]][j]
            num_vb1 = vax_schedule[[k]]$vb1[[i]][j]
            num_vb2 = vax_schedule[[k]]$vb2[[i]][j]
            df = data.table(date = rep(day,4), 
                            region = rep(k,4), 
                            age.group = rep(j,4),
                            vaccine = c('AZ', 'AZ', 'Pfizer', 'Pfizer'),
                            dose = c('First', 'Second', 'First', 'Second'),
                            num = c(num_va1, num_va2, num_vb1, num_vb2))
            master_df = rbind(master_df, df)
            age_sum[1] = age_sum[1] + num_va1
            age_sum[2] = age_sum[2] + num_va2
            age_sum[3] = age_sum[3] + num_vb1
            age_sum[4] = age_sum[4] + num_vb2
        }
    }
}

# plot region 1 only, use facet wrap for vaccine types, plot dose as linetype
ggplot(master_df[master_df$region == 1,]) + 
    geom_line(aes(date, num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) +
    facet_wrap(~vaccine, scales = "free")

# plot all regions using facet wrap
ggplot(master_df[master_df$vaccine == 'Pfizer',]) + 
    geom_line(aes(date, num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) + 
    facet_wrap(~region, scales = "free")

setDT(master_df)
master_df[,cum_doses := cumsum(num), by =.(age.group, region, vaccine, dose)]
master_df[,overall_doses := cumsum(num), by =.(region, dose)]

ggplot(master_df[master_df$vaccine == 'Pfizer',]) + 
    geom_line(aes(date, overall_doses, group = dose, linetype = dose), lwd = 0.5) + 
    facet_wrap(~region, scales = "free")

ggplot(master_df[master_df$vaccine == 'Pfizer',]) + 
    geom_line(aes(date, cum_doses, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) + 
    facet_wrap(~region, scales = "free")




# compare old and new vaccine schedules to look for differences


vax_schedule_old = readRDS('./fitting_data/vax-covidm20211219052551.rds')
vax_schedule_new = readRDS('./fitting_data/vax-covidm20220107085654.rds')

num_to_loop_old = length(vax_schedule_old[[1]]$va1)
num_to_loop_new = length(vax_schedule_new[[1]]$va1)

region_idxs = c(1,3,4,5,6,9,10)

master_df_old = data.table(date = NULL, 
                       region = NULL, 
                       age.group = NULL, 
                       vaccine = NULL,
                       dose = NULL,
                       num = NULL)

for (i in 1:num_to_loop_old){
    print(i/num_to_loop_old)
    day = vax_schedule_old[[1]]$vt[i]
    for (j in 1:16){
        age_sum = rep(0,4)
        for (k in region_idxs){
            num_va1 = vax_schedule_old[[k]]$va1[[i]][j]
            num_va2 = vax_schedule_old[[k]]$va2[[i]][j]
            num_vb1 = vax_schedule_old[[k]]$vb1[[i]][j]
            num_vb2 = vax_schedule_old[[k]]$vb2[[i]][j]
            df = data.table(date = rep(day,4), 
                            region = rep(k,4), 
                            age.group = rep(j,4),
                            vaccine = c('AZ', 'AZ', 'Pfizer', 'Pfizer'),
                            dose = c('First', 'Second', 'First', 'Second'),
                            num = c(num_va1, num_va2, num_vb1, num_vb2))
            master_df_old = rbind(master_df_old, df)
            age_sum[1] = age_sum[1] + num_va1
            age_sum[2] = age_sum[2] + num_va2
            age_sum[3] = age_sum[3] + num_vb1
            age_sum[4] = age_sum[4] + num_vb2
        }
    }
}


master_df_new = data.table(date = NULL, 
                           region = NULL, 
                           age.group = NULL, 
                           vaccine = NULL,
                           dose = NULL,
                           num = NULL)

for (i in 1:num_to_loop_new){
    print(i/num_to_loop_new)
    day = vax_schedule_new[[1]]$vt[i]
    for (j in 1:16){
        age_sum = rep(0,4)
        for (k in region_idxs){
            num_va1 = vax_schedule_new[[k]]$va1[[i]][j]
            num_va2 = vax_schedule_new[[k]]$va2[[i]][j]
            num_vb1 = vax_schedule_new[[k]]$vb1[[i]][j]
            num_vb2 = vax_schedule_new[[k]]$vb2[[i]][j]
            df = data.table(date = rep(day,4), 
                            region = rep(k,4), 
                            age.group = rep(j,4),
                            vaccine = c('AZ', 'AZ', 'Pfizer', 'Pfizer'),
                            dose = c('First', 'Second', 'First', 'Second'),
                            num = c(num_va1, num_va2, num_vb1, num_vb2))
            master_df_new = rbind(master_df_new, df)
            age_sum[1] = age_sum[1] + num_va1
            age_sum[2] = age_sum[2] + num_va2
            age_sum[3] = age_sum[3] + num_vb1
            age_sum[4] = age_sum[4] + num_vb2
        }
    }
}

master_df_old$type = rep('Old', dim(master_df_old)[1])
master_df_new$type = rep('New', dim(master_df_new)[1])

master_df_new[,cum_doses := cumsum(num), by =.(age.group, region, vaccine, dose)]
master_df_new[,overall_doses := cumsum(num), by =.(region, dose)]

master_df_old[,cum_doses := cumsum(num), by =.(age.group, region, vaccine, dose)]
master_df_old[,overall_doses := cumsum(num), by =.(region, dose)]

overall_df = rbind(master_df_old, master_df_new)

ggplot(overall_df[overall_df$region == 1 & dose == "First",]) + 
    geom_line(aes(date, num, colour=factor(age.group), group = interaction(age.group, type), linetype = type), lwd = 0.5) +
    facet_wrap(~vaccine, scales = "free")

ggplot(overall_df[overall_df$region == 1 & dose == "First",]) + 
    +     geom_line(aes(date, num, colour=factor(age.group), group = interaction(age.group, type), linetype = type), lwd = 0.5) +
    +     facet_wrap(~vaccine, scales = "free")


ggplot(overall_df[overall_df$vaccine == 'Pfizer' && overall_df$region == 1,]) + 
    geom_line(aes(date, cum_doses, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) + 
    facet_wrap(~as.factor(overall_df$type), scales = "free")


######

master_df = overall_df

# plot region 1 only, use facet wrap for vaccine types, plot dose as linetype
ggplot(master_df[master_df$region == 1,]) + 
    geom_line(aes(date, num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) +
    facet_wrap(~vaccine, scales = "free")

# plot all regions using facet wrap
ggplot(master_df[master_df$vaccine == 'Pfizer',]) + 
    geom_line(aes(date, num, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) + 
    facet_wrap(~region, scales = "free")

setDT(master_df)
master_df[,cum_doses := cumsum(num), by =.(age.group, region, vaccine, dose)]
master_df[,overall_doses := cumsum(num), by =.(region, dose)]

ggplot(master_df[master_df$vaccine == 'Pfizer',]) + 
    geom_line(aes(date, overall_doses, group = dose, linetype = dose), lwd = 0.5) + 
    facet_wrap(~region, scales = "free")

ggplot(master_df[master_df$vaccine == 'Pfizer',]) + 
    geom_line(aes(date, cum_doses, colour=factor(age.group), group = interaction(age.group, dose), linetype = dose), lwd = 0.5) + 
    facet_wrap(~region, scales = "free")

