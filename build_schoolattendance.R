library(data.table)
library(ggplot2)
library(lubridate)
library(zoo)
library(cowplot)
library(stringr)

# load table 1b from attendance data available at: 
# https://explore-education-statistics.service.gov.uk/find-statistics/attendance-in-education-and-early-years-settings-during-the-coronavirus-covid-19-outbreak
# table 1b corresponds to attendance data at England level
dt = fread('./data/attendance-in-education-and-early-years-settings-during-the-coronavirus-covid-19-outbreak_2022-week-18/data/table_1b_daily_attendance_in_state_schools_during_covid_19_.csv')

# reformat dates
dt$date = as.Date(dt$date)

# combine unadjusted with adjusted estimates: adjusted estimates correspond to 
# estimates which are adjusted for the absence of students in years 11-13 (Summer term 2021)

# get unadjusted attendance data (not adjusted for Y11-13 absences in June/July 2021)
uadj <- data.table(date = as.Date(dt$date), phase = dt$phase, 
                   value = as.numeric(dt$proportion_students_attending, na.rm = T), 
                   datatype = rep('Unadjusted', dim(dt)[1]),
                   number_of_schools = (as.numeric(dt$number_open_including_inset, na.rm = T) / as.numeric(dt$proportion_open_including_inset, na.rm = T)) * 100,
                   number_of_students = (as.numeric(dt$number_students_attending, na.rm = T) / as.numeric(dt$proportion_students_attending, na.rm = T)) * 100)
# plot to check unadjusted attendance data
p_uadj = ggplot(data = uadj, aes(x = as.Date(date), y = as.numeric(value, na.rm = TRUE)), group = phase, colour = phase) +
    geom_rect(xmin=as.Date('2020-10-25'), xmax=as.Date('2020-10-31'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2020-12-20'), xmax=as.Date('2021-01-04'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-01-05'), xmax=as.Date('2021-03-08'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='slategray2') +
    geom_rect(xmin=as.Date('2021-04-01'), xmax=as.Date('2021-04-19'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-05-29'), xmax=as.Date('2021-06-07'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-07-24'), xmax=as.Date('2021-09-01'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_point(aes(color = phase), size = 1.5) +
    labs(x = 'Date', y = 'Unadjusted proportion of students attending (England)', color = 'Institution type')
p_uadj    
# get adjusted attendance data (adjusted for Y11-13 absences in June/July 2021)
adj <- data.table(date  = as.Date(dt$date), phase = dt$phase, 
                  value = as.numeric(dt$overall_attendance_rate_excl_year_11_13_students_not_in_attendance_adjusted, na.rm = T), 
                  datatype = rep('Adjusted for Y11-13 absence', dim(dt)[1]),
                  number_of_schools = (as.numeric(dt$number_open_including_inset) / as.numeric(dt$proportion_open_including_inset)) * 100,
                  number_of_students = (as.numeric(dt$number_students_attending) / as.numeric(dt$overall_attendance_rate_excl_year_11_13_students_not_in_attendance_adjusted)) * 100)
# plot to check adjusted attendance data
p_adj = ggplot(data = adj, aes(x = as.Date(date), y = as.numeric(value, na.rm = TRUE)), group = phase, colour = phase) +
    geom_rect(xmin=as.Date('2020-10-25'), xmax=as.Date('2020-10-31'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2020-12-20'), xmax=as.Date('2021-01-04'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-01-05'), xmax=as.Date('2021-03-08'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='slategray2') +
    geom_rect(xmin=as.Date('2021-04-01'), xmax=as.Date('2021-04-19'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-05-29'), xmax=as.Date('2021-06-07'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-07-24'), xmax=as.Date('2021-09-01'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_point(aes(color = phase), size = 1.5) +
    labs(x = 'Date', y = 'Adjusted proportion of students attending (England)', color = 'Institution type')
p_adj
# look at unadjusted and adjusted data together
plot_grid(p_uadj, p_adj, nrow = 2, align = 'v')

# we use unadjusted data from *outside* of 7th June 2021 to 16th July 2021 and
# prior to the 28th of April 2022, and adjusted data from those dates inclusive

# remove unadjusted data from 7th June 2021 to 16th July 2021 inclusive
uadj <- uadj[!(uadj$date %in% seq(as.Date('2021-06-07'), as.Date('2021-07-16'), by=1)),]
# remove unadjusted data from 28th April 2022 onwards
uadj <- uadj[!(uadj$date >= as.Date('2022-04-28')),]

# remove adjusted data from before 7th June 2021
adj <- adj[!(adj$date < as.Date('2021-06-07')),]
# remove adjusted data from after 16th July 2021 and before 28th April 2022
adj <- adj[!(adj$date %in% seq(as.Date('2021-07-17'), as.Date('2022-04-27'), by=1)),]

# combine the two datasets together
combo <- rbind(uadj, adj)

# compute the average values across institution type

# calculate total number of students across type of institution by date and by datatype
combo[, total_students_across_phase := sum(number_of_students, na.rm = TRUE), by = .(datatype, date)]
# calculate weighted average of attendance rates using numbers of students as calculated above (across type of institution)
combo[, weight := as.numeric(number_of_students, na.rm = TRUE) / as.numeric(total_students_across_phase, na.rm = TRUE)]
# remove all entries with NA weight
combo <- combo[!(is.na(combo$weight) == T),]
# calculate weighted mean across institution types
combo[, weighted_mean := weighted.mean(as.numeric(value, na.rm = TRUE), as.numeric(weight, na.rm = TRUE), na.rm = TRUE), by = .(datatype, date)]

# add average value across institutions
subset_to_keep <- combo[combo$phase == 'All state-funded schools',]
subset_to_keep$phase <- rep('Weighted mean across institutions', dim(subset_to_keep)[1])
subset_to_keep$value <- subset_to_keep$weighted_mean
subset_to_keep$datatype <- rep('Weighted mean', dim(subset_to_keep)[1])
another_combo <- rbind(combo, subset_to_keep)
# remove averaged values which are from before 8th March 2021
another_combo <- another_combo[!(another_combo$datatype == 'Weighted mean' & another_combo$date < as.Date('2021-03-08')),]

# remove unadjusted estimates from 7th June 2021 to 16th July 2021 inclusive
another_combo <- another_combo[!(another_combo$datatype == 'Unadjusted' & another_combo$date %in% seq(as.Date('2021-06-07'), as.Date('2021-07-16'), by=1)),]
# remove weighted mean column
another_combo$weighted_mean <- NULL
another_combo$datatype <- factor(another_combo$datatype, levels = c('Unadjusted', 'Adjusted for Y11-13 absence', 'Weighted mean'))

# plotting

# data points only (full time series)
ggplot(data = another_combo, aes(x = as.Date(date), y = as.numeric(value))) +
    geom_rect(xmin=as.Date('2020-10-25'), xmax=as.Date('2020-10-31'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2020-12-20'), xmax=as.Date('2021-01-04'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-01-05'), xmax=as.Date('2021-03-08'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='slategray2') +
    geom_rect(xmin=as.Date('2021-04-01'), xmax=as.Date('2021-04-19'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-05-29'), xmax=as.Date('2021-06-07'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-07-24'), xmax=as.Date('2021-09-01'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    scale_x_date(date_breaks = '2 months', date_labels = "%b %Y") +
    geom_point(aes(color = phase, shape = datatype), size = 1.5) +
    scale_shape_manual(values = c('Unadjusted' = 2, 'Adjusted for Y11-13 absence' = 6, 'Weighted mean' = 16)) +
    labs(x = 'Date', y = 'Proportion of students attending (England)', color = 'Institution type', shape = 'Data type')

# data to use for England-level school attendance

# pre 8th March 2021, use rolling average of All state-funded schools attendance
pre__8March <- another_combo[another_combo$date <  as.Date('2021-03-08'),]
first       <- data.table(date = pre__8March$date, attendance = pre__8March$value)
# post 8th March 2021, use weighted mean rolling average over all types of institution
post_8March <- another_combo[(another_combo$date >= as.Date('2021-03-08') & another_combo$datatype == 'Weighted mean'),]
second      <- data.table(date = post_8March$date, attendance = post_8March$value)

data_to_use <- rbind(first, second)

# plot to check
ggplot(data = data_to_use, aes(x = as.Date(date))) +
    geom_rect(xmin=as.Date('2020-10-25'), xmax=as.Date('2020-10-31'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2020-12-20'), xmax=as.Date('2021-01-04'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-01-05'), xmax=as.Date('2021-03-08'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='slategray2') +
    geom_rect(xmin=as.Date('2021-04-01'), xmax=as.Date('2021-04-19'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-05-29'), xmax=as.Date('2021-06-07'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-07-24'), xmax=as.Date('2021-09-01'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    scale_x_date(date_breaks = '2 months', date_labels = "%b %Y") +
    geom_point(aes(y = as.numeric(attendance, na.rm = TRUE)), size = 1.5) +
    labs(x = 'Date', y = 'Proportion of students attending (England)') +
    scale_y_continuous(breaks = c(5, 15, 25, 35, 45, 55, 65, 75, 85, 95), 
                       labels = c('5%', '15%', '25%', '35%', '45%', '55%', '65%', '75%', '85%', '95%'),
                       limits = c(5, 95))

# table 3 - Daily attendance in education settings during the COVID-19 outbreak Pre17 July 2020
t3 = fread('./data/attendance-in-education-and-early-years-settings-during-the-coronavirus-covid-19-outbreak_2022-week-18/data/table_3_daily_attendance_in_education_settings_during_covid19_pre17july2020.csv') 
t3$date = as.Date(t3$date, format = '%d/%m/%Y')

ggplot(data = t3, aes(x = date, y = as.numeric(proportion_children_attending, na.rm = T))) +
    geom_rect(xmin=as.Date('2020-10-25'), xmax=as.Date('2020-10-31'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2020-12-20'), xmax=as.Date('2021-01-04'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-01-05'), xmax=as.Date('2021-03-08'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='slategray2') +
    geom_rect(xmin=as.Date('2021-04-01'), xmax=as.Date('2021-04-19'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-05-29'), xmax=as.Date('2021-06-07'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-07-24'), xmax=as.Date('2021-09-01'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    scale_x_date(date_breaks = '2 months', date_labels = "%b %Y") +
    geom_point() +
    labs(x = 'Date', y = 'Proportion of students attending (England)')

# combine attendance from table 3 (Spring 2020) with attendance from table 1b (Autumn 2020 ->)
third = data.table(date = t3$date, attendance = t3$proportion_children_attending)
final_attendance = rbind(data_to_use, third)
# plot full time series to check
ggplot(data = final_attendance, aes(x = as.Date(date))) +
    geom_rect(xmin=as.Date('2020-03-23'), xmax=as.Date('2020-07-17'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='slategray2') +
    geom_rect(xmin=as.Date('2020-07-18'), xmax=as.Date('2020-09-01'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2020-10-25'), xmax=as.Date('2020-10-31'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2020-12-20'), xmax=as.Date('2021-01-04'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-01-05'), xmax=as.Date('2021-03-08'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='slategray2') +
    geom_rect(xmin=as.Date('2021-04-01'), xmax=as.Date('2021-04-19'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-05-29'), xmax=as.Date('2021-06-07'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    geom_rect(xmin=as.Date('2021-07-24'), xmax=as.Date('2021-09-01'), ymin=-Inf, ymax=Inf, alpha=0.5, fill='grey80') +
    #geom_line(aes(y = as.numeric(rolling_average, na.rm = TRUE), color = phase, linetype = datatype)) +
    #scale_linetype_manual(values=c("longdash", "twodash", 'solid')) +
    #scale_size_manual(values = c(0.6, 0.6, 1)) +
    scale_x_date(date_breaks = '2 months', date_labels = "%b %Y", limits = c(as.Date(min(final_attendance$date)), as.Date(max(final_attendance$date)))) +
    geom_point(aes(y = as.numeric(final_attendance$attendance, na.rm = TRUE)), size = 1.5) +
    labs(x = 'Date', y = 'Proportion of students attending (England)') +
    scale_y_continuous(breaks = c(5, 15, 25, 35, 45, 55, 65, 75, 85, 95), 
                       labels = c('5%', '15%', '25%', '35%', '45%', '55%', '65%', '75%', '85%', '95%'),
                       limits = c(0, 95))
final_dt <- data.table(date = final_attendance$date, 
                       attendance_percentage = final_attendance$attendance)
final_dt = final_dt[order(as.Date(final_dt$date)),]

datetime <- str_replace_all(Sys.time(), "[ :GMTBST-]", "")
write.csv(final_dt, file = paste0("./fitting_data/schoolattendancedata-England-", datetime, ".csv"), 
          row.names = FALSE)
