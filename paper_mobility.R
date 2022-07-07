source("build_schedule_google_mobility.R")

# Load google mobility data, calculate rolling 7-day averages and plot to check
xs = make_xs("~/Desktop/Global_Mobility_Report-2022-05-04.csv", TRUE, NA)

# Load school attendance data for England
schools_data = fread('./fitting_data/schoolattendancedata-England-20220506120901.csv')
schools_data[, date := as.Date(date)]

# Process mobility and school attendance data, generating future schedule
y_schools = make_y(xs, schools_data);

# Check output
mob_plot_3a = plot_y_england(y_schools, end_date = '2022-12-31')
mob_plot_3a
mob_plot_schools = plot_y_england_schools(y_schools, end_date = '2022-09-30')
mob_plot_schools

####
# BUILD SCENARIOS
####

# Mobility scenarios for paper (resubmission May 2022)

HOLD_FROM = max(xs$date) 

paper_aw_3wk = steps_y(
    xs,
    list(c(1.0, 1.0, 1.0, 1.0)), 
    hold_from = HOLD_FROM, 
    ramp_start= c(HOLD_FROM), 
    ramp_end =  c(HOLD_FROM+21),
    SMOOTH_ADJ = 7)

paper_aw_3mo = steps_y(
    xs,
    list(c(1.0, 1.0, 1.0, 1.0)), 
    hold_from = HOLD_FROM, 
    ramp_start= c(HOLD_FROM), 
    ramp_end =  c(HOLD_FROM+91), # HOLD_FROM + 91 falls on the same date 3 months later
    SMOOTH_ADJ = 7) 

paper_aw_6mo = steps_y(
    xs,
    list(c(1.0, 1.0, 1.0, 1.0)), 
    hold_from = HOLD_FROM, 
    ramp_start= c(HOLD_FROM), 
    ramp_end =  c(HOLD_FROM+183), # HOLD_FROM + 183 falls on the same day 6 months later
    SMOOTH_ADJ = 7) 

paper_aw_flt = steps_y(
    xs,
    list(c(1.0, 1.0, 1.0, 1.0)), 
    hold_from = HOLD_FROM, 
    ramp_start= c(as.Date('2023-02-01')), 
    ramp_end =  c(as.Date('2023-02-08')), 
    SMOOTH_ADJ = 7) 

# we don't need to generate an additional flat schedule - that file already 
# exists and was used for the model fitting (schedule3-MTPs-20220506121302.rds)

# Bring in correct schools data
paper_aw_3wk[, school := NULL]
paper_aw_3mo[, school := NULL]
paper_aw_6mo[, school := NULL]
paper_aw_flt[, school := NULL]

paper_aw_3wk = merge(paper_aw_3wk, y_schools[, .(date, region_name, school)], by = c("date", "region_name"), all.x = TRUE)
paper_aw_3mo = merge(paper_aw_3mo, y_schools[, .(date, region_name, school)], by = c("date", "region_name"), all.x = TRUE)
paper_aw_6mo = merge(paper_aw_6mo, y_schools[, .(date, region_name, school)], by = c("date", "region_name"), all.x = TRUE)
paper_aw_flt = merge(paper_aw_flt, y_schools[, .(date, region_name, school)], by = c("date", "region_name"), all.x = TRUE)

setcolorder(paper_aw_3wk, names(y_schools))
setcolorder(paper_aw_3mo, names(y_schools))
setcolorder(paper_aw_6mo, names(y_schools))
setcolorder(paper_aw_flt, names(y_schools))

setorder(paper_aw_3wk, pop, date)
setorder(paper_aw_3mo, pop, date)
setorder(paper_aw_6mo, pop, date)
setorder(paper_aw_flt, pop, date)

# Save all schedules
datetime <- str_replace_all(Sys.time(), "[- :BSTGMT]", "")
saveRDS(make_schedule(paper_aw_3wk), paste0("./fitting_data/paper_aw_3wk_", datetime, ".rds"))
saveRDS(make_schedule(paper_aw_3mo), paste0("./fitting_data/paper_aw_3mo_", datetime, ".rds"))
saveRDS(make_schedule(paper_aw_6mo), paste0("./fitting_data/paper_aw_6mo_", datetime, ".rds"))
saveRDS(make_schedule(paper_aw_flt), paste0("./fitting_data/paper_aw_flt_", datetime, ".rds"))

# # get flat mobility schedule
# paper_aw_flt = y_schools

# Plot all schedules
paw = plot_y_fancy(list(
    `3 weeks`  = paper_aw_3wk,
    `3 months` = paper_aw_3mo,
    `6 months` = paper_aw_6mo,
    `flat`     = paper_aw_flt
),  proj_start_date = max(xs$date), 
xmark_label = c("S4", "PBA", "PBE"), xmark_date = c('2021-07-19', '2021-12-08', '2022-01-27'),
start_date = "2021-09-01", end_date = "2022-12-31", 
dbreaks = "1 month", dlabels = "%b %Y", ybreaks = seq(0.2, 1.2, by = 0.2), ylim = c(0.2, 1.2),
colours_list = c('#1b9e77','#d95f02','#7570b3','#e7298a')) + 
    labs(title = "Mobility scenarios")
paw

paw_full = plot_y_fancy(list(
    `3 weeks`  = paper_aw_3wk,
    `3 months` = paper_aw_3mo,
    `6 months` = paper_aw_6mo,
    `None`     = paper_aw_flt
),  proj_start_date = max(xs$date), 
xmark_label = c("S4", "PBA", "PBE"), xmark_date = c('2021-07-19', '2021-12-08', '2022-01-27'),
start_date = "2020-03-01", end_date = "2022-09-30", 
dbreaks = "2 months", dlabels = "%m/%y", ybreaks = seq(0.2, 1.4, by = 0.2), ylim = c(0.2, 1.4),
colours_list = c('#1b9e77','#d95f02','#7570b3','#e7298a')) + 
    labs(title = "Mobility", colour = "Return to\nbaseline", x = NULL)
paw_full

# Save figures
datetime <- str_replace_all(Sys.time(), "[- :BSTGMT]", "")
ggsave(paste0("./output/paperfigs/may22", "/mobility_scenarios_", datetime, ".png"), paw, width = 14, height = 18, units = "cm")
ggsave(paste0("./output/paperfigs/may22", "/mobility_scenarios_full_", datetime, ".png"), paw_full, width = 18, height = 18, units = "cm")
