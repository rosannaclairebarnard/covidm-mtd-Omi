library(stringr)
library(zoo)
library(data.table)
library(lubridate)
library(ggplot2)

make_xs = function(mobility_file, plot_check = FALSE, date_cutoff = NA)
{
    # Needed components
    cm_structure_UK = readRDS("./covidm_for_fitting/data/structure_UK.rds");
    uk_covid_data_path = "./fitting_data/";
    datapath = function(x) paste0(uk_covid_data_path, x)
    
    # Google visits data
    # We don't upload full Google Mobility Report files to the repo because they are 200+ MB, 
    # but if you download the "global CSV" from https://www.google.com/covid19/mobility/, that's it.
    googmo = fread(mobility_file);
    # between 2021-02-12 and 2021-02-18, Google have added a new column `place_id` (need to remove for this script)
    googmo$place_id <- NULL
    googmo = googmo[country_region_code == "GB" & sub_region_1 != "" & sub_region_2 == "" & metro_area == ""]
    
    # Apply date cutoff if requested
    if (!is.na(date_cutoff)) {
        googmo = googmo[date <= date_cutoff]
    }
    
    # Melt mobility data
    googmo = melt(googmo, id.vars = 1:8, variable.name = "GPL_TYPE", value.name = "value")
    googmo[, GPL_TYPE := str_replace(GPL_TYPE, "_percent_change_from_baseline", "")]
    
    # Match GSS with Google placenames
    process_name = function(x)
    {
        x = str_replace_all(x, "Council", "");
        x = str_replace_all(x, "Principle Area", "");
        x = str_replace_all(x, "Principal Area", "");
        x = str_replace_all(x, "County of", "");
        x = str_replace_all(x, "City of", "");
        x = str_replace_all(x, "City", "");
        x = str_replace_all(x, "County Borough", "");
        x = str_replace_all(x, "Borough of", "");
        x = str_replace_all(x, "Greater", "");
        x = str_replace_all(x, "Islands", "");
        x = str_replace_all(x, ",", "");
        x = str_replace_all(x, "  ", " ");
        x = str_trim(x);
        x
    }
    
    # Create table of all GSS names relevant for us
    wards = fread("~/Dropbox/uk_covid_data/fitting/data/Ward_to_Local_Authority_District_to_County_to_Region_to_Country_December_2019_Lookup_in_United_Kingdom.csv")
    
    placematch = rbind(
        unique(data.table(gss_code = wards$LAD19CD, gss_name = wards$LAD19NM)),
        unique(data.table(gss_code = wards$CTY19CD, gss_name = wards$CTY19NM)),
        unique(data.table(gss_code = wards$RGN19CD, gss_name = wards$RGN19NM))
    )[gss_name != ""];
    
    placematch[, match_name := process_name(gss_name)];
    
    # Create table of Google placenames
    googmo[sub_region_1 %like% "Eileanan", sub_region_1 := "Na h-Eileanan Siar"] # Rename places with spelling differences
    googmo[sub_region_1 %like% "Rhondda",  sub_region_1 := "Rhondda Cynon Taf"]  # in Google mobility data resp. ONS data
    gplacematch = data.table(goog_name = googmo[, unique(sub_region_1)]);
    gplacematch[, match_name := process_name(goog_name)];
    
    # Merge
    match_table = merge(placematch, gplacematch, by = "match_name", all.y = T);
    
    # Check for duplicates
    resolve_duplicates = function(match_table, goog_nm, gss_nm, gss_cd = "")
    {
        if (is.na(gss_nm)) {
            match_table[!(goog_name == goog_nm & (gss_code != gss_cd))]
        } else {
            match_table[!(goog_name == goog_nm & (gss_name != gss_nm))]
        }
    }
    match_table = resolve_duplicates(match_table, "Greater London", "London") # London (Google) refers to Greater London (ONS), not City of London (ONS)
    match_table = resolve_duplicates(match_table, "Greater Manchester", "Greater Manchester") # Greater Manchester (Google) refers to Greater Manchester (ONS), not Manchester (ONS)
    match_table = resolve_duplicates(match_table, "West Midlands", NA, "E11000005") # West Midlands (Google) refers to the metropolitan county, not the metropolitan district
    
    if (nrow(match_table[duplicated(goog_name)]) > 0) {
        stop("Duplicates in match_table: ", paste(match_table[, goog_name[duplicated(goog_name)]], collapse = ", "))
    }
    
    # Amalgamate by NHS England region / DA
    lads = fread("~/Dropbox/uk_covid_data/fitting/data/Lower_Layer_Super_Output_Area_2011_to_Clinical_Commissioning_Group_to_Local_Authority_District_April_2020_Lookup_in_England.csv");
    nhsers = fread("~/Dropbox/uk_covid_data/fitting/data/Clinical_Commissioning_Group_to_STP_and_NHS_England_Region_April_2020_Lookup_in_England.csv");
    
    # Match LAD to CCG to NHS region
    lad_to_ccg = unique(lads[, .(gss_code = LAD20CD, CCG20CD)]);
    ccg_to_nhs = unique(nhsers[, .(CCG20CD, region_name = NHSER20NM)]);
    
    match_table = merge(match_table, lad_to_ccg, by = "gss_code", all.x = T)
    match_table = merge(match_table, ccg_to_nhs, by = "CCG20CD", all.x = T)
    
    # Fill in DAs and some regions
    match_table[gss_code %like% "^N", region_name := "Northern Ireland"]
    match_table[gss_code %like% "^S", region_name := "Scotland"]
    match_table[gss_code %like% "^W", region_name := "Wales"]
    match_table[gss_name == "London", region_name := "London"]
    match_table[gss_name == "Buckinghamshire", region_name := "South East"];
    
    # Fill in missing counties E10 and metropolitan counties E11
    ladcounty = unique(wards[, .(LAD19CD, LAD19NM, CTY19CD, CTY19NM)])
    ladcounty = merge(ladcounty, unique(lads[, .(LAD19CD = LAD20CD, CCG20CD)]), all.x = T)
    ladcounty = merge(ladcounty, unique(nhsers[, .(CCG20CD, region_name = NHSER20NM)]), by = "CCG20CD", all.x = T)
    county_to_nhs = unique(ladcounty[CTY19CD != "", .(gss_code = CTY19CD, region_name)])
    match_table = merge(match_table, unique(ladcounty[, .(gss_code = CTY19CD, region_name)]), by = "gss_code", all.x = T);
    
    # Take region name x or y
    if (nrow(match_table[!is.na(region_name.x) & !is.na(region_name.y)]) > 0) {
        stop("Match table has both region_name.x and region_name.y.");
    }
    region_names = unique(c(match_table[!is.na(region_name.x), region_name.x], match_table[!is.na(region_name.y), region_name.y]));
    either = function(x, y) ifelse(is.na(x) & is.na(y), NULL, ifelse(is.na(x), y, x))
    match_table[, region_name := either(region_name.x, region_name.y)];
    match_table[, region_name.x := NULL];
    match_table[, region_name.y := NULL];
    
    # Remove entries that are duplicates because of two CCGs but which are in the same region
    match_table[, CCG20CD := NULL];
    match_table = unique(match_table);
    
    # Count number of regions each google name matches to
    match_table[, mappings := .N, by = goog_name];
    
    # Get population info
    match_table = merge(match_table, cm_structure_UK[, .(gss_code = Code, population = `All ages`)], by = "gss_code", all.x = T);
    if (nrow(match_table[is.na(population)]) > 0) {
        stop("Missing population information.")
    }
    
    # Bring together and amalgamate by NHS region
    googen2 = merge(googmo, match_table, by.x = "sub_region_1", by.y = "goog_name", allow.cartesian = T);
    
    # Fill missing values in data sets, from median for date within region
    googen2[, value := as.numeric(value)]
    googen2[, value := na.aggregate(value, FUN = median), by = .(date, GPL_TYPE, region_name)]
    
    # Amalgamate by NHS region
    x = googen2[, .(value = weighted.mean(value, population/mappings, na.rm = T)), keyby = .(region_name, GPL_TYPE, date)]
    
    # Add England and United Kingdom
    x = rbind(x,
        googen2[, .(region_name = "United Kingdom", value = weighted.mean(value, population/mappings, na.rm = T)), keyby = .(GPL_TYPE, date)],
        googen2[!region_name %in% c("Northern Ireland", "Scotland", "Wales"), 
            .(region_name = "England", value = weighted.mean(value, population/mappings, na.rm = T)), keyby = .(GPL_TYPE, date)],
        use.names = TRUE
    )
    
    # Transform data and add date
    x[, date := ymd(as.character(date))]
    x[, value := 1 + value/100]
    x[, t := as.numeric(date - ymd("2020-01-01"))]
    
    # Fill in NAs within time series
    x[, value := na.approx(value, na.rm = FALSE), by = .(region_name, GPL_TYPE)]
    
    # Ensmoothen
    xs = x[, .(date, t, value = rollmean(value, 7, fill = NA)), by = .(region_name, GPL_TYPE)]
    xs = xs[!is.na(value)]
    
    if (plot_check) {
        # Plot to check
        theplot = ggplot(xs) +
            geom_line(aes(x = date, y = value, colour = region_name, group = region_name)) +
            geom_hline(aes(yintercept = 1), linetype = "22") +
            facet_wrap(~GPL_TYPE) +
            ylim(0, NA) +
            labs(x = NULL, y = "Mobility index", colour = "Region")
        print(theplot)
    }
    
    if (any(is.na(xs$value))) {
        stop("NAs remain.")
    }
    
    return (xs)
}

make_y = function(xs, schools)
{
    # Reformat
    y = cbind(
        xs[GPL_TYPE == "grocery_and_pharmacy",           .(grocpharm = value), keyby = .(region_name, t)],
        parks = xs[GPL_TYPE == "parks",                  .(parks = value), keyby = .(region_name, t)][, parks],
        residential = xs[GPL_TYPE == "residential",      .(residential = value), keyby = .(region_name, t)][, residential],
        retrec = xs[GPL_TYPE == "retail_and_recreation", .(retrec = value), keyby = .(region_name, t)][, retrec],
        transit = xs[GPL_TYPE == "transit_stations",     .(transit = value), keyby = .(region_name, t)][, transit],
        workplace = xs[GPL_TYPE == "workplaces",         .(workplaces = value), keyby = .(region_name, t)][, workplaces]
    );
    
    # Fill in from Jan 1 2020
    y = rbind(y, 
        data.table(region_name = rep(unique(y$region_name), each = 48), t = rep(0:47, length(unique(y$region_name))),
            grocpharm = 1, parks = 1, residential = 1, retrec = 1, transit = 1, workplace = 1)
    )
    y = y[order(region_name, t)]
    
    # Carry forward last day
    last_day = y[t == max(t)]
    last_fortnight = rbindlist(rep(list(last_day), 14))
    last_fortnight[, t := t + 0:13, by = region_name]
    
    for (added_fortnights in 1:45) {
        y = rbind(y, last_fortnight);
        last_fortnight[, t := t + 14];
    }
    
    # Add school terms: England
    # opened for 4th Jan 2021, closed again on 5th Jan 2021, closed until Monday 8th March 2021
    school_close =  c("2020-2-16", "2020-3-22", "2020-10-25", "2020-12-20", "2021-01-05", "2021-04-01", "2021-05-29", "2021-07-24", "2021-10-23", "2021-12-18", "2022-02-12", "2022-04-09", "2022-05-28", "2022-07-23", "2022-10-24", "2022-12-17")
    school_reopen = c("2020-2-22", "2020-9-01", "2020-10-31", "2021-01-04", "2021-03-08", "2021-04-19", "2021-06-07", "2021-09-01", "2021-11-01", "2022-01-04", "2022-02-21", "2022-04-25", "2022-06-06", "2022-09-01", "2022-10-31", "2023-01-02")
    school_c = as.numeric(ymd(school_close) - ymd("2020-01-01"))
    school_r = as.numeric(ymd(school_reopen) - ymd("2020-01-01")) - 1
    
    days = y[, unique(t)]
    
    n_closures = rowSums(matrix(days, ncol = length(school_c), nrow = length(days), byrow = FALSE) >= 
            matrix(school_c, ncol = length(school_c), nrow = length(days), byrow = TRUE));
    n_reopenings = rowSums(matrix(days, ncol = length(school_r), nrow = length(days), byrow = FALSE) > 
            matrix(school_r, ncol = length(school_r), nrow = length(days), byrow = TRUE));
    y_school = data.table(t = days, school = 1 - (n_closures - n_reopenings))
    y_school[, date := ymd("2020-01-01") + t]
    
    # now overwrite elements in y_school that exist in schools data
    y_school = merge(y_school, schools, by = "date", all.x = TRUE)
    y_school[!is.na(attendance_percentage), school := attendance_percentage / 100]

    # for each weekend throughout the schools schedule, calculate school 
    # contacts as the average across the preceding week (since I think we are 
    # using POLYMOD contact matrices averaged over days of the week, not adjusted for 
    # weekdays, see https://rdrr.io/cran/socialmixr/man/contact_matrix.html
    # and makecontact.R script in covidm repo)
    
    y_school$weekday <- weekdays(y_school$date)
    y_school$weekday <- factor(y_school$weekday, 
                               levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"), 
                               ordered = TRUE)
    y_school$weekdaynum <- as.integer(y_school$weekday)
    y_school$weeknum <- rep(1, dim(y_school)[1])
    for (i in 1:(dim(y_school)[1]-1)){
        if (y_school$weekdaynum[i] == 7){
            y_school$weeknum[(i+1):dim(y_school)[1]] = y_school$weeknum[(i+1):dim(y_school)[1]] + 1
        }
    }
    
    # loop through all calendar weeks
    for (i in 1:max(y_school$weeknum)){
        # if there are only two elements in each week, skip to next iteration
        if (sum(y_school$weeknum == i) < 3){
            next 
        } else { # otherwise, update weekend attendance using average weekday attendance
            weekday_idxs <- which(y_school$weeknum == i & y_school$weekday %in% c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'))
            weekend_idxs <- which(y_school$weeknum == i & y_school$weekday %in% c('Saturday', 'Sunday'))
            this_attendance <- mean(y_school$school[weekday_idxs])
            y_school$school[weekend_idxs] <- this_attendance
        }
    }
    
    # now that we've done the above calculations, we need to remove the extra
    # columns we added to y_school before merging with y
    y_school$weekday <- NULL
    y_school$weekdaynum <- NULL
    y_school$weeknum <- NULL
    y_school$date <- NULL
    y_school$attendance_percentage <- NULL
    
    y = merge(y, y_school, by = "t")
    y[, date := ymd("2020-01-01") + t]
    
    # Adjust for mass testing in schools
    y[(date >= "2021-03-08" & date <= '2021-09-01'), school := school * 0.7]
    
    region_names = y[, sort(unique(region_name))]
    y[, pop := match(region_name, region_names) - 1]
    y = y[order(pop, t)]
    
    return (y)
}

make_schedule = function(y)
{
    # Create schedule from Google data and school closures
    schedule = list()
    schedule[[1]]  = list(parameter = "contact", pops = 0,  mode = "assign", values = list(), times = numeric(0))
    schedule[[2]]  = list(parameter = "contact", pops = 1,  mode = "assign", values = list(), times = numeric(0))
    schedule[[3]]  = list(parameter = "contact", pops = 2,  mode = "assign", values = list(), times = numeric(0))
    schedule[[4]]  = list(parameter = "contact", pops = 3,  mode = "assign", values = list(), times = numeric(0))
    schedule[[5]]  = list(parameter = "contact", pops = 4,  mode = "assign", values = list(), times = numeric(0))
    schedule[[6]]  = list(parameter = "contact", pops = 5,  mode = "assign", values = list(), times = numeric(0))
    schedule[[7]]  = list(parameter = "contact", pops = 6,  mode = "assign", values = list(), times = numeric(0))
    schedule[[8]]  = list(parameter = "contact", pops = 7,  mode = "assign", values = list(), times = numeric(0))
    schedule[[9]]  = list(parameter = "contact", pops = 8,  mode = "assign", values = list(), times = numeric(0))
    schedule[[10]] = list(parameter = "contact", pops = 9,  mode = "assign", values = list(), times = numeric(0))
    schedule[[11]] = list(parameter = "contact", pops = 10, mode = "assign", values = list(), times = numeric(0))
    schedule[[12]] = list(parameter = "contact", pops = 11, mode = "assign", values = list(), times = numeric(0))

    for (r in 1:nrow(y))
    {
        schedule[[y[r, pop]+1]]$values = c(schedule[[y[r, pop]+1]]$values,
            list(y[r, c(residential, workplace, grocpharm, retrec, transit, school)]));
        schedule[[y[r, pop]+1]]$times = c(schedule[[y[r, pop]+1]]$times, y[r, t])
    }

    return (schedule)
}

plot_y_fancy = function(roadmap_list, proj_start_date, xmark_label = "S4", 
                        xmark_date = "2021-07-19", start_date = "2020-01-01", 
                        end_date = "2021-10-01", dbreaks = "2 months", 
                        dlabels = "%b '%y", ybreaks = seq(0, 1.4, by = 0.2), 
                        ylim = c(0, 1.4), colours_list = NULL)
{
    rm0 = rbindlist(roadmap_list, idcol = "scenario")
    rm0 = rm0[region_name == "England"]
    rm0 = melt(rm0, id.vars = c("t", "region_name", "date", "pop", "scenario"))
    rm0 = rm0[variable %in% c("grocpharm", "retrec", "transit", "workplace")]

    rm0[variable == "grocpharm", DisplayName := "Grocery and pharmacy"]
    rm0[variable == "retrec", DisplayName := "Retail and recreation"]
    rm0[variable == "transit", DisplayName := "Transit stations"]
    rm0[variable == "workplace", DisplayName := "Workplaces"]
    rm0[, scenario := factor(scenario, names(roadmap_list))]

    #rm0[date >= GM_CUTOFF, value := zoo::rollmean(value, 7, fill = "extend"), by = .(region_name, scenario, variable)]

    marks = fread(
    "label,date
    L1,2020-03-26
    L2,2020-11-05
    L3,2021-01-05
    S1,2021-03-08
    S2,2021-04-12
    S3,2021-05-17")
    marks[, date := ymd(date)]
    marks = rbind(marks, data.table(label = xmark_label, date = ymd(xmark_date)))
    marks[, y := 1.3]
    
    colgrey = "#606060"
    if (is.null(colours_list)==TRUE){
        col3 = rev(c("#33c5ff", "#8a03fa", "#8f0c0c"))
    } else {
        col3 = colours_list
    }
    
    names(col3) = names(roadmap_list)[1:length(roadmap_list)]

    ggplot(rm0) +
        geom_hline(aes(yintercept = 1), colour = "#d0d0d0") +
        geom_line(data = rm0[date >= start_date & date <= proj_start_date & scenario == names(roadmap_list)[1]], aes(x = date, y = value), colour = colgrey) +
        geom_line(data = rm0[date >= proj_start_date & date <= end_date], aes(x = date, y = value, colour = scenario)) +
        geom_vline(data = marks[date >= start_date], aes(xintercept = date), linetype = "22", size = 0.25) +
        geom_label(data = marks[date >= start_date], aes(x = date, y = y, label = label), size = 2.5, label.padding = unit(0.15, "lines")) +
        facet_wrap(~DisplayName, ncol = 1) +
        scale_colour_manual(values = col3, aesthetics = c("colour", "fill")) +
        scale_x_date(date_breaks = dbreaks, date_labels = dlabels, expand = c(0, 0)) +
        scale_y_continuous(breaks = ybreaks, limits = ylim) +
        labs(x = "Date", y = "Mobility relative to baseline", colour = "Scenario") +
        cowplot::theme_cowplot(font_size = 9) +
        theme(panel.background = element_rect(fill = "#f4f4f4"), strip.background = element_blank(),
            panel.grid = element_line(colour = "#ffffff"))
}

plot_y_england = function(y, 
    start_date = "2020-01-01", end_date = "2021-12-01", 
    dbreaks = "2 months", dlabels = "%b '%y", 
    ybreaks = seq(0, 1.4, by = 0.2), ylim = c(0, 1.4))
{
    rm0 = copy(y)
    rm0 = rm0[region_name == "England"]
    rm0 = melt(rm0, id.vars = c("t", "region_name", "date", "pop"))
    rm0 = rm0[variable %in% c("grocpharm", "retrec", "transit", "workplace")]

    rm0[variable == "grocpharm", DisplayName := "Grocery and pharmacy"]
    rm0[variable == "retrec", DisplayName := "Retail and recreation"]
    rm0[variable == "transit", DisplayName := "Transit stations"]
    rm0[variable == "workplace", DisplayName := "Workplaces"]

    marks = fread(
    "label,date
    L1,2020-03-26
    L2,2020-11-05
    L3/S0,2021-01-05
    S1,2021-03-08
    S2,2021-04-12
    S3,2021-05-17
    S4,2021-07-19")
    marks[, date := ymd(date)]
    marks[, y := 1.3]
    
    colgrey = "#606060"

    ggplot(rm0) +
        geom_hline(aes(yintercept = 1), colour = "#d0d0d0") +
        geom_line(data = rm0[date >= start_date & date <= end_date], aes(x = date, y = value, colour = DisplayName)) +
        geom_vline(data = marks[date >= start_date], aes(xintercept = date), linetype = "22", size = 0.25) +
        geom_label(data = marks[date >= start_date], aes(x = date, y = y, label = label), size = 2.5, label.padding = unit(0.15, "lines")) +
        scale_x_date(date_breaks = dbreaks, date_labels = dlabels, expand = c(0, 0)) +
        scale_y_continuous(breaks = ybreaks, limits = ylim) +
        labs(x = NULL, y = "Mobility relative\nto baseline", colour = "Mobility index") +
        cowplot::theme_cowplot(font_size = 9) +
        theme(panel.background = element_rect(fill = "#f4f4f4"), strip.background = element_blank(),
            panel.grid = element_line(colour = "#ffffff"), legend.position = "bottom")
}


plot_y_england_schools = function(y, 
                          start_date = "2020-01-01", end_date = "2021-12-01", 
                          dbreaks = "2 months", dlabels = "%b '%y", 
                          ybreaks = seq(0, 1.4, by = 0.2), ylim = c(0, 1.4))
{
    rm0 = copy(y)
    rm0 = rm0[region_name == "England"]
    rm0 = melt(rm0, id.vars = c("t", "region_name", "date", "pop"))
    rm0 = rm0[variable %in% c("school")]
    
    rm0[variable == "school", DisplayName := "Schools"]
    
    marks = fread(
        "label,date
    L1,2020-03-26
    L2,2020-11-05
    L3/S0,2021-01-05
    S1,2021-03-08
    S2,2021-04-12
    S3,2021-05-17
    S4,2021-07-19")
    marks[, date := ymd(date)]
    marks[, y := 1.3]
    
    colgrey = "#606060"
    
    ggplot(rm0) +
        geom_line(data = rm0[date >= start_date & date <= end_date], aes(x = date, y = value, colour = DisplayName)) +
        geom_vline(data = marks[date >= start_date], aes(xintercept = date), linetype = "22", size = 0.25) +
        geom_label(data = marks[date >= start_date], aes(x = date, y = y, label = label), size = 2.5, label.padding = unit(0.15, "lines")) +
        scale_x_date(date_breaks = dbreaks, date_labels = dlabels, expand = c(0, 0)) +
        scale_y_continuous(breaks = ybreaks, limits = ylim) +
        labs(x = NULL, y = "School attendance", colour = "Mobility index") +
        cowplot::theme_cowplot(font_size = 9) +
        theme(panel.background = element_rect(fill = "#f4f4f4"), strip.background = element_blank(),
              panel.grid = element_line(colour = "#ffffff"), legend.position = "bottom")
}


# Alternative schedules for projections
# Version based on absolute rather than relative changes
steps_y = function(xs, grtw, hold_from, ramp_start, ramp_end, schools, SMOOTH_ADJ)
{
    # Reformat
    y = cbind(
        xs[GPL_TYPE == "grocery_and_pharmacy",           .(grocpharm = value), keyby = .(region_name, t)],
        parks = xs[GPL_TYPE == "parks",                  .(parks = value), keyby = .(region_name, t)][, parks],
        residential = xs[GPL_TYPE == "residential",      .(residential = value), keyby = .(region_name, t)][, residential],
        retrec = xs[GPL_TYPE == "retail_and_recreation", .(retrec = value), keyby = .(region_name, t)][, retrec],
        transit = xs[GPL_TYPE == "transit_stations",     .(transit = value), keyby = .(region_name, t)][, transit],
        workplace = xs[GPL_TYPE == "workplaces",         .(workplaces = value), keyby = .(region_name, t)][, workplaces]
    );
    
    # Convert grtw into diffs
    grtw_ref = y[region_name == "England" & t == max(t), c(grocpharm, retrec, transit, workplace)]
    grtw_mat = diff(matrix(c(grtw_ref, unlist(grtw)), byrow = TRUE, ncol = 4))
    grtw = unname(split(grtw_mat, row(grtw_mat)))

    # Fill in from Jan 1 2020
    y = rbind(y, 
        data.table(region_name = rep(unique(y$region_name), each = 48), t = rep(0:47, length(unique(y$region_name))),
            grocpharm = 1, parks = 1, residential = 1, retrec = 1, transit = 1, workplace = 1)
    )
    y = y[order(region_name, t)]
    
    # Carry forward last day
    last_day = y[t == max(t)]
    last_fortnight = rbindlist(rep(list(last_day), 14))
    last_fortnight[, t := t + 0:13, by = region_name]
    
    for (added_fortnights in 1:45) {
        y = rbind(y, last_fortnight);
        last_fortnight[, t := t + 14];
    }
    
    y[, date := ymd("2020-01-01") + t]
    yy = split(y, by = "region_name")

    do_ramp = function(yr, date_start, date_end, grtw)
    {
        yr[date %between% ymd(c(date_start, date_end)), grocpharm := grocpharm[1] + seq(0, grtw[1], length.out = .N)]
        yr[date %between% ymd(c(date_start, date_end)), retrec    := retrec[1]    + seq(0, grtw[2], length.out = .N)]
        yr[date %between% ymd(c(date_start, date_end)), transit   := transit[1]   + seq(0, grtw[3], length.out = .N)]
        yr[date %between% ymd(c(date_start, date_end)), workplace := workplace[1] + seq(0, grtw[4], length.out = .N)]
        return (yr)
    }

    for (reg in names(yy)) {
        # Holding period
        yy[[reg]] = do_ramp(yy[[reg]], hold_from, ramp_start[1], c(0, 0, 0, 0));

        # Each ramp period
        for (i in seq_along(ramp_start)) {
            # Ramp part
            yy[[reg]] = do_ramp(yy[[reg]], ramp_start[i], ramp_end[i], grtw[[i]]);
            final_date = if (i == length(ramp_start)) as.character(max(yy[[reg]]$date)) else ramp_start[i + 1]
            yy[[reg]] = do_ramp(yy[[reg]], ramp_end[i], final_date, c(0, 0, 0, 0));
        }
    }
    
    y = rbindlist(yy)
    
    # Add school terms: England
    school_close =  c("2020-2-16", "2020-3-22", "2020-10-25", "2020-12-20", "2021-01-05", "2021-04-01", "2021-05-29", "2021-07-24", "2021-10-23", "2021-12-18", "2022-02-12", "2022-04-09", "2022-05-28", "2022-07-23", "2022-10-24", "2022-12-17")
    school_reopen = c("2020-2-22", "2020-9-01", "2020-10-31", "2021-01-04", "2021-03-08", "2021-04-19", "2021-06-07", "2021-09-01", "2021-11-01", "2022-01-04", "2022-02-21", "2022-04-25", "2022-06-06", "2022-09-01", "2022-10-31", "2023-01-02")
    school_c = as.numeric(ymd(school_close) - ymd("2020-01-01"))
    school_r = as.numeric(ymd(school_reopen) - ymd("2020-01-01")) - 1

    days = y[, unique(t)]

    n_closures = rowSums(matrix(days, ncol = length(school_c), nrow = length(days), byrow = FALSE) >= 
            matrix(school_c, ncol = length(school_c), nrow = length(days), byrow = TRUE));
    n_reopenings = rowSums(matrix(days, ncol = length(school_r), nrow = length(days), byrow = FALSE) > 
            matrix(school_r, ncol = length(school_r), nrow = length(days), byrow = TRUE));
    y_school = data.table(t = days, school = 1 - (n_closures - n_reopenings))
    y = merge(y, y_school, by = "t")
    
    region_names = y[, sort(unique(region_name))]
    y[, pop := match(region_name, region_names) - 1]
    y = y[order(pop, t)]
    
    # Adjust for mass testing in schools
    y[(date >= "2021-03-08" & date <= '2021-09-01'), school := school * 0.7]
    
    # Smooth
    y[date >= hold_from-SMOOTH_ADJ, grocpharm   := zoo::rollmean(grocpharm  , 7, fill = "extend"), by = .(region_name)]
    y[date >= hold_from-SMOOTH_ADJ, parks       := zoo::rollmean(parks      , 7, fill = "extend"), by = .(region_name)]
    y[date >= hold_from-SMOOTH_ADJ, residential := zoo::rollmean(residential, 7, fill = "extend"), by = .(region_name)]
    y[date >= hold_from-SMOOTH_ADJ, retrec      := zoo::rollmean(retrec     , 7, fill = "extend"), by = .(region_name)]
    y[date >= hold_from-SMOOTH_ADJ, transit     := zoo::rollmean(transit    , 7, fill = "extend"), by = .(region_name)]
    y[date >= hold_from-SMOOTH_ADJ, workplace   := zoo::rollmean(workplace  , 7, fill = "extend"), by = .(region_name)]


    return (y)
}






####
# CODE TO ASSEMBLE MAIN SCHEDULE FOR FITTING
####

if (0)
{
    # NOTE 4 June 2021 -- cutting off on Wednesday 26th May because of unusual patterns due to May bank holiday weekend.
    # xs = make_xs("~/Dropbox/uk_covid_data/fitting/data/Global_Mobility_Report-2021-06-03.csv", TRUE, "2021-05-26")
    xs = make_xs("~/Dropbox/uk_covid_data/fitting/data/Global_Mobility_Report-2022-05-04.csv", TRUE)
    schools_data = fread('./fitting_data/schoolattendancedata-England-20220506120901.csv')
    schools_data[, date := as.Date(date)]
    y = make_y(xs, schools_data);
    schedule = make_schedule(y);
    datetime <- str_replace_all(Sys.time(), "[ :GMTBST-]", "")
    saveRDS(schedule, paste0("./fitting_data/schedule3-MTPs-", datetime, ".rds"))
    
    mob_plot_3a = plot_y_england(y, end_date = '2022-12-31')
    mob_plot_3a
    mob_plot_schools = plot_y_england_schools(y, end_date = '2022-12-31')
    mob_plot_schools
    # near = function(date, ds) date >= ymd(ds) - 3 & date <= ymd(ds) + 3
    # y[region_name == "England" & near(date, "2021-03-18")]
}



####
# SCRATCH BELOW
####

if (0) {

# use the following to plot school attendance by day of week
y$weekday <- weekdays(y$date)
y$weekday <- factor(y$weekday, 
                           levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"), 
                           ordered = TRUE)
# plot by day of week
ggplot(y, aes(date, school, colour=weekday)) + 
    geom_point() +
    labs(x = "Date", y = "Attendance (%)", colour = "Weekday") +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y")

y$weekend <- rep(NA, dim(y)[1])
for (i in 1:dim(y)[1]){
    if (y$weekday[i] %in% c('Saturday', 'Sunday')){
        y$weekend[i] <- 'Weekend'
    } else {
        y$weekend[i] <- 'Weekday'
    }
}

# plot weekday vs. weekend data
ggplot(y, aes(date, school, colour=weekend)) + 
    geom_point() +
    labs(x = "Date", y = "Attendance (%)", colour = "Weekday") +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y")

}

