# Process Delta variant sequencing data @ NHS England regional scale for covidm

library(data.table)
library(lubridate)
library(ggplot2)
library(qs)

# Load files - these are sensitive data and are not provided within the Git repo
folder_path = '~/Documents/uk_covid_data_sensitive/phe/2021-12-03/'
var  = fread('~/Documents/uk_covid_data_sensitive/phe/2021-12-02/Modellers_rapid_travel_20211201.csv')
pos  = fread(paste0(folder_path, 'Anonymised Combined Line List 20211203.csv'))
sgtf = fread('~/Documents/uk_covid_data_sensitive/phe/2021-12-03/SGTF_linelist_20211203.csv')

# Clean dates and variable names
pos[, specimen_date := dmy(specimen_date)]
pos[, Onsetdate := dmy(Onsetdate)]
pos[, lab_report_date := dmy(lab_report_date)]
var[, specimen_date_sk := ymd(specimen_date_sk)]
sgtf[, specimen_date := ymd(specimen_date)]
names(var) = paste0("v_", names(var))
names(sgtf) = paste0("s_", names(sgtf))

# Merge to master data
total = merge(pos, var, by.x = "finalid", by.y = "v_finalid", all = TRUE)
total = merge(total, sgtf, by.x = "finalid", by.y = "s_FINALID", all = TRUE)

# Get a subset of data from October 2020 onwards and plot to check
w = total[v_specimen_date_sk >= "2020-10-01" & specimen_date >= "2020-10-01"]
ggplot(w) +
    geom_bar(aes(x = v_specimen_date_sk, fill = v_variant, 
                 group = interaction(v_specimen_date_sk, v_variant)), 
             position = position_fill(), width = 1) +
    facet_wrap(~NHSER_name)

# Save final subset of (sensitive) data
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
if (0){
    fwrite(w, paste0(folder_path, "merged-", datetime, '.csv'))
}

# Categorise Delta as a combination of VOC-21APR-02 and AY.4.2 sublineage VUI-21OCT-01,
# see https://www.gov.uk/government/news/covid-19-variants-identified-in-the-uk

w$delta = rep(FALSE, dim(w)[1])
w$delta[w$v_variant %in% c("VOC-21APR-02", "VUI-21OCT-01")] = TRUE

# Process data for covidm fitting
delta = w[pillar == "Pillar 2" & 
              v_overall_travel == 'Unknown' & 
              NHSER_name != "", 
          .(delta = sum(delta == TRUE), 
            other = sum(delta == FALSE)), 
            keyby = .(date = v_specimen_date_sk, nhs_name = NHSER_name)]
delta = rbind(delta, 
              delta[!nhs_name %in% c("Northern Ireland", "Scotland", "Wales"),
                    .(delta = sum(delta, na.rm = T), 
                      other = sum(other, na.rm = T), 
                      nhs_name = "England"),
                    by = date], fill = TRUE)

# Plot delta data to check
ggplot(delta) +
    geom_line(aes(x = date, y = delta / (delta + other)), size = 0.25) +
    facet_wrap(~nhs_name) +
    labs(x = NULL, y = "Relative frequency of\nDelta B.1.617.2 VOC") +
    scale_x_date(date_breaks = "2 months", date_labels = "%b")

# Save aggregated delta data for model fitting in Git repo
qsave(delta, paste0("./fitting_data/delta-", datetime, '.qs'))
