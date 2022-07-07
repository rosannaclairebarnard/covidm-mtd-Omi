# calculate Table S1 estimates for relative transmissibility of VOCs

library(data.table)
library(qs)
library(stringr)

# load particle filter fit file
pffit = qread('./fits/pf_relu_yeswane_sev2.0_22050605_20220511163014.qs')

# calculate v2_relu, v3_relu and v4_relu from each fit file

which_pops = c(1,3,4,5,6,9,10)

table = data.table(region=NULL, v2_relu=NULL, v3_relu=NULL, v4_relu=NULL)

for (i in which_pops){

    pf_v2 = pffit$posteriors[[i]]$v2_relu
    pf_v3 = pffit$posteriors[[i]]$v3_relu
    pf_v4 = pffit$posteriors[[i]]$v4_relu
    
    region_name = pffit$parameters[[i]]$pop[[1]]$name
    
    row = data.table(region = region_name, v2_relu = pf_v2, 
                     v3_relu = pf_v3, v4_relu = pf_v4)
    table = rbind(table, row)
}

table$v3_relative_to_WT = table$v2_relu * table$v3_relu
table$v4_relative_to_WT = table$v3_relative_to_WT * table$v4_relu
table
# we want to save the output as a .csv file
datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
fwrite(table, paste0("./output/paper/may22/tableS1_", datetime, ".csv"))

