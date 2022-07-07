library(ggplot2)
theme_set(theme_classic())

scenario_sheet = fread(
    "codenm,          mobility,                    wane_yn, vac_eff, wane_scen, bdur, bage, pby,  pbo,  pbst_all, seas, cert, mask,  wfh_date,   title,                               vax_file,                                   result_type,    basecaseYN
mob_flt_220511,       paper_aw_flt_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         No relaxation,                      vax-covidm20220505205235.rds,                Behaviour,        N
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         6-month relaxation*,                vax-covidm20220505205235.rds,                Behaviour,        Y
mob_3mo_220511,       paper_aw_3mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         3-month relaxation,                 vax-covidm20220505205235.rds,                Behaviour,        N
mob_3wk_220511,       paper_aw_3wk_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         3-week relaxation,                  vax-covidm20220505205235.rds,                Behaviour,        N
boostNO_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,    0.0,  0.0, NULL,     0.1,  0,    0,     NA,         No booster uptake,                  vax-covidm20220505205235.rds,                Boosters,         N
boost50_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,    0.0, 0.95, NULL,     0.1,  0,    0,     NA,         Boosters for 50+,                   vax-covidm20220505205235.rds,                Boosters,         N
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Actual booster uptake*,             vax-covidm20220505205235.rds,                Boosters,         Y
boostHI_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.90, 0.98, NULL,     0.1,  0,    0,     NA,         High booster uptake,                vax-covidm20220505205235.rds,                Boosters,         N
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Central waning*,                    vax-covidm20220505205235.rds,                Waning,           Y
waneHI_220511,        paper_aw_6mo_20220509105131, yeswane, central, hiwane,    180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         High waning,                        vax-covidm20220505205235.rds,                Waning,           N
waneVHI_220511,       paper_aw_6mo_20220509105131, yeswane, central, vhiwane,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Very high waning,                   vax-covidm20220505205235.rds,                Waning,           N
seas_10_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.05,  0,    0,     NA,         10% seasonality,                    vax-covidm20220505205235.rds,                Seasonality,      N
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         20% seasonality*,                   vax-covidm20220505205235.rds,                Seasonality,      Y
seas_30_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.15,  0,    0,     NA,         30% seasonality,                    vax-covidm20220505205235.rds,                Seasonality,      N
seas_40_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals, 0.20,  0,    0,     NA,         40% seasonality,                    vax-covidm20220505205235.rds,                Seasonality,      N
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Vax 5+ 80% uptake*,                vax-covidm20220505205235.rds,                Vaccinations,     Y
vax0550_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Vax 5+ 50% uptake,                  vax-covidm202205052222505plus_50percent.rds, Vaccinations,     N
shrtbst_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,    90,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Short booster duration,             vax-covidm20220505205235.rds,                Booster duration, N
basecase_220511,      paper_aw_6mo_20220509105131, yeswane, central, central,   180,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Central booster duration*,          vax-covidm20220505205235.rds,                Booster duration, Y
longbst_220511,       paper_aw_6mo_20220509105131, yeswane, central, central,   270,  50,   0.85, 0.95, actuals,  0.1,  0,    0,     NA,         Long booster duration,              vax-covidm20220505205235.rds,                Booster duration, N
    ")

# Note that 'High booster uptake' scenario corresponds to 90% in under 50s and 96% in 50+

list = list()
for (ROW in 1:nrow(scenario_sheet)) {
    print(ROW)
    cat(scenario_sheet[ROW, codenm], "\n");
    list[[ROW]] = fread(paste0('./output/paper/may22/extratotals_', scenario_sheet[ROW, codenm], '.csv'))
}

# function to summarise required data
box_plot_scenarios <- function(variables_to_plot, times_to_plot){
    all_d = data.table(variable = NULL, when = NULL, median = NULL, q05 = NULL, q25 = NULL, q75 = NULL, q95 = NULL, scenario = NULL, result_type = NULL, basecaseYN = NULL, title = NULL)
    for (i in 1:length(list)){
        this_d = list[[i]]
        this_d = this_d[this_d$variable %in% variables_to_plot & when %in% times_to_plot,]
        this_d$scenario = rep(scenario_sheet[i, codenm], dim(this_d)[1])
        this_d$result_type = rep(scenario_sheet[i, result_type], dim(this_d)[1])
        this_d$basecaseYN = rep(scenario_sheet[i, basecaseYN], dim(this_d)[1])
        this_d$title = rep(scenario_sheet[i, title], dim(this_d)[1])
        print(scenario_sheet[i, title])
        all_d = rbind(all_d, this_d)
    }
    return(all_d)
}

d = box_plot_scenarios(c('pcr_positive_i','deaths','hosp_adm'), c('from_oct_21', 'from_jan_22', 'from_may_22_to_dec_22'))

# add units for plotting
d$unit = rep(1, dim(d)[1])
d[variable == 'pcr_positive_i', unit := 1000000]
d[variable == 'hosp_adm', unit := 1000]
d[variable == 'deaths', unit := 1000]

# force ordering of results
d$title = factor(d$title, levels = c("No relaxation",                                       
                                     "6-month relaxation*",                                 
                                     "3-month relaxation",
                                     "3-week relaxation",
                                     "No booster uptake",
                                     "Boosters for 50+",
                                     "Actual booster uptake*",                              
                                     "High booster uptake",
                                     "Central waning*",                                     
                                     "High waning",                                         
                                     "Very high waning",                                    
                                     "10% seasonality",                                     
                                     "20% seasonality*",                                    
                                     "30% seasonality",                                     
                                     "40% seasonality",                                     
                                     "Vax 5+ 80% uptake*",                                  
                                     "Vax 5+ 50% uptake",                                  
                                     "Short booster duration",  
                                     "Central booster duration*",                 
                                     "Long booster duration"))

d$result_type = factor(d$result_type, levels = c('Behaviour', 
                                                 'Booster uptake', 
                                                 'Waning', 
                                                 'Seasonality', 
                                                 'Vaccinations',
                                                 'Booster duration'))

d$variable[d$variable == 'deaths'] = 'Deaths\n(thousands)'
d$variable[d$variable == 'hosp_adm'] = 'Hospitalisations\n(thousands)'
d$variable[d$variable == 'pcr_positive_i'] = 'Infections\n(millions)'

six_colours = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')

OG_jan22 = ggplot(d[d$when == 'from_may_22_to_dec_22' & d$result_type %in% c('Behaviour', 'Waning', 'Seasonality', 'Vaccinations')], aes(x=title)) +
    geom_boxplot(aes(ymin=q05/unit, lower=q25/unit, middle=median/unit, upper=q75/unit, ymax=q95/unit, group = interaction(title, variable), colour = result_type, fill = result_type), alpha = 0.5,
                 stat = "identity") +
    labs(title='May to December 2022', x = NULL, color = 'Result type', fill = 'Result type') + 
    facet_wrap(~ variable, scales = 'free_y', strip.position = 'left', ncol = 1) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.placement = 'outside', strip.background = element_blank(), 
          axis.line.x = element_line(), text = element_text(size = 7, family = 'sans'),
          legend.position = 'bottom') + 
    geom_hline(aes(yintercept=-Inf)) + 
    geom_vline(aes(xintercept=-Inf)) + 
    coord_cartesian(clip="off") +
    scale_colour_manual(values = six_colours, aesthetics = c("fill", "colour")) 

OG_jan22

datetime <- str_replace_all(Sys.time(), "[- :BST]", "")
ggsave(paste0("./output/paperfigs/july22/fig5_", datetime, "_2column.pdf"), OG_jan22, width = 180, height = 225, units = "mm", useDingbats = FALSE)
ggsave(paste0("./output/paperfigs/july22/fig5_", datetime, "_2column.png"), OG_jan22, width = 180, height = 225, units = "mm")
