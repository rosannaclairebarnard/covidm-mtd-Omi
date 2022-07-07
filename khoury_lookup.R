library(data.table)
library(ggplot2)

# Khoury et al, Nature Medicine, Fig. 3a. Efficacy against severe infection versus efficacy against any infection.
# https://www.nature.com/articles/s41591-021-01377-8
# Data extracted from GitHub: https://github.com/InfectionAnalytics/COVID19-ProtectiveThreshold

plot_khoury_3a = function()
{
    khoury_nat_med_3a = fread("./data/Khoury_et_al_Nat_Med_fig_3a.csv")
    ggplot(khoury_nat_med_3a) +
        geom_ribbon(aes(mild, ymin = severe_lb, ymax = severe_ub), fill = "lightblue") +
        geom_line(aes(mild, severe_estimate), size = 1, colour = "#666666") +
        labs(x = "Any infection (%)", y = "Severe infection (%)")
}

khoury_get_3a = function(any_infection_pct)
{
    khoury_nat_med_3a = fread("./data/Khoury_et_al_Nat_Med_fig_3a.csv")

    low  = approx(khoury_nat_med_3a$mild, khoury_nat_med_3a$severe_lb, any_infection_pct)
    mid  = approx(khoury_nat_med_3a$mild, khoury_nat_med_3a$severe_estimate, any_infection_pct)
    high = approx(khoury_nat_med_3a$mild, khoury_nat_med_3a$severe_ub, any_infection_pct)
    
    return (list(any_infection = mid$x, low = low$y, mid = mid$y, high = high$y))
}

khoury_get_1a = function(any_infection_pct, fold_reduction)
{
    khoury_nat_med_1a = fread("./data/Khoury_et_al_Nat_Med_fig_1a.csv")
    
    neut_titre = approx(khoury_nat_med_1a$Efficacy, khoury_nat_med_1a$NeutTitreRelConvPlasma, any_infection_pct)$y
    neut_titre_new = neut_titre / fold_reduction
    efficacy_new = approx(khoury_nat_med_1a$NeutTitreRelConvPlasma, khoury_nat_med_1a$Efficacy, neut_titre_new)$y
    
    efficacy_new
}

if (0) {
    plot_khoury_3a()
    khoury_get_3a(50)
    
    points = fread(
"variant,   imm, VE_inf, VE_hosp, VE_death
   delta,  AZ 1,     43,      84,       95
   delta,  AZ 2,     63,      93,       95
   delta,  AZ WH,     48,      76,       79
   delta,  Pf 1,     62,      92,       92
   delta,  Pf 2,     80,      96,       96
   delta,  Pf WH,     60,      89,       88"
    )
    
    plot_khoury_3a() +
        geom_point(data = points, aes(x = VE_inf, y = VE_hosp, colour = variant, shape = "Hospitalisation")) +
        geom_point(data = points, aes(x = VE_inf, y = VE_death, colour = variant, shape = "Death")) +
        scale_shape_manual(values = c(0, 4)) +
        scale_colour_manual(values = "black") +
        geom_text(data = points, aes(x = VE_inf, y = VE_hosp, label = imm), hjust = 0, vjust = 0.5, nudge_x = 1)
    
    khoury_1a = fread("./data/Khoury_et_al_Nat_Med_fig_1a.csv")
    ggplot(khoury_1a[NeutTitreRelConvPlasma < 1]) +
        geom_ribbon(aes(NeutTitreRelConvPlasma, ymin = Lower, ymax = Upper), alpha = 0.4) +
        geom_line(aes(NeutTitreRelConvPlasma, Efficacy))
    
    khoury_1a[, efficacy_floor := floor(Efficacy)]
    khoury_1a[, diff := c(1, diff(efficacy_floor))]
    khoury_1a[.N, diff := 1]
    
    khoury_1a[diff == 1, paste0(NeutTitreRelConvPlasma, collapse = ", ")]
}


