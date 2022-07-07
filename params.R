source("./khoury_lookup.R")

reinfection_defaults = function(params, i)
{
    params$pop[[i]]$pd_ri   = rep(0, 16);
    params$pop[[i]]$pd_ri2  = rep(0, 16);
    params$pop[[i]]$pd_ri3  = rep(0, 16);
    params$pop[[i]]$pd_r2i  = rep(0, 16);
    params$pop[[i]]$pd_r2i2 = rep(0, 16);
    params$pop[[i]]$pd_r2i3 = rep(0, 16);
    params$pop[[i]]$pd_r3i  = rep(0, 16);
    params$pop[[i]]$pd_r3i2 = rep(0, 16);
    params$pop[[i]]$pd_r3i3 = rep(0, 16);
    params$pop[[i]]$ph_rd   = rep(0, 16);
    params$pop[[i]]$ph_rd2  = rep(0, 16);
    params$pop[[i]]$ph_rd3  = rep(0, 16);
    params$pop[[i]]$ph_r2d  = rep(0, 16);
    params$pop[[i]]$ph_r2d2 = rep(0, 16);
    params$pop[[i]]$ph_r2d3 = rep(0, 16);
    params$pop[[i]]$ph_r3d  = rep(0, 16);
    params$pop[[i]]$ph_r3d2 = rep(0, 16);
    params$pop[[i]]$ph_r3d3 = rep(0, 16);
    params$pop[[i]]$pm_rd   = rep(0, 16);
    params$pop[[i]]$pm_rd2  = rep(0, 16);
    params$pop[[i]]$pm_rd3  = rep(0, 16);
    params$pop[[i]]$pm_r2d  = rep(0, 16);
    params$pop[[i]]$pm_r2d2 = rep(0, 16);
    params$pop[[i]]$pm_r2d3 = rep(0, 16);
    params$pop[[i]]$pm_r3d  = rep(0, 16);
    params$pop[[i]]$pm_r3d2 = rep(0, 16);
    params$pop[[i]]$pm_r3d3 = rep(0, 16);
    params$pop[[i]]$pt_ri   = rep(0, 16);
    params$pop[[i]]$pt_ri2  = rep(0, 16);
    params$pop[[i]]$pt_ri3  = rep(0, 16);
    params$pop[[i]]$pt_r2i  = rep(0, 16);
    params$pop[[i]]$pt_r2i2 = rep(0, 16);
    params$pop[[i]]$pt_r2i3 = rep(0, 16);
    params$pop[[i]]$pt_r3i  = rep(0, 16);
    params$pop[[i]]$pt_r3i2 = rep(0, 16);
    params$pop[[i]]$pt_r3i3 = rep(0, 16);

    return (params);
}

waning_scenario = function(WANE_YN, params, i)
{
    if (WANE_YN == "yeswane") {

        # Natural immunity - using assumption of 15% of individuals with
        # immunity having a reinfection occurring within 6 months, translated to
        # 15% of people with a complete loss of protection in 1 year

        params$pop[[i]]$wn  = rep(log(0.85)/-365,16)
        params$pop[[i]]$wn2 = rep(log(0.85)/-365,16)
        params$pop[[i]]$wn3 = rep(log(0.85)/-365,16)

        # Vaccine induced immunity - updated assumptions on 11th February 2022
        # using Andrews et al. NEJM article on "Duration of Protection against
        # Mild and Severe Disease by Covid-19 Vaccines", available online at:
        # https://www.nejm.org/doi/full/10.1056/NEJMoa2115481

        # We calculate the % change in protection against hospitalisation
        # between the measured level at 1 week post second vaccine dose and the
        # measured level at 20+ weeks post second vaccine dose, for the
        # AstraZeneca and Pfizer vaccines

        # AZ dose 2 + waned overall waning back to susceptible rate
        # log(0.851)/-140, corresponding to loss of 14.9% protection after 20 weeks
        # Pfizer dose 2 waning
        # log(0.923)/-140, corresponding to loss of 7.7% protection in 20 weeks

        # AZ waning
        params$pop[[i]]$wva1 = rep(0, 16) # zero since all first doses get followed up with second doses
        params$pop[[i]]$wva2 = rep(0, 16) # zero since after dose 2 duration, individuals either wane (move to va3) or are boosted and return to start of va2
        params$pop[[i]]$wva3 = rep(log(0.851)/-140, 16) # may need to revisit this, currently setup related to waning of VE against severe outcomes
        params$pop[[i]]$wva21 = rep(0, 16)

        # Pfizer waning
        params$pop[[i]]$wvb1 = rep(0, 16) # zero since all first doses get followed up with second doses
        params$pop[[i]]$wvb2 = rep(0, 16) # zero since after dose 2 duration, individuals either wane (move to vb3) or are boosted and return to start of vb2
        params$pop[[i]]$wvb3 = rep(log(0.923)/-140, 16) # may need to revisit this, currently setup related to waning of VE against severe outcomes
        params$pop[[i]]$wvb21 = rep(0, 16)

    } else if (WANE_YN == "nowane") {
        params$pop[[i]]$wn  = rep(0, 16)
        params$pop[[i]]$wn2 = rep(0, 16)
        params$pop[[i]]$wn3 = rep(0, 16)
        params$pop[[i]]$wva1 = rep(0, 16)
        params$pop[[i]]$wva2 = rep(0, 16)
        params$pop[[i]]$wva21 = rep(0, 16)
        params$pop[[i]]$wvb1 = rep(0, 16)
        params$pop[[i]]$wvb2 = rep(0, 16)
        params$pop[[i]]$wvb21 = rep(0, 16)
    } else {
        stop("WANE_YN must be yeswane or nowane.")
    }

    return (params)
}

VE_scenario = function(VAC_EFF, P_BOOST_VEC, params, i, delta_fold_reduction)
{
    # set up vaccine efficacy parameters
    VE = function(x, given = 0)
    {
        ifelse(rep_len(given, 16) == 1, 1, (rep_len(x, 16) - given) / (1 - given))
    }

    if (VAC_EFF == "central") {
        L3 = 0.75
    } else if (VAC_EFF == "centralHI") {
        L3 = 1.00
    } else if (VAC_EFF == "centralLO") {
        L3 = 0.50
    } else {
        stop("VAC_EFF must be central, centralHI, or centralLO")
    }

    # Probability of having a booster dose
    params$pop[[i]]$p_boost_va2 = P_BOOST_VEC
    params$pop[[i]]$p_boost_vb2 = P_BOOST_VEC

    # Delays between vaccine compartments

    # Dose 1 to dose 2
    # initially the first to second dose delay was 21.5 days on average up
    # to and including second doses on the 26th January 2021 when JCVI
    # issued guidance on dosing schedules. dVa1 and dVa2 are updated in
    # cpp_funcs.R code to a longer delay later in the rollout

    # note that 7.5 = 21.5 + 14 - 28
    #               = 21.5 + d2_delay_to_efficacy - d1_delay_to_efficacy
    params$pop[[i]]$dVa1 = cm_delay_gamma(7.5, 1000, t_max = 10, t_step = 0.25)$p
    params$pop[[i]]$dVb1 = cm_delay_gamma(7.5, 1000, t_max = 10, t_step = 0.25)$p

    # Dose 2 to dose 3 / dose 2 + waned
    # JCVI recommends 6+ month delay between second dose and booster dose
    # Vaccine rollout (first doses) began 8th December 2020
    # Booster rollout began 24th September 2021
    # We assume dose 1 efficacy kicks in after 28 days
    # We assume dose 2 efficacy kicks in after 14 days
    # Average first:second dose delay is 71 days (using data from 3rd October 2021)
    # 290 days between 8th December 2020 and 24th September 2021
    # 28 days delay to first dose efficacy
    # + (71 - 14) = average d1:d2 delay - delay to second dose efficacy
    # 28 + (57) = 85
    # 290 - 85 = 205 days in dose 2
    params$pop[[i]]$dVa2 = cm_delay_gamma(205, 2000, t_max = 225, t_step = 0.25)$p
    params$pop[[i]]$dVb2 = cm_delay_gamma(205, 2000, t_max = 225, t_step = 0.25)$p

    # AstraZeneca infection
    # Dose 1
    params$pop[[i]]$ei_va1  = VE(0.70)
    params$pop[[i]]$ei2_va1 = VE(0.70)
    params$pop[[i]]$ei3_va1 = VE(0.43)
    # Dose 2
    params$pop[[i]]$ei_va2  = VE(0.75)
    params$pop[[i]]$ei2_va2 = VE(0.75)
    params$pop[[i]]$ei3_va2 = VE(0.63)
    # Dose 3
    params$pop[[i]]$ei_va3  = VE(L3 * 0.513)  ## These remain the same as before for Khoury update
    params$pop[[i]]$ei2_va3 = VE(L3 * 0.513)  ## These remain the same as before for Khoury update
    params$pop[[i]]$ei3_va3 = VE(L3 * 0.431)  ## These remain the same as before for Khoury update

    # Pfizer / Moderna infection
    # Dose 1
    params$pop[[i]]$ei_vb1  = VE(0.70)
    params$pop[[i]]$ei2_vb1 = VE(0.70)
    params$pop[[i]]$ei3_vb1 = VE(0.62)
    # Dose 2
    params$pop[[i]]$ei_vb2  = VE(0.85)
    params$pop[[i]]$ei2_vb2 = VE(0.85)
    params$pop[[i]]$ei3_vb2 = VE(0.80)
    # Dose 3
    params$pop[[i]]$ei_vb3  = VE(L3 * 0.613)  ## These remain the same as before
    params$pop[[i]]$ei2_vb3 = VE(L3 * 0.613)  ## These remain the same as before
    params$pop[[i]]$ei3_vb3 = VE(L3 * 0.574)   ## These remain the same as before

    # AstraZeneca symptomatic disease
    # Dose 1
    params$pop[[i]]$ed_va1i  = VE(0.70,   params$pop[[i]]$ei_va1)
    params$pop[[i]]$ed_va1i2 = VE(0.70,   params$pop[[i]]$ei2_va1)
    params$pop[[i]]$ed_va1i3 = VE(0.52,   params$pop[[i]]$ei3_va1)
    # Dose 2
    params$pop[[i]]$ed_va2i  = VE(0.80,   params$pop[[i]]$ei_va2)
    params$pop[[i]]$ed_va2i2 = VE(0.80,   params$pop[[i]]$ei2_va2)
    params$pop[[i]]$ed_va2i3 = VE(0.65,   params$pop[[i]]$ei3_va2)
    # Dose 3
    ## params$pop[[i]]$ed_va3i  = VE(L3 * 0.6,   params$pop[[i]]$ei_va3)
    ## params$pop[[i]]$ed_va3i2 = VE(L3 * 0.6,   params$pop[[i]]$ei2_va3)
    ## params$pop[[i]]$ed_va3i3 = VE(L3 * 0.49,  params$pop[[i]]$ei3_va3)
    params$pop[[i]]$ed_va3i  = VE(0.547,   params$pop[[i]]$ei_va3)  ## Removed L3 multiplier
    params$pop[[i]]$ed_va3i2 = VE(0.547,   params$pop[[i]]$ei2_va3) ## Removed L3 multiplier
    params$pop[[i]]$ed_va3i3 = VE(0.445,  params$pop[[i]]$ei3_va3) ## Removed L3 multiplier

    # Pfizer / Moderna symptomatic disease
    # Dose 1
    params$pop[[i]]$ed_vb1i  = VE(0.70,   params$pop[[i]]$ei_vb1)
    params$pop[[i]]$ed_vb1i2 = VE(0.70,   params$pop[[i]]$ei2_vb1)
    params$pop[[i]]$ed_vb1i3 = VE(0.62,   params$pop[[i]]$ei3_vb1)
    # Dose 2
    params$pop[[i]]$ed_vb2i  = VE(0.90,   params$pop[[i]]$ei_vb2)
    params$pop[[i]]$ed_vb2i2 = VE(0.90,   params$pop[[i]]$ei2_vb2)
    params$pop[[i]]$ed_vb2i3 = VE(0.81,   params$pop[[i]]$ei3_vb2)
    # Dose 3
    ## params$pop[[i]]$ed_vb3i  = VE(L3 * 0.68,   params$pop[[i]]$ei_vb3)
    ## params$pop[[i]]$ed_vb3i2 = VE(L3 * 0.68,   params$pop[[i]]$ei2_vb3)
    ## params$pop[[i]]$ed_vb3i3 = VE(L3 * 0.61,   params$pop[[i]]$ei3_vb3)
    params$pop[[i]]$ed_vb3i  = VE(0.646,   params$pop[[i]]$ei_vb3)  ## Removed L3 multiplier
    params$pop[[i]]$ed_vb3i2 = VE(0.646,   params$pop[[i]]$ei2_vb3) ## Removed L3 multiplier
    params$pop[[i]]$ed_vb3i3 = VE(0.582,   params$pop[[i]]$ei3_vb3) ## Removed L3 multiplier

    kh_inf_to_sev = function(ve_inf) {
        khoury_get_3a(ve_inf * 100)$mid / 100
    }

    # Checks on above function
    if (0) {
        khoury_get_3a(33)$mid # 75.7
        kh_inf_to_sev(0.33)   # 0.757
        VE(kh_inf_to_sev(0.33), 1 - (1 - 0.33) * (1 - 0.0)) # 0.638
        1 - (1 - 0.33) * (1 - 0.0) * (1 - 0.638) # 0.757
    }

    # AstraZeneca hospitalisation
    # Dose 1
    params$pop[[i]]$eh_va1d  = VE(0.85, 1 - (1 - params$pop[[i]]$ei_va1)  * (1 - params$pop[[i]]$ed_va1i))
    params$pop[[i]]$eh_va1d2 = VE(0.85, 1 - (1 - params$pop[[i]]$ei2_va1) * (1 - params$pop[[i]]$ed_va1i2))
    params$pop[[i]]$eh_va1d3 = VE(0.84, 1 - (1 - params$pop[[i]]$ei3_va1) * (1 - params$pop[[i]]$ed_va1i3))
    # Dose 2
    params$pop[[i]]$eh_va2d  = VE(0.90, 1 - (1 - params$pop[[i]]$ei_va2)  * (1 - params$pop[[i]]$ed_va2i))
    params$pop[[i]]$eh_va2d2 = VE(0.90, 1 - (1 - params$pop[[i]]$ei2_va2) * (1 - params$pop[[i]]$ed_va2i2))
    params$pop[[i]]$eh_va2d3 = VE(0.93, 1 - (1 - params$pop[[i]]$ei3_va2) * (1 - params$pop[[i]]$ed_va2i3))
    # Dose 3
    ## params$pop[[i]]$eh_va3d  = VE(L3 * 0.74, 1 - (1 - params$pop[[i]]$ei_va3) * (1 - params$pop[[i]]$ed_va3i))
    ## params$pop[[i]]$eh_va3d2 = VE(L3 * 0.74, 1 - (1 - params$pop[[i]]$ei2_va3) * (1 - params$pop[[i]]$ed_va3i2))
    ## params$pop[[i]]$eh_va3d3 = VE(L3 * 0.76, 1 - (1 - params$pop[[i]]$ei3_va3) * (1 - params$pop[[i]]$ed_va3i3))
    params$pop[[i]]$eh_va3d  = VE(kh_inf_to_sev(params$pop[[i]]$ei_va3),  1 - (1 - params$pop[[i]]$ei_va3)  * (1 - params$pop[[i]]$ed_va3i))  ## added kh_inf_to_sev
    params$pop[[i]]$eh_va3d2 = VE(kh_inf_to_sev(params$pop[[i]]$ei2_va3), 1 - (1 - params$pop[[i]]$ei2_va3) * (1 - params$pop[[i]]$ed_va3i2)) ## added kh_inf_to_sev
    params$pop[[i]]$eh_va3d3 = VE(kh_inf_to_sev(params$pop[[i]]$ei3_va3), 1 - (1 - params$pop[[i]]$ei3_va3) * (1 - params$pop[[i]]$ed_va3i3)) ## added kh_inf_to_sev

    # Pfizer / Moderna hospitalisation
    # Dose 1
    params$pop[[i]]$eh_vb1d  = VE(0.85, 1 - (1 - params$pop[[i]]$ei_vb1)  * (1 - params$pop[[i]]$ed_vb1i))
    params$pop[[i]]$eh_vb1d2 = VE(0.85, 1 - (1 - params$pop[[i]]$ei2_vb1) * (1 - params$pop[[i]]$ed_vb1i2))
    params$pop[[i]]$eh_vb1d3 = VE(0.92, 1 - (1 - params$pop[[i]]$ei3_vb1) * (1 - params$pop[[i]]$ed_vb1i3))
    # Dose 2
    params$pop[[i]]$eh_vb2d  = VE(0.95, 1 - (1 - params$pop[[i]]$ei_vb2)  * (1 - params$pop[[i]]$ed_vb2i))
    params$pop[[i]]$eh_vb2d2 = VE(0.95, 1 - (1 - params$pop[[i]]$ei2_vb2) * (1 - params$pop[[i]]$ed_vb2i2))
    params$pop[[i]]$eh_vb2d3 = VE(0.96, 1 - (1 - params$pop[[i]]$ei3_vb2) * (1 - params$pop[[i]]$ed_vb2i3))
    # Dose 3 (currently same as dose 2)
    ## params$pop[[i]]$eh_vb3d  = VE(L3 * 0.88, 1 - (1 - params$pop[[i]]$ei_vb3)  * (1 - params$pop[[i]]$ed_vb3i))
    ## params$pop[[i]]$eh_vb3d2 = VE(L3 * 0.88, 1 - (1 - params$pop[[i]]$ei2_vb3)  * (1 - params$pop[[i]]$ed_vb3i2))
    ## params$pop[[i]]$eh_vb3d3 = VE(L3 * 0.89, 1 - (1 - params$pop[[i]]$ei3_vb3)  * (1 - params$pop[[i]]$ed_vb3i3))
    params$pop[[i]]$eh_vb3d  = VE(kh_inf_to_sev(params$pop[[i]]$ei_vb3),  1 - (1 - params$pop[[i]]$ei_vb3 ) * (1 - params$pop[[i]]$ed_vb3i))  ## ditto
    params$pop[[i]]$eh_vb3d2 = VE(kh_inf_to_sev(params$pop[[i]]$ei2_vb3), 1 - (1 - params$pop[[i]]$ei2_vb3) * (1 - params$pop[[i]]$ed_vb3i2)) ## ditto
    params$pop[[i]]$eh_vb3d3 = VE(kh_inf_to_sev(params$pop[[i]]$ei3_vb3), 1 - (1 - params$pop[[i]]$ei3_vb3) * (1 - params$pop[[i]]$ed_vb3i3)) ## ditto

    # AstraZeneca mortality
    # Dose 1
    params$pop[[i]]$em_va1d  = VE(0.85, 1 - (1 - params$pop[[i]]$ei_va1)  * (1 - params$pop[[i]]$ed_va1i))
    params$pop[[i]]$em_va1d2 = VE(0.85, 1 - (1 - params$pop[[i]]$ei2_va1) * (1 - params$pop[[i]]$ed_va1i2))
    params$pop[[i]]$em_va1d3 = VE(0.95, 1 - (1 - params$pop[[i]]$ei3_va1) * (1 - params$pop[[i]]$ed_va1i3))
    # Dose 2
    params$pop[[i]]$em_va2d  = VE(0.95, 1 - (1 - params$pop[[i]]$ei_va2)  * (1 - params$pop[[i]]$ed_va2i))
    params$pop[[i]]$em_va2d2 = VE(0.95, 1 - (1 - params$pop[[i]]$ei2_va2) * (1 - params$pop[[i]]$ed_va2i2))
    params$pop[[i]]$em_va2d3 = VE(0.95, 1 - (1 - params$pop[[i]]$ei3_va2) * (1 - params$pop[[i]]$ed_va2i3))
    # Dose 3
    ## params$pop[[i]]$em_va3d  = VE(L3 * 0.79, 1 - (1 - params$pop[[i]]$ei_va3)  * (1 - params$pop[[i]]$ed_va3i))
    ## params$pop[[i]]$em_va3d2 = VE(L3 * 0.79, 1 - (1 - params$pop[[i]]$ei2_va3) * (1 - params$pop[[i]]$ed_va3i2))
    ## params$pop[[i]]$em_va3d3 = VE(L3 * 0.79, 1 - (1 - params$pop[[i]]$ei3_va3) * (1 - params$pop[[i]]$ed_va3i3))
    params$pop[[i]]$em_va3d  = VE(kh_inf_to_sev(params$pop[[i]]$ei_va3),  1 - (1 - params$pop[[i]]$ei_va3 ) * (1 - params$pop[[i]]$ed_va3i))  ## ditto
    params$pop[[i]]$em_va3d2 = VE(kh_inf_to_sev(params$pop[[i]]$ei2_va3), 1 - (1 - params$pop[[i]]$ei2_va3) * (1 - params$pop[[i]]$ed_va3i2)) ## ditto
    params$pop[[i]]$em_va3d3 = VE(kh_inf_to_sev(params$pop[[i]]$ei3_va3), 1 - (1 - params$pop[[i]]$ei3_va3) * (1 - params$pop[[i]]$ed_va3i3)) ## ditto

    # Pfizer / Moderna mortality
    # Dose 1
    params$pop[[i]]$em_vb1d  = VE(0.85, 1 - (1 - params$pop[[i]]$ei_vb1)  * (1 - params$pop[[i]]$ed_vb1i))
    params$pop[[i]]$em_vb1d2 = VE(0.85, 1 - (1 - params$pop[[i]]$ei2_vb1) * (1 - params$pop[[i]]$ed_vb1i2))
    params$pop[[i]]$em_vb1d3 = VE(0.92, 1 - (1 - params$pop[[i]]$ei3_vb1) * (1 - params$pop[[i]]$ed_vb1i3))
    # Dose 2
    params$pop[[i]]$em_vb2d  = VE(0.95, 1 - (1 - params$pop[[i]]$ei_vb2)  * (1 - params$pop[[i]]$ed_vb2i))
    params$pop[[i]]$em_vb2d2 = VE(0.95, 1 - (1 - params$pop[[i]]$ei2_vb2) * (1 - params$pop[[i]]$ed_vb2i2))
    params$pop[[i]]$em_vb2d3 = VE(0.96, 1 - (1 - params$pop[[i]]$ei3_vb2) * (1 - params$pop[[i]]$ed_vb2i3))
    # Dose 3
    ## params$pop[[i]]$em_vb3d  = VE(L3 * 0.87, 1 - (1 - params$pop[[i]]$ei_vb3 ) * (1 - params$pop[[i]]$ed_vb3i))
    ## params$pop[[i]]$em_vb3d2 = VE(L3 * 0.87, 1 - (1 - params$pop[[i]]$ei2_vb3) * (1 - params$pop[[i]]$ed_vb3i2))
    ## params$pop[[i]]$em_vb3d3 = VE(L3 * 0.88, 1 - (1 - params$pop[[i]]$ei3_vb3) * (1 - params$pop[[i]]$ed_vb3i3))
    params$pop[[i]]$em_vb3d  = VE(kh_inf_to_sev(params$pop[[i]]$ei_vb3),  1 - (1 - params$pop[[i]]$ei_vb3 ) * (1 - params$pop[[i]]$ed_vb3i))  ## ditto
    params$pop[[i]]$em_vb3d2 = VE(kh_inf_to_sev(params$pop[[i]]$ei2_vb3), 1 - (1 - params$pop[[i]]$ei2_vb3) * (1 - params$pop[[i]]$ed_vb3i2)) ## ditto
    params$pop[[i]]$em_vb3d3 = VE(kh_inf_to_sev(params$pop[[i]]$ei3_vb3), 1 - (1 - params$pop[[i]]$ei3_vb3) * (1 - params$pop[[i]]$ed_vb3i3)) ## ditto

    # AstraZeneca onward transmission
    # Dose 1
    params$pop[[i]]$et_va1i  = VE(0.47)
    params$pop[[i]]$et_va1i2 = VE(0.47)
    params$pop[[i]]$et_va1i3 = VE(0.05)
    # Dose 2
    params$pop[[i]]$et_va2i  = VE(0.47)
    params$pop[[i]]$et_va2i2 = VE(0.47)
    params$pop[[i]]$et_va2i3 = VE(0.27)
    # Dose 3
    params$pop[[i]]$et_va3i  = VE(L3 * 0.380)
    params$pop[[i]]$et_va3i2 = VE(L3 * 0.380)
    params$pop[[i]]$et_va3i3 = VE(L3 * 0.218)

    # Pfizer / Moderna onward transmission
    # Dose 1
    params$pop[[i]]$et_vb1i  = VE(0.47)
    params$pop[[i]]$et_vb1i2 = VE(0.47)
    params$pop[[i]]$et_vb1i3 = VE(0.24)
    # Dose 2
    params$pop[[i]]$et_vb2i  = VE(0.47)
    params$pop[[i]]$et_vb2i2 = VE(0.47)
    params$pop[[i]]$et_vb2i3 = VE(0.37)
    # Dose 3
    params$pop[[i]]$et_vb3i  = VE(L3 * 0.403)
    params$pop[[i]]$et_vb3i2 = VE(L3 * 0.403)
    params$pop[[i]]$et_vb3i3 = VE(L3 * 0.317)

    # Cross protection from reinfection (3 virus variants)
    params$pop[[i]]$pi_r   = VE(1.00);
    params$pop[[i]]$pi2_r  = VE(1.00);
    params$pop[[i]]$pi3_r  = VE(khoury_get_1a(85, delta_fold_reduction) / 85);
    params$pop[[i]]$pi_r2  = VE(1.00);
    params$pop[[i]]$pi2_r2 = VE(1.00);
    params$pop[[i]]$pi3_r2 = VE(khoury_get_1a(85, delta_fold_reduction) / 85);
    params$pop[[i]]$pi_r3  = VE(1.00);
    params$pop[[i]]$pi2_r3 = VE(1.00);
    params$pop[[i]]$pi3_r3 = VE(1.00);

    params$pop[[i]]$pd_ri   = params$pop[[i]]$ed_vb2i;
    params$pop[[i]]$pd_ri2  = params$pop[[i]]$ed_vb2i2;
    params$pop[[i]]$pd_ri3  = params$pop[[i]]$ed_vb2i3;
    params$pop[[i]]$pd_r2i  = params$pop[[i]]$ed_vb2i;
    params$pop[[i]]$pd_r2i2 = params$pop[[i]]$ed_vb2i2;
    params$pop[[i]]$pd_r2i3 = params$pop[[i]]$ed_vb2i3;
    params$pop[[i]]$pd_r3i  = params$pop[[i]]$ed_vb2i;
    params$pop[[i]]$pd_r3i2 = params$pop[[i]]$ed_vb2i2;
    params$pop[[i]]$pd_r3i3 = params$pop[[i]]$ed_vb2i3;

    params$pop[[i]]$ph_rd   = VE(kh_inf_to_sev(params$pop[[i]]$pi_r  ), 1 - (1 - params$pop[[i]]$pi_r ) * (1 - params$pop[[i]]$pd_ri))
    params$pop[[i]]$ph_rd2  = VE(kh_inf_to_sev(params$pop[[i]]$pi2_r ), 1 - (1 - params$pop[[i]]$pi2_r) * (1 - params$pop[[i]]$pd_ri2))
    params$pop[[i]]$ph_rd3  = VE(kh_inf_to_sev(params$pop[[i]]$pi3_r ), 1 - (1 - params$pop[[i]]$pi3_r) * (1 - params$pop[[i]]$pd_ri3))
    params$pop[[i]]$ph_r2d  = VE(kh_inf_to_sev(params$pop[[i]]$pi_r2 ), 1 - (1 - params$pop[[i]]$pi_r2 ) * (1 - params$pop[[i]]$pd_r2i))
    params$pop[[i]]$ph_r2d2 = VE(kh_inf_to_sev(params$pop[[i]]$pi2_r2), 1 - (1 - params$pop[[i]]$pi2_r2) * (1 - params$pop[[i]]$pd_r2i2))
    params$pop[[i]]$ph_r2d3 = VE(kh_inf_to_sev(params$pop[[i]]$pi3_r2), 1 - (1 - params$pop[[i]]$pi3_r2) * (1 - params$pop[[i]]$pd_r2i3))
    params$pop[[i]]$ph_r3d  = VE(kh_inf_to_sev(params$pop[[i]]$pi_r3 ), 1 - (1 - params$pop[[i]]$pi_r3 ) * (1 - params$pop[[i]]$pd_r3i))
    params$pop[[i]]$ph_r3d2 = VE(kh_inf_to_sev(params$pop[[i]]$pi2_r3), 1 - (1 - params$pop[[i]]$pi2_r3) * (1 - params$pop[[i]]$pd_r3i2))
    params$pop[[i]]$ph_r3d3 = VE(kh_inf_to_sev(params$pop[[i]]$pi3_r3), 1 - (1 - params$pop[[i]]$pi3_r3) * (1 - params$pop[[i]]$pd_r3i3))
    params$pop[[i]]$pm_rd   = VE(kh_inf_to_sev(params$pop[[i]]$pi_r  ), 1 - (1 - params$pop[[i]]$pi_r ) * (1 - params$pop[[i]]$pd_ri))
    params$pop[[i]]$pm_rd2  = VE(kh_inf_to_sev(params$pop[[i]]$pi2_r ), 1 - (1 - params$pop[[i]]$pi2_r) * (1 - params$pop[[i]]$pd_ri2))
    params$pop[[i]]$pm_rd3  = VE(kh_inf_to_sev(params$pop[[i]]$pi3_r ), 1 - (1 - params$pop[[i]]$pi3_r) * (1 - params$pop[[i]]$pd_ri3))
    params$pop[[i]]$pm_r2d  = VE(kh_inf_to_sev(params$pop[[i]]$pi_r2 ), 1 - (1 - params$pop[[i]]$pi_r2 ) * (1 - params$pop[[i]]$pd_r2i))
    params$pop[[i]]$pm_r2d2 = VE(kh_inf_to_sev(params$pop[[i]]$pi2_r2), 1 - (1 - params$pop[[i]]$pi2_r2) * (1 - params$pop[[i]]$pd_r2i2))
    params$pop[[i]]$pm_r2d3 = VE(kh_inf_to_sev(params$pop[[i]]$pi3_r2), 1 - (1 - params$pop[[i]]$pi3_r2) * (1 - params$pop[[i]]$pd_r2i3))
    params$pop[[i]]$pm_r3d  = VE(kh_inf_to_sev(params$pop[[i]]$pi_r3 ), 1 - (1 - params$pop[[i]]$pi_r3 ) * (1 - params$pop[[i]]$pd_r3i))
    params$pop[[i]]$pm_r3d2 = VE(kh_inf_to_sev(params$pop[[i]]$pi2_r3), 1 - (1 - params$pop[[i]]$pi2_r3) * (1 - params$pop[[i]]$pd_r3i2))
    params$pop[[i]]$pm_r3d3 = VE(kh_inf_to_sev(params$pop[[i]]$pi3_r3), 1 - (1 - params$pop[[i]]$pi3_r3) * (1 - params$pop[[i]]$pd_r3i3))

    params$pop[[i]]$pt_ri   = params$pop[[i]]$et_vb2i;
    params$pop[[i]]$pt_ri2  = params$pop[[i]]$et_vb2i2;
    params$pop[[i]]$pt_ri3  = params$pop[[i]]$et_vb2i3;
    params$pop[[i]]$pt_r2i  = params$pop[[i]]$et_vb2i;
    params$pop[[i]]$pt_r2i2 = params$pop[[i]]$et_vb2i2;
    params$pop[[i]]$pt_r2i3 = params$pop[[i]]$et_vb2i3;
    params$pop[[i]]$pt_r3i  = params$pop[[i]]$et_vb2i;
    params$pop[[i]]$pt_r3i2 = params$pop[[i]]$et_vb2i2;
    params$pop[[i]]$pt_r3i3 = params$pop[[i]]$et_vb2i3;

    return (params)
}


extra_waning_scenario = function(wane_scen, params, i){

    if (wane_scen == "hiwane"){

        # set dose 2 waning equal to dose 3 waning
        params$pop[[i]]$wva2 = rep(log(0.851)/-140, 16)
        params$pop[[i]]$wvb2 = rep(log(0.923)/-140, 16)

    } else if(wane_scen == "vhiwane"){

        # set dose 2 waning equal to dose 3 waning and then double both rates
        params$pop[[i]]$wva2 = rep(log(0.851)/-70, 16)
        params$pop[[i]]$wvb2 = rep(log(0.923)/-70, 16)
        params$pop[[i]]$wva3 = rep(log(0.851)/-70, 16)
        params$pop[[i]]$wvb3 = rep(log(0.923)/-70, 16)

        # double the rate of 'natural' waning
        params$pop[[i]]$wn  = rep(log(0.85)/-182.5, 16)
        params$pop[[i]]$wn2 = rep(log(0.85)/-182.5, 16)
        params$pop[[i]]$wn3 = rep(log(0.85)/-182.5, 16)

    } else {
        stop("wane_scen must be hiwane or vhiwane.")
    }

    return(params)
}
