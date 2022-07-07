# New plan
# 24th of September (day 632) to 30th of November: phase 1 (15,580,000)
# 1st December onwards: phase 2

get_booster_schedule = function(target_phase1 = 229000, target_phase2 = 229000, proportion_booster = rep(1, 16),
    phase1_start = 632, phase2_start = 714, last_day = 1003)
{
    # England population
    boost_plan = fread(
    "group pop
     1  3239447
     2  3539458
     3  3435579
     4  3115871
     5  3472522
     6  3771493
     7  3824652
     8  3738209
     9  3476303
    10  3638639
    11  3875351
    12  3761782
    13  3196813
    14  2784300
    15  2814128
    16  4865591")

    boost_plan[, proportion_booster := proportion_booster]
    boost_plan[, to_boost := pop * proportion_booster]
    
    working = copy(boost_plan)
    booster_schedule = NULL
    
    for (day in phase1_start:(phase2_start - 1))
    {
        booster_schedule = rbind(booster_schedule, data.table(t = day, group = working[to_boost > 0, tail(group, 1)], n = target_phase1))
        working[to_boost > 0, to_boost := pmax(0, c(head(to_boost, -1), tail(to_boost, 1) - target_phase1))]
    }
    
    for (day in phase2_start:last_day)
    {
        booster_schedule = rbind(booster_schedule, data.table(t = day, group = working[to_boost > 0, tail(group, 1)], n = target_phase2))
        working[to_boost > 0, to_boost := pmax(0, c(head(to_boost, -1), tail(to_boost, 1) - target_phase2))]
    }
    
    bsched = booster_schedule[!is.na(group), .(boost_start = min(t)), by = group]
    if (length(setdiff(1:16, bsched$group)) > 0){
        bsched = rbind(bsched, data.table(group = setdiff(1:16, bsched$group), boost_start = last_day + 1))
    }
    bsched[order(group)]
}

# get_booster_schedule(target_phase1 = 229000, target_phase2 = 229000, proportion_booster = rep(0.9, 16))
# get_booster_schedule(target_phase1 = 229000, target_phase2 = 350000, proportion_booster = rep(0.9, 16))
# get_booster_schedule(target_phase1 = 229000, target_phase2 = 500000, proportion_booster = rep(0.9, 16))


