# Probability of ICU given hospitalisation (derived from CO-CIN)
picu_cocin_func = function(age)
{
    model = qread(datapath("icu_model.qs"))
    pred = data.table(patient_age = age)
    p = predict(model, pred)
    exp(p) / (1 + exp(p))
}
picu_cocin = picu_cocin_func(0:85)

# Infection fatality rate (derived from Levin et al., preprint)
ifr_levin = 100 * exp(-7.56 + 0.121 * 0:85) / (100 + exp(-7.56 + 0.121 * 0:85)) / 100;
ifr_manual = 100 * exp(-12.75 + 0.19 * 0:85) / (100 + exp(-12.75 + 0.19 * 0:85)) / 100;

# Infection hospitalisation rate (derived from Salje et al., Science)
ihr_salje = exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85));
ihr_logit = c(0.025 * (10 - 0:10)^2 - 6.98, rep(-6.98, 20), -10.63 + 0.118 * 31:85)
ihr_manual = exp(ihr_logit) / (1 + exp(ihr_logit));

# plot(ihr_logit)
# plot(ihr_salje, type = "l")
# plot(c(ihr_manual,0), type = "l")

# Amalgamate probabilities
probabilities = data.table(age = 0:85, ihr = ihr_manual, ifr = ifr_manual, picu = picu_cocin)
probabilities[, age_group := pmin(15, age %/% 5)]
probabilities = probabilities[, lapply(.SD, mean), by = age_group, .SDcols = 2:4]

# Create model burden processes
P.hosp     = probabilities[, ihr];
P.critical = probabilities[, ihr * picu];
P.severe   = probabilities[, ihr * (1 - picu)];
P.death    = probabilities[, ifr];
# paste(signif(P.death, 4), collapse = ", ")

burden_processes = list(
    # Sero and PCR positivity
    # note -- delay should actually be gamma (symptom onset) plus normal -- to be fixed.
    # from Borremans et al
    # 0
    list(source = "newII2I3", type = "multinomial", names = "to_lfia_positive", report = c(""),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_normal(11.9 + 2.5, 5.3, 60, 0.25)$p, nrow = 1, byrow = T)),
    # 1
    list(source = "to_lfia_positive", type = "multinomial", names = "lfia_positive", report = c("ip"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_erlang(550, 10)$p, nrow = 1, byrow = T)),
    # 2
    list(source = "newEL123", type = "multinomial", names = "to_pcr_positive", report = c(""),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(2.76, 4.79, 60, 0.25)$p), nrow = 1, byrow = T)),
    # 3
    list(source = "to_pcr_positive", type = "multinomial", names = "pcr_positive", report = c("ip"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_erlang(8.47, 2)$p, nrow = 1, byrow = T)),
    
    # Pillar 2
    # 4
    list(source = "newI", type = "multinomial", names = "test", report = c("o"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(2.5 + 3.9, 2, 60, 0.25)$p, nrow = 1, byrow = T)),
    # 5
    list(source = "newI2", type = "multinomial", names = "test2", report = c("o"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(2.5 + 3.9, 2, 60, 0.25)$p, nrow = 1, byrow = T)),
    # 6
    list(source = "newI3", type = "multinomial", names = "test3", report = c("o"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(2.5 + 3.9, 2, 60, 0.25)$p, nrow = 1, byrow = T)),
    
    # Vaccine derived sero positivity (using same duration as natural infection)
    # 7
    list(source = "newVa1", type = "multinomial", names = "vaccsero_a", report = c("ip"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_erlang(550, 10)$p, nrow = 1, byrow = T)),
    # 8
    list(source = "newVb1", type = "multinomial", names = "vaccsero_b", report = c("ip"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_erlang(550, 10)$p, nrow = 1, byrow = T)),
    # 9
    list(source = "HospAdm", type = "multinomial", names = "hosp_undetected", report = c("po"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(6, 1, 60, 0.25)$p, nrow = 1, byrow = T))
)
