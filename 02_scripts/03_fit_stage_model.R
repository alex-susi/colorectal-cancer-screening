## 03 - Fit three-window reminder model ========================================
##
## Input : 01_data/clean/kits.rds
## Model : 02_scripts/stan/kit_returns_stage.stan
## Output: 03_output/models/stage_model/        posterior: CmdStan CSVs, summaries
##         03_output/models/stage_model_prior/  prior-only draws (prior predictive check)
##
## Run from the project root (the .Rproj directory).
##
## Variants (first command-line argument):
##   (none)  main model
##   no_q4   sensitivity: drops the Oct-Dec (Q4) term; output in
##           03_output/models/sens_no_q4/, no prior-only run

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(cmdstanr)
  library(posterior)
})

variant <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(variant)) variant <- "main"
stopifnot(variant %in% c("main", "no_q4"))
run_prior <- variant == "main"

stan_file <- "02_scripts/stan/kit_returns_stage.stan"
out_dir   <- switch(variant, main = "03_output/models/stage_model",
                             no_q4 = "03_output/models/sens_no_q4")
prior_dir <- "03_output/models/stage_model_prior"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (run_prior) dir.create(prior_dir, showWarnings = FALSE, recursive = TRUE)
message("Variant: ", variant, " -> ", out_dir)

covid_start <- as.Date("2020-03-15")
covid_end   <- as.Date("2021-06-30")
remind_start <- 14   # first text on day 14
remind_end   <- 90   # texts every 2 weeks until day 90, then stop




## 1 - Model data ===============================================================

kits <- read_rds("01_data/clean/kits.rds")

d <- kits %>%
  mutate(outcome   = case_when(returned == 0L                 ~ 4L,
                               is.na(days_to_return)          ~ 5L,   # returned, date unknown
                               days_to_return <= remind_start ~ 1L,
                               days_to_return <= remind_end   ~ 2L,
                               TRUE                           ~ 3L),
         covid     = as.integer(dispense_date >= covid_start & dispense_date <= covid_end),
         q4        = as.integer(dispense_month >= 10),
         first_kit = as.integer(n_prior_kits == 0),
         pat       = as.integer(factor(patient_id)))


# Shared covariates: treatment contrasts against the largest / reference level
x_formula <- ~ age_c10 + family_size_c + first_kit + q4 + covid + time_yrs +
  sex + relationship + employment + race_eth + language + text_consent
if (variant == "no_q4") x_formula <- update(x_formula, ~ . - q4)
X <- model.matrix(x_formula, data = d)[, -1]
C_con <- model.matrix(~ text_consent, data = d)[, -1, drop = FALSE]


# Reporting groups for subgroup reminder effects: one global id per factor level
rep_factors <- c("sex", "race_eth", "language", "employment", "relationship", "text_consent")
grp_levels  <- bind_rows(lapply(rep_factors, function(f)
  tibble(factor = f, level = levels(d[[f]]))))
grp_levels$id <- seq_len(nrow(grp_levels))
grp <- sapply(rep_factors, function(f)
  match(paste(f, as.character(d[[f]])), paste(grp_levels$factor, grp_levels$level)))

stan_data <- list(N = nrow(d), 
                  J = max(d$pat), 
                  pat = d$pat, 
                  outcome = d$outcome,
                  post = d$post, 
                  covid = d$covid,
                  K = ncol(X), 
                  X = X, 
                  beta_sd = rep(1, ncol(X)),
                  n_sex  = nlevels(d$sex),          sex  = as.integer(d$sex),
                  n_race = nlevels(d$race_eth),     race = as.integer(d$race_eth),
                  n_lang = nlevels(d$language),     lang = as.integer(d$language),
                  n_emp  = nlevels(d$employment),   emp  = as.integer(d$employment),
                  n_rel  = nlevels(d$relationship), rel  = as.integer(d$relationship),
                  age_c10 = d$age_c10, 
                  family_size_c = d$family_size_c,
                  K_con = ncol(C_con), C_con = C_con,
                  n_fac = length(rep_factors), 
                  n_grp = nrow(grp_levels), grp = grp,
                  save_log_lik = 0L, prior_only = 0L)

message("Kits: ", stan_data$N, " | patients: ", stan_data$J,
        " | outcomes (1 day 0-14, 2 day 15-90, 3 day 91+, 4 never, 5 date unknown): ",
        paste(tabulate(d$outcome, 5), collapse = " / "))

write_rds(list(X_cols = colnames(X), 
               con_cols = colnames(C_con), 
               grp_levels = grp_levels,
               patients = d %>% distinct(pat, patient_id)),
          file.path(out_dir, "design_labels.rds"))





## 2 - Compile & sample ========================================================

mod <- cmdstan_model(stan_file)

fit <- mod$sample(data            = stan_data,
                  chains          = 4,
                  parallel_chains = 4,
                  iter_warmup     = 1000,
                  iter_sampling   = 1500,
                  adapt_delta     = 0.9,
                  seed            = 2022,
                  refresh         = 250,
                  output_dir      = out_dir,
                  output_basename = "stage_model")
write_rds(fit$time(), file.path(out_dir, "timing.rds"))
message("Sampling time (s): total ", round(fit$time()$total, 1))


# Prior-only run: same model, likelihood switched off (main model only)
if (run_prior) {
  fit_prior <- mod$sample(data            = modifyList(stan_data, list(prior_only = 1L)),
                          chains          = 4,
                          parallel_chains = 4,
                          iter_warmup     = 500,
                          iter_sampling   = 500,
                          seed            = 2022,
                          refresh         = 0,
                          output_dir      = prior_dir,
                          output_basename = "stage_model_prior")
}





## 3 - Diagnostics =============================================================

diag <- fit$diagnostic_summary()
print(diag)

core_pars <- c("alpha", "beta", "beta_post", "beta_covid_late", "beta_remind", "beta_carry",
               "tau", "gamma_age", "gamma_fam", "gamma_con", 
               "gamma_carry_con", "sigma_u", "rho_u",
               "d_sex", "d_race", "d_lang", "d_emp", "d_rel")
summ <- fit$summary(variables = core_pars, "mean", "sd",
                    ~ quantile(.x, probs = c(0.025, 0.5, 0.975)), 
                    "rhat", "ess_bulk", "ess_tail")
message("Max R-hat: ", round(max(summ$rhat, na.rm = TRUE), 3),
        " | min bulk ESS: ", round(min(summ$ess_bulk, na.rm = TRUE)),
        " | min tail ESS: ", round(min(summ$ess_tail, na.rm = TRUE)))





## 4 - Label & save summaries ===================================================

lab <- c(setNames(paste0("beta[", seq_along(colnames(X)), "]"), colnames(X)),
         setNames(paste0("gamma_con[", seq_along(colnames(C_con)), "]"), 
                  paste0("remind x ", colnames(C_con))),
         setNames(paste0("gamma_carry_con[", seq_along(colnames(C_con)), "]"), 
                  paste0("carry x ", colnames(C_con))))
summ <- summ %>%
  mutate(label = coalesce(names(lab)[match(variable, lab)], variable), .after = variable)
write_csv(summ, file.path(out_dir, "summary_parameters.csv"))

effects <- fit$summary(variables = c("or_remind", "or_carry", "or_post", 
                                     "p_on_post", "p_off_post", "rd_post"),
                       "mean", ~ quantile(.x, probs = c(0.025, 0.5, 0.975)))
write_csv(effects, file.path(out_dir, "summary_reminder_effect.csv"))

sub <- fit$summary(variables = "rd_grp", "mean", 
                   ~ quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  bind_cols(grp_levels %>% 
              select(factor, level)) %>%
  relocate(factor, level)
write_csv(sub, file.path(out_dir, "summary_subgroup_rd.csv"))


# Observed counterparts of the predictive-check statistics
# (kits with unknown return dates are left out of the conditional window shares)
obs_stats <- function(o) c(
  ret = mean(o %in% c(1, 2, 3, 5)),
  w1  = mean(o == 1),
  w2  = sum(o == 2) / sum(o %in% c(2, 3, 4)),
  w3  = sum(o == 3) / sum(o %in% c(3, 4)))
pre  <- obs_stats(d$outcome[d$post == 0])
post <- obs_stats(d$outcome[d$post == 1])
pat_obs <- d %>% 
  group_by(pat) %>% 
  summarise(n = n(), 
            r = sum(returned)) %>% 
  filter(n >= 2)
obs <- c(rep_ret_pre = pre[["ret"]], rep_ret_post = post[["ret"]],
         rep_w1_pre = pre[["w1"]], rep_w1_post = post[["w1"]],
         rep_w2_pre = pre[["w2"]], rep_w2_post = post[["w2"]],
         rep_w3_pre = pre[["w3"]], rep_w3_post = post[["w3"]],
         rep_all_multi = mean(pat_obs$r == pat_obs$n), rep_none_multi = mean(pat_obs$r == 0))

ppc <- fit$summary(variables = names(obs), "mean", ~ quantile(.x, probs = c(0.025, 0.975))) %>%
  mutate(observed = unname(obs[variable]))
write_csv(ppc, file.path(out_dir, "ppc_summary.csv"))


# Patient random effects: u_j1 (window 1) and u_j2 (windows 2-3), rebuilt from
# the non-centred draws: u_j = diag(sigma_u) * L_u * z_j
zd <- as_draws_matrix(fit$draws(c("sigma_u", "L_u", "z_u")))
J  <- stan_data$J
z1 <- zd[, paste0("z_u[1,", 1:J, "]")]
z2 <- zd[, paste0("z_u[2,", 1:J, "]")]
u1 <- z1 * as.numeric(zd[, "sigma_u[1]"])
u2 <- (z1 * as.numeric(zd[, "L_u[2,1]"]) + 
         z2 * as.numeric(zd[, "L_u[2,2]"])) * as.numeric(zd[, "sigma_u[2]"])
pat_rec <- d %>% 
  group_by(pat) %>%
  summarise(n_kits = n(), 
            n_returned = sum(returned), 
            n_w1 = sum(outcome == 1),
            n_pre = sum(post == 0), 
            n_post = sum(post == 1), 
            .groups = "drop")
re_sum <- function(m, w) tibble(pat = 1:J, window = w,
                                mean = colMeans(m),
                                q2.5 = apply(m, 2, quantile, 0.025),
                                q50  = apply(m, 2, quantile, 0.5),
                                q97.5 = apply(m, 2, quantile, 0.975))
pat_re <- bind_rows(re_sum(u1, "prompt"), re_sum(u2, "later")) %>% 
  left_join(pat_rec, by = "pat")
write_csv(pat_re, file.path(out_dir, "patient_effects.csv"))

capture.output(print(diag), file = file.path(out_dir, "diagnostics.txt"))
message("Done. Outputs in ", out_dir, if (run_prior) paste(" and", prior_dir) else "")
