## 05 - Sensitivity: main model vs. model without the Q4 (Oct-Dec) term ======
##
## Input : 03_output/models/stage_model/   (main fit)
##         03_output/models/sens_no_q4/    (refit without q4: 03_fit_stage_model.R no_q4)
## Output: 03_output/models/sens_no_q4/comparison.csv
##         03_output/figures/sens_no_q4_comparison.png
##
## Run from the project root (the .Rproj directory).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggdist)
  library(cmdstanr)
  library(posterior)
})
options(cmdstanr_no_ver_check = TRUE)

models <- c(`With Q4 (main)` = "03_output/models/stage_model",
            `Without Q4`     = "03_output/models/sens_no_q4")

# Quantities compared, built from each model's own draws. Covariate effects are
# looked up by column name, because removing q4 shifts the beta[] indices.
quantities <- function(dir) {
  labs  <- read_rds(file.path(dir, "design_labels.rds"))
  csvs  <- list.files(dir, pattern = "^stage_model-\\d+\\.csv$", full.names = TRUE)
  dr    <- as_draws_df(read_cmdstan_csv(csvs, variables = c("beta", "beta_post", "beta_remind",
                                                            "beta_carry", "gamma_age", "sigma_u",
                                                            "rd_post", "p_on_post", "p_off_post", "rd_grp"),
                                        sampler_diagnostics = "")$post_warmup_draws)
  b     <- function(col) dr[[paste0("beta[", match(col, labs$X_cols), "]")]]
  g     <- function(f, l) dr[[paste0("rd_grp[", labs$grp_levels$id[labs$grp_levels$factor == f &
                                                                    labs$grp_levels$level == l], "]")]]
  tibble(
    `Reminder effect (points)`                = 100 * dr$rd_post,
    `Return rate, reminders on (%)`           = 100 * dr$p_on_post,
    `Return rate, reminders off (%)`          = 100 * dr$p_off_post,
    `Reminder odds ratio, days 15–90`         = exp(dr$beta_remind),
    `Carry-over odds ratio, after day 90`     = exp(dr$beta_carry),
    `Other post-2022 change (odds ratio)`     = exp(dr$beta_post),
    `COVID period (odds ratio)`               = exp(b("covid")),
    `Calendar trend per year (odds ratio)`    = exp(b("time_yrs")),
    `Effect, declined texts (points)`         = 100 * g("text_consent", "No"),
    `Effect, consented (points)`              = 100 * g("text_consent", "Yes"),
    `Effect, women minus men (points)`        = 100 * (g("sex", "Female") - g("sex", "Male")),
    `Age moderation (log-odds per decade)`    = dr$gamma_age,
    `Patient SD, later returns`               = dr$`sigma_u[2]`
  )
}

draws <- bind_rows(lapply(names(models), function(m)
  quantities(models[[m]]) %>% mutate(model = m))) %>%
  pivot_longer(-model, names_to = "quantity") %>%
  mutate(quantity = factor(quantity, levels = unique(quantity)),
         model    = factor(model, levels = names(models)))

comparison <- draws %>%
  group_by(quantity, model) %>%
  summarise(median = median(value), lo = quantile(value, 0.025), hi = quantile(value, 0.975),
            .groups = "drop") %>%
  mutate(cell = sprintf("%.2f, 95%% CrI [%.2f, %.2f]", median, lo, hi)) %>%
  select(quantity, model, cell) %>%
  pivot_wider(names_from = model, values_from = cell)
write_csv(comparison, file.path(models[["Without Q4"]], "comparison.csv"))
print(comparison, n = Inf, width = 200)

# Q4's own estimate in the main model, for context
main_labs <- read_rds(file.path(models[[1]], "design_labels.rds"))
q4 <- as_draws_df(read_cmdstan_csv(list.files(models[[1]], "^stage_model-\\d+\\.csv$", full.names = TRUE),
                                   variables = "beta", sampler_diagnostics = "")$post_warmup_draws)[[
  paste0("beta[", match("q4", main_labs$X_cols), "]")]]
message(sprintf("Q4 odds ratio in main model: OR = %.2f, 95%% CrI [%.2f, %.2f]",
                median(exp(q4)), quantile(exp(q4), 0.025), quantile(exp(q4), 0.975)))

# Side-by-side intervals; each quantity on its own scale
p <- ggplot(draws, aes(x = value, y = model, colour = model)) +
  stat_pointinterval(.width = c(0.5, 0.95), point_size = 1.8) +
  facet_wrap(~ quantity, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = c("#2a78d6", "#eb6834")) +
  labs(title = "Sensitivity: Removing the Oct–Dec (Q4) Term",
       subtitle = "Posterior median, 50% and 95% credible intervals",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5), plot.title.position = "plot",
        plot.subtitle = element_text(hjust = 0.5),
        strip.text = element_text(face = "bold", hjust = 0))
ggsave("03_output/figures/sens_no_q4_comparison.png", p, width = 10, height = 7.5, dpi = 150)
