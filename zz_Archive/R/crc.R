
## 0 - Configure =================================================================

library(dplyr)
library(tidyr)
library(bayesplot)
library(ggplot2)  
library(loo)
library(rstan)
library(posterior)
library(patchwork)
library(cmdstanr)
library(ggrepel)
library(scales)
library(purrr)
library(stringr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())





## 1 - Load Data =================================================================

df <- read.csv("data/crc.csv", stringsAsFactors = FALSE)


# Feature engineering
df_long <- df %>%
  # Center age and family size
  mutate(age_c     = patientage - mean(patientage),
         famsize_c = familysize - mean(familysize)) %>%
  pivot_longer(cols = c(n_sent_pre, n_sent_post),
               names_to = "period",
               values_to = "n_sent") %>%
  pivot_longer(cols = c(n_returned_pre, n_returned_post),
               names_to = "period2",
               values_to = "n_ret") %>%
  filter(sub(".*_(pre|post)$", "\\1", period) == 
           sub(".*_(pre|post)$", "\\1", period2)) %>%
  mutate(post    = if_else(period=="n_sent_post",  1L, 0L),
         sex01   = if_else(patientsex=="M", 1L, 0L),
         lang_id = as.integer(factor(lang)),
         emp_id  = as.integer(factor(employment)),
         rel_id  = as.integer(factor(relationship)),
         race_id = as.integer(factor(race)),
         rate = n_ret / n_sent,
         is_native_lang = ifelse(lang %in% c("ENG", "POL", "SPN"), 1, 0)) %>%
  filter(n_sent > 0) %>%
  as.data.frame()





## 2 - Random Intercept Checks ===================================================

# 1a) Boxplots of pre‐reminder return‐rates by group
df_long %>%
  filter(post == 0) %>%                            # pre-reminder only
  # coerce all grouping vars to character so pivot_longer can stack them
  mutate(across(c(is_native_lang, race, employment, relationship), as.character)) %>%
  pivot_longer(cols       = c(is_native_lang, race, employment, relationship),
               names_to   = "variable",
               values_to  = "group") %>%
  mutate(variable = factor(variable,
                           levels = c("is_native_lang", "race",
                                      "employment", "relationship")),
         group    = factor(group)) %>%
  ggplot(aes(x = group, y = rate)) +
  geom_boxplot() +
  coord_flip() +
  facet_wrap(~ variable, scales = "free_x", ncol = 2) +
  labs(title = "Pre-reminder Return Rates by Group",
       x     = NULL,
       y     = "Return rate") +
  theme_minimal(base_size = 13) +
  theme(strip.text       = element_text(face = "bold"),
        panel.grid.major = element_line(size = 0.1))

# If those boxes sit at very different vertical locations,
# that’s strong evidence you need random intercepts for that grp.

# 1b) Mean & 95% CI by group (pre only)
df_long %>%
  filter(post == 0) %>%
  group_by(is_native_lang) %>%
  summarise(mu = sum(n_ret) / sum(n_sent),
            se = sqrt(mu * (1 - mu) / sum(n_sent)),
            lo = mu - 1.96 * se,
            hi = mu + 1.96 * se) %>%
  ggplot(aes(x = is_native_lang, y = mu)) +
  geom_point() +
  geom_errorbar(aes(ymin = lo, ymax = hi)) +
  coord_flip() +
  labs(title = "Pre‐reminder mean ±95% CI by language",
       y = "return rate")





## 2 - Random Slope Checks =======================================================

# 2a) Pre vs Post “spaghetti” lines by group
group_summary_all <- df_long %>%
  # Make sure all grouping vars are the same type
  mutate(across(c(patientsex, is_native_lang, employment, relationship, race),
                as.character)) %>%
  # Pivot into long form
  pivot_longer(cols      = c(patientsex, is_native_lang, employment, 
                             relationship, race),
               names_to  = "variable",
               values_to = "group") %>%
  # Compute mean & SE of the return rate per (variable, group, post)
  group_by(variable, group, post) %>%
  summarise(n_sent_sum = sum(n_sent),
            n_ret_sum  = sum(n_ret),
            mu = sum(n_ret) / sum(n_sent),
            se = sqrt(mu * (1 - mu) / sum(n_sent)),
            .groups = "drop") %>%
  # Slopegraph-friendly x & labels
  mutate(x = ifelse(post == 0, 1, 2),
         period = factor(ifelse(post == 0, "Pre", "Post"), 
                         levels = c("Pre","Post")))


# Nice facet labels
pretty_var <- c(employment     = "Employment",
                is_native_lang = "Native language",
                patientsex     = "Sex",
                relationship   = "Relationship",
                race           = "Race")

# SLOPEGRAPH: facet-specific y scales 
okabe_ito <- c("#000000","#E69F00","#56B4E9","#009E73",
               "#F0E442","#0072B2","#D55E00","#CC79A7")


labs_post <- group_summary_all %>%
  mutate(x = ifelse(post == 0, 1, 2),
         period = factor(ifelse(post == 0, "Pre", "Post"), 
                         levels = c("Pre", "Post"))) %>%
  dplyr::filter(period == "Post")

gg_slope_freeY <- ggplot(
  group_summary_all %>%
    mutate(x = ifelse(post == 0, 1, 2),
           period = factor(ifelse(post == 0, "Pre", "Post"), 
                           levels = c("Pre", "Post"))),
  aes(x = x, y = mu, group = group, color = group)) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = mu - 1.96*se, ymax = mu + 1.96*se),
                width = 0.02, linewidth = 0.6) +
  ggrepel::geom_text_repel(data = labs_post,
                           aes(label = group),
                           nudge_x = 0.12, 
                           direction = "y", 
                           hjust = 0,
                           size = 3.4, box.padding = 0.2, segment.size = 0.3,
                           show.legend = FALSE, seed = 123) +
  scale_x_continuous(breaks = c(1,2), labels = c("Pre","Post"),
                     limits = c(0.95, 2.45), expand = c(0,0)) +
  # facet-specific y scales (the key change)
  facet_wrap(~ variable, ncol = 2, scales = "free_y",
             labeller = as_labeller(pretty_var)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     breaks = scales::breaks_pretty(n = 4),
                     expand = expansion(mult = c(0.05, 0.10))) +
  scale_color_manual(values = rep_len(okabe_ito, 
                                      length(unique(group_summary_all$group)))) +
  labs(
    title = "Return rate by group: Pre vs Post reminder",
    subtitle = "Points = mean return rate; error bars = 95% binomial SE",
    x = NULL, y = "Return rate"
  ) +
  theme_minimal(base_size = 15) +
  theme(plot.title       = element_text(face = "bold"),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3),
        legend.position  = "none",
        plot.margin      = margin(10, 30, 10, 10))

gg_slope_freeY
ggsave("pre_post_slopegraph.png", gg_slope_freeY, width = 12, height = 7, dpi = 300)



# 2b) DELTA FOREST: end caps + numeric labels 
fmt_pct <- scales::label_percent(accuracy = 0.1)

deltas <- group_summary_all %>%
  select(variable, group, post, mu, se) %>%
  mutate(period = ifelse(post == 0, "Pre", "Post")) %>%
  select(-post) %>%
  tidyr::pivot_wider(names_from = period, values_from = c(mu, se)) %>%
  mutate(delta    = mu_Post - mu_Pre,
         se_delta = sqrt(se_Post^2 + se_Pre^2),
         lo       = delta - 1.96 * se_delta,
         hi       = delta + 1.96 * se_delta) %>%
  group_by(variable) %>%
  arrange(delta, .by_group = TRUE) %>%
  mutate(group_ord = factor(group, levels = unique(group))) %>%
  ungroup()

gg_delta_caps <- ggplot(deltas,
                        aes(y = group_ord)) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             linewidth = 0.5, 
             color = "grey50") +
  # end-capped horizontal intervals
  geom_errorbarh(aes(xmin = lo, xmax = hi, x = delta, color = variable),
                 height = 0.2, linewidth = 0.7) +
  geom_point(aes(x = delta, color = variable), size = 2.3) +
  # numeric labels: lo (above), mean (to the right), hi (below)
  geom_text(aes(x = lo,    label = fmt_pct(lo)),  
            nudge_y =  0.175, size = 3.2, color = "grey30") +
  geom_text(aes(x = delta, label = fmt_pct(delta)), 
            nudge_y = -0.175, nudge_x = 0.0, size = 3.2, fontface = "bold") +
  geom_text(aes(x = hi,    label = fmt_pct(hi)),  
            nudge_y = 0.175, size = 3.2, color = "grey30") +
  facet_wrap(~ variable, scales = "free_y", ncol = 2,
             labeller = as_labeller(pretty_var)) +
  scale_x_continuous(labels = fmt_pct) +
  scale_color_manual(values = rep_len(okabe_ito, length(pretty_var))) +
  labs(
    title = "Change in return rate from reminder (Post − Pre)",
    subtitle = "Points = estimated change; bars = 95% interval with end caps",
    x = "Change in return rate", y = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(plot.title      = element_text(face = "bold"),
        strip.text      = element_text(face = "bold"),
        legend.position = "none")

gg_delta_caps
ggsave("delta_slopegraph.png", gg_delta_caps, width = 12, height = 7, dpi = 300)





## 3 - Interaction Checks ========================================================

# 3a) Does the reminder effect vary with age?
age_summary <- df_long %>%
  mutate(age_bin = cut(age_c, breaks=10)) %>%
  group_by(age_bin, post) %>%
  summarise(mu = sum(n_ret)/sum(n_sent),
            .groups = "drop")

ggplot(age_summary, 
       aes(x = age_bin, y = mu, color = factor(post), group = post)) +
  geom_line() +
  geom_point() +
  labs(title = "Return‐rate by Age Bin: pre vs post",
       x = "Age bin", y = "return rate", color = "Period") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# If the distance between the pre vs post curves changes across age,
# you might consider an age×post interaction (i.e. random slope over age).

# 3b) Same for family size
fams_summary <- df_long %>%
  group_by(family_size = familysize, post) %>%
  summarise(mu = sum(n_ret)/sum(n_sent), .groups="drop")

ggplot(fams_summary, 
       aes(x = family_size, y = mu, color = factor(post), group = post)) +
  geom_line() +
  geom_point() +
  labs(title = "Return‐rate by Family Size: pre vs post",
       x = "family size", y = "return rate", color = "Period")





## 4 - Simple Coefficient Estimates ==============================================

# Fit separate glm’s per language to get their own reminder slopes:
lang_coefs <- df_long %>%
  nest(data = -is_native_lang) %>%
  mutate(fit = map(data, ~ glm(cbind(n_ret, n_sent - n_ret) ~ post,
                               data = .x, 
                               family = binomial)),
         beta_post = map_dbl(fit, ~ coef(.x)["post"]),
         se_post   = map_dbl(fit, ~ sqrt(vcov(.x)["post","post"])))

ggplot(lang_coefs, aes(x = is_native_lang, y = beta_post)) +
  geom_point() +
  geom_errorbar(aes(ymin = beta_post - 1.96*se_post,
                    ymax = beta_post + 1.96*se_post),
                width = 0.2) +
  coord_flip() +
  labs(title = "Per‐language reminder slope (log‐odds scale)",
       y = "β_post estimate", x = "language")

# If these confidence intervals overlap heavily and are all similar,
# a common β_post is fine; if some languages have much larger slopes,
# consider a random slope for language.





## 5 - Stan Model ================================================================

# Stan data
stan_data <- df_long %>%
  mutate(concord   = as.integer(is_native_lang),          # MUST be 0/1 integer
         sex_idx  = as.integer(factor(patientsex)),     # 1..S
         rel_idx  = as.integer(factor(relationship)),   # 1..J_rel
         race_idx = as.integer(factor(race)),           # 1..J_race
         lang_idx = as.integer(factor(is_native_lang)), # 1..J_lang
         emp_idx  = as.integer(factor(employment))) %>% # 1..J_emp
  {
    list(N         = nrow(.),
         n_ret     = .$n_ret,
         n_sent    = .$n_sent,
         post      = .$post,
         concord   = .$concord,
         age_c     = .$age_c,
         famsize_c = .$famsize_c,
         sex       = .$sex_idx,
         lang      = .$lang_idx,
         emp       = .$emp_idx,
         rel       = .$rel_idx,
         race      = .$race_idx,
         S         = length(unique(.$sex_idx)),
         J_lang    = length(unique(.$lang_idx)),
         J_emp     = length(unique(.$emp_idx)),
         J_rel     = length(unique(.$rel_idx)),
         J_race    = length(unique(.$race_idx)))
  }


# add grainsize for reduce_sum
threads_per_chain <- 3L                 # 4 chains × 3 threads = 12 total
stan_data$grainsize <- max(1L, stan_data$N %/% (threads_per_chain * 2L))

# write + compile with threading
mod <- cmdstan_model(stan_file = 'stan/kit_returns.stan', 
                     cpp_options = list(stan_threads = TRUE))

# sample ~3000 total iterations (1500 warmup + 1500 sampling)
fit_cmd <- mod$sample(data = stan_data,
                      chains = 4,
                      parallel_chains = 4,
                      threads_per_chain = threads_per_chain,
                      iter_warmup = 1500,
                      iter_sampling = 1500,
                      seed = 123,
                      refresh = 150,    # progress every ~5%
                      max_treedepth = 10,
                      adapt_delta = 0.95)

# Diagnose and Summarize
fit_cmd$cmdstan_diagnose()
fit_cmd$summary(variables = c("alpha", "beta_concord_post", 
                              "beta_post", "beta_age", 
                              "beta_famsize", "beta_sex", "beta_rel",
                              "sigma_emp", "sigma_race",
                              "z_emp", "z_race",
                              "z_sex_slope","sigma_sex_slope", 
                              "z_rel_slope","sigma_rel_slope", 
                              "z_race_slope","sigma_race_slope"), 
                ~quantile(.x, probs = c(0.05, 0.5, 0.95))) %>% as.data.frame()


# Observed rates
obs_rate <- df_long$n_ret / df_long$n_sent

# Matrix of p (draws x N)
p_mat <- fit_cmd$draws(variables = "p", format = "draws_matrix")
p_sub <- as.matrix(p_mat)[1:200, , drop = FALSE]    # first 200 reps

# Dens overlay of rates
ppc_dens_overlay(y = obs_rate, yrep = p_sub) +
  ggtitle("Marginal PPC (return rate)")

# Matrix of replicated counts y_rep (draws x N)
yrep_mat <- fit_cmd$draws(variables = "y_rep", format = "draws_matrix")
yrep_sub <- as.matrix(yrep_mat)[1:200, , drop = FALSE]

# Bars: counts
ppc_bars(y = df_long$n_ret, yrep = yrep_sub) +
  ggtitle("Counts: observed vs posterior replicates")

# ECDF overlay for rates from replicated counts
rate_rep <- sweep(as.matrix(yrep_mat), 2, df_long$n_sent, FUN = "/")
ppc_ecdf_overlay(y = obs_rate, yrep = rate_rep[1:200, , drop = FALSE]) +
  ggtitle("ECDF: observed rate vs replicated")


ppc_bars_grouped(y = df_long$n_ret, 
                 yrep = yrep_sub, 
                 group = df_long$race) + ggtitle("PPC by Race")
ppc_bars_grouped(y = df_long$n_ret, 
                 yrep = yrep_sub, 
                 group = df_long$patientsex) + ggtitle("PPC by Sex")
ppc_bars_grouped(y = df_long$n_ret, 
                 yrep = yrep_sub, 
                 group = df_long$relationship) + ggtitle("PPC by Relationship")
ppc_bars_grouped(y = df_long$n_ret, 
                 yrep = yrep_sub, 
                 group = df_long$post) + ggtitle("PPC by Period (pre/post)")


library(ggplot2)

# helper to extract matrices (draws x K) for vector params
grab_mat <- function(prefix) {
  dm <- fit_cmd$draws(variables = prefix, format = "draws_matrix")
  as.matrix(dm)  # draws x K
}

beta_post_vec <- as.numeric(fit_cmd$draws("beta_post", format = "draws_matrix")[,1])

b_sex_mat  <- grab_mat("b_sex")    # draws x S
b_rel_mat  <- grab_mat("b_rel")    # draws x J_rel
b_race_mat <- grab_mat("b_race")   # draws x J_race
b_lang_mat <- grab_mat("beta_concord_post")   # draws x J_lang

# add the global beta_post to each column (per draw)
sex_samps  <- sweep(b_sex_mat, 1, beta_post_vec, `+`)
rel_samps  <- sweep(b_rel_mat, 1, beta_post_vec, `+`)
race_samps <- sweep(b_race_mat, 1, beta_post_vec, `+`)
lang_samps <- sweep(b_lang_mat, 1, beta_post_vec, `+`)

# labels from data
sex_levels  <- levels(factor(df_long$patientsex))
rel_levels  <- levels(factor(df_long$relationship))
race_levels <- levels(factor(df_long$race))
lang_levels <- levels(factor(df_long$is_native_lang))

summarize_effects <- function(mat, levels_vec) {
  data.frame(
    level = levels_vec,
    mean  = colMeans(mat),
    low   = apply(mat, 2, quantile, 0.05),
    high  = apply(mat, 2, quantile, 0.95)
  )
}

sex_df  <- summarize_effects(sex_samps,  sex_levels)
rel_df  <- summarize_effects(rel_samps,  rel_levels)
race_df <- summarize_effects(race_samps, race_levels)
lang_df <- summarize_effects(lang_samps, lang_levels)

sex_plot <- ggplot(sex_df, aes(x = level, y = mean)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  labs(title = "Sex-specific reminder effects (log-odds)",
       x = "Sex", y = expression(beta[post] + b[sex])) +
  theme_minimal()

rel_plot <- ggplot(rel_df, aes(x = level, y = mean)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  labs(title = "Relationship-specific reminder effects (log-odds)",
       x = "Relationship", y = expression(beta[post] + b[rel])) +
  theme_minimal()

race_plot <- ggplot(race_df, aes(x = level, y = mean)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  labs(title = "Race-specific reminder effects (log-odds)",
       x = "Race", y = expression(beta[post] + b[race])) +
  coord_flip() +
  theme_minimal()

lang_plot <- ggplot(lang_df, aes(x = level, y = mean)) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  labs(title = "Language-specific reminder effects (log-odds)",
       x = "Language", y = expression(beta[post] + b[lang])) +
  theme_minimal()

# arrange (if you use patchwork)
# library(patchwork)
sex_plot + rel_plot + race_plot + lang_plot + plot_layout(ncol = 2)

