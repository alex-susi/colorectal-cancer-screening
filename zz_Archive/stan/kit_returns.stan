// ---------- model_kit_returns_concord_fixed.stan ----------
//
// Bayesian hierarchical Binomial-logit model for FIT kit returns,
// revised to estimate a *fixed* "same-language reminder" effect via
// a post × concordance interaction.
//
// DATA CONTEXT (from the report):
// - Outcome: FIT kit return counts.
// - QI intervention: biweekly text reminders implemented in the post period (2022–2023).
// - Goal: estimate the overall post-intervention reminder effect AND quantify whether
//         reminders being in the patient's preferred language ("same-language") changes
//         the post-period response, while controlling for demographics.
//
// WHAT CHANGED RELATIVE TO THE PRIOR MODEL:
// - Previously: "language" entered as a random slope on post: b_lang[lang[i]] * post[i].
// - Now: language concordance is represented directly as a binary indicator concord[i] ∈ {0,1}
//        (1 = message language matches patient's preferred language as defined in your data).
//        The model includes a single fixed interaction term:
//            beta_concord_post * post[i] * concord[i]
//   Interpretation:
//     - beta_concord_post is the additional log-odds boost in the post period when reminders
//       are same-language (concord=1), after adjusting for covariates and other random effects.
//     - exp(beta_concord_post) is the *odds ratio* for same-language vs not, specifically for
//       the reminder effect in the post period.
//
// MODEL SUMMARY:
//   n_ret[i] ~ Binomial(n_sent[i], p[i])
//   logit(p[i]) = alpha
//    + u_race[race[i]] + u_emp[emp[i]]              // baseline differences
//    + X_i * beta                                   // fixed demographics
//    + post[i] * (beta_post + subgroup deviations)  // overall QI effect + heterogeneity
//    + post[i] * concord[i] * beta_concord_post     // NEW: same-language increment
//
// Where subgroup deviations (random slopes on post) are included for:
//   sex, relationship status, race
// and baseline random intercepts are included for:
//   race, employment status
//
// COMPUTATION:
// - Uses reduce_sum for within-chain parallelization; grainsize tunes chunk size.
//
// POSTERIOR PREDICTIVE:
// - Generates p[i] and replicated counts y_rep[i] for PPCs.
//

functions {
  // binom_slice_lpmf:
  //   Sliced log-likelihood used by reduce_sum.
  //   Computes the binomial-logit log PMF for observations i in [start, end]
  //   and returns the summed log-likelihood for that slice.
  //
  // NOTE: The first argument `slice_n` is required by reduce_sum's signature
  // but is unused (common Stan pattern).
  real binom_slice_lpmf(
      array[] int  slice_n, int start, int end,
      array[] int  n_ret,
      array[] int  n_sent,
      array[] int  post,
      array[] int  concord,   // NEW: 0/1 same-language indicator
      vector       age_c,
      vector       famsize_c,
      array[] int  sex,
      array[] int  emp,
      array[] int  rel,
      array[] int  race,
      // random intercepts:
      // u_emp[j]  = baseline log-odds deviation for employment group j
      // u_race[k] = baseline log-odds deviation for race group k
      vector u_emp,
      vector u_race,
      // random slopes for post:
      // subgroup-specific deviations to the post-period reminder effect
      // (language-specific slope removed in this revision)
      vector b_sex,
      vector b_rel,
      vector b_race,
      // fixed parameters (population-level effects)
      real alpha,
      real beta_post,
      real beta_concord_post,  // NEW: fixed post × concord effect
      real beta_age,
      real beta_famsize,
      real beta_sex,
      real beta_rel
  ) {
    real lp = 0;

    for (i in start:end) {

      // Linear predictor for observation i:
      //
      // alpha = baseline log-odds for an "average" patient in reference groups,
      //         because continuous covariates are centered (age_c, famsize_c).
      //
      // Baseline subgroup differences (random intercepts):
      //   + u_race[race[i]]  : race-specific baseline deviation
      //   + u_emp[emp[i]]    : employment-specific baseline deviation
      //
      // Fixed (non-interaction) covariates:
      //   + beta_age     * age_c[i]
      //   + beta_famsize * famsize_c[i]
      //
      // Sex and relationship are coded here as binary indicators:
      //   (sex[i] == 2) means "sex level 2 vs everyone else"
      //   (rel[i] == 2) means "relationship level 2 vs everyone else"
      // With data: S=2 and J_rel=2, effectively a standard binary coding.
      //
      // Post/intervention effect:
      //   post[i] * (beta_post + subgroup deviations)
      // where subgroup deviations allow the reminder effect to differ by:
      //   sex, relationship, race
      //
      // NEW: Same-language effect in the post period:
      //   + beta_concord_post * post[i] * concord[i]
      //
      // Interpretation of beta_concord_post:
      //   - When post=0 (pre period), this term is 0 regardless of concord.
      //   - When post=1 (post period):
      //       concord=0 => no additional boost beyond beta_post + subgroup deviations
      //       concord=1 => adds beta_concord_post to the post-period log-odds
      //   Therefore, exp(beta_concord_post) is the odds ratio for same-language vs not
      //   *among post-period reminders*, controlling for other covariates and random effects.
      real nu = alpha
              + u_race[race[i]]
              + u_emp[emp[i]]
              + beta_age     * age_c[i]
              + beta_famsize * famsize_c[i]
              + beta_sex     * (sex[i] == 2)
              + beta_rel     * (rel[i] == 2)
              + (beta_post
                 + b_sex[sex[i]]
                 + b_rel[rel[i]]
                 + b_race[race[i]]) * post[i]
              + beta_concord_post * post[i] * concord[i];

      // Binomial-logit likelihood:
      //   n_ret[i] successes out of n_sent[i] trials with logit(p)=nu.
      lp += binomial_logit_lpmf(n_ret[i] | n_sent[i], nu);
    }
    return lp;
  }
}

data {
  // Number of observations.
  int<lower=1> N;

  // Outcome:
  // - n_sent[i] is the number of kits distributed (trials).
  // - n_ret[i]  is the number returned (successes).
  // This supports both:
  //   (a) patient-level rows: n_sent[i]=1, n_ret[i]∈{0,1}
  //   (b) aggregated rows:    n_sent[i]>1, n_ret[i] can be 0..n_sent[i]
  array[N] int<lower=0> n_ret;
  array[N] int<lower=1> n_sent;

  // Intervention period indicator:
  // post=0: pre-intervention (2019–2021) - no reminder program
  // post=1: post-intervention (2022–2023) - biweekly text reminders
  array[N] int<lower=0, upper=1> post;

  // NEW: same-language indicator for post-period reminders
  // concord=1: message language matches patient's preferred language (per your data definition)
  // concord=0: message not in preferred language (or otherwise not concordant)
  //
  // NOTE: This is modeled as only affecting the *post* period via an interaction,
  // because only post-period includes reminders. If concord is defined for all rows,
  // it still has no effect in pre because it's multiplied by post.
  array[N] int<lower=0, upper=1> concord;

  // Centered continuous covariates (centered in R/Python before passing to Stan):
  // age_c[i]     = age - mean(age)
  // famsize_c[i] = family_size - mean(family_size)
  vector[N] age_c;
  vector[N] famsize_c;

  // Grouping indices (all 1-based, as required in Stan):
  array[N] int<lower=1> sex;     // 1..S
  array[N] int<lower=1> emp;     // 1..J_emp
  array[N] int<lower=1> rel;     // 1..J_rel
  array[N] int<lower=1> race;    // 1..J_race

  // Number of levels for each grouping factor:
  int<lower=1> S;
  int<lower=1> J_emp;
  int<lower=1> J_rel;
  int<lower=1> J_race;

  // reduce_sum tuning:
  // grainsize controls how many observations per slice when threading.
  // Larger grainsize => fewer slices (less overhead), but less load balancing.
  int<lower=1> grainsize;
}


parameters {
  // -----------------------
  // FIXED EFFECTS (population-level)
  // -----------------------
  real alpha;               // baseline log-odds for average patient in reference groups
  real beta_post;           // overall (average) post-period reminder effect on log-odds
  real beta_concord_post;   // NEW: additional post effect if reminders are same-language
  real beta_age;            // effect per 1 unit increase in centered age (often 1 year)
  real beta_famsize;        // effect per 1 unit increase in centered family size
  real beta_sex;            // main effect for (sex == 2) vs baseline
  real beta_rel;            // main effect for (rel == 2) vs baseline

  // -----------------------
  // RANDOM INTERCEPTS (baseline subgroup differences)
  // Non-centered parameterization: u = z * sigma
  // -----------------------
  vector[J_emp]  z_emp;         // standard normal latent for employment intercepts
  real<lower=0>  sigma_emp;     // SD of employment intercepts (half-normal prior)

  vector[J_race] z_race;        // standard normal latent for race intercepts
  real<lower=0>  sigma_race;    // SD of race intercepts (half-normal prior)

  // -----------------------
  // RANDOM SLOPES on post (heterogeneous reminder effects)
  // Non-centered parameterization: b = z * sigma
  //
  // In this revision, language heterogeneity is not random; it's captured by the
  // fixed interaction beta_concord_post * post * concord.
  // -----------------------
  vector[S]      z_sex_slope;         // sex-specific post deviations
  real<lower=0>  sigma_sex_slope;     // SD of sex post deviations

  vector[J_rel]  z_rel_slope;         // relationship-specific post deviations
  real<lower=0>  sigma_rel_slope;     // SD of relationship post deviations

  vector[J_race] z_race_slope;        // race-specific post deviations
  real<lower=0>  sigma_race_slope;    // SD of race post deviations
}


transformed parameters {
  // Convert non-centered latent z's into actual random effects.

  // Random intercepts:
  // u_emp[j]  = baseline deviation for employment group j
  // u_race[k] = baseline deviation for race group k
  vector[J_emp]  u_emp   = z_emp  * sigma_emp;
  vector[J_race] u_race  = z_race * sigma_race;

  // Random slopes on post:
  // These are subgroup-specific deviations from the overall beta_post.
  // The post effect for observation i is:
  //   beta_post + b_sex[sex[i]] + b_rel[rel[i]] + b_race[race[i]]
  vector[S]      b_sex   = z_sex_slope  * sigma_sex_slope;
  vector[J_rel]  b_rel   = z_rel_slope  * sigma_rel_slope;
  vector[J_race] b_race  = z_race_slope * sigma_race_slope;
}

model {
  // -----------------------
  // PRIORS (weakly informative)
  // -----------------------

  // Fixed effects priors:
  // alpha is wider because it represents baseline log-odds.
  alpha             ~ normal(0, 5);

  // Regularizing priors on log-odds scale.
  beta_post         ~ normal(0, 1);
  beta_concord_post ~ normal(0, 1);  // NEW
  beta_age          ~ normal(0, 1);
  beta_famsize      ~ normal(0, 1);
  beta_sex          ~ normal(0, 1);
  beta_rel          ~ normal(0, 1);

  // SD priors: with <lower=0>, normal(0,1) implies half-normal(0,1).
  // These control shrinkage of group deviations toward 0.
  sigma_emp        ~ normal(0, 1);
  sigma_race       ~ normal(0, 1);
  sigma_sex_slope  ~ normal(0, 1);
  sigma_rel_slope  ~ normal(0, 1);
  sigma_race_slope ~ normal(0, 1);

  // Latent z's are standard normal in non-centered parameterization.
  z_emp        ~ normal(0, 1);
  z_race       ~ normal(0, 1);
  z_sex_slope  ~ normal(0, 1);
  z_rel_slope  ~ normal(0, 1);
  z_race_slope ~ normal(0, 1);

  // -----------------------
  // LIKELIHOOD (parallelized)
  // -----------------------
  //
  // reduce_sum splits the likelihood into slices and sums them in parallel.
  // rep_array(1, N) is a dummy array needed as the "sliced" argument.
  // The actual data vectors/arrays are passed as shared arguments.
  target += reduce_sum(
    binom_slice_lpmf,
    rep_array(1, N), grainsize,
    n_ret, n_sent, post, concord, age_c, famsize_c,
    sex, emp, rel, race,
    u_emp, u_race, b_sex, b_rel, b_race,
    alpha, beta_post, beta_concord_post, beta_age, beta_famsize, beta_sex, beta_rel
  );
}


generated quantities {
  // Posterior predictive / fitted quantities:
  //
  // p[i]     = fitted return probability for observation i
  // y_rep[i] = posterior predictive replicated returned counts,
  //            drawn from Binomial(n_sent[i], p[i])
  //
  // These support model checks like:
  // - predicted vs observed return rates
  // - subgroup PPCs
  // - calibration across bins
  vector[N] p;
  array[N] int y_rep;

  for (i in 1:N) {
    // Recompute linear predictor for each observation (same as in likelihood).
    real nu = alpha
            + u_race[race[i]]
            + u_emp[emp[i]]
            + beta_age     * age_c[i]
            + beta_famsize * famsize_c[i]
            + beta_sex     * (sex[i] == 2)
            + beta_rel     * (rel[i] == 2)
            + (beta_post
               + b_sex[sex[i]]
               + b_rel[rel[i]]
               + b_race[race[i]]) * post[i]
            + beta_concord_post * post[i] * concord[i];

    // Convert logit to probability.
    p[i] = inv_logit(nu);

    // Simulate replicated outcomes for PPC.
    y_rep[i] = binomial_rng(n_sent[i], p[i]);
  }
}
