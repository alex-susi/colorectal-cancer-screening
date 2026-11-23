// =============================================================================
// kit_returns_stage.stan
//
// Three-window (discrete-time) hierarchical logistic model for FIT kit returns.
//
// Reminder schedule (post period, 2022+): a text every 2 weeks from day 14 to
// day 90 after the kit was sent, stopping early if the kit is returned. Each kit
// therefore passes through three windows:
//   Window 1 (days 0-14):  returned by day 14?                 no reminders yet
//   Window 2 (days 15-90): returned, given still out at day 14? reminders active
//   Window 3 (day 91+):    returned, given still out at day 90? reminders stopped
//
//   P(kit returned) = p1 + (1 - p1) p2 + (1 - p1)(1 - p2) p3
//
// The reminder effect R_i enters window 2 of post-period kits. Window 3 gets a
// separate post-period "carry-over" term C_i (late responses to the last texts,
// or a changed pool of kits still out at day 90); like R_i, it has its own
// unshrunk terms for patients who declined texts. Changes after 2022 that are not
// reminders (beta_post) act on all windows and are identified mainly from
// window 1. Identifying assumption: absent reminders, post-2022 changes would
// have shifted returns in every window equally on the log-odds scale.
//
// Patients get correlated random intercepts: one for prompt returns (window 1)
// and one for later returns (windows 2 and 3 share it).
//
// Reminder-effect heterogeneity: sum-to-zero deviations by sex, race/ethnicity,
// language, employment and relationship, all sharing one learned shrinkage
// scale tau; linear moderation by age and family size; unshrunk terms for text
// consent (non-consenting patients presumably got no texts -> placebo check).
//
// prior_only = 1 skips the likelihood, giving draws from the prior (for prior
// predictive checks).
// =============================================================================

functions {
  // Linear predictors for every kit. Columns 1-3 = window 1-3 logits under the
  // regime the kit actually had; column 4 = that kit's reminder effect R_i
  // (window 2); column 5 = its carry-over effect C_i (window 3).
  matrix linear_predictors(vector alpha, vector beta, real beta_post,
                           real beta_covid_late, real beta_remind, real beta_carry,
                           vector d_sex, vector d_race, vector d_lang,
                           vector d_emp, vector d_rel,
                           real gamma_age, real gamma_fam, vector gamma_con,
                           vector gamma_carry_con,
                           matrix u,
                           matrix X, vector post, vector covid, array[] int pat,
                           array[] int sex, array[] int race, array[] int lang,
                           array[] int emp, array[] int rel,
                           vector age_c10, vector family_size_c, matrix C_con) {
    int N = rows(X);
    matrix[N, 5] out;
    vector[N] reff = beta_remind
                     + d_sex[sex] + d_race[race] + d_lang[lang]
                     + d_emp[emp] + d_rel[rel]
                     + gamma_age * age_c10 + gamma_fam * family_size_c
                     + C_con * gamma_con;
    vector[N] carry  = beta_carry + C_con * gamma_carry_con;
    vector[N] shared = X * beta + beta_post * post;
    vector[N] late   = shared + beta_covid_late * covid + u[pat, 2];
    out[, 1] = alpha[1] + shared + u[pat, 1];
    out[, 2] = alpha[2] + late + post .* reff;
    out[, 3] = alpha[3] + late + post .* carry;
    out[, 4] = reff;
    out[, 5] = carry;
    return out;
  }

  // Log-likelihood of one kit's observed outcome
  real kit_ll(int outcome, real e1, real e2, real e3) {
    if (outcome == 1) return log_inv_logit(e1);
    if (outcome == 2) return log1m_inv_logit(e1) + log_inv_logit(e2);
    if (outcome == 3) return log1m_inv_logit(e1) + log1m_inv_logit(e2) + log_inv_logit(e3);
    if (outcome == 4) return log1m_inv_logit(e1) + log1m_inv_logit(e2) + log1m_inv_logit(e3);
    // returned, timing unknown: marginalise over the three windows
    return log_sum_exp({log_inv_logit(e1),
                        log1m_inv_logit(e1) + log_inv_logit(e2),
                        log1m_inv_logit(e1) + log1m_inv_logit(e2) + log_inv_logit(e3)});
  }

  // P(returned in any window)
  real p_return(real e1, real e2, real e3) {
    real p1 = inv_logit(e1);
    real p2 = inv_logit(e2);
    return p1 + (1 - p1) * (p2 + (1 - p2) * inv_logit(e3));
  }
}

data {
  int<lower=1> N;                              // kits
  int<lower=1> J;                              // patients
  array[N] int<lower=1, upper=J> pat;          // patient index
  array[N] int<lower=1, upper=5> outcome;      // 1 = returned day 0-14
                                               // 2 = returned day 15-90
                                               // 3 = returned after day 90
                                               // 4 = not returned
                                               // 5 = returned, date unknown

  vector<lower=0, upper=1>[N] post;            // dispensed on/after 2022-01-01
  vector<lower=0, upper=1>[N] covid;           // dispensed in COVID disruption

  // Covariates shared by all windows (demographics, age, time trend, season)
  int<lower=1> K;
  matrix[N, K] X;
  vector<lower=0>[K] beta_sd;                  // prior SD per column of X

  // Reminder-effect moderators
  int<lower=2> n_sex;  array[N] int<lower=1, upper=n_sex>  sex;
  int<lower=2> n_race; array[N] int<lower=1, upper=n_race> race;
  int<lower=2> n_lang; array[N] int<lower=1, upper=n_lang> lang;
  int<lower=2> n_emp;  array[N] int<lower=1, upper=n_emp>  emp;
  int<lower=2> n_rel;  array[N] int<lower=1, upper=n_rel>  rel;
  vector[N] age_c10;
  vector[N] family_size_c;
  int<lower=0> K_con;
  matrix[N, K_con] C_con;                      // text-consent dummies

  // Reporting groups for subgroup effects (ids unique across all factors)
  int<lower=1> n_fac;
  int<lower=1> n_grp;
  array[N, n_fac] int<lower=1, upper=n_grp> grp;

  int<lower=0, upper=1> save_log_lik;
  int<lower=0, upper=1> prior_only;
}

transformed data {
  // Index sets by outcome type, for a vectorised likelihood
  array[5] int n_out = rep_array(0, 5);
  for (i in 1:N) n_out[outcome[i]] += 1;
  array[n_out[1]] int idx1;
  array[n_out[2]] int idx2;
  array[n_out[3]] int idx3;
  array[n_out[4]] int idx4;
  array[n_out[5]] int idx5;
  {
    array[5] int c = rep_array(1, 5);
    for (i in 1:N) {
      int o = outcome[i];
      if (o == 1) idx1[c[1]] = i;
      else if (o == 2) idx2[c[2]] = i;
      else if (o == 3) idx3[c[3]] = i;
      else if (o == 4) idx4[c[4]] = i;
      else idx5[c[5]] = i;
      c[o] += 1;
    }
  }

  // Post-period kits per reporting group (denominators for subgroup effects)
  real n_post = sum(post);
  vector[n_grp] n_post_grp = rep_vector(0, n_grp);
  for (i in 1:N) {
    if (post[i] == 1) {
      for (f in 1:n_fac) n_post_grp[grp[i, f]] += 1;
    }
  }

  // Kits per patient, for the all-or-none posterior predictive check
  array[J] int n_kits_pat = rep_array(0, J);
  for (i in 1:N) n_kits_pat[pat[i]] += 1;
  int n_multi = 0;
  for (j in 1:J) n_multi += n_kits_pat[j] >= 2;

  // Sum-to-zero vectors: this scale gives each element a marginal SD of 1
  real s_sex  = sqrt(n_sex  / (n_sex  - 1.0));
  real s_race = sqrt(n_race / (n_race - 1.0));
  real s_lang = sqrt(n_lang / (n_lang - 1.0));
  real s_emp  = sqrt(n_emp  / (n_emp  - 1.0));
  real s_rel  = sqrt(n_rel  / (n_rel  - 1.0));
}

parameters {
  vector[3] alpha;                    // window 1 / 2 / 3 baseline log-odds
  vector[K] beta;                     // shared covariate effects
  real beta_post;                     // non-reminder change from 2022 (all windows)
  real beta_covid_late;               // extra COVID effect on later returns (windows 2-3)
  real beta_remind;                   // REMINDER EFFECT (window 2, post), log-odds
  real beta_carry;                    // post-period change in window 3 (after texts stop), consenters

  // Heterogeneity in the reminder effect (non-centred, sum-to-zero)
  real<lower=0> tau;                  // shared shrinkage scale for all factors
  sum_to_zero_vector[n_sex]  z_sex;
  sum_to_zero_vector[n_race] z_race;
  sum_to_zero_vector[n_lang] z_lang;
  sum_to_zero_vector[n_emp]  z_emp;
  sum_to_zero_vector[n_rel]  z_rel;
  real gamma_age;                     // reminder effect change per decade of age
  real gamma_fam;                     // ... per additional household member
  vector[K_con] gamma_con;            // ... for consent = No / Unknown (unshrunk)
  vector[K_con] gamma_carry_con;      // carry-over difference for consent = No / Unknown

  // Patient random intercepts (prompt, later), correlated (non-centred)
  vector<lower=0>[2] sigma_u;
  cholesky_factor_corr[2] L_u;
  matrix[2, J] z_u;
}

transformed parameters {
  vector[n_sex]  d_sex  = tau * z_sex;
  vector[n_race] d_race = tau * z_race;
  vector[n_lang] d_lang = tau * z_lang;
  vector[n_emp]  d_emp  = tau * z_emp;
  vector[n_rel]  d_rel  = tau * z_rel;
}

model {
  // Priors (log-odds scale)
  alpha           ~ normal(0, 1.5);
  beta            ~ normal(0, beta_sd);
  beta_post       ~ normal(0, 1);
  beta_covid_late ~ normal(0, 1);
  beta_remind     ~ normal(0, 1);
  beta_carry      ~ normal(0, 1);

  tau    ~ normal(0, 0.5);            // half-normal: most subgroup deviations small
  z_sex  ~ normal(0, s_sex);
  z_race ~ normal(0, s_race);
  z_lang ~ normal(0, s_lang);
  z_emp  ~ normal(0, s_emp);
  z_rel  ~ normal(0, s_rel);
  gamma_age ~ normal(0, 0.5);
  gamma_fam ~ normal(0, 0.5);
  gamma_con ~ normal(0, 1);
  gamma_carry_con ~ normal(0, 1);

  sigma_u ~ normal(0, 2);             // half-normal
  L_u     ~ lkj_corr_cholesky(2);
  to_vector(z_u) ~ std_normal();

  // Likelihood
  if (!prior_only) {
    matrix[J, 2] u = (diag_pre_multiply(sigma_u, L_u) * z_u)';
    matrix[N, 5] lp = linear_predictors(alpha, beta, beta_post, beta_covid_late,
                                        beta_remind, beta_carry, d_sex, d_race,
                                        d_lang, d_emp, d_rel, gamma_age, gamma_fam,
                                        gamma_con, gamma_carry_con, u, X, post, covid, pat,
                                        sex, race, lang, emp, rel,
                                        age_c10, family_size_c, C_con);
    vector[N] e1 = lp[, 1];
    vector[N] e2 = lp[, 2];
    vector[N] e3 = lp[, 3];

    target += log_inv_logit(e1[idx1]);
    target += log1m_inv_logit(e1[idx2]) + log_inv_logit(e2[idx2]);
    target += log1m_inv_logit(e1[idx3]) + log1m_inv_logit(e2[idx3]) + log_inv_logit(e3[idx3]);
    target += log1m_inv_logit(e1[idx4]) + log1m_inv_logit(e2[idx4]) + log1m_inv_logit(e3[idx4]);
    for (k in 1:n_out[5]) target += kit_ll(5, e1[idx5[k]], e2[idx5[k]], e3[idx5[k]]);
  }
}

generated quantities {
  real or_remind = exp(beta_remind);  // conditional odds ratio, window 2
  real or_carry  = exp(beta_carry);
  real or_post   = exp(beta_post);
  real rho_u     = multiply_lower_tri_self_transpose(L_u)[1, 2];

  // Reminder program effect on the probability scale (g-computation over post
  // kits): each post-period kit's P(returned) with the reminder program vs.
  // with it switched off (window-2 effect and window-3 carry-over removed)
  real p_on_post  = 0;
  real p_off_post = 0;
  real rd_post    = 0;
  vector[n_grp] rd_grp = rep_vector(0, n_grp);

  // Posterior (or prior) predictive checks
  real rep_ret_pre;          // share of kits returned
  real rep_ret_post;
  real rep_w1_pre;           // share returned by day 14
  real rep_w1_post;
  real rep_w2_pre;           // share returned day 15-90 | still out at day 14
  real rep_w2_post;
  real rep_w3_pre;           // share returned after day 90 | still out at day 90
  real rep_w3_post;
  real rep_all_multi;        // share of multi-kit patients returning every kit
  real rep_none_multi;       // ... returning no kit

  vector[save_log_lik ? N : 0] log_lik;

  {
    matrix[J, 2] u = (diag_pre_multiply(sigma_u, L_u) * z_u)';
    matrix[N, 5] lp = linear_predictors(alpha, beta, beta_post, beta_covid_late,
                                        beta_remind, beta_carry, d_sex, d_race,
                                        d_lang, d_emp, d_rel, gamma_age, gamma_fam,
                                        gamma_con, gamma_carry_con, u, X, post, covid, pat,
                                        sex, race, lang, emp, rel,
                                        age_c10, family_size_c, C_con);
    array[J] int ret_pat = rep_array(0, J);
    // [period, statistic]: kits, returned, w1, at risk w2, w2, at risk w3, w3
    array[2, 7] real s = rep_array(0.0, 2, 7);

    for (i in 1:N) {
      int per = post[i] == 1 ? 2 : 1;

      if (post[i] == 1) {
        real e2_off = lp[i, 2] - lp[i, 4];
        real e3_off = lp[i, 3] - lp[i, 5];
        real on  = p_return(lp[i, 1], lp[i, 2], lp[i, 3]);
        real off = p_return(lp[i, 1], e2_off, e3_off);
        p_on_post  += on;
        p_off_post += off;
        for (f in 1:n_fac) rd_grp[grp[i, f]] += on - off;
      }

      // Replicate the kit under the regime it actually had
      int w1 = bernoulli_logit_rng(lp[i, 1]);
      int w2 = w1 ? 0 : bernoulli_logit_rng(lp[i, 2]);
      int w3 = (w1 || w2) ? 0 : bernoulli_logit_rng(lp[i, 3]);
      ret_pat[pat[i]] += w1 + w2 + w3;
      s[per, 1] += 1;
      s[per, 2] += w1 + w2 + w3;
      s[per, 3] += w1;
      s[per, 4] += 1 - w1;
      s[per, 5] += w2;
      s[per, 6] += 1 - w1 - w2;
      s[per, 7] += w3;

      if (save_log_lik) log_lik[i] = kit_ll(outcome[i], lp[i, 1], lp[i, 2], lp[i, 3]);
    }

    p_on_post  /= n_post;
    p_off_post /= n_post;
    rd_post = p_on_post - p_off_post;
    for (g in 1:n_grp) {
      rd_grp[g] = n_post_grp[g] > 0 ? rd_grp[g] / n_post_grp[g] : not_a_number();
    }

    rep_ret_pre  = s[1, 2] / s[1, 1];
    rep_ret_post = s[2, 2] / s[2, 1];
    rep_w1_pre   = s[1, 3] / s[1, 1];
    rep_w1_post  = s[2, 3] / s[2, 1];
    rep_w2_pre   = s[1, 5] / s[1, 4];
    rep_w2_post  = s[2, 5] / s[2, 4];
    rep_w3_pre   = s[1, 6] > 0 ? s[1, 7] / s[1, 6] : not_a_number();
    rep_w3_post  = s[2, 6] > 0 ? s[2, 7] / s[2, 6] : not_a_number();

    {
      real all_ret = 0;
      real none_ret = 0;
      for (j in 1:J) {
        if (n_kits_pat[j] >= 2) {
          all_ret  += ret_pat[j] == n_kits_pat[j];
          none_ret += ret_pat[j] == 0;
        }
      }
      rep_all_multi  = all_ret / n_multi;
      rep_none_multi = none_ret / n_multi;
    }
  }
}
