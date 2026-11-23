## 01 - Audit & Clean =========================================================
##
## Input : 01_data/raw/crc_raw.xlsx  (one row per patient, kits spread wide)
## Output: 01_data/clean/kits.{rds,csv}            one row per kit  (primary modeling table)
##         01_data/clean/patients.{rds,csv}        one row per patient (demographics)
##         01_data/clean/patient_period.{rds,csv}  one row per patient x period (binomial counts)
##         01_data/clean/data_dictionary.csv
##         03_output/audit/*.csv                   audit log, corrections, cross-checks
##
## Run from the project root (the .Rproj directory).
##
## Key structural fact (verified below): kit slots 1-8 hold the COMPLETE kit history
## (2019-2023). The "2022/23" columns are a partial COPY of the post-2022 kits from
## slots 1-8, not additional kits. Slots 1-8 are the source of truth; the 2022/23
## block is used only to cross-check and to resolve date typos.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
})

raw_path   <- "01_data/raw/crc_raw.xlsx"
clean_dir  <- "01_data/clean"
audit_dir  <- "03_output/audit"
dir.create(clean_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(audit_dir, showWarnings = FALSE, recursive = TRUE)

reminder_start <- as.Date("2022-01-01")   # text reminders began in 2022
excel_origin   <- "1899-12-30"

num       <- function(x) suppressWarnings(as.numeric(x))
xl_date   <- function(x) as.Date(num(x), origin = excel_origin)

# Audit log: every check appends one row
audit <- list()
log_check <- function(area, check, n, action) {
  audit[[length(audit) + 1]] <<- tibble(area = area, check = check, n = as.character(n), action = action)
  message(sprintf("[%s] %s: %s -> %s", area, check, n, action))
}




## 1 - Load =====================================================================

# Read everything as text so nothing is silently coerced
raw <- read_excel(raw_path, sheet = "Sheet1", col_types = "text")

log_check("structure", "rows in raw file", nrow(raw), "one row per patient")
log_check("structure", "duplicate patientid", sum(duplicated(raw$patientid)), "none expected")
stopifnot(!any(duplicated(raw$patientid)))




## 2 - Kit-level reshape (slots 1-8) ==========================================

slot_cols <- function(s) paste(c("Dispense date", "Returned", "Returned date", "Result"), s)

kits_raw <- map_dfr(1:8, function(s) {
  cols <- slot_cols(s)
  tibble(patient_id      = raw$patientid,
         slot            = s,
         dispense_raw    = raw[[cols[1]]],
         returned_raw    = raw[[cols[2]]],
         return_date_raw = raw[[cols[3]]],
         result_raw      = raw[[cols[4]]])
}) %>%
  filter(!(is.na(dispense_raw) & is.na(returned_raw) & is.na(return_date_raw) & is.na(result_raw)))

log_check("kits", "kit records in slots 1-8", nrow(kits_raw), "source of truth")
log_check("kits", "kit records missing dispense date", sum(is.na(kits_raw$dispense_raw)), "none expected")
stopifnot(!any(is.na(kits_raw$dispense_raw)))

# The 2022/23 duplicate block, for cross-checking
b2223_slots <- c("2022/23", "2022/23 2nd", "2022/23 3rd", "2022/23 4th")
b2223 <- map_dfr(b2223_slots, function(s) {
  cols <- slot_cols(s)
  tibble(patient_id        = raw$patientid,
         b_slot            = s,
         dispense_date     = xl_date(raw[[cols[1]]]),
         b_returned_raw    = raw[[cols[2]]],
         b_return_date_raw = raw[[cols[3]]],
         b_result_raw      = raw[[cols[4]]])
}) %>%
  filter(!is.na(dispense_date))




## 3 - Kit field cleaning =======================================================

kits <- kits_raw %>%
  mutate(dispense_date = xl_date(dispense_raw),
         returned_raw  = str_to_upper(str_trim(returned_raw)),
         return_date   = xl_date(return_date_raw))

log_check("kits", "lowercase 'y' in Returned", sum(kits_raw$returned_raw == "y", na.rm = TRUE), "uppercased to Y")

# -- 3a. Return dates: non-numeric strings and implausible values -------------
# Every correction is listed explicitly so it can be reviewed.
date_fixes <- tribble(
  ~patient_id, ~slot, ~fixed_return_date, ~reason,
  "23460", 1, "2021-11-22", "typed as '11/22/221'; dispensed 2021-11-09",
  "23083", 1, "2021-12-02", "typed as '12/2/221'; dispensed 2021-11-17",
  "15120", 2, "2021-11-16", "typed as '11/16/221'; dispensed 2021-10-06",
  "17633", 4, "2022-06-24", "stored as text '6/24/2022'",
  "14467", 1, "2021-12-09", "recorded 1921-12-09 (century typo); dispensed 2021-09-14",
  "21702", 2, "2022-07-13", "recorded 2002-07-13 (year typo); dispensed 2022-06-29",
  "11984", 2, "2024-01-09", "recorded 2023-01-09, before 2023-12-29 dispense; 2022/23 block has 2024-01-09",
  "23225", 1, NA,           "recorded 2021-07-14, before 2021-08-12 dispense; no recoverable date",
  "25967", 1, NA,           "'UNK'; 2022/23 copy holds 1964-04-03 (invalid)",
  "19132", 4, NA,           "'UNK'",
  "15458", 4, NA,           "'UNK'"
) %>%
  mutate(fixed_return_date = as.Date(fixed_return_date))

kits <- kits %>%
  left_join(date_fixes, by = c("patient_id", "slot")) %>%
  mutate(return_date_fixed = !is.na(reason),
         return_date       = if_else(return_date_fixed, fixed_return_date, return_date)) %>%
  select(-fixed_return_date, -reason)

# Anything still unparseable or before dispense means a new problem the table does not cover
unhandled <- kits %>%
  filter((!is.na(return_date_raw) & is.na(return_date) & !return_date_fixed) |
           (!is.na(return_date) & return_date < dispense_date))
if (nrow(unhandled) > 0) {
  print(unhandled)
  stop("Unhandled return-date problems; extend `date_fixes`.")
}
log_check("kits", "return dates corrected (typo / text)", sum(!is.na(date_fixes$fixed_return_date)), "corrected; see corrections.csv")
log_check("kits", "return dates unrecoverable (UNK / before dispense)", sum(is.na(date_fixes$fixed_return_date)), "set NA; returned status kept as Y")

# -- 3b. Returned flag consistency --------------------------------------------
# Patient 24647 slot 1: flagged N but has a return date (166 d later) and a result,
# and the 2022/23 copy says Y. Treat as returned.
n_N_with_return <- sum(kits$returned_raw == "N" & !is.na(kits$return_date_raw), na.rm = TRUE)
log_check("kits", "Returned = N but return date present", n_N_with_return, "recoded to Y (2022/23 copy agrees)")

kits <- kits %>%
  mutate(returned = case_when(returned_raw == "Y" ~ 1L,
                              returned_raw == "N" & !is.na(return_date_raw) ~ 1L,
                              returned_raw == "N" ~ 0L,
                              TRUE ~ NA_integer_))

log_check("kits", "Returned flag missing", sum(is.na(kits$returned)), "none expected")
stopifnot(!any(is.na(kits$returned)))
log_check("kits", "Returned = Y but no return date (after fixes)",
          sum(kits$returned == 1 & is.na(kits$return_date)), "kept as returned; days_to_return NA")

# -- 3c. Results ---------------------------------------------------------------
log_check("kits", "distinct raw Result spellings", n_distinct(kits$result_raw, na.rm = TRUE), "normalised to 4 levels")

kits <- kits %>%
  mutate(result_key = str_to_lower(str_remove_all(result_raw, "[^A-Za-z ]")) %>% str_squish(),
         result = case_when(is.na(result_raw)                                   ~ NA_character_,
                            result_key %in% c("negative", "negaitve", "neagtive", "negtive") ~ "Negative",
                            result_key %in% c("positive", "postive")            ~ "Positive",
                            result_key %in% c("tnp", "no result")               ~ "Inadequate",   # test not performed
                            result_key %in% c("unk", "ukn")                     ~ "Unknown",
                            TRUE ~ "UNMAPPED"))
stopifnot(!any(kits$result == "UNMAPPED", na.rm = TRUE))
kits <- kits %>% mutate(result = if_else(returned == 0L, NA_character_, result)) %>% select(-result_key)

# -- 3d. Same-day duplicate kits ----------------------------------------------
# Three patients have two records with the same dispense date where one is a
# returned kit and the other is a blank "N" record: a double entry of one kit.
# Two other patients have two RETURNED records on one dispense date (distinct
# return dates/results): genuinely two kits, kept.
same_day <- kits %>%
  group_by(patient_id, dispense_date) %>%
  filter(n() > 1) %>%
  ungroup()
dup_drop <- same_day %>%
  group_by(patient_id, dispense_date) %>%
  filter(any(returned == 1L) & any(returned == 0L)) %>%
  filter(returned == 0L) %>%
  ungroup() %>%
  select(patient_id, slot)
log_check("kits", "same-day kit pairs", nrow(same_day) / 2, "inspected")
log_check("kits", "same-day blank duplicate removed", nrow(dup_drop), "dropped (one returned record kept)")

write_csv(same_day %>% select(patient_id, slot, dispense_date, returned, return_date, result) %>%
            mutate(dropped = paste(patient_id, slot) %in% paste(dup_drop$patient_id, dup_drop$slot)),
          file.path(audit_dir, "same_day_duplicates.csv"))

kits <- kits %>% anti_join(dup_drop, by = c("patient_id", "slot"))

# -- 3e. Derived kit fields ---------------------------------------------------
data_max_return <- max(kits$return_date, na.rm = TRUE)   # proxy for extraction date

kits <- kits %>%
  arrange(patient_id, dispense_date, slot) %>%
  group_by(patient_id) %>%
  mutate(kit_seq             = row_number(),                       # chronological order
         n_prior_kits        = kit_seq - 1L,
         n_prior_returned    = lag(cumsum(returned), default = 0L),
         prior_return_rate   = if_else(n_prior_kits > 0, n_prior_returned / n_prior_kits, NA_real_),
         days_since_prev_kit = as.numeric(dispense_date - lag(dispense_date)),
         # a non-returned kit followed by another within 30 days looks like a re-issue
         reissued_next_30d   = returned == 0L & !is.na(lead(dispense_date)) &
                               as.numeric(lead(dispense_date) - dispense_date) <= 30) %>%
  ungroup() %>%
  mutate(days_to_return  = as.numeric(return_date - dispense_date),
         returned_90d    = as.integer(returned == 1L & !is.na(days_to_return) & days_to_return <= 90),
         returned_90d    = if_else(returned == 1L & is.na(days_to_return), NA_integer_, returned_90d),
         days_followup   = as.numeric(data_max_return - dispense_date),
         post            = as.integer(dispense_date >= reminder_start),
         period          = factor(if_else(post == 1L, "Post", "Pre"), levels = c("Pre", "Post")),
         dispense_year   = as.integer(format(dispense_date, "%Y")),
         dispense_month  = as.integer(format(dispense_date, "%m")),
         # calendar time in years relative to reminder start (for secular trend / ITS terms)
         time_yrs        = as.numeric(dispense_date - reminder_start) / 365.25)

log_check("kits", "kits with days_to_return > 365", sum(kits$days_to_return > 365, na.rm = TRUE), "kept; flagged via returned_90d = 0")
log_check("kits", "non-returned kits re-issued within 30 d", sum(kits$reissued_next_30d), "kept; flagged reissued_next_30d")
log_check("kits", "latest return date in file (extraction proxy)", format(data_max_return), "used for days_followup")




## 4 - Cross-check against the 2022/23 copy block ============================

xchk <- kits %>%
  filter(post == 1L) %>%
  select(patient_id, slot, dispense_date, returned, return_date, result) %>%
  full_join(b2223, by = c("patient_id", "dispense_date")) %>%
  mutate(b_returned    = as.integer(str_to_upper(b_returned_raw) == "Y"),
         b_return_date = xl_date(b_return_date_raw),
         status = case_when(is.na(slot)   ~ "only in 2022/23 block",
                            is.na(b_slot) ~ "only in slots 1-8",
                            returned != b_returned ~ "returned flag differs",
                            is.na(return_date) != is.na(b_return_date) |
                              coalesce(return_date != b_return_date, FALSE) ~ "return date differs",
                            TRUE ~ "match"))
log_check("crosscheck", "post kits matched in 2022/23 block", sum(xchk$status == "match"), "agree")
log_check("crosscheck", "post kits only in slots 1-8", sum(xchk$status == "only in slots 1-8"), "2022/23 copy is incomplete")
log_check("crosscheck", "kits only in 2022/23 block", sum(xchk$status == "only in 2022/23 block"), "none expected (copy, not extra kits)")
log_check("crosscheck", "post kits with disagreeing fields", sum(xchk$status %in% c("returned flag differs", "return date differs")),
          "slots 1-8 kept unless shown wrong (see corrections)")
write_csv(xchk %>% filter(status != "match"), file.path(audit_dir, "crosscheck_2223_block.csv"))




## 5 - Patient-level demographics =============================================

hisp_eth <- c("Central American", "Chilean", "Colombian", "Dominican", "Ecuadorian",
              "Guatemalan", "Hispanic or Latino/Spanish", "Honduran", "Latin American/Latin, Latino",
              "Mexican", "Mexicano", "Peruvian", "Puerto Rican", "Salvadoran", "South American",
              "Spaniard", "Uruguayan", "Venezuelan")

patients <- raw %>%
  transmute(
    patient_id     = patientid,
    dob            = xl_date(patientdob),
    sex            = factor(patientsex, levels = c("F", "M"), labels = c("Female", "Male")),

    # Text-message consent: reminders can only reach consenting patients
    text_consent   = factor(coalesce(`ptnt cnsnt t txt ysn`, "Unknown"),
                            levels = c("Y", "N", "Unknown"), labels = c("Yes", "No", "Unknown")),

    # Race/ethnicity re-derived with one consistent rule (Hispanic ethnicity takes precedence).
    # The spreadsheet's 'Race/ethnicity adjusted' applied this inconsistently (e.g. Mexican + White
    # was coded White for 6 patients and Hispanic for 59).
    ethnicity_raw  = ethnicity,
    race_raw       = `multiple races`,
    race_eth       = case_when(
      ethnicity %in% hisp_eth                                               ~ "Hispanic/Latino",
      `multiple races` %in% c("White", "Polish")                            ~ "White",
      `multiple races` %in% c("Black or African American", "African")       ~ "Black",
      `multiple races` %in% c("Asian", "Chinese")                           ~ "Asian",
      `multiple races` %in% c("Other Race", "American Indian or Alaska Native") ~ "Other",
      TRUE                                                                  ~ "Unknown"),
    race_eth       = factor(race_eth, levels = c("Hispanic/Latino", "White", "Black", "Asian", "Other", "Unknown")),
    race_eth_sheet = `Race/ethnicity adjusted`,

    # Preferred language. ENG/SPN/POL received reminders in that language;
    # OTH did not have a language-concordant option.
    language       = case_when(Lang %in% c("ENG", "SPN", "POL") ~ Lang,
                               Lang %in% c("OTH", "SPN/OTH")    ~ "OTH",
                               TRUE                             ~ "Unknown"),
    language       = factor(language, levels = c("SPN", "POL", "ENG", "OTH", "Unknown"),
                            labels = c("Spanish", "Polish", "English", "Other", "Unknown")),
    lang_concordant = case_when(language %in% c("Spanish", "Polish", "English") ~ 1L,
                                language == "Other" ~ 0L,
                                TRUE ~ NA_integer_),

    # Employment. Raw NA was coded as 'unemployed' (0) in the sheet's derived columns; kept as Unknown here.
    employment_raw = `ptnt crrnt mplymnt stts`,
    employment     = case_when(str_detect(employment_raw, "^Full Time")     ~ "Full-time",
                               str_detect(employment_raw, "^Part Time")     ~ "Part-time",
                               employment_raw %in% c("Seasonal", "Temporary") ~ "Seasonal/Temporary",
                               employment_raw == "Other"                    ~ "Other",
                               employment_raw == "Unemployed"               ~ "Unemployed",
                               TRUE                                         ~ "Unknown"),
    employment     = factor(employment, levels = c("Full-time", "Part-time", "Seasonal/Temporary",
                                                   "Other", "Unemployed", "Unknown")),
    employed       = case_when(employment == "Unemployed" ~ 0L,
                               employment == "Unknown"    ~ NA_integer_,
                               TRUE                       ~ 1L),

    marital_raw    = `patient marital status`,
    relationship   = case_when(marital_raw %in% c("MARRIED", "PARTNER") ~ "Partnered",
                               marital_raw %in% c("SINGLE", "DIVORCED", "SEPARATED", "WIDOWED") ~ "Not partnered",
                               TRUE ~ "Unknown"),
    relationship   = factor(relationship, levels = c("Not partnered", "Partnered", "Unknown")),

    family_size    = as.integer(familysize),
    order_chartdate = xl_date(`order chartdate`),
    order_name     = `order name`
  )

log_check("patients", "sex missing (patientsex)", sum(is.na(patients$sex)), "complete; gender-identity columns not used")
log_check("patients", "text consent missing", sum(patients$text_consent == "Unknown"), "level 'Unknown'")
log_check("patients", "text consent = No", sum(patients$text_consent == "No"), "potential comparison group (no reminders)")
log_check("patients", "employment missing", sum(patients$employment == "Unknown"),
          "level 'Unknown' (sheet had coded these as unemployed)")
log_check("patients", "marital status missing/UNKNOWN", sum(patients$relationship == "Unknown"), "level 'Unknown'")
log_check("patients", "language missing", sum(patients$language == "Unknown"), "level 'Unknown'")
log_check("patients", "language 'SPN/OTH'", sum(raw$Lang == "SPN/OTH", na.rm = TRUE), "mapped to Other")
log_check("patients", "race/ethnicity differs from sheet's adjusted column",
          sum(as.character(patients$race_eth) != recode(patients$race_eth_sheet, "African American" = "Black",
                                                        "Declined" = "Unknown", "Native American" = "Other")),
          "re-derived consistently; see race_eth_recode.csv")
log_check("patients", "family size outside 1-8", sum(!between(patients$family_size, 1, 8)), "none")

write_csv(patients %>% count(ethnicity_raw, race_raw, race_eth_sheet, race_eth),
          file.path(audit_dir, "race_eth_recode.csv"))

# Verify the sheet's own derived summary columns against what we rebuild
sheet_totals <- raw %>%
  transmute(patient_id = patientid,
            sheet_n_kits     = num(`Total number of kits`),
            sheet_n_returned = coalesce(num(`Number of returned kits`), 0))   # sheet leaves 0 blank
our_totals <- kits %>%
  group_by(patient_id) %>%
  summarise(n_kits = n(), n_returned = sum(returned), .groups = "drop")
tot_chk <- full_join(sheet_totals, our_totals, by = "patient_id") %>%
  filter(sheet_n_kits != n_kits | sheet_n_returned != n_returned)
log_check("crosscheck", "patients where sheet's kit totals differ from rebuilt", nrow(tot_chk),
          "rebuilt counts used; sheet summary columns dropped")
write_csv(tot_chk, file.path(audit_dir, "sheet_totals_mismatch.csv"))




## 6 - Join, age, eligibility flags ============================================

kits <- kits %>%
  left_join(patients %>% select(patient_id, dob), by = "patient_id") %>%
  mutate(age_at_dispense = as.numeric(dispense_date - dob) / 365.25,
         # Screening guideline age band (USPSTF: 45-75 from May 2021; 50-75 before)
         age_in_45_75    = as.integer(age_at_dispense >= 45 & age_at_dispense < 76)) %>%
  select(-dob)

log_check("eligibility", "kits dispensed at age < 45", sum(kits$age_at_dispense < 45), "kept; flagged age_in_45_75 = 0")
log_check("eligibility", "kits dispensed at age >= 76", sum(kits$age_at_dispense >= 76), "kept; flagged age_in_45_75 = 0")

first_kit <- kits %>% group_by(patient_id) %>%
  summarise(first_dispense = min(dispense_date), last_dispense = max(dispense_date),
            n_kits = n(), n_returned = sum(returned),
            n_kits_pre = sum(post == 0L), n_kits_post = sum(post == 1L),
            ever_positive = as.integer(any(result == "Positive", na.rm = TRUE)), .groups = "drop") %>%
  mutate(periods_observed = factor(case_when(n_kits_pre > 0 & n_kits_post > 0 ~ "Both",
                                             n_kits_pre > 0 ~ "Pre only",
                                             TRUE ~ "Post only"),
                                   levels = c("Pre only", "Post only", "Both")))

patients <- patients %>%
  left_join(first_kit, by = "patient_id") %>%
  mutate(age_first_kit = as.numeric(first_dispense - dob) / 365.25)

log_check("patients", "patients with kits in both periods", sum(patients$periods_observed == "Both"), "within-patient pre/post contrast")
log_check("patients", "patients with pre-period kits only", sum(patients$periods_observed == "Pre only"), "")
log_check("patients", "patients with post-period kits only", sum(patients$periods_observed == "Post only"), "")

# Centering/scaling constants: from the kit-level analysis data so they match the model rows
age_center <- 55
kits <- kits %>%
  left_join(patients %>% select(patient_id, sex, text_consent, race_eth, language, lang_concordant,
                                employment, employed, relationship, family_size),
            by = "patient_id") %>%
  mutate(age_c10        = (age_at_dispense - age_center) / 10,     # per decade, centred at 55
         family_size_c  = family_size - mean(patients$family_size)) %>%
  select(patient_id, kit_seq, slot, dispense_date, dispense_year, dispense_month, time_yrs,
         period, post, returned, returned_90d, return_date, days_to_return, days_followup, result,
         age_at_dispense, age_c10, age_in_45_75,
         sex, text_consent, race_eth, language, lang_concordant, employment, employed,
         relationship, family_size, family_size_c,
         n_prior_kits, n_prior_returned, prior_return_rate, days_since_prev_kit, reissued_next_30d,
         return_date_fixed)




## 7 - Aggregated binomial table (patient x period) ===========================

patient_period <- kits %>%
  group_by(patient_id, period, post) %>%
  summarise(n_sent        = n(),
            n_returned    = sum(returned),
            return_rate   = n_returned / n_sent,
            mean_age      = mean(age_at_dispense),
            first_dispense = min(dispense_date),
            .groups = "drop") %>%
  left_join(patients %>% select(patient_id, sex, text_consent, race_eth, language, lang_concordant,
                                employment, employed, relationship, family_size, periods_observed),
            by = "patient_id") %>%
  mutate(age_c10       = (mean_age - age_center) / 10,
         family_size_c = family_size - mean(patients$family_size))




## 8 - Final checks & write ====================================================

stopifnot(
  nrow(kits) == nrow(kits_raw) - nrow(dup_drop),
  n_distinct(kits$patient_id) == nrow(raw),
  all(kits$returned %in% c(0L, 1L)),
  all(kits$days_to_return >= 0, na.rm = TRUE),
  sum(patient_period$n_sent) == nrow(kits),
  !any(is.na(kits$sex)), !any(is.na(kits$age_at_dispense))
)

log_check("final", "kits in analysis table", nrow(kits), "")
log_check("final", "patients", n_distinct(kits$patient_id), "")
log_check("final", "pre-period kits", sum(kits$post == 0L), "")
log_check("final", "post-period kits", sum(kits$post == 1L), "")

patients <- patients %>% select(-dob)   # DOB not needed downstream; age is derived

write_rds(kits,           file.path(clean_dir, "kits.rds"))
write_rds(patients,       file.path(clean_dir, "patients.rds"))
write_rds(patient_period, file.path(clean_dir, "patient_period.rds"))
write_csv(kits,           file.path(clean_dir, "kits.csv"))
write_csv(patients,       file.path(clean_dir, "patients.csv"))
write_csv(patient_period, file.path(clean_dir, "patient_period.csv"))

write_csv(bind_rows(audit), file.path(audit_dir, "audit_log.csv"))
write_csv(date_fixes, file.path(audit_dir, "corrections.csv"))




## 9 - Data dictionary =========================================================

dictionary <- tribble(
  ~table, ~variable, ~description,
  "kits", "patient_id", "Patient identifier (grouping factor for random effects)",
  "kits", "kit_seq", "Chronological kit number within patient (1 = first kit)",
  "kits", "slot", "Source slot (1-8) in raw sheet",
  "kits", "dispense_date", "Date kit was dispensed",
  "kits", "dispense_year / dispense_month", "Calendar year / month of dispense",
  "kits", "time_yrs", "Years from 2022-01-01 to dispense (negative = pre); for secular-trend terms",
  "kits", "period / post", "Pre (dispensed before 2022-01-01) vs Post (reminders); post = 1 for Post",
  "kits", "returned", "OUTCOME: 1 = kit returned, 0 = not returned",
  "kits", "returned_90d", "1 = returned within 90 days; NA if returned with unknown date",
  "kits", "return_date / days_to_return", "Return date and days from dispense (NA if not returned or date unknown)",
  "kits", "days_followup", "Days from dispense to latest return date in file (censoring check)",
  "kits", "result", "Negative / Positive / Inadequate / Unknown (NA if not returned)",
  "kits", "age_at_dispense", "Age in years at dispense (time-varying)",
  "kits", "age_c10", "(age_at_dispense - 55) / 10",
  "kits", "age_in_45_75", "1 = age 45-75 at dispense (screening guideline band)",
  "kits", "sex", "Female / Male (patientsex column)",
  "kits", "text_consent", "Consented to text messages: Yes / No / Unknown",
  "kits", "race_eth", "Hispanic/Latino, White, Black, Asian, Other, Unknown (re-derived)",
  "kits", "language", "Preferred language: Spanish, Polish, English, Other, Unknown",
  "kits", "lang_concordant", "1 = preferred language had a reminder version (ENG/SPN/POL), 0 = Other",
  "kits", "employment / employed", "6-level employment status; employed 1/0 (NA if unknown)",
  "kits", "relationship", "Partnered (married/partner), Not partnered, Unknown",
  "kits", "family_size / family_size_c", "Household size (1-8); mean-centred",
  "kits", "n_prior_kits / n_prior_returned / prior_return_rate", "Patient's kit history before this kit",
  "kits", "days_since_prev_kit", "Days since the patient's previous dispense",
  "kits", "reissued_next_30d", "TRUE = not returned and a new kit was dispensed within 30 days",
  "kits", "return_date_fixed", "TRUE = return date was manually corrected (see audit/corrections.csv)",
  "patient_period", "n_sent / n_returned / return_rate", "Kits dispensed / returned per patient-period (binomial form)",
  "patient_period", "mean_age / age_c10", "Mean age at dispense across that period's kits",
  "patients", "periods_observed", "Pre only / Post only / Both",
  "patients", "ever_positive", "1 = any kit had a positive result",
  "patients", "race_eth_sheet", "Original 'Race/ethnicity adjusted' column (for comparison only)"
)
write_csv(dictionary, file.path(clean_dir, "data_dictionary.csv"))

message("\nDone. Kits: ", nrow(kits), " | Patients: ", nrow(patients), " | Patient-periods: ", nrow(patient_period))
