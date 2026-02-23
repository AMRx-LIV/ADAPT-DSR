#FEATURE ENGINEERING

##Functions

###Assigning previous event feature variable
prev_event_assign <- function(
  df,
  B_var,
  event_df,
  event_var,
  no_days,
  no_events
) {
  df <- df %>%
    mutate(charttime = as.POSIXct(charttime, format = '%Y-%m-%d %H:%M:%S'))

  event_df %>%
    mutate(event = {{ event_var }}) %>%
    select('subject_id', "event", charttime = 'admittime') %>%
    mutate(charttime = as.POSIXct(charttime, format = '%Y-%m-%d %H:%M:%S')) %>%
    filter(!is.na(event)) %>%
    bind_rows(df) %>%
    mutate(event = case_when(!is.na(event) ~ "Yes", TRUE ~ "No")) %>%
    MIMER::check_previous_events(
      cols = "event",
      sort_by_col = 'charttime',
      patient_id_col = 'subject_id',
      event_indi_value = 'Yes',
      new_col_prefix = "pr_",
      time_period_in_days = no_days,
      minimum_prev_events = no_events,
      default_na_date = '9999-12-31 00:00:00'
    ) %>%
    mutate({{ B_var }} := case_when(pr_event == TRUE ~ TRUE, TRUE ~ FALSE)) %>%
    mutate(event = NULL, pr_event = NULL) %>%
    filter(grepl('URINE', spec_type_desc))
}

###Assigning previous event type feature variable
prev_event_type_assign <- function(
  df,
  B_var,
  event_df,
  event_var,
  event_type,
  no_days,
  no_events
) {
  df <- df %>%
    mutate(charttime = as.POSIXct(charttime, format = '%Y-%m-%d %H:%M:%S'))

  event_df %>%
    mutate(event = {{ event_var }}) %>%
    select('subject_id', "event", charttime = 'admittime') %>%
    mutate(charttime = as.POSIXct(charttime, format = '%Y-%m-%d %H:%M:%S')) %>%
    filter(grepl(event_type, event)) %>%
    bind_rows(df) %>%
    mutate(event = case_when(!is.na(event) ~ "Yes", TRUE ~ "No")) %>%
    MIMER::check_previous_events(
      cols = "event",
      sort_by_col = 'charttime',
      patient_id_col = 'subject_id',
      event_indi_value = 'Yes',
      new_col_prefix = "pr_",
      time_period_in_days = no_days,
      minimum_prev_events = no_events,
      default_na_date = '9999-12-31 00:00:00'
    ) %>%
    mutate({{ B_var }} := case_when(pr_event == TRUE ~ TRUE, TRUE ~ FALSE)) %>%
    mutate(event = NULL, pr_event = NULL) %>%
    filter(grepl('URINE', spec_type_desc))
}

###Assigning "NT" variable to relevant NAs in microbiology dataframe
NT_assigner <- function(df) {
  micaborgs <- df %>% filter(!is.na(org_name))
  micabnas <- df %>% filter(is.na(org_name))
  micaborgab <- micaborgs %>% select(PEN:MTR)
  micaborgab[is.na(micaborgab)] <- "NT"
  micaborgs[, 17:81] <- micaborgab
  df2 <- tibble(rbind(micaborgs, micabnas))
  df2 %>% rename(admittime = "charttime")
}

###Applying previous AST result search across multiple result types
prev_AST_applier <- function(
  df1,
  micro_data,
  suffix,
  result,
  timeframe = 365,
  n_events = 1
) {
  params <- paste0("p", antibiotics, suffix)

  apply_prev_event <- function(df, param, antibiotic) {
    df %>%
      prev_event_type_assign(
        !!sym(param),
        micro_data,
        !!sym(antibiotic),
        result,
        timeframe,
        n_events
      )
  }
  df1 <- reduce(
    seq_along(antibiotics),
    function(df, i) {
      apply_prev_event(df, params[i], antibiotics[i])
    },
    .init = df1
  ) %>%
    ungroup()
}

###Assigning previous antimicrobial treatment variable
prev_rx_assign <- function(
  df,
  B_var,
  drug_df,
  abx,
  abx_groupvar,
  no_days,
  no_events
) {
  ur_df <- df %>%
    mutate(charttime = as.POSIXct(charttime, format = '%Y-%m-%d %H:%M:%S'))

  abx_groupvar <- enquo(abx_groupvar)

  drug_df %>%
    select('subject_id', ab_name, charttime = 'starttime') %>%
    mutate(charttime = as.POSIXct(charttime, format = '%Y-%m-%d %H:%M:%S')) %>%
    filter(grepl(glue("{abx}"), !!abx_groupvar)) %>%
    bind_rows(ur_df) %>%
    mutate(abx_treatment = case_when(!is.na(ab_name) ~ "Yes", TRUE ~ "No")) %>%
    MIMER::check_previous_events(
      cols = "abx_treatment",
      sort_by_col = 'charttime',
      patient_id_col = 'subject_id',
      event_indi_value = 'Yes',
      new_col_prefix = "pr_rx_",
      time_period_in_days = no_days,
      minimum_prev_events = no_events,
      default_na_date = '9999-12-31 00:00:00'
    ) %>%
    mutate(
      {{ B_var }} := case_when(pr_rx_abx_treatment == TRUE ~ TRUE, TRUE ~ FALSE)
    ) %>%
    mutate(abx_treatment = NULL, pr_rx_abx_treatment = NULL) %>%
    filter(grepl('URINE', spec_type_desc))
}

###Finding abnormal inflammatory markers on day of urine test
labevent_search <- function(df, search_term, feature_name) {
  feature_name <- enquo(feature_name)

  filter_term <- labitems %>%
    filter(grepl(search_term, label, ignore.case = T)) %>%
    count(itemid) %>%
    arrange(n) %>%
    slice(1) %>%
    select(itemid) %>%
    unlist()
  filtered_df <- labevents %>%
    filter(itemid == filter_term) %>%
    filter(!is.na(valuenum)) %>%
    rename(admittime = "charttime")
  df %>%
    prev_event_type_assign(
      !!feature_name,
      filtered_df,
      flag,
      "abnormal",
      1,
      1
    ) %>%
    ungroup()
}

###Assigning gender feature variable
gender_assign <- function(df, B_var, gender_df) {
  gender_df %>%
    select('subject_id', 'gender') %>%
    right_join(df) %>%
    mutate({{ B_var }} := case_when(gender == "M" ~ TRUE, TRUE ~ FALSE)) %>%
    mutate(gender = NULL)
}

###Finding patient demographic characeristics
demographic_assign <- function(df, demographic) {
  demographic <- enquo(demographic)

  hadm_demographic <- hadm %>%
    select(subject_id, !!demographic) %>%
    distinct(subject_id, .keep_all = T)
  df %>%
    left_join(hadm_demographic, by = "subject_id") %>%
    mutate(
      !!demographic := case_when(
        is.na(!!demographic) ~ "UNKNOWN",
        TRUE ~ !!demographic
      )
    )
}

###Applying ICD-1O code search across multiple ICD-10 code prefixes
prev_ICD_applier <- function(df, icd_df, prefix, codes) {
  apply_prev_event_assignments <- function(df, code) {
    param_name <- paste0(prefix, code)
    df %>%
      prev_event_type_assign(!!sym(param_name), icd_df, icd_group, code, 365, 1)
  }

  urine_dsr <- reduce(
    codes,
    function(df, code) {
      apply_prev_event_assignments(df, code)
    },
    .init = urine_dsr
  ) %>%
    mutate(pDIAG_U = FALSE) %>%
    ungroup()
}

###Checking for previous care events
care_event_assigner <- function(
  df,
  search_df,
  search_term,
  search_column,
  feature_name,
  event_date_col,
  timeframe,
  n_events = 1
) {
  feature_name <- enquo(feature_name)
  search_column <- enquo(search_column)

  care_event <- search_df %>%
    filter(grepl(search_term, !!search_column, ignore.case = T)) %>%
    mutate(
      !!search_column := search_term
    ) %>%
    rename(admittime = event_date_col)
  df %>%
    prev_event_type_assign(
      !!feature_name,
      care_event,
      !!search_column,
      search_term,
      timeframe,
      n_events
    ) %>%
    ungroup()
}

###Applying BMI category search across multiple categories
assign_bmi_events <- function(df, bmi_df, categories, days, min_events) {
  reduce(
    categories,
    function(acc, category) {
      param <- paste0("p", category)
      prev_event_type_assign(
        acc,
        !!sym(param),
        bmi_df,
        BMI_cat,
        category,
        days,
        min_events
      )
    },
    .init = df
  )
}

###Attaching previous diagnosis label
diag_label <- function(df, filter_term, label, timeframe, newvar) {
  filter_term <- enquo(filter_term)
  newvar <- enquo(newvar)

  df <- df %>%
    prev_event_type_assign(!!newvar, hadm, description, label, timeframe, 1) %>%
    ungroup() %>%
    filter(!is.na(!!filter_term))
}

###Attaching cytotoxins
cytotoxic_join <- function(df) {
  df %>%
    left_join(cytotoxics_key) %>%
    mutate(
      Cytotoxic_agent = case_when(
        is.na(Cytotoxic_agent) ~ FALSE,
        TRUE ~ Cytotoxic_agent
      )
    )
}

###Filter to those with triage obs available
obs_filter <- function(df) {
  df %>%
    filter(
      !is.na(temperature) &
        !is.na(heartrate) &
        !is.na(resprate) &
        !is.na(o2sat) &
        !is.na(sbp) &
        !is.na(dbp) &
        !is.na(acuity) &
        !is.na(chiefcomplaint)
    )
}

###Standardising observation values
standardize_obs <- function(df) {
  df$heartrate <- standardize(df$heartrate)
  df$resprate <- standardize(df$resprate)
  df$o2sat <- standardize(df$o2sat)
  df$sbp <- standardize(df$sbp)
  df$dbp <- standardize(df$dbp)

  df
}

###Dummy variables for main presenting complains
pc_dummies <- function(df) {
  df %>%
    mutate(
      pc_dyspnea = case_when(
        grepl("Dyspnea", chiefcomplaint, ignore.case = T) |
          grepl("shortness", chiefcomplaint, ignore.case = T) |
          grepl("sob", chiefcomplaint, ignore.case = T) |
          grepl("hypoxia", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_abdopain = case_when(
        grepl("abd", chiefcomplaint, ignore.case = T) &
          grepl("pain", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_confusion = case_when(
        (grepl("altered", chiefcomplaint, ignore.case = T) &
          grepl("mental", chiefcomplaint, ignore.case = T) |
          grepl("confus", chiefcomplaint, ignore.case = T)) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_chestpain = case_when(
        grepl("chest", chiefcomplaint, ignore.case = T) &
          grepl("pain", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_weakness = case_when(
        grepl("weakness", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_dyspnea = case_when(
        grepl("fever", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_wound = case_when(
        grepl("wound", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_fall = case_when(
        grepl("fall", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_prbleed = case_when(
        grepl("brbpr", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_vomiting = case_when(
        grepl("N/V", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_backpain = case_when(
        grepl("back", chiefcomplaint, ignore.case = T) &
          grepl("pain", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_lethargy = case_when(
        grepl("lethargy", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_diarvom = case_when(
        grepl("N/V/D", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_diarrhea = case_when(
        grepl("diarrhea", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_headache = case_when(
        grepl("headache", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_syncope = case_when(
        grepl("syncope", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_seizure = case_when(
        grepl("seizure", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_flankpain = case_when(
        grepl("flank", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_lowbp = case_when(
        grepl("hypotension", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_anemia = case_when(
        grepl("anemia", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_pain = case_when(
        grepl("pain", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_swelling = case_when(
        grepl("swelling", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_cough = case_when(
        grepl("cough", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_sepsis = case_when(
        grepl("cough", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      ),
      pc_fever = case_when(
        grepl("fever", chiefcomplaint, ignore.case = T) |
          grepl("pyrexia", chiefcomplaint, ignore.case = T) ~ TRUE,
        TRUE ~ FALSE
      )
    )
}

###Data upload (CSV files accessible at https://physionet.org/content/mimiciv/2.2/)
setwd("/Users/alexhoward/Documents/Projects/UDAST_code")
omr <- read_csv("omr.csv") #Measurements e.g., height, weight
hadm <- read_csv("admissions.csv") #Admission data
labevents <- read_csv("labevents.csv") #Laboratory tests (non-micro)
labitems <- read_csv("d_labitems.csv") #Laboratory test codes
pats <- read_csv("patients.csv") #Patient demographics
services <- read_csv("services.csv") #Service providers
d_icd_diagnoses <- read_csv("d_icd_diagnoses.csv") #icd codes
diagnoses_raw <- read_csv("diagnoses_icd.csv") #icd epi
diagnoses <- read_csv("diagnoses_clean.csv")
procedures <- read_csv("procedures_clean.csv")
poe <- read_csv("poe_clean.csv")
micro <- read_csv("micro_clean2.csv")
drugs <- read_csv("drugs_clean.csv")
drgcodes <- read_csv("drgcodes.csv")
vitalsign <- read_csv("vitalsign.csv")
edstays <- read_csv("edstays.csv")
triage <- read_csv("triage.csv")

urine_dsr <- read_csv("urine_pre_features.csv")
blood_dsr <- read_csv("blood_pre_features.csv")

admitkey <- hadm %>% select(hadm_id, admittime)
blood_dsr <- blood_dsr %>% left_join(admitkey)
urine_dsr <- urine_dsr %>% left_join(admitkey)

admit_chart <- function(df) {
  df %>% select(-charttime) %>% rename(charttime = "admittime")
}

urine_dsr <- admit_chart(urine_dsr)
blood_dsr <- admit_chart(blood_dsr)


##Finding previous AST results

###Assigning modified microbiology dataframes to enable prev_event_type_assign
micro3 <- micro %>% rename(admittime = "charttime")
micro2 <- micro %>% NT_assigner()

###At least one resistant isolate in the last year
antibiotics <- c(
  "AMP",
  "SAM",
  "TZP",
  "CZO",
  "CRO",
  "CAZ",
  "FEP",
  "MEM",
  "CIP",
  "GEN",
  "SXT",
  "NIT",
  "VAN",
  "AMPC",
  "TCY",
  "PEN",
  "CLI",
  "LVX",
  "AMK",
  "TOB"
)
urine_dsr <- prev_AST_applier(urine_dsr, micro3, "r", "R")

###At least one previous resistant isolate in the last week
urine_dsr <- prev_AST_applier(urine_dsr, micro2, "7dr", "R", 7, 1)

###At least one susceptible isolate in the last year
urine_dsr <- prev_AST_applier(urine_dsr, micro3, "s", "S")

###At least one 'I' isolate in the last year
antibiotics <- antibiotics[
  antibiotics != "AMPC" &
    antibiotics != "SXT" &
    antibiotics != "VAN"
] #No 'I' results
urine_dsr <- prev_AST_applier(urine_dsr, micro3, "i", "I")

###At least one isolate with the antimicrobial not tested for in the last year
antibiotics <- c(
  "AMP",
  "SAM",
  "TZP",
  "CZO",
  "CRO",
  "CAZ",
  "FEP",
  "MEM",
  "CIP",
  "SXT",
  "VAN",
  "PEN"
)
urine_dsr <- prev_AST_applier(urine_dsr, micro2, "nt", "NT")

###Finding previous complicated uti
urine_dsr <- urine_dsr %>%
  MIMER::check_previous_events(
    cols = "complicated_vs_uncomplicated_uti",
    sort_by_col = 'charttime',
    patient_id_col = 'subject_id',
    event_indi_value = '1',
    new_col_prefix = "pr_",
    time_period_in_days = 1e4,
    minimum_prev_events = 1,
    default_na_date = '9999-12-31 00:00:00'
  ) %>%
  ungroup()

##Finding previous antimicrobial treatment

###Modifying prescriptions dataframes to enable prev_event_type_assign
drugs <- drugs %>% rename(ab_name = "abx_name")

###Assigning reference lists of antimicrobials and new feature suffixes
antibiotics <- c(
  "Ampicillin",
  "Amoxicillin",
  "Amoxicillin/clavulanic acid",
  "Ampicillin/sulbactam",
  "Piperacillin/tazobactam",
  "Cefazolin",
  "Cefalexin",
  "Cefpodoxime proxetil",
  "Ceftriaxone",
  "Ceftazidime",
  "Cefepime",
  "Meropenem",
  "Ertapenem",
  "Aztreonam",
  "Ciprofloxacin",
  "Levofloxacin",
  "Gentamicin",
  "Tobramycin",
  "Amikacin",
  "Rifampicin",
  "Trimethoprim/sulfamethoxazole",
  "Nitrofurantoin",
  "Erythromycin",
  "Clarithromycin",
  "Azithromycin",
  "Clindamycin",
  "Vancomycin",
  "Metronidazole",
  "Linezolid",
  "Daptomycin",
  "Doxycycline"
)
suffixes <- c(
  "AMPrx",
  "AMXrx",
  "AMCrx",
  "SAMrx",
  "TZPrx",
  "CZOrx",
  "CZOrx",
  "CZOrx",
  "CROrx",
  "CAZrx",
  "FEPrx",
  "MEMrx",
  "ETPrx",
  "ATMrx",
  "CIPrx",
  "CIPrx",
  "GENrx",
  "TOBrx",
  "AMKrx",
  "RIFrx",
  "SXTrx",
  "NITrx",
  "ERYrx",
  "CLRrx",
  "AZMrx",
  "CLIrx",
  "VANrx",
  "MTRrx",
  "LNZrx",
  "DAPrx",
  "DOXrx"
)

###At least one inpatient antimicrobial prescription in the last year
apply_prev_rx <- function(df, suffix, antibiotic) {
  param_name <- paste0("p", suffix)
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, 365, 1)
}
urine_dsr <- reduce(
  seq_along(antibiotics),
  function(df, i) {
    apply_prev_rx(df, suffixes[i], antibiotics[i])
  },
  .init = urine_dsr
) %>%
  ungroup()

###At least one inpatient antimicrobial prescription in the last week
apply_prev_rx <- function(df, suffix, antibiotic) {
  param_name <- paste0("d7", suffix)
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, 7, 1)
}
urine_dsr <- reduce(
  seq_along(antibiotics),
  function(df, i) {
    apply_prev_rx(df, suffixes[i], antibiotics[i])
  },
  .init = urine_dsr
) %>%
  ungroup()

##Find inflammatory marker results on admission

###Elevated C-reactive protein on the same day as the urine specimen
urine_dsr <- urine_dsr %>% labevent_search("reactive", highCRP)

###Abnormal total peripheral white cell count on the same day as the urine specimen
urine_dsr <- urine_dsr %>% labevent_search("White", abnormalWCC)

##Find patient characteristic and history variables

###At least one hospital admission in the last year
urine_dsr <- urine_dsr %>%
  prev_event_assign(pHADM, hadm, hadm_id, 365, 1) %>%
  ungroup()

###At least one discharge to a nursing home
urine_dsr <- urine_dsr %>%
  prev_event_type_assign(pNH, hadm, discharge_location, "NURSING", 1e4, 1) %>%
  ungroup()

###Hospital admission from outpatient location
outpatient_check <- function(df) {
  df %>%
    mutate(
      admission_location = case_when(
        is.na(admission_location) ~ "OUTPATIENT",
        TRUE ~ admission_location
      )
    )
}
hadm_admission <- hadm %>%
  select(hadm_id, admission_location) %>%
  outpatient_check() %>%
  mutate(hadm_id = case_when(is.na(hadm_id) ~ 0, TRUE ~ hadm_id)) %>%
  distinct(hadm_id, .keep_all = T)
urine_dsr <- left_join(urine_dsr, hadm_admission, by = "hadm_id") %>%
  outpatient_check()

###Coded ICD-10 diagnosis groups for admissions in the last year
diag_codes <- c(
  "A",
  "B",
  "C",
  "D",
  "E",
  "F",
  "G",
  "H",
  "I",
  "J",
  "K",
  "L",
  "M",
  "N",
  "O",
  "P",
  "Q",
  "R",
  "S",
  "T",
  "V",
  "W",
  "X",
  "Y",
  "Z"
)
urine_dsr <- urine_dsr %>% prev_ICD_applier(diagnoses, "pDIAG_", diag_codes)

###Coded ICD-10 procedure groups for admissions in the last year
proc_codes <- c(
  "0",
  "3",
  "8",
  "5",
  "T",
  "4",
  "S",
  "A",
  "9",
  "H",
  "I",
  "B",
  "7",
  "G",
  "1",
  "R",
  "J",
  "Q",
  "K",
  "6",
  "M",
  "P",
  "L",
  "D",
  "F",
  "2",
  "N",
  "C",
  "E",
  "X",
  "O"
)
urine_dsr <- urine_dsr %>% prev_ICD_applier(procedures, "pPROC_", proc_codes)

###At least one coded previous UTI diagnosis in the last year
uti_key <- d_icd_diagnoses %>%
  filter(
    grepl("urinary tract infection", long_title, ignore.case = T) |
      grepl("acute pyelon", long_title, ignore.case = T) |
      (grepl("urinary catheter", long_title, ignore.case = T) &
        grepl("infec", long_title, ignore.case = T))
  )
hadm_key <- hadm %>% select(hadm_id, admittime)
uti_df <- diagnoses_raw %>%
  left_join(uti_key, by = c("icd_code", "icd_version")) %>%
  filter(!is.na(long_title)) %>%
  left_join(hadm_key, by = "hadm_id")
urine_dsr <- urine_dsr %>%
  care_event_assigner(
    uti_df,
    "(urin|pyelo|cath)",
    long_title,
    pUTI,
    "admittime",
    365,
    1
  )

###Presence of an outpatient provider ID
urine_dsr <- urine_dsr %>%
  mutate(provider_id = case_when(order_provider_id != "" ~ TRUE, TRUE ~ FALSE))

###Current inpatient specialty
serv_key <- services %>%
  select(hadm_id, curr_service) %>%
  distinct(hadm_id, .keep_all = T)
urine_dsr <- urine_dsr %>%
  left_join(serv_key, by = "hadm_id") %>%
  mutate(
    curr_service = case_when(
      is.na(curr_service) ~ "UNKNOWN",
      TRUE ~ curr_service
    )
  )

###At least one measured obese, underweight, or overweight BMI category in the last 3 years
categorise_bmi <- function(df) {
  df %>%
    filter(grepl("BMI", result_name)) %>%
    mutate(
      BMI_cat = case_when(
        as.numeric(result_value) >= 30 ~ "Obese",
        as.numeric(result_value) >= 25 &
          as.numeric(result_value) < 30 ~ "Overweight",
        as.numeric(result_value) >= 18.5 &
          as.numeric(result_value) < 25 ~ "Normal weight",
        as.numeric(result_value) < 18.5 ~ "Underweight"
      ),
      admittime = as.POSIXct(chartdate, format = '%Y-%m-%d %H:%M:%S')
    )
}
bmi <- categorise_bmi(omr)
bmi_categories <- c("Obese", "Underweight", "Overweight")
urine_dsr <- assign_bmi_events(urine_dsr, bmi, bmi_categories, 1095, 1) %>%
  ungroup()

###Observation frequency on day of test
obs <- poe %>%
  filter(order_subtype == "Vitals/Monitoring") %>%
  mutate(ordertime = as.Date(ordertime)) %>%
  group_by(subject_id, ordertime) %>%
  count(order_subtype) %>%
  arrange(desc(n)) %>%
  select(-order_subtype)
urine_dsr <- urine_dsr %>%
  mutate(ordertime = chartdate) %>%
  left_join(obs, by = c("subject_id", "ordertime")) %>%
  rename(ob_freq = "n") %>%
  mutate(
    ob_freq = case_when(is.na(ob_freq) ~ 0, TRUE ~ ob_freq),
    ob_freq = standardize(ob_freq)
  ) %>%
  select(-ordertime)

###Other specific previous care events
urine_dsr <- urine_dsr %>%
  care_event_assigner(poe, "cath", field_value, pCATH, "ordertime", 28) %>% ###At least one urinary catheter insertion in the last 28 days
  care_event_assigner(poe, "DNR", field_value, pDNR, "ordertime", 365) %>% ###At least one 'do not resuscitate' order in the last year
  care_event_assigner(poe, "Discharge", field_value, pDISC, "ordertime", 28) %>% ###At least one discharge from hospital in the last 28 days
  care_event_assigner(poe, "ICU", field_value, pICU, "ordertime", 28) %>% ###At least one intensive care admission in the last 28 days
  care_event_assigner(
    poe,
    "Psychiatry",
    field_value,
    pPsych,
    "ordertime",
    365
  ) %>% ###At least one psychiatry review in the last year
  care_event_assigner(
    poe,
    "Nephrostomy",
    field_value,
    pNeph,
    "ordertime",
    365
  ) %>% ###At least one nephrostomy insertion in the last year
  care_event_assigner(poe, "Surgery", field_value, pSURG, "ordertime", 365) %>% ###At least one surgical procedure in the last year
  care_event_assigner(poe, "Hydration", field_value, pHyd, "ordertime", 28) %>% ###At least one hydration order in the last 28 days
  care_event_assigner(poe, "NGT", field_value, pNGT, "ordertime", 28) %>% ###At least one nasogastric tube insertion in the last 28 days
  care_event_assigner(poe, "Chemo", field_value, pChemo, "ordertime", 28) %>% ###At least one administration of cancer chemotherapy in the last 28 days
  care_event_assigner(
    poe,
    "Nutrition consult",
    order_subtype,
    pNUTR,
    "ordertime",
    365
  ) %>% ###At least one nutrition consultation in the last year
  care_event_assigner(
    poe,
    "Physical Therapy",
    order_subtype,
    pPhysio,
    "ordertime",
    365
  ) %>% ###At least one physiotherapy consultation in the last year
  care_event_assigner(
    poe,
    "Restraints",
    order_subtype,
    pRestr,
    "ordertime",
    365
  ) %>% ###At least one requirement for restraints in the last year
  care_event_assigner(
    poe,
    "Occupational Therapy",
    order_subtype,
    pOT,
    "ordertime",
    365
  ) %>% ###At least one occupational therapy consultation in the last year
  care_event_assigner(poe, "Central TPN", order_subtype, pTPN, "ordertime", 365) ###At least one administration of total parenteral nutrition in the last year

###Revert charttime to admittime
urine_dsr <- urine_dsr %>% rename(admittime = "charttime")

##Write CSVs
write_csv(urine_dsr, "urine_dsr_w_features.csv")
urine_dsr <- urine_dsr %>% distinct(hadm_id, .keep_all = T)
write_csv(urine_dsr, "preprocessed_urines.csv")


#Blood
blood_dsr <- blood_dsr %>% mutate(spec_type_desc = "URINE")
blood_dsr <- blood_dsr %>% filter(!is.na(org_name))

prev_ICD_applier <- function(df, icd_df, prefix, codes) {
  apply_prev_event_assignments <- function(df, code) {
    param_name <- paste0(prefix, code)
    df %>%
      prev_event_type_assign(!!sym(param_name), icd_df, icd_group, code, 365, 1)
  }

  blood_dsr <- reduce(
    codes,
    function(df, code) {
      apply_prev_event_assignments(df, code)
    },
    .init = blood_dsr
  ) %>%
    mutate(pDIAG_U = FALSE) %>%
    ungroup()
}

##Finding previous AST results

###Assigning modified microbiology dataframes to enable prev_event_type_assign
micro3 <- micro %>% rename(admittime = "charttime")
micro2 <- micro %>% NT_assigner()

###At least one resistant isolate in the last year
antibiotics <- c(
  "AMP",
  "SAM",
  "TZP",
  "CZO",
  "CRO",
  "CAZ",
  "FEP",
  "MEM",
  "CIP",
  "GEN",
  "SXT",
  "NIT",
  "VAN",
  "AMPC",
  "TCY",
  "PEN",
  "CLI",
  "LVX",
  "AMK",
  "TOB"
)
blood_dsr <- prev_AST_applier(blood_dsr, micro3, "r", "R")

###At least one previous resistant isolate in the last week
blood_dsr <- prev_AST_applier(blood_dsr, micro2, "7dr", "R", 7, 1)

###At least one susceptible isolate in the last year
blood_dsr <- prev_AST_applier(blood_dsr, micro3, "s", "S")

###At least one 'I' isolate in the last year
antibiotics <- antibiotics[
  antibiotics != "AMPC" &
    antibiotics != "SXT" &
    antibiotics != "VAN"
] #No 'I' results
blood_dsr <- prev_AST_applier(blood_dsr, micro3, "i", "I")

###At least one isolate with the antimicrobial not tested for in the last year
antibiotics <- c(
  "AMP",
  "SAM",
  "TZP",
  "CZO",
  "CRO",
  "CAZ",
  "FEP",
  "MEM",
  "CIP",
  "SXT",
  "VAN",
  "PEN"
)
blood_dsr <- prev_AST_applier(blood_dsr, micro2, "nt", "NT")

###Finding previous complicated uti
blood_dsr <- blood_dsr %>%
  MIMER::check_previous_events(
    cols = "uti_vs_not_uti",
    sort_by_col = 'charttime',
    patient_id_col = 'subject_id',
    event_indi_value = '1',
    new_col_prefix = "pr_",
    time_period_in_days = 1e4,
    minimum_prev_events = 1,
    default_na_date = '9999-12-31 00:00:00'
  ) %>%
  ungroup()

##Finding previous antimicrobial treatment

###Assigning reference lists of antimicrobials and new feature suffixes
antibiotics <- c(
  "Ampicillin",
  "Amoxicillin",
  "Amoxicillin/clavulanic acid",
  "Ampicillin/sulbactam",
  "Piperacillin/tazobactam",
  "Cefazolin",
  "Cefalexin",
  "Cefpodoxime proxetil",
  "Ceftriaxone",
  "Ceftazidime",
  "Cefepime",
  "Meropenem",
  "Ertapenem",
  "Aztreonam",
  "Ciprofloxacin",
  "Levofloxacin",
  "Gentamicin",
  "Tobramycin",
  "Amikacin",
  "Rifampicin",
  "Trimethoprim/sulfamethoxazole",
  "Nitrofurantoin",
  "Erythromycin",
  "Clarithromycin",
  "Azithromycin",
  "Clindamycin",
  "Vancomycin",
  "Metronidazole",
  "Linezolid",
  "Daptomycin",
  "Doxycycline"
)
suffixes <- c(
  "AMPrx",
  "AMXrx",
  "AMCrx",
  "SAMrx",
  "TZPrx",
  "CZOrx",
  "CZOrx",
  "CZOrx",
  "CROrx",
  "CAZrx",
  "FEPrx",
  "MEMrx",
  "ETPrx",
  "ATMrx",
  "CIPrx",
  "CIPrx",
  "GENrx",
  "TOBrx",
  "AMKrx",
  "RIFrx",
  "SXTrx",
  "NITrx",
  "ERYrx",
  "CLRrx",
  "AZMrx",
  "CLIrx",
  "VANrx",
  "MTRrx",
  "LNZrx",
  "DAPrx",
  "DOXrx"
)

###At least one inpatient antimicrobial prescription in the last year
apply_prev_rx <- function(df, suffix, antibiotic) {
  param_name <- paste0("p", suffix)
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, 365, 1)
}
blood_dsr <- reduce(
  seq_along(antibiotics),
  function(df, i) {
    apply_prev_rx(df, suffixes[i], antibiotics[i])
  },
  .init = blood_dsr
) %>%
  ungroup()

###At least one inpatient antimicrobial prescription in the last week
apply_prev_rx <- function(df, suffix, antibiotic) {
  param_name <- paste0("d7", suffix)
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, 7, 1)
}
blood_dsr <- reduce(
  seq_along(antibiotics),
  function(df, i) {
    apply_prev_rx(df, suffixes[i], antibiotics[i])
  },
  .init = blood_dsr
) %>%
  ungroup()

##Find inflammatory marker results on admission

###Elevated C-reactive protein on the same day as the urine specimen
blood_dsr <- blood_dsr %>% labevent_search("reactive", highCRP)

###Abnormal total peripheral white cell count on the same day as the urine specimen
blood_dsr <- blood_dsr %>% labevent_search("White", abnormalWCC)

##Find patient characteristic and history variables

###At least one hospital admission in the last year
blood_dsr <- blood_dsr %>%
  prev_event_assign(pHADM, hadm, hadm_id, 365, 1) %>%
  ungroup()

###At least one discharge to a nursing home
blood_dsr <- blood_dsr %>%
  prev_event_type_assign(pNH, hadm, discharge_location, "NURSING", 1e4, 1) %>%
  ungroup()

###Hospital admission from outpatient location
outpatient_check <- function(df) {
  df %>%
    mutate(
      admission_location = case_when(
        is.na(admission_location) ~ "OUTPATIENT",
        TRUE ~ admission_location
      )
    )
}
hadm_admission <- hadm %>%
  select(hadm_id, admission_location) %>%
  outpatient_check() %>%
  mutate(hadm_id = case_when(is.na(hadm_id) ~ 0, TRUE ~ hadm_id)) %>%
  distinct(hadm_id, .keep_all = T)
blood_dsr <- left_join(blood_dsr, hadm_admission, by = "hadm_id") %>%
  outpatient_check()

###Coded ICD-10 diagnosis groups for admissions in the last year
diag_codes <- c(
  "A",
  "B",
  "C",
  "D",
  "E",
  "F",
  "G",
  "H",
  "I",
  "J",
  "K",
  "L",
  "M",
  "N",
  "O",
  "P",
  "Q",
  "R",
  "S",
  "T",
  "V",
  "W",
  "X",
  "Y",
  "Z"
)
blood_dsr <- blood_dsr %>% prev_ICD_applier(diagnoses, "pDIAG_", diag_codes)

###Coded ICD-10 procedure groups for admissions in the last year
proc_codes <- c(
  "0",
  "3",
  "8",
  "5",
  "T",
  "4",
  "S",
  "A",
  "9",
  "H",
  "I",
  "B",
  "7",
  "G",
  "1",
  "R",
  "J",
  "Q",
  "K",
  "6",
  "M",
  "P",
  "L",
  "D",
  "F",
  "2",
  "N",
  "C",
  "E",
  "X",
  "O"
)
blood_dsr <- blood_dsr %>% prev_ICD_applier(procedures, "pPROC_", proc_codes)

###At least one coded previous UTI diagnosis in the last year
uti_key <- d_icd_diagnoses %>%
  filter(
    grepl("urinary tract infection", long_title, ignore.case = T) |
      grepl("acute pyelon", long_title, ignore.case = T) |
      (grepl("urinary catheter", long_title, ignore.case = T) &
        grepl("infec", long_title, ignore.case = T))
  )
hadm_key <- hadm %>% select(hadm_id, admittime)
uti_df <- diagnoses_raw %>%
  left_join(uti_key, by = c("icd_code", "icd_version")) %>%
  filter(!is.na(long_title)) %>%
  left_join(hadm_key, by = "hadm_id")
blood_dsr <- blood_dsr %>%
  care_event_assigner(
    uti_df,
    "(urin|pyelo|cath)",
    long_title,
    pUTI,
    "admittime",
    365,
    1
  )

###Presence of an outpatient provider ID
blood_dsr <- blood_dsr %>%
  mutate(provider_id = case_when(order_provider_id != "" ~ TRUE, TRUE ~ FALSE))

###Current inpatient specialty
serv_key <- services %>%
  select(hadm_id, curr_service) %>%
  distinct(hadm_id, .keep_all = T)
blood_dsr <- blood_dsr %>%
  left_join(serv_key, by = "hadm_id") %>%
  mutate(
    curr_service = case_when(
      is.na(curr_service) ~ "UNKNOWN",
      TRUE ~ curr_service
    )
  )
colnames(blood_dsr)
###At least one measured obese, underweight, or overweight BMI category in the last 3 years
categorise_bmi <- function(df) {
  df %>%
    filter(grepl("BMI", result_name)) %>%
    mutate(
      BMI_cat = case_when(
        as.numeric(result_value) >= 30 ~ "Obese",
        as.numeric(result_value) >= 25 &
          as.numeric(result_value) < 30 ~ "Overweight",
        as.numeric(result_value) >= 18.5 &
          as.numeric(result_value) < 25 ~ "Normal weight",
        as.numeric(result_value) < 18.5 ~ "Underweight"
      ),
      admittime = as.POSIXct(chartdate, format = '%Y-%m-%d %H:%M:%S')
    )
}
bmi <- categorise_bmi(omr)
bmi_categories <- c("Obese", "Underweight", "Overweight")
blood_dsr <- assign_bmi_events(blood_dsr, bmi, bmi_categories, 1095, 1) %>%
  ungroup()

###Observation frequency on day of test
obs <- poe %>%
  filter(order_subtype == "Vitals/Monitoring") %>%
  mutate(ordertime = as.Date(ordertime)) %>%
  group_by(subject_id, ordertime) %>%
  count(order_subtype) %>%
  arrange(desc(n)) %>%
  select(-order_subtype)
blood_dsr <- blood_dsr %>%
  mutate(ordertime = chartdate) %>%
  left_join(obs, by = c("subject_id", "ordertime")) %>%
  rename(ob_freq = "n") %>%
  mutate(
    ob_freq = case_when(is.na(ob_freq) ~ 0, TRUE ~ ob_freq),
    ob_freq = standardize(ob_freq)
  ) %>%
  select(-ordertime)

###Other specific previous care events
blood_dsr <- blood_dsr %>%
  care_event_assigner(poe, "cath", field_value, pCATH, "ordertime", 28) %>% ###At least one urinary catheter insertion in the last 28 days
  care_event_assigner(poe, "DNR", field_value, pDNR, "ordertime", 365) %>% ###At least one 'do not resuscitate' order in the last year
  care_event_assigner(poe, "Discharge", field_value, pDISC, "ordertime", 28) %>% ###At least one discharge from hospital in the last 28 days
  care_event_assigner(poe, "ICU", field_value, pICU, "ordertime", 28) %>% ###At least one intensive care admission in the last 28 days
  care_event_assigner(
    poe,
    "Psychiatry",
    field_value,
    pPsych,
    "ordertime",
    365
  ) %>% ###At least one psychiatry review in the last year
  care_event_assigner(
    poe,
    "Nephrostomy",
    field_value,
    pNeph,
    "ordertime",
    365
  ) %>% ###At least one nephrostomy insertion in the last year
  care_event_assigner(poe, "Surgery", field_value, pSURG, "ordertime", 365) %>% ###At least one surgical procedure in the last year
  care_event_assigner(poe, "Hydration", field_value, pHyd, "ordertime", 28) %>% ###At least one hydration order in the last 28 days
  care_event_assigner(poe, "NGT", field_value, pNGT, "ordertime", 28) %>% ###At least one nasogastric tube insertion in the last 28 days
  care_event_assigner(poe, "Chemo", field_value, pChemo, "ordertime", 28) %>% ###At least one administration of cancer chemotherapy in the last 28 days
  care_event_assigner(
    poe,
    "Nutrition consult",
    order_subtype,
    pNUTR,
    "ordertime",
    365
  ) %>% ###At least one nutrition consultation in the last year
  care_event_assigner(
    poe,
    "Physical Therapy",
    order_subtype,
    pPhysio,
    "ordertime",
    365
  ) %>% ###At least one physiotherapy consultation in the last year
  care_event_assigner(
    poe,
    "Restraints",
    order_subtype,
    pRestr,
    "ordertime",
    365
  ) %>% ###At least one requirement for restraints in the last year
  care_event_assigner(
    poe,
    "Occupational Therapy",
    order_subtype,
    pOT,
    "ordertime",
    365
  ) %>% ###At least one occupational therapy consultation in the last year
  care_event_assigner(poe, "Central TPN", order_subtype, pTPN, "ordertime", 365) ###At least one administration of total parenteral nutrition in the last year

###Revert charttime to admittime
blood_dsr <- blood_dsr %>% rename(admittime = "charttime")
blood_dsr <- blood_dsr %>% mutate(spec_type_desc = "BLOOD")

##Write CSVs
write_csv(blood_dsr, "blood_dsr_w_features.csv")
blood_dsr <- blood_dsr %>% distinct(hadm_id, .keep_all = T)
write_csv(blood_dsr, "preprocessed_blood.csv")

##Filter out urine cases with deterministic relationship with complicated UTI

###Diabetes
diabkey <- drgcodes %>%
  filter(grepl(
    "(DIAB)",
    description
  )) %>%
  select(hadm_id, description) %>%
  distinct(hadm_id)
hadm <- hadm %>% left_join(diabkey, by = "hadm_id")
urine_dsr <- urine_dsr %>%
  mutate(charttime = admittime) %>%
  diag_label(spec_type_desc, "DIAB", 1e4, pDIAB) %>%
  filter(!pDIAB)

###Male
urine_dsr <- urine_dsr %>% gender_assign(male, pats) %>% filter(!male)

###Immunosuppressants
cytotoxins <- c(
  "adalimumab",
  "aldesleukin",
  "alemtuzumab",
  "amsacrine",
  "arsenic trioxide",
  "asparaginase",
  "axitinib",
  "azacitidine",
  "azathioprine",
  "belatacept",
  "bendamustine",
  "bevacizumab",
  "bexarotene",
  "bleomycin",
  "blinatumomab",
  "bortezomib",
  "bosutinib",
  "brentuximab vedotin",
  "busulfan",
  "cabazitaxel",
  "cabozantinib",
  "canakinumab",
  "capecitabine",
  "carboplatin",
  "carfilzomib",
  "carmustine",
  "ceritinib",
  "certolizumab pegol",
  "chlorambucil",
  "cisplatin",
  "cladribine",
  "clofarabine",
  "crisantaspase",
  "cyclophosphamide",
  "cytarabine",
  "dacarbazine",
  "dactinomycin",
  "daratumumab",
  "dasatinib",
  "daunorubicin",
  "decitabine",
  "dexrazoxane",
  "dinutuximab",
  "docetaxel",
  "doxorubicin",
  "epirubicin",
  "eribulin",
  "estramustine",
  "etoposide",
  "fludarabine",
  "fluorouracil",
  "ganciclovir",
  "gemcitabine",
  "gemtuzumab ozogamicin",
  "golimumab",
  "hydroxycarbamide",
  "ibrutinib",
  "idarubicin",
  "ifosfamide",
  "imatinib",
  "infliximab",
  "inotuzumab ozogamicin",
  "ipilimumab",
  "irinotecan",
  "leflunomide",
  "lenalidomide",
  "lomustine",
  "melphalan",
  "mercaptopurine",
  "methotrexate",
  "mifamurtide",
  "mitomycin",
  "mitotane",
  "mitoxantrone",
  "mogamulizumab",
  "nelarabine",
  "nilotinib",
  "niraparib",
  "nivolumab",
  "obinutuzumab",
  "olaparib",
  "oxaliplatin",
  "paclitaxel",
  "palbociclib",
  "panobinostat",
  "pegaspargase",
  "peginterferon alfa",
  "pembrolizumab",
  "pemetrexed",
  "pentostatin",
  "pixantrone",
  "pomalidomide",
  "procarbazine",
  "raltitrexed",
  "ramucirumab",
  "regorafenib",
  "ribociclib",
  "rituximab",
  "ropeginterferon alfa",
  "rucaparib",
  "ruxolitinib",
  "sorafenib",
  "streptozocin",
  "sulfasalazine",
  "sunitinib",
  "talazoparib",
  "tegafur",
  "temozolomide",
  "temsirolimus",
  "thalidomide",
  "thiotepa",
  "tioguanine",
  "topotecan",
  "trabectedin",
  "trastuzumab",
  "trastuzumab",
  "deruxtecan",
  "trastuzumab",
  "emtansine",
  "treosulfan",
  "valganciclovir",
  "vinblastine",
  "vincristine",
  "vindesine",
  "vinorelbine"
)
cytotoxins <- str_to_title(str_to_lower(cytotoxins))
cytotoxics_key <- drugs %>%
  filter(drug %in% cytotoxins) %>%
  distinct(hadm_id) %>%
  mutate(Cytotoxic_agent = TRUE)
urine_dsr <- urine_dsr %>% cytotoxic_join() %>% filter(!Cytotoxic_agent)

##Standardised observations

###Observation key from ED datasets
staykey <- edstays %>%
  select(hadm_id, stay_id) %>%
  distinct(hadm_id, .keep_all = T) %>%
  filter(!is.na(hadm_id))
triage <- triage %>% left_join(staykey)
triagekey <- triage %>%
  select(
    hadm_id,
    temperature,
    heartrate,
    resprate,
    o2sat,
    sbp,
    dbp,
    acuity,
    chiefcomplaint
  ) %>%
  distinct(hadm_id, .keep_all = T)
urine_dsr <- urine_dsr %>% left_join(triagekey)
blood_dsr <- blood_dsr %>% left_join(triagekey)

###Filter to only those with triage obs available
urine_dsr <- urine_dsr %>% obs_filter()
blood_dsr <- blood_dsr %>% obs_filter()

urine_dsr <- urine_dsr %>% standardize_obs()
blood_dsr <- blood_dsr %>% standardize_obs()

###Adding presenting complaint variables
urine_dsr %>% count(chiefcomplaint) %>% arrange(desc(n)) %>% print(n = 200)

urine_dsr <- urine_dsr %>% pc_dummies() %>% select(-chiefcomplaint)
blood_dsr <- blood_dsr %>% pc_dummies() %>% select(-chiefcomplaint)

write_csv(urine_dsr, "preprocessed_urines8.csv")
write_csv(blood_dsr, "preprocessed_blood8.csv")
