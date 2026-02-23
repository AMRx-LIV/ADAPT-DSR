#ADAPT-DSR MICROSIMULATION

##Load-in
blood_dsr_best <- read_csv("preprocessed_blood8_XGBoost_holdout_prob_best 2.csv")
urine_dsr_best <- read_csv("preprocessed_urines8_XGBoost_holdout_prob_best 3.csv")
drugs <- read_csv("drugs_clean.csv")
pyxis <- read_csv("pyxis.csv") %>% rename(drug="name")
edstays <- read_csv("edstays.csv")
micro <- read_csv("microbiologyevents.csv")

##Convert AMP I to R

###Classifications
urine_assign <- function(df,cutoff) {
  
  df %>% mutate(DSR_AMP = case_when(org_order=="Enterobacterales"&AMP=="S"  &
                                      probabilities>cutoff ~ "I",
                                    org_order=="Enterobacterales"&AMP=="S" &
                                      probabilities<=cutoff ~ "S",
                                    org_order=="Enterobacterales"&AMP=="R" ~ "R",
                                    TRUE ~ AMP),
                NA_AMP = case_when(org_order=="Enterobacterales"&AMP=="S" ~ NA,
                                   org_order=="Enterobacterales"&AMP=="I" ~ "I",
                                   org_order=="Enterobacterales"&AMP=="R" ~ "R",
                                   TRUE ~ AMP),
                I_AMP = case_when(org_order=="Enterobacterales"&AMP=="S" ~ "I",
                                  org_order=="Enterobacterales"&AMP=="I" ~ "I",
                                  org_order=="Enterobacterales"&AMP=="R" ~ "R",
                                   TRUE ~ AMP),
                S_AMP = case_when(org_order=="Enterobacterales"&AMP=="S" ~ "S",
                                  org_order=="Enterobacterales"&AMP=="I" ~ "I",
                                  org_order=="Enterobacterales"&AMP=="R" ~ "R",
                                   TRUE ~ AMP),
                true_pos=case_when(org_order=="Enterobacterales"&AMP=="S"&probabilities>=cutoff&
                                     complicated_vs_uncomplicated_uti==1 ~ TRUE,
                                   TRUE~FALSE),
                true_neg=case_when(org_order=="Enterobacterales"&AMP=="S"&probabilities<cutoff&
                                    complicated_vs_uncomplicated_uti==0 ~ TRUE,
                                  TRUE~FALSE),
                false_pos=case_when(org_order=="Enterobacterales"&AMP=="S"&probabilities>=cutoff&
                                      complicated_vs_uncomplicated_uti==0 ~ TRUE,
                                    TRUE~FALSE),
                false_neg=case_when(org_order=="Enterobacterales"&AMP=="S"&probabilities<cutoff&
                                      complicated_vs_uncomplicated_uti==1 ~ TRUE,
                                    TRUE~FALSE),
                amp_res=case_when(org_order=="Enterobacterales"&(AMP=="R"|AMP=="I") ~ TRUE,
                          TRUE~FALSE))
  
}

blood_assign <- function(df,cutoff) {
  
  df %>% mutate(DSR_AMP = case_when(org_order=="Enterobacterales"&AMP=="S"  &
                                      probabilities>cutoff ~ "I",
                                    org_order=="Enterobacterales"&AMP=="S" &
                                      probabilities<=cutoff ~ NA,
                                    org_order=="Enterobacterales"&AMP=="R" ~ "R",
                                    TRUE ~ AMP),
                DSR_comment = case_when(org_order=="Enterobacterales"&AMP=="S"  &
                                      probabilities>cutoff ~ NA,
                                      org_order=="Enterobacterales"&AMP=="S" &
                                      probabilities<=cutoff ~ "Absence of ampicillin resistance, but amoxicillin should be used at high dose with an adjunctive agent",
                                      org_order=="Enterobacterales"&AMP=="R" ~ "R",
                                    TRUE ~ AMP),
                NA_AMP = case_when(org_order=="Enterobacterales"&AMP=="S" ~ NA,
                                   org_order=="Enterobacterales"&AMP=="I" ~ "I",
                                   org_order=="Enterobacterales"&AMP=="R" ~ "R",
                                   TRUE ~ AMP),
                I_AMP = case_when(org_order=="Enterobacterales"&AMP=="S" ~ "I",
                                  org_order=="Enterobacterales"&AMP=="I" ~ "I",
                                  org_order=="Enterobacterales"&AMP=="R" ~ "R",
                                  TRUE ~ AMP),
                S_AMP = case_when(org_order=="Enterobacterales"&AMP=="S" ~ "S",
                                  org_order=="Enterobacterales"&AMP=="I" ~ "I",
                                  org_order=="Enterobacterales"&AMP=="R" ~ "R",
                                  TRUE ~ AMP),
                true_pos=case_when(org_order=="Enterobacterales"&AMP=="S"&probabilities>=cutoff&
                                     uti_vs_not_uti==1 ~ TRUE,
                                   TRUE~FALSE),
                true_neg=case_when(org_order=="Enterobacterales"&AMP=="S"&probabilities<cutoff&
                                     uti_vs_not_uti==0 ~ TRUE,
                                   TRUE~FALSE),
                false_pos=case_when(org_order=="Enterobacterales"&AMP=="S"&probabilities>=cutoff&
                                      uti_vs_not_uti==0 ~ TRUE,
                                    TRUE~FALSE),
                false_neg=case_when(org_order=="Enterobacterales"&AMP=="S"&probabilities<cutoff&
                                      uti_vs_not_uti==1 ~ TRUE,
                                    TRUE~FALSE),
                amp_res=case_when(org_order=="Enterobacterales"&(AMP=="R"|AMP=="I") ~ TRUE,
                                  TRUE~FALSE))
  
}

urine_dsr <- urine_dsr_best %>% urine_assign(0.5)
blood_dsr <- blood_dsr_best %>% blood_assign(0.5)

##Under and over-dosed patients in UTI
overunder <- function(df,func,dftype,neg_cons,pos_cons,right_cons,consequence) {
threshseq <- seq(0,1,0.025)
underdosed <- data.frame(matrix(nrow=0,ncol=4))
overdosed <- data.frame(matrix(nrow=0,ncol=4))
rightdose <- data.frame(matrix(nrow=0,ncol=4))
resistant <- data.frame(matrix(nrow=0,ncol=4))
unknown <- data.frame(matrix(nrow=0,ncol=4))
for (i in seq_along(threshseq)) {
  dsr <- df %>% func(threshseq[i])
  underdosed[i,1] <- threshseq[i]
  overdosed[i,1] <- threshseq[i]
  rightdose[i,1] <- threshseq[i]
  resistant[i,1] <- threshseq[i]
  unknown[i,1] <- threshseq[i]
  underdosed[i,2] <- (nrow(dsr %>% filter(false_neg))/nrow(dsr %>% filter(org_order=="Enterobacterales")))*100
  overdosed[i,2] <- (nrow(dsr %>% filter(false_pos))/nrow(dsr %>% filter(org_order=="Enterobacterales")))*100
  rightdose[i,2] <- (nrow(dsr %>% filter(true_pos|true_neg))/nrow(dsr %>% filter(org_order=="Enterobacterales")))*100
  resistant[i,2] <- (nrow(dsr %>% filter(amp_res))/nrow(dsr %>% filter(org_order=="Enterobacterales")))*100
  unknown[i,2] <- (nrow(dsr %>% filter(org_order=="Enterobacterales"&is.na(AMP)))/nrow(dsr %>% filter(org_order=="Enterobacterales")))*100
  underdosed[i,3] <- nrow(dsr %>% filter(false_neg))
  overdosed[i,3] <- nrow(dsr %>% filter(false_pos))
  rightdose[i,3] <- nrow(dsr %>% filter(true_pos|true_neg))
  resistant[i,3] <- nrow(dsr %>% filter(amp_res))
  unknown[i,3] <- nrow(dsr %>% filter(org_order=="Enterobacterales"&is.na(AMP)))
}
underdosed[,4] <- neg_cons
overdosed[,4] <- pos_cons
rightdose[,4] <- right_cons
resistant[,4] <- "Resistant"
unknown[,4] <- "Unknown"
dosedf <- data.frame(rbind(underdosed,overdosed,rightdose,resistant,unknown))
colnames(dosedf) <- c("Probability decision threshold","Percentage of Enterobacterales","N","Outcome")
dosedf$Outcome <- factor(dosedf$Outcome,levels=c("Unknown","Resistant",pos_cons,right_cons,neg_cons))
ggplot(dosedf,aes(x=`Probability decision threshold`,y=`Percentage of Enterobacterales`,group=Outcome,color=Outcome)) +
  geom_line() +
  ggtitle(glue("Enterobacterales {consequence} for {dftype} infection by\nADAPT-DSR in microsimulation"))+
  ylim(0,100)

dosedf <<- dosedf

}
overunder_2 <- function(df,func,dftype,neg_cons,pos_cons,right_cons,consequence) {
  threshseq <- seq(0,1,0.025)
  underdosed <- data.frame(matrix(nrow=0,ncol=4))
  overdosed <- data.frame(matrix(nrow=0,ncol=4))
  rightdose <- data.frame(matrix(nrow=0,ncol=4))
  for (i in seq_along(threshseq)) {
    dsr <- df %>% func(threshseq[i])
    underdosed[i,1] <- threshseq[i]
    overdosed[i,1] <- threshseq[i]
    rightdose[i,1] <- threshseq[i]
    underdosed[i,2] <- (nrow(dsr %>% filter(false_neg))/nrow(dsr %>% filter(org_order=="Enterobacterales"&(AMP=="S"|AMP=="I"))))*100
    overdosed[i,2] <- (nrow(dsr %>% filter(false_pos))/nrow(dsr %>% filter(org_order=="Enterobacterales"&(AMP=="S"|AMP=="I"))))*100
    rightdose[i,2] <- (nrow(dsr %>% filter(true_pos|true_neg))/nrow(dsr %>% filter(org_order=="Enterobacterales"&(AMP=="S"|AMP=="I"))))*100
    underdosed[i,3] <- nrow(dsr %>% filter(false_neg))
    overdosed[i,3] <- nrow(dsr %>% filter(false_pos))
    rightdose[i,3] <- nrow(dsr %>% filter(true_pos|true_neg))
  }
  underdosed[,4] <- neg_cons
  overdosed[,4] <- pos_cons
  rightdose[,4] <- right_cons
  dosedf <- data.frame(rbind(underdosed,overdosed,rightdose))
  colnames(dosedf) <- c("Probability decision threshold","Percentage of Enterobacterales","N","Outcome")
  dosedf$Outcome <- factor(dosedf$Outcome,levels=c(pos_cons,right_cons,neg_cons))
  ggplot(dosedf,aes(x=`Probability decision threshold`,y=`Percentage of Enterobacterales`,group=Outcome,color=Outcome)) +
    geom_line() +
    ggtitle(glue("Enterobacterales {consequence} for {dftype} infection by\nADAPT-DSR in microsimulation"))+
    ylim(0,100)
  
  dosedf2 <<- dosedf
  
}

urine_dsr_best %>% overunder(urine_assign,"urine","Under-dosed","Over-dosed","Correctly dosed","under- and over-dosed with amoxicillin")
urinedf <- dosedf
blood_dsr_best %>% overunder(blood_assign,"bloodstream","Over-treated","Under-treated","Correctly treated","under- and over-treated with amoxicillin +/- another agent")
blooddf <- dosedf

urine_dsr_best %>% overunder_2(urine_assign,"urine","Under-dosed","Over-dosed","Correctly dosed","under- and over-dosed with amoxicillin")
urinedf2 <- dosedf2
blood_dsr_best %>% overunder_2(blood_assign,"bloodstream","Over-treated","Under-treated","Correctly treated","under- and over-treated with amoxicillin +/- another agent")
blooddf2 <- dosedf2


##Plot
overplot <- function(df,bloodorurine) {

ggplot(df,aes(x=`Probability decision threshold`,y=`Percentage of Enterobacterales`,
                     group=Outcome,fill=Outcome)) +
  geom_area(alpha=0.6) +
  ylim(0,100)+
  ggtitle(glue("Treatment recommendation impact of ADAPT-DSR for {bloodorurine} cultures"))+
  theme_minimal()+
  scale_fill_viridis(option="D",discrete = T)
  
}

urinedf %>% overplot("urine")
blooddf %>% overplot("blood")
urinedf %>% filter(`Probability decision threshold`==0.5) %>% 
  select(-`Probability decision threshold`)
blooddf %>% filter(`Probability decision threshold`==0.5) %>% 
  select(-`Probability decision threshold`)
view(urinedf %>% filter(grepl("Correct",Outcome)))
urinedf %>% filter(Outcome=="Correctly dosed") %>% 
  arrange(desc(N))

##Nitrofurantoin additional analysis


ur_nit_sens <- urine_dsr_best %>% mutate(
  NIT = case_when(NIT=="I"~"S",TRUE~NIT),
  nit_pred = case_when(
    NIT=="S"&probabilities<=0.5~TRUE,
    is.na(NIT)~NA,
    TRUE~FALSE
  ),
  nit_usable = case_when(
    NIT=="S"&complicated_vs_uncomplicated_uti==0~TRUE,
    is.na(NIT)~NA,
    TRUE~FALSE
  ),
  nit_match = 
    case_when(
      NIT=="R"~"Resistant",
      nit_pred&nit_usable~"Correctly usable",
      nit_pred==FALSE&nit_usable==FALSE~"Correctly non-usable",
      is.na(NIT)~"Unknown",
      nit_pred==TRUE&nit_usable==FALSE~"Potential undertreatment",
      nit_pred==FALSE&nit_usable==TRUE~"Potential overtreatment"
    )
)

ur_nit_sens %>% count(nit_match) %>% 
  arrange(desc(n)) %>% 
  mutate(
    perc_undertreat=round(((n/(nrow(urine_dsr_best)))*100),1)
  )

##Descriptive AST data

preproc_ur <- read_csv("preprocessed_urines8.csv")
preproc_blood <- read_csv("preprocessed_blood8.csv")

abslist <- preproc_ur %>% select(PEN:MTR) %>% colnames()
abslist <- abslist[c(3,9,10,11,14,15,17,18,20,21,22,23,28,29,31,33,35,36,40,42,
          41,46,58)]
abslist
abcounter <- function(df,antibiotic,name) {

df %>% count(!!sym(antibiotic)) %>% 
  mutate(Antibiotic=ab_name(antibiotic)) %>% 
  relocate(Antibiotic,.before=1) %>% 
  mutate(n=glue("{n} ({round((n/nrow(df))*100,2)}%)")) %>% 
  arrange(!!sym(antibiotic)) %>% 
  mutate(Antibiotic = ifelse(row_number() == 1, Antibiotic, "··")) %>% 
  rename("Result"=2, !!sym(name):=3) 
  
}

sensdf <- data.frame(matrix(nrow=0,ncol=4))
colnames(sensdf) <- c("Antibiotic","Result","Urine","Blood")

for (i in seq_along(abslist)) {

urine <- preproc_ur %>% abcounter(abslist[i],"Urine")
blood <- preproc_blood %>% abcounter(abslist[i],"Blood")

urineblood <- left_join(urine,blood)

sensdf <- sensdf %>% rbind(urineblood)

}

write_csv(sensdf,"dsr_sensdf.csv")

##Prescription data
pyxis <- MIMER::clean_antibiotics(pyxis,drug_col=drug)
pyxis <- pyxis %>% filter(grepl("Amoxicillin",abx_name))

drugs <- drugs %>% rename(abx_name="ab_name") %>% filter(grepl("Amoxicillin",abx_name)) %>% 
  mutate(charttime=starttime) %>% 
  bind_rows(pyxis)

urine_dsr2 <- urine_dsr %>% filter(org_order=="Enterobacterales")
blood_dsr2 <- blood_dsr %>% filter(org_order=="Enterobacterales")

ur_drugs <- drugs %>% semi_join(urine_dsr2,by="subject_id")
bl_drugs <- drugs %>% semi_join(blood_dsr2,by="subject_id")
charttime_key <- micro %>% select(micro_specimen_id,charttime)

urdrugs <- urine_dsr2 %>% left_join(charttime_key,by="micro_specimen_id")
bldrugs <- blood_dsr2 %>% left_join(charttime_key,by="micro_specimen_id")
urdrugs <- urdrugs %>% distinct(micro_specimen_id,.keep_all = T)
bldrugs <- bldrugs %>% distinct(micro_specimen_id,.keep_all = T)

filterer <- function(df) {
  
  df %>% select(subject_id,charttime,DSR_AMP,NA_AMP,I_AMP,S_AMP,
                true_pos,true_neg,false_pos,false_neg,amp_res)
  
}

urdrugs <- urdrugs %>% filterer()
bldrugs <- bldrugs %>% filterer()

urdrugs <- urdrugs %>% bind_rows(ur_drugs)
bldrugs <- bldrugs %>% bind_rows(bl_drugs)

urdrugs <- urdrugs %>% semi_join(ur_drugs,by="subject_id")
bldrugs <- bldrugs %>% semi_join(bl_drugs,by="subject_id")

drug_beforeafter <- function(df){
  
  df %>% arrange(subject_id,charttime) %>%
    mutate(
      on_drug = 
        case_when(!is.na(DSR_AMP)&!is.na(lag(abx_name))&
                    subject_id==lag(subject_id)&
                    lag(stoptime)>charttime~lag(abx_name),
                  TRUE~NA),
      drug_7d_after=
        case_when(!is.na(DSR_AMP)&!is.na(lead(abx_name))&
                    subject_id==lead(subject_id)&
                    difftime(lead(starttime),charttime,units="days")<7~lead(abx_name),
                  TRUE~NA)
    )
  
}

urdrugs <- urdrugs %>% drug_beforeafter()
bldrugs <- bldrugs %>% drug_beforeafter()

ur_amox <- urdrugs %>% filter(!is.na(on_drug)|!is.na(drug_7d_after))
bl_amox <- bldrugs %>% filter(!is.na(on_drug)|!is.na(drug_7d_after))

amox_outcomes <- function(df) {
  
  df %>% mutate(outcome = 
                  case_when(
                    true_pos&on_drug=="Amoxicillin"~"tp_onamox",
                    true_neg&on_drug=="Amoxicillin"~"tn_onamox",
                    false_pos&on_drug=="Amoxicillin"~"fp_onamox",
                    false_neg&on_drug=="Amoxicillin"~"fn_onamox",
                    amp_res&on_drug=="Amoxicillin"~"res_onamox",
                    true_pos&drug_7d_after=="Amoxicillin"~"tp_7damox",
                    true_neg&drug_7d_after=="Amoxicillin"~"tn_7damox",
                    false_pos&drug_7d_after=="Amoxicillin"~"fp_7damox",
                    false_neg&drug_7d_after=="Amoxicillin"~"fn_7damox",
                    amp_res&drug_7d_after=="Amoxicillin"~"res_7damox",
                    true_pos&on_drug=="Amoxicillin/clavulanic acid"~"tp_onamc",
                    true_neg&on_drug=="Amoxicillin/clavulanic acid"~"tn_onamc",
                    false_pos&on_drug=="Amoxicillin/clavulanic acid"~"fp_onamc",
                    false_neg&on_drug=="Amoxicillin/clavulanic acid"~"fn_onamc",
                    amp_res&on_drug=="Amoxicillin/clavulanic acid"~"res_onamc",
                    true_pos&drug_7d_after=="Amoxicillin/clavulanic acid"~"tp_7damc",
                    true_neg&drug_7d_after=="Amoxicillin/clavulanic acid"~"tn_7damc",
                    false_pos&drug_7d_after=="Amoxicillin/clavulanic acid"~"fp_7damc",
                    false_neg&drug_7d_after=="Amoxicillin/clavulanic acid"~"fn_7damc",
                    amp_res&drug_7d_after=="Amoxicillin/clavulanic acid"~"res_7damc"
                  ))
  
}

ur_amox <- ur_amox %>% amox_outcomes()
bl_amox <- bl_amox %>% amox_outcomes()

ur_amox %>% count(outcome)
bl_amox %>% count(outcome)

