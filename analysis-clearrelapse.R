library(reshape2)
library(lmerTest)
library(survival)
library(ggplot2)
library(dplyr)
library(survminer)
library(sjPlot)
library(EValue)
library(tableone)
library(perm)
library(rms)
setwd("~/Developer/canbind-depression-anxiety-coupling")
set.seed(235)

################################################################################
# LOAD DATA
################################################################################
# Load the clinical data, convert date formats, and calculate time to event 
load_clinical_data <- function() {
    data <- read.csv("data/clinical_data.csv", na.strings=c(""))
    data$eventtime <- data$ADY
    data$baselinedate <- as.Date(data$STARTDT)
    data$relapseenddate <- as.Date(data$ADT)
    data$relapseenddate - data$baselinedate
    data <- data %>% distinct(USUBJID, .keep_all = TRUE)

    # Merge the eating, anxiety and substance use disorders 
    data$eating_disorder <- apply(
      data[,c("anorexia", "bulimia", "binge_eating")], 
      MARGIN = 1, 
      FUN = function (x) { any(x == "Yes", na.rm=TRUE) })

    data$anxiety_disorder <- apply(data[,c(
      "panic_disorder_curr", "agoraphobia", 
      "social_phobia", "ocd", "ptsd", "gad")], 
      MARGIN = 1, 
      FUN = function (x) { any(x == "Yes", na.rm=TRUE) })

    data$substance_use <- apply(
      data[,c("drug", "etoh")], 
      MARGIN = 1, 
      FUN = function (x) { any(x == "Yes", na.rm=TRUE) })

    # Order the income levels 
    data$INCMLVL <- factor(
      data$INCMLVL, 
      levels=c(
        "LESS THAN $10, 000", 
        "$10, 000 - $24, 999", 
        "$25, 000 - $49, 999", 
        "$50, 000 - $74, 999", 
        "$75, 000 - $99, 999", 
        "$100,000 - $149,999", 
        "$150,000 - $199,999", 
        "$200, 000 OR MORE", 
        "PREFER NOT TO ANSWER", 
        "DON'T KNOW"
      ))
    
    # Order the work levels 
    data$EMPSTAT <- factor(
      data$EMPSTAT, 
      level = c(
        "RETIRED",
      "WORKING NOW", 
      "SELF-EMPLOYED", 
      "STUDENT",
      "KEEPING HOUSE",
      "CASUAL WORK", 
      "WORKING PART-TIME DUE TO DISABILITY",
      "ONLY TEMPORARILY LAID OFF, SICK LEAVE, OR MATERNITY LEAVE", 
      "LOOKING FOR WORK, UNEMPLOYED", 
      "NOT WORKING, NOT LOOKING FOR WORK", 
      "DISABLED, PERMANENTLY OR TEMPORARILY"
      )
    )

    # Order the education levels 
    data$EDULEVEL <- factor(
      data$EDULEVEL, 
      level = c(
        "10TH GRADE",
        "11TH GRADE",
        "HIGH SCHOOL GRADUATE",
        "SOME COLLEGE, NO DEGREE", 
        "ASSOCIATE DEGREE: OCCUPATIONAL, TECHNICAL, OR VOCATIONAL PROGRAM", 
        "ASSOCIATE DEGREE: ACADEMIC PROGRAM", 
        "BACHELOR'S DEGREE (E.G., BA, AB, BS, BBA)",
        "MASTER'S DEGREE (E.G., MA, MS, MENG, MED,MBA)", 
        "PROFESSIONAL SCHOOL DEGREE (E.G., MD, DDS, DVM, JD)"
      )
    )

    # Order job classes 
    data$JOBCLAS <- factor(
      data$JOBCLAS, 
      level=c(
        "SERVICE WORKER",
        "LABORER/HELPER", 
        "CRAFT WORKER", 
        "TECHNICIAN", 
        "PROFESSIONAL", 
        "SALES WORKER", 
        "ADMINISTRATIVE SUPPORT WORKER",
        "OPERATIVE",
        "OFFICIAL/MANAGER", 
        "NONE", 
        "UNKNOWN"
      ))

    return(data)
}

data <- load_clinical_data() 

# Load the qids and GAD7 data 
qids <- merge(
    read.csv("data/qids.csv"), 
    read.csv("data/gad7.csv"),
    by=c("USUBJID", "QSDTC")
)

qids$timestamp <- as.Date(qids$QSDTC)

# Remove the ambiguous relapsers 
ambig <- filter(
    data, 
    subject_id %in% c(
        "CBN_00004_010",
        "CBN_00003_015",
        "CBN_00004_029",
        "CBN_00003_020",
        "CBN_00002_014",
        "CBN_00001_013",
        "CBN_00001_016"
    )
)[,c("subject_id", "USUBJID")]

data <- filter(data, !(subject_id %in% ambig$subject_id))
qids <- filter(qids, !(USUBJID %in% ambig$USUBJID))

# Compute total qids score
qids$QIDS <- pmax(qids$QIDS0201, qids$QIDS0202, qids$QIDS0203, qids$QIDS0204) + 
  qids$QIDS0205 + pmax(qids$QIDS0206, qids$QIDS0207, qids$QIDS0208, qids$QIDS0209) + 
  qids$QIDS0210 + qids$QIDS0211 + qids$QIDS0212 + qids$QIDS0213 + qids$QIDS0214 + 
  pmax(qids$QIDS0215, qids$QIDS0216)
qids$GAD7 <- qids$GAD0101 + qids$GAD0102 + qids$GAD0103 + 
    qids$GAD0104 + qids$GAD0105 + qids$GAD0106 + qids$GAD0107

# Filter out qids/gad7 data that are not between baseline or relapse
# and calculate mean gad7 level
qids$include  = 0
data$gad7mean = NA
data$qidsmean = NA
for (s in unique(data$USUBJID)) {
    start_date <- data[data$USUBJID == s, "baselinedate"]
    end_date <- data[data$USUBJID == s, "relapseenddate"] - 30
    qids[
        (qids$USUBJID == s) & (qids$timestamp >= start_date) & (qids$timestamp < end_date), 
        "include"] <- 1
    
    data[data$USUBJID == s, "gad7mean"] <- mean(qids[
      (qids$USUBJID == s) & (qids$timestamp >= start_date) & (qids$timestamp < end_date), 
      "GAD7"], na.rm = TRUE)
    
    data[data$USUBJID == s, "qidsmean"] <- mean(qids[
      (qids$USUBJID == s) & (qids$timestamp >= start_date) & (qids$timestamp < end_date), 
      "QIDS"], na.rm = TRUE)
}

# Identify excluded subjects
excluded_subjects <- filter(aggregate(include ~ USUBJID, qids, FUN=sum), include == 0)
excluded_subjects_table <- merge(excluded_subjects, data, by="USUBJID")


# Filter the subjects from dataset
qids <- filter(qids, !(USUBJID %in% excluded_subjects$USUBJID))
qids <- filter(qids, !is.na(QIDS))
data <- filter(data, !(USUBJID %in% excluded_subjects$USUBJID))

################################################################################
# ESTIMATE THE DEPRESSION-ANXIETY COUPLING
################################################################################
set.seed(27243)

# Fit model to estimate effect of GAD7 on QIDS
m <- lmer(QIDS ~ GAD7 + (1 + GAD7|USUBJID), data=qids, 
  control=lmerControl(calc.derivs=FALSE, optCtrl = list(maxfun=100000)))
tab_model(m, file="tables/lmm-table-clearrelapse.doc")

# Plot individual-level random slopes 
set_theme(base = theme_bw())
pdf("figures/DACSre-clearrelapse.pdf", height=3, width=3)
p = plot_model(m, type="pred", pred.type="re", 
    terms=c("GAD7", "USUBJID"), ci.lvl=FALSE, line.size=0.1,
    dot.size=0.5, axis.lim = list(c(0, 21), c(0, 27)),
    show.data=TRUE, show.legend=FALSE, colors="black", 
    title="Patient-Level DACS")
p + theme(
  axis.title.x=element_text(colour="black"), 
  axis.text.x=element_text(colour="black"), 
  axis.title.y=element_text(colour="black"), 
  axis.text.y=element_text(colour="black"))
dev.off()

# Extract random effects
rand_eff <- dcast(as.data.frame(ranef(m)), grp ~ term, value.var="condval")
names(rand_eff) <- c("USUBJID", "Intercept", "GAD7")

write.csv(rand_eff, "results/dacs-random-effects-clearrelapse.csv")

# Merge random effects with relapse data 
data <- merge(data, rand_eff, by="USUBJID")

# Plot overall survival curve 
sfit0 <- survfit(Surv(eventtime, relapse) ~ 1, data = data)
pdf("figures/survival-all-clearrelapse.pdf", width=6, height=6)
p = plot(sfit0, conf.int=TRUE, 
         xlab="Days", ylab="Proportion Free of Relapse", 
         lwd=2, col="black", mark.time=TRUE, 
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
)
print(p)
dev.off()

# Create survival curves for DACS Groups
data$DACS <- cut(data$GAD7, quantile(data$GAD7), labels=c("Low", "Low", "Low", "High"))
sfit <- survfit(Surv(eventtime, relapse) ~ DACS, data = data)
pdf("figures/survival-dacs-clearrelapse.pdf", width=6, height=6)
p = plot(sfit, conf.int=TRUE, 
         xlab="Days", ylab="Proportion Free of Relapse", 
         lwd=2, col=c("gray", "black"), mark.time=TRUE, 
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
)
legend(1, 0.1, legend=c("Top 25% of Coupling Strength", "Bottom 75% of Coupling Strength"),
       col=c("black","gray"), lty=1, cex=0.8)
print(p)
dev.off()

################################################################################
# RUN THE COX PROPORTIONAL HAZARDS MODEL
################################################################################
# Fit Cox model
cm = coxph(Surv(eventtime, relapse) ~ 
    scale(GAD7) + 
    scale(gad7bl) + 
    scale(madrsbl) + 
    scale(gad7mean) + 
    scale(qidsmean), data = data)
summary(cm)
tab_model(cm, file="tables/coxmodel-clearrelapse.doc")
vif(cm)

cm_unscaled = coxph(Surv(eventtime, relapse) ~ 
    GAD7 + 
    gad7bl + 
    madrsbl + 
    gad7mean + 
    qidsmean, data = data)
summary(cm_unscaled)
tab_model(cm, cm_unscaled, file="tables/coxmodel-clearrelapse.doc")

# Cox model swapping the qids mean for intercept
cm2 = coxph(Surv(eventtime, relapse) ~ 
    scale(GAD7) + 
    scale(gad7bl) + 
    scale(madrsbl) + 
    scale(gad7mean) + 
    scale(Intercept), data = data)
summary(cm2)
tab_model(cm, cm2, file="tables/coxmodel-intercept-clearrelapse.doc")

# Examine Schoenfeld residuals
cm_z = cox.zph(cm)

pdf("figures/schoenfeld-clearrelapse.pdf", height=10, width=10)
p = ggcoxzph(cm_z)
print(p)
dev.off()

# Compute E-Value
eval = evalue(HR(exp(cm$coefficients)[1], rare=FALSE), 
       lo = exp(confint(cm))[1,1], 
       hi = exp(confint(cm))[1,2], true = 1)
eval

point_eval = eval[2,1]
ci_eval = eval[2,2]

