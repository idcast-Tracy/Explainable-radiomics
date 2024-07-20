# -------------- Kaplan Meier Curves ---------------- #
cat("\014"); rm(list = ls()); options(warn = -1); options(digits=3);library(survival); library(survminer)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); cont=0; cut.off=''

K = 3; filename = c('risk-train.csv','risk-test.csv','risk-val.csv')[K]
mydata <- read.csv(filename, encoding = 'gbk'); names(mydata)[1]='p_ID'; colnames(mydata)
if (max(mydata$futime)> 365) {mydata$futime =mydata$futime/365} # Day -> Year
if (max(mydata$futime)> 30)  {mydata$futime =mydata$futime/ 12} # Mon -> Year


if(K==3){maxfutime=4}else {maxfutime=8}
Break.time = ifelse(max(mydata$futime) < 7, 1, 2)
variate = 'riskScore'
edata = data.frame(p_ID=mydata$p_ID, futime=mydata$futime, fustat=mydata$fustat, Group=mydata[, variate],ICT=mydata$ICT)

# ------- Cox: HR 95%CI -------------
Coxph <- coxph(Surv(as.numeric(edata$futime), as.numeric(edata$fustat)) ~ edata$Group)
HR = paste('HR =', round(exp(coef(Coxph)),3), '(',round(exp(confint(Coxph))[1],3), '-', round(exp(confint(Coxph))[2],3),')'); HR
# ------- Cox: C-Index -------------
library(Hmisc); cdex0 <- round(rcorr.cens(as.matrix(mydata[variate]), Surv(mydata$futime, mydata$fustat)),3); c_train <- 1-cdex0; vn_train_upper <- round(cdex0[3]/2*1.96+c_train, 3); vn_train_lower <- round(c_train-cdex0[3]/2*1.96, 3)
Ci = paste0(' C-index = ',c_train[1], ' ',vn_train_lower[1], '~', vn_train_upper[1]); Ci

if (nrow(unique(edata['Group']))>5){cont=1
res.cut = surv_cutpoint(edata, time = "futime",  event = "fustat",  variables = "Group"); summary(res.cut) 
cut.off = as.numeric(summary(res.cut)[1])
if (K!=1){cut.off = 0.527}
print(paste0('cutoff = ', round(cut.off, 3)))
edata['Group']=as.vector(ifelse(edata['Group']>cut.off,"High-Risk"," Low-Risk"))
}

# ---------------------------survival plot------------------------------
fit = survfit(Surv(edata$futime, edata$fustat) ~ edata$Group)
pdf(paste0(variate, ' ', strsplit(filename, "\\.")[[1]][1], '  cutoff = ', cut.off, '  ', HR, Ci, ".pdf"), width=4.5, height=5.35, onefile=FALSE)
ggsurvplot(fit, data = edata, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
           palette = c( '#2274C6','#FF3232'), 
           xlab = "Time (Year)",            # customize X axis label.
           ggtheme = theme_survminer(),      # customize plot and risk table with a theme.
           risk.table.y.text.col = T,        # colour risk table text annotations.
           risk.table.height = 0.32,         # the height of the risk table
           risk.table.y.text = FALSE,        # show bars instead of names in text annotations in legend of risk table.
           # risk.table.col="Group",
           ncensor.plot = FALSE,             # plot the number of censored subjects at futime t
           pval.method = TRUE,
           conf.int.style = "ribbon",          # ribbon, customize style of confidence intervals
           legend = 'top', 
           legend.title = "",
           xlim = c(0, maxfutime),
           break.x.by = Break.time) 
dev.off()


print(paste0('variable: ', variate)); if (cont==1){print(paste0('Continuous value, cutoff = ', cut.off))}; HR; Ci
print('finish')
"ALL_people"; quantile(mydata$futime)[c(3,2,4)]
unique(edata$Group)[1]; quantile(edata[which(edata$Group==unique(edata$Group)[1]), 'futime'])[c(3,2,4)]
unique(edata$Group)[2]; quantile(edata[which(edata$Group==unique(edata$Group)[2]), 'futime'])[c(3,2,4)]
