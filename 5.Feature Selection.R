# -------------------- Radiomics feature selection ------------------
cat("\014"); rm(list = ls()); options(warn = -1)
library(survival);library(caret);library(glmnet);library(survminer);library(timeROC); pacman::p_load("rms")

coxPfilter=0.05
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

hospital1_file = "ComBat_Hospital_I.csv"
hospital2_file  = "ComBat_Hospital_II.csv"; ex=1
hospital3_file = "ComBat_Hospital_III.csv"; ex=2
data1=read.csv(hospital1_file)
data2 =read.csv(hospital2_file)  
data3=read.csv(hospital3_file)
Event = 1
sum(data1$fustat)/nrow(data1)
sum(data2$fustat) /nrow(data2)
sum(data3$fustat)/nrow(data3)
cc = 19
data1_clin = data1[,c(1:cc)]; data1 = data1[,-c(4:cc)]
data2_clin  =  data2[,c(1:cc)]; data2  =  data2[,-c(4:cc)]
two_center_clin = rbind(data1_clin,data2_clin)
data3_clin = data3[,c(1:cc)]; data3 = data3[,-c(4:cc)]

for (j in 4:ncol(data1)) {
  index05 = which(data1[,j] < quantile(data1[,j],0.025));  index95 = which(data1[,j] > quantile(data1[,j],0.975)); data1[index05,j] = quantile(data1[,j],0.025);  data1[index95,j] = quantile(data1[,j],0.975)
  index05 = which(data2[,j]  < quantile( data2[,j],0.025));  index95 = which( data2[,j] > quantile(data2[,j],0.975));  data2[index05,j] = quantile(data2[,j],0.025);  data2[index95,j] = quantile(data2[,j],0.975)
  index05 = which(data3[,j] < quantile(data3[,j],0.025));  index95 = which(data3[,j] > quantile(data3[,j],0.975)); data3[index05,j] = quantile(data3[,j],0.025);  data3[index95,j] = quantile(data3[,j],0.975)}


for (v in 4:ncol(data1)) {
  Mean = mean(data1[,v]); SD = sd(data1[,v])
  data1[,v] = (data1[,v]-min(data1[,v]))/(max(data1[,v])-min(data1[,v]));   data2[,v]  = (data2[,v]-min(data2[,v]))/(max(data2[,v])-min(data2[,v]));   data3[,v] = (data3[,v]-min(data3[,v]))/(max(data3[,v])-min(data3[,v]))

}


# ----------------------- Uni_cox + LASSO-cox -----------------------
j=0; P=data.frame(train=numeric(), test=numeric(), val1=numeric(), val2=numeric(), seed=numeric())
AUC=data.frame(train=numeric(), test=numeric(), val1=numeric(), val2=numeric(), seed=numeric())
Score = 0
j = j+1; seed=1
set.seed(seed);  tow_center = rbind(data1,data2)
split = sample(nrow(tow_center), 0.7*nrow(tow_center))
train_data = tow_center[split, ];   test_data = tow_center[-split, ]

# ---------- 1. Uni_cox ------------
outUniTab=data.frame();   sigGenes=c("futime","fustat")
for(i in colnames(train_data[,4:ncol(train_data)])){cox <- coxph(Surv(futime, fustat) ~ train_data[,i], data = train_data); coxSummary = summary(cox); coxP=coxSummary$coefficients[,"Pr(>|z|)"]
if(coxP<coxPfilter){sigGenes=c(sigGenes,i);outUniTab=rbind(outUniTab,cbind(id=i,HR=round(coxSummary$conf.int[,"exp(coef)"],3),HR.95L=round(coxSummary$conf.int[,"lower .95"],3),HR.95H=round(coxSummary$conf.int[,"upper .95"],3),pvalue=round(coxSummary$coefficients[,"Pr(>|z|)"],3)))}};uniSigExp=train_data[,sigGenes]; uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)
length(sigGenes); if(length(sigGenes)<5){next}

# ---------- 2. lasso_cox ----------
x=as.matrix(train_data[,sigGenes[3:length(sigGenes)]])
y=data.matrix(Surv(train_data$futime, train_data$fustat))
set.seed(100); fit <- glmnet(x, y, family = "cox", maxit = 1000)
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000, type.measure = "C",nfolds = 10)
coef <- coef(fit, s = cvfit$lambda.min);index <- which(coef != 0); actCoef <- coef[index]; lassoGene=row.names(coef)[index]; lassoSigExp=uniSigExp[,c("futime", "fustat", lassoGene)]; lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp);geneCoef=cbind(Gene=lassoGene, Coef=actCoef); print(paste0('LASSO-COX feature num：', nrow(geneCoef)));if(nrow(geneCoef)<2){next}
multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
if(nrow(geneCoef)>20){multiCox=step(multiCox,direction = "both", trace = FALSE)}
multiCoxSum=summary(multiCox); multiCoxSum; AIC = min(multiCox[["anova"]][["AIC"]]); AIC
if (is.na(sum(multiCoxSum$coefficients[, "coef"]))){next}
print(paste0('Step_Cox feature num:', length(multiCoxSum$coefficients[, "coef"])))
outMultiTab=data.frame(cbind(coef=multiCoxSum$coefficients[,"coef"],HR=round(multiCoxSum$conf.int[,"exp(coef)"],3),HR.95L=round(multiCoxSum$conf.int[,"lower .95"],3),HR.95H=round(multiCoxSum$conf.int[,"upper .95"],3),pvalue=round(multiCoxSum$coefficients[,"Pr(>|z|)"],3)))
print('finish')
