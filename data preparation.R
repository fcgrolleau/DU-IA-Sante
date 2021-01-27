data_raw1 <- read.csv('mimic_1.csv')
data_raw2 <- read.csv('mimic_2.csv')
data_raw3 <- read.csv('mimic_3.csv')
prepdata <- data_raw3[,-1]

extra_var <- c("lactate_max", "creatinine_max", "ptt_max", "first_icu_stay", "los_icu")

prepdata2 <- merge(prepdata, data_raw2[, c("subject_id", extra_var)], by="subject_id")

round(sapply(prepdata2, function(x) mean(is.na(x)))*100,2)
round(sapply(data_raw2, function(x) mean(is.na(x)))*100,2)

par(mfrow=c(5,6))
for(i in which(sapply(prepdata2, is.numeric)==TRUE)) {
hist(prepdata2[,i], main=names(prepdata2)[i], xlab="")
        }        

prepdata2$gender <- factor(prepdata2$gender)
prepdata2$CURR_CAREUNIT_transfers <- factor(prepdata2$CURR_CAREUNIT_transfers)

library(dplyr)
prepdata2 <- filter(prepdata2, age <299 & weight<500 & creatinine_max<40)
prepdata2 <- prepdata2[,-which(names(prepdata2) %in% "creatinine_min")]

set.seed(1234657)
prepdata2 <- prepdata2[sample(nrow(prepdata2), 7500),]

library(mice)
library(pROC)
library(CalibrationCurves)
library(gam)
library(glmnet)

imp <- mice(prepdata2, maxit=0)
predM <- imp$predictorMatrix
nonpredictors <- c("subject_id","hadm_id","sapsii","sapsii_prob")

predM[,nonpredictors]<-0
#predM[rownames(predM) %in% nonpredictors]<-0

set.seed(123)
impdata<-mice(prepdata2, m=4, seed = 1, predictorMatrix = predM)
imp1 <- complete(impdata, 1)

set.seed(12352)
rand_n <- sample(nrow(prepdata2), 7500)
trainset <- imp1[rand_n[1:5000],]
testset <- imp1[rand_n[5001:7500],]

#write.csv(trainset, "icu_training.csv")
#write.csv(testset, "icu_testing.csv")

sapply(trainset, function(x) mean(is.na(x)))

par(mfrow=c(5,6))
for(i in which(sapply(trainset, is.numeric)==TRUE)) {
        hist(trainset[,i], main=names(trainset)[i], xlab="")
}     

dev.off()

plot(los_icu~lactate_max, col=rgb(0,0,1, alpha=.2), pch=19, cex=.1, data=trainset)
fit <- lm(los_icu~poly(lactate_max,5), data=trainset)
summary(fit)
lines(sort(trainset$lactate_max), predict(fit, newdata = data.frame(lactate_max=sort(trainset$lactate_max))), col="red", type="l", lwd=2)
lines(sort(trainset$lactate_max), predict(fit, newdata = data.frame(lactate_max=sort(trainset$lactate_max)), interval ="confidence")[,2], col="red", type="l", lty=2, lwd=2)
lines(sort(trainset$lactate_max), predict(fit, newdata = data.frame(lactate_max=sort(trainset$lactate_max)), interval ="confidence")[,3], col="red", type="l", lty=2, lwd=2)
lines(sort(trainset$lactate_max), predict(fit, newdata = data.frame(lactate_max=sort(trainset$lactate_max)), interval ="prediction")[,2], col="red", type="l", lty=2, lwd=1)
lines(sort(trainset$lactate_max), predict(fit, newdata = data.frame(lactate_max=sort(trainset$lactate_max)), interval ="prediction")[,3], col="red", type="l", lty=2, lwd=1)

par(mfrow=c(2,2))
plot(gam(los_icu~s(lactate_max,3), data=trainset), main="3 dl")
plot(gam(los_icu~s(lactate_max,5), data=trainset), main="5 dl")
plot(gam(los_icu~s(lactate_max,7), data=trainset), main="7 dl")
plot(gam(los_icu~s(lactate_max), data=trainset), main="dl déterminés automantiquement par leave one out crossvalidation ")

dev.off()
plot(gam(los_icu~s(weight), data=trainset))
plot(los_icu~weight, col=rgb(0,0,1, alpha=.2), pch=19, cex=.1, data=trainset)
fit <- lm(los_icu~poly(weight,2), data=trainset)
lines(sort(trainset$weight), predict(fit, newdata = data.frame(weight=sort(trainset$weight))), col="red", type="l", lwd=2)
lines(sort(trainset$weight), predict(fit, newdata = data.frame(weight=sort(trainset$weight)), interval ="confidence")[,2], col="red", type="l", lty=2, lwd=2)
lines(sort(trainset$weight), predict(fit, newdata = data.frame(weight=sort(trainset$weight)), interval ="confidence")[,3], col="red", type="l", lty=2, lwd=2)
lines(sort(trainset$weight), predict(fit, newdata = data.frame(weight=sort(trainset$weight)), interval ="prediction")[,2], col="red", type="l", lty=2, lwd=1)
lines(sort(trainset$weight), predict(fit, newdata = data.frame(weight=sort(trainset$weight)), interval ="prediction")[,3], col="red", type="l", lty=2, lwd=1)
summary(fit)

fit <- lm(los_icu~weight+lactate_max, data=trainset)
summary(fit)

fit <- lm(los_icu~poly(weight,5)+poly(lactate_max,5), data=trainset)
summary(fit)

plot(los_icu~age, col=rgb(0,0,1, alpha=.2), pch=19, cex=.1, data=trainset)
par(mfrow=c(2,2))
fit <- lm(los_icu~age, data = trainset)
plot(fit)
dev.off()

plot(gam(hospital_mortality~s(age), data=trainset))
fit <- glm(hospital_mortality~age, family = "binomial", data=trainset)
summary(fit)

plot(hospital_mortality~age, data=trainset, cex=.5, pch=3, col=rgb(0,0,1, alpha=.05))
lines(sort(trainset$age),
      predict(fit, newdata=data.frame(age=sort(trainset$age)), type="response"),
      type="l", col="green")
legend("right", legend="Predicted Probability", lty=1, col="green", cex=.7, box.lty=0)

plot(gam(hospital_mortality~s(creatinine_max), data=trainset))
plot(hospital_mortality~creatinine_max, data=trainset, cex=.5, pch=3, col=rgb(0,0,1, alpha=.01))
fit <- glm(hospital_mortality~poly(creatinine_max,2), family = "binomial", data=trainset)
lines(sort(trainset$creatinine_max),
      predict(fit, newdata=data.frame(creatinine_max=sort(trainset$creatinine_max)), type="response"),
      type="l", col="green")
legend("right", legend="Predicted Probability", cex=.7, lty=1, col="green", box.lty=0)

fit <- glm(hospital_mortality~age+poly(creatinine_max,2), family = "binomial", data=trainset)
summary(fit)

roc2 <- roc(fit$data$hospital_mortality, fit$fitted.values)
plot(roc2, print.auc=TRUE)
roc2

val.prob.ci.2(fit$fitted.values, fit$data$hospital_mortality, g = 5, col.ideal="blue", smooth = "loess")

numvar <- c("age_score",            "hr_score",                  "sysbp_score",            
"temp_score",             "pao2fio2_score",            "uo_score",
"bun_score",              "wbc_score",                 "potassium_score",
"sodium_score",           "bicarbonate_score",         "bilirubin_score",
"gcs_score",              "comorbidity_score",         "admissiontype_score",
"age",                    "weight",                    "bun_min",
"hemoglobin_max",         "height",                    "lactate_max", 
"creatinine_max",         "ptt_max") 


xvars <- as.matrix(trainset[, numvar])

continuous_outcome <- "los_icu"
ycont <- trainset[, continuous_outcome]
fit <- glmnet(xvars, ycont, alpha=1, standardize = FALSE)
plot(fit, xvar = "lambda", label = TRUE)

xvars_mean <- apply(xvars, 2, mean)
xvars_sd <- apply(xvars, 2, sd)

xvars_prep <- list()
for (i in 1:ncol(xvars)){
xvars_prep[[i]] <-(xvars[,i]-xvars_mean[i])/xvars_sd[i]
}
xvars_scaled <- sapply(xvars_prep, c)
#scale(as.matrix(trainset[, numvar]))

fit <- glmnet(xvars_scaled, ycont, alpha=1)
plot(fit, xvar = "lambda", label = TRUE)

colnames(xvars_scaled)[5]
colnames(xvars_scaled)[18]

cvfit <- cv.glmnet(xvars_scaled, ycont, type.measure = "mse", alpha=1, nfold=10)
plot(cvfit)

cvfit$lambda.min
cvfit$lambda.1se
coef(cvfit, s = "lambda.min")
predict(fit, newx = xvars_scaled[1:5,], type = "response", s = c(cvfit$lambda.min, cvfit$lambda.1se))


binary_outcome <- "hospital_mortality"
ybin <- trainset[, binary_outcome]

fit <- glmnet(xvars_scaled, ybin, family = "binomial", alpha=0)
plot(fit, label = TRUE)

cvfit <- cv.glmnet(xvars_scaled, ybin, family = "binomial", type.measure = "auc", alpha=0, nfold=10)
plot(cvfit)

coef(cvfit, s = "lambda.min")

predict(cvfit, newx = xvars_scaled[1:5,], s = "lambda.min", type = "response")

val.prob.ci.2(
predict(cvfit, newx = xvars_scaled, s = "lambda.min", type = "response"),
ybin,
g = 5, col.ideal="blue", smooth = "loess")

###
xvars_test <- as.matrix(testset[, numvar])

xvars_prep <- list()
for (i in 1:ncol(xvars_test)){
xvars_prep[[i]] <-(xvars_test[,i]-xvars_mean[i])/xvars_sd[i]
}
xvars_test_scaled_identically <- sapply(xvars_prep, c)

ybin_test <- testset[, binary_outcome]

preds_test <- predict(cvfit, newx = xvars_test_scaled_identically, s = "lambda.min", type = "response")

val.prob.ci.2(preds_test, ybin_test, g = 5, col.ideal="blue", smooth = "loess")

##
preimp <- prepdata2[1:1000,]

apply(preimp, 2, function(x) mean(is.na(x)))

set.seed(123)
imp <- mice(prepdata2, maxit=0)
predM <- imp$predictorMatrix

predM

nonpredictors <- c("subject_id","hadm_id","sapsii","sapsii_prob")

predM[,nonpredictors]<-0
set.seed(123)
impdata1<-mice(preimp, m=4, predictorMatrix = predM)
imp1 <- complete(impdata1, 1)

pool(with(impdata1,
                  glm(hospital_mortality~age+poly(creatinine_max,2), family = "binomial")))


