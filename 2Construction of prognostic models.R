
train_his <- read.csv("train_his.csv", row.names = 1)
outTab = data.frame()
pFilter = 0.05  
sighis = c("futime", "fustat")
for(i in colnames(train_his[,3:ncol(train_his)])){
  cox <- coxph(Surv(futime, fustat) ~ train_his[,i], data = train_his)
  coxSummary <- summary(cox)
  coxP <- coxSummary$coefficients[,"Pr(>|z|)"]
  outTab <- rbind(outTab,
                  cbind(id = i,
                        coef = coxSummary$coefficients[,"coef"],
                        HR = coxSummary$conf.int[,"exp(coef)"],
                        HR.95L = coxSummary$conf.int[,"lower .95"],
                        HR.95H = coxSummary$conf.int[,"upper .95"],
                        pvalue = coxSummary$coefficients[,"Pr(>|z|)"])
  )
  if(coxP < pFilter){
    sighis = c(sighis,i)
  }
}
write.csv(outTab, file = "outTab1.csv", row.names = T)
rownames(outTab) <- outTab[,1]
sig_his <- sighis[3:length(sighis)]
outTab_sig_his <- outTab[sig_his,]
write.csv(outTab_sig_his, file = "outTab_sig_his.csv", row.names = T)
uniSigExp_his = train_his[,sighis]
write.csv(uniSigExp_his,file = "train_his_unicox2.csv",row.names = T)
train_lasso <- read.csv("train_his_unicox2.csv")
rownames(train_lasso) <- train_lasso[,1]
train_lasso <- train_lasso[,2:ncol(train_lasso)]
x = as.matrix(train_lasso[,c(3:ncol(train_lasso))])
y = data.matrix(Surv(train_lasso$futime, train_lasso$fustat))
fit <- glmnet(x, y, family = "cox")
cvfit <- cv.glmnet(x, y, family = "cox")
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassohis = row.names(coef)[index]
lassohis = c("futime", "fustat", lassohis)
lassoSigExp_his = train_lasso[,lassohis]
write.csv(lassoSigExp_his,file = "lassoSigExp_his.csv",row.names = T)


lassoSigExp_xx <- read.csv("lassoSigExp_his.csv")
rownames(lassoSigExp_xx) <- lassoSigExp_xx[,1]
lassoSigExp_xx <- lassoSigExp_xx[,2:ncol(lassoSigExp_xx)]
multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp_xx)
multiCoxSum <- summary(multiCox)
C_index <- multiCoxSum$concordance
Se_C <- C_index[2]
C_up95 <- C_index[1] + 1.96 * Se_C
C_down95 <- C_index[1] - 1.96 * Se_C
C_index <- cbind(C_index[1], C_down95, C_up95)
C_index <- as.data.frame(C_index)
colnames(C_index) <- c("C_value", "lower .95", "upper .95")
write.csv(C_index, file = "C_index.csv", row.names = F)

outMultiTab <- data.frame()
outMultiTab <- cbind(
  coef = multiCoxSum$coefficients[,"coef"],
  HR = multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L = multiCoxSum$conf.int[,"lower .95"],
  HR.95H = multiCoxSum$conf.int[,"upper .95"],
  pvalue = multiCoxSum$coefficients[,"Pr(>|z|)"]
)
outMultiTab <- cbind(id = row.names(outMultiTab), outMultiTab)
write.csv(outMultiTab,file = "outMultiTab_his1.csv",row.names = T)
riskScore = predict(multiCox, type = "risk", newdata = lassoSigExp_xx)
medianTrainRisk = median(riskScore)
risk = as.vector(ifelse(riskScore>medianTrainRisk, "high", "low"))
trainRiskOut = cbind(id=rownames(lassoSigExp_xx, riskScore, risk),
                     cbind(lassoSigExp_xx, riskScore, risk))
write.csv(trainRiskOut,file = "trainRiskOut.csv",row.names = F)

test_his <- read.csv("test_his.csv")
rownames(test_his) <- test_his[,1]
lassoSig <- colnames(lassoSigExp_xx)
lassoSigExp_xxx <- test_his[,lassoSig]
write.csv(lassoSigExp_xxx, file = "lassoSigExp_xxx.csv", row.names = T)
c_index_test <- survConcordance(Surv(lassoSigExp_xxx$futime, lassoSigExp_xxx$fustat) ~ 
                                  predict(multiCox,lassoSigExp_xxx)) 
C_index_test <- c_index_test$concordance
Se_C <- c_index_test$std.err
C_up95 <- C_index_test + 1.96 * Se_C
C_down95 <- C_index_test - 1.96 * Se_C
C_index_test <- cbind(C_index_test, C_down95, C_up95)
C_index_test <- as.data.frame(C_index_test)
colnames(C_index_test) <- c("C_value_test", "lower .95", "upper .95")
write.csv(C_index_test, file = "C_index_test.csv", row.names = F)
riskScoreTest = predict(multiCox, type = "risk", newdata = lassoSigExp_xxx)
riskTest = as.vector(ifelse(riskScoreTest > medianTrainRisk, "high", "low"))
testRiskOut = cbind(id=rownames(cbind(lassoSigExp_xxx, riskScoreTest, riskTest)),
                    cbind(lassoSigExp_xxx, riskScoreTest, riskTest))
write.csv(testRiskOut,file = "testRiskOut1.csv",row.names = F)
diff = survdiff(Surv(futime, fustat) ~ risk, data = lassoSigExp_xx)
pValue = 1 - pchisq(diff$chisq, df = 1)
roc = survivalROC(Stime = lassoSigExp_xx$futime, status = lassoSigExp_xx$fustat, marker = riskScore,
                  predict.time = 3, method = "KM")
roc_TNR <- 1 - roc$FP
roc_TPR_TNR <- cbind(roc$TP, roc_TNR)
roc_TPR_TNR <- as.data.frame(roc_TPR_TNR)
colnames(roc_TPR_TNR) <- c("Sensitivity", "Specificity")
write.csv(roc_TPR_TNR, file = "roc_TRP_TNR.csv", row.names = F)
diffTest = survdiff(Surv(futime, fustat) ~ riskTest, data = lassoSigExp_xxx)
pValueTest = 1 - pchisq(diffTest$chisq, df = 1)
rocTest = survivalROC(Stime = lassoSigExp_xxx$futime, status = lassoSigExp_xxx$fustat, marker = riskScoreTest,
                      predict.time = 3, method = "KM")
rocTest_TNR <- 1 - rocTest$FP
rocTest_TPR_TNR <- cbind(rocTest$TP, rocTest_TNR)
rocTest_TPR_TNR <- as.data.frame(rocTest_TPR_TNR)
colnames(rocTest_TPR_TNR) <- c("Sensitivity", "Specificity")
write.csv(rocTest_TPR_TNR, file = "rocTest_TRP_TNR.csv", row.names = F)
P_value <- cbind(pValue, pValueTest)
P_value <- as.data.frame(P_value)
colnames(P_value) <- c("pValue_train", "pValue_test")
write.csv(P_value, file = "P_value.csv", row.names = F)
AUC_value<- cbind(roc$AUC, rocTest$AUC)
AUC_value <- as.data.frame(AUC_value)
colnames(AUC_value) <- c("roc$AUC", "rocTest$AUC")
write.csv(AUC_value, file = "AUC_value.csv", row.names = T)

rt = read.csv("trainRiskOut.csv")
diff = survdiff(Surv(futime, fustat) ~ risk, data = rt)
print(diff)
pValue = 1 - pchisq(diff$chisq, df = 1)
pValue = signif(pValue, 4)
pValue = format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)
pdf(file = "survivalTrain.pdf")
par(oma=c(0.5, 1, 0, 1), mar=c(5,5.3,4,2), font.lab=2, font.axis=2)
plot(fit,
     col = c("red", "blue"),
     xlab = "Time(year)",
     ylab = "Survival rate",
     lwd = 4, cex.main = 2.5, cex.lab = 2, cex.axis = 2, font = 2,
     main = paste("K-M curve of training set", sep = ""),
     mark.time = TRUE)
legend("topright",
       c("high risk", "low risk"),
       lwd = 4,
       text.font=2, text.width=2.3,cex=1.5,
       col = c("red", "blue"))
dev.off()

rt = read.csv("testRiskOut.csv")
diff = survdiff(Surv(futime, fustat) ~ riskTest, data = rt)
print(diff)
pValue = 1 - pchisq(diff$chisq, df = 1)
pValue = signif(pValue, 4)
pValue = format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ riskTest, data = rt)
summary(fit)
pdf(file = "survivalTest.pdf")
par(oma=c(0.5, 1, 0, 1), mar=c(5,5.3,4,2), font.lab=2, font.axis=2)
plot(fit,
     
     col = c("red", "blue"),
     xlab = "Time(year)",
     ylab = "Survival rate",
     lwd = 4, cex.main = 2.5, cex.lab = 2, cex.axis = 2, font = 2,
     main = paste("K-M curve of validation set", sep = ""),
     mark.time = TRUE)
legend("topright",
       c("high risk", "low risk"),
       lwd = 4,
       text.font=2, text.width=2.3,cex=1.5,
       col = c("red", "blue"))
dev.off()


pdf(file = "rocTrian.pdf")
par(oma=c(0.5, 1, 0, 1), mar=c(5,5.3,4,2), font.lab=2, font.axis=2)
roc = survivalROC(Stime = rt$futime, status = rt$fustat, marker = rt$riskScore,
                  predict.time = 3, method = "KM")
plot(roc$FP, roc$TP, type = "l", xlim = c(0,1), ylim = c(0,1), col = "red", 
     xlab = "False positive rate", ylab = "True positive rate",
     main = paste("ROC curve of training set"),
     lwd = 4, cex.main = 2.5, cex.lab = 2, cex.axis = 2, font = 2)
abline(0,1)
dev.off()
rt = read.csv("testRiskOut1.csv")
pdf(file = "rocTest.pdf")
par(oma=c(0.5, 1, 0, 1), mar=c(5,5.3,4,2), font.lab=2, font.axis=2)
roc = survivalROC(Stime = rt$futime, status = rt$fustat, marker = rt$riskScore,
                  predict.time = 3, method = "KM")
plot(roc$FP, roc$TP, type = "l", xlim = c(0,1), ylim = c(0,1), col = "red", 
     xlab = "False positive rate", ylab = "True positive rate",
     main = paste("ROC curve of validation set"),
     #main = paste("ROC curve (", "AUC = ", round(roc$AUC, 3),")"),
     lwd = 4, cex.main = 2.5, cex.lab = 2, cex.axis = 2, font = 2)
abline(0,1)
dev.off()

nomo_train <- read.csv("nomo_train.csv", row.names = 1)
dd <- datadist(nomo_train)
options(datadist="dd")
f2 <- psm(Surv(futime, fustat) ~ Age + Grade + FIGO_stage + RiskScore_CT + RiskScore_his, data = nomo_train, dist='lognormal')
med <- Quantile(f2)
surv <- Survival(f2) 

nom <- nomogram(f2, fun=function(x) med(lp=x),
                funlabel="Median Survival Time")
plot(nom)

nom <- nomogram(f2, fun=list(function(x) surv(365, x),
                             function(x) surv(1095, x),
                             function(x) surv(1825, x)),
                funlabel=c("1-year Survival Probability",
                           "3-year Survival Probability",
                           "5-year Survival Probability"))
par(cex = 1.5,lwd=2,font.lab=2, font.axis=2, font=2)
plot(nom, xfrac=.25,col.grid =gray(c(0.8,0.95)))
