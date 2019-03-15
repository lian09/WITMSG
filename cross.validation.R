library(randomForest)
library("caret")
library(pROC)
cross_validation <- function(data,dir,filename){

  #train
  n=length(names(data))
  folds<-createFolds(y=data$Target_label,k=10) #根据training的laber-Species把数据集切分成10等份
  TPR=c(length=10)
  FPR=c(length=10)
  MCC=c(length=10)
  F1=c(length=10)
  ACC=c(length=10)
  precision=c(length=10)
  recall=c(length=10)
  ROCArea=c(length=10)
  PRCArea=c(length=10)
  prob=c(length=nrow(data))
  for(i in 1:10){
    print(paste(i,"-fold"))
    index=1:nrow(data)
    index_train=sample(index[-folds[[i]]])
    train<-data[index_train,]
    index_test=sample(index[folds[[i]]])
    test<-data[index_test,]
    rf <- randomForest(Target_label ~ ., data=train)
    classification=predict(rf,newdata=test)
    pro=predict(rf,newdata=test,type="prob")
    TP <- as.numeric(sum(classification=="m6A_positive" & test$Target_label=="m6A_positive"))
    FP <- as.numeric(sum(classification=="m6A_positive" & test$Target_label=="m6A_negative"))
    TN <- as.numeric(sum(classification=="m6A_negative" & test$Target_label=="m6A_negative"))
    FN <- as.numeric(sum(classification=="m6A_negative" & test$Target_label=="m6A_positive"))
    TPR[i]=TP/(TP+FN)
    FPR[i]=TN/(TN+FP)
    MCC[i]=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    F1[i]=(2*TP)/(2*TP+FP+FN)
    ACC[i]=(TP+TN)/(TP+TN+FP+FN)
    precision[i]=TP/(TP+FP)
    recall[i]=TP/(TP+FN)
    ROCArea[i]=auc(roc(test$Target_label,pro[,2]))
    prob[index_test]=pro[,2]
  }
  result=data.frame(mean(TPR),mean(FPR),mean(F1),mean(precision),mean(ACC),mean(recall),mean(MCC),mean(ROCArea))
  colnames(result)<-c("TPR","FPR","F-Measure","Precision","ACC","Recall","MCC","AUC")
  file1=paste(dir,paste(filename,".result.csv"),sep = '/')
  write.csv(result,paste(dir,paste(filename,".result.csv"),sep = '/'),row.names = F)
  write.csv(pro[,2],paste(dir,paste(filename,".probability.csv"),sep = '/'),row.names = F)
  
  return(result)
}