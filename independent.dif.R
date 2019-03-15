library(randomForest)
library("caret")
library(pROC)
independent_dif<-function(train,test,dir,filename){
  exclude=c("UTR5","UTR3","cds","exon_stop","constitutive_exon","internal_exon","long_exon","last_exon","last_exon_400bp","last_exon_sc400",
            "intron","pos_UTR5","pos_UTR3","pos_cds","pos_exons","length_UTR3","length_UTR5","length_cds","miR_targeted_genes")
  # #get feature and label
  print("Data processing...")
  name=colnames(train)
  index=which(!is.na(match(name,exclude)))
  if(length(index)>0){
    train=train[,-index]
  }
  
  #get train and test
  name=colnames(test)
  index=which(!is.na(match(name,exclude)))
  if(length(index)>0){
    test=test[,-index]
  }
  
  
  train.name=colnames(train)
  test.name=colnames(test)
  index=which(!is.na(match(train.name,test.name)))
  if(length(index)>0){
    train=train[,index]
  }
  
  
  index=which(!is.na(match(test.name,train.name)))
  if(length(index)>0){
    test=test[,index]
  }
  
  print("Predicting...")
  #train
  rf <- randomForest(Target_label ~ ., data=train)
  classification=predict(rf,newdata=test)
  pro=predict(rf,newdata=test,type="prob")
  TP <- as.numeric(sum(classification=="m6A_positive" & test$Target_label=="m6A_positive"))
  FP <- as.numeric(sum(classification=="m6A_positive" & test$Target_label=="m6A_negative"))
  TN <- as.numeric(sum(classification=="m6A_negative" & test$Target_label=="m6A_negative"))
  FN <- as.numeric(sum(classification=="m6A_negative" & test$Target_label=="m6A_positive"))
  TPR=TP/(TP+FN)
  FPR=TN/(TN+FP)
  MCC=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  F1=(2*TP)/(2*TP+FP+FN)
  ACC=(TP+TN)/(TP+TN+FP+FN)
  precision=TP/(TP+FP)
  recall=TP/(TP+FN)
  ROCArea=auc(roc(test$label,as.numeric(pro)))
  #prob=pro[,2]
  
  result=data.frame(TPR,FPR,F1,precision,ACC,recall,MCC,ROCArea)
  colnames(result)<-c("TPR","FPR","F-Measure","Precision","ACC","Recall","MCC","AUC")
  
  write.csv(result,paste(dir,paste(filename,".result.csv"),sep = '/'),row.names = F)
  write.csv(pro[,2],paste(dir,paste(filename,".probability.csv"),sep = '/'),row.names = F)
  
  return(result)
  
}