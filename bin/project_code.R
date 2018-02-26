data1=read.table(file.choose(),header=T,sep=",")
data2<-as.matrix(data1)
colnames(data2)[1]<-"age"
data3<-data.frame(data2)
#################################hot.deck
h1<-hot.deck(data3,method='p.draw')
write.table(h1[1],"C:\\Users\\zhenm\\Documents\\STAT 504\\Stats504_Project\\Stats504_Project\\Data\\imp1.csv", sep =",")
#################################
data4<-data3[,1:10]
for(i in 1:10){
data4[,i]<-as.numeric(as.character(data4[,i]))
}

m1<-glm(act_demented~.,family=binomial(link='logit'),data=data4)
r1=rep(0,length())
predict(m1)
plot(m1)
summary(m1)

m2<-lm(act_demented~.,data=data4)
summary(m2)
attach(data4)
library(glmnet)
names(data4)
y=as.factor(act_demented)
x=data4[,1:12]
x$age=as.factor(x$sex)
x$race=as.factor(x$race)
glmmod <- glmnet(data4[,1:12],act_demented,alpha=1, family="binomial")

# Plot variable coefficients vs. shrinkage parameter lambda.
plot(glmmod, xvar="lambda")

x1=data4$braak
y1=data4$act_demented
prob1=1:7
for(i in 1:7){
  prob1[i]=mean(y1[x1==(i-1)])
  prob1[i]=log(prob1[i]/(1-prob1[i]))
}
plot(0:6,prob1)




