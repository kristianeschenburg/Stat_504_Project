data1=read.table(file.choose(),header=T,sep=",")
data2<-as.matrix(data1)
colnames(data2)[1]<-"age"
data3<-data.frame(data2)
h1<-hot.deck(data3,method='p.draw')























