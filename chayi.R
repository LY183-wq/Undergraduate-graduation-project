

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#OR值
#首先对保守性数据，进行差异表达分析
setwd("C:/Users/jiang/Desktop/chayitongji/conservations")
data <- read.table('Reviewed_provean_score.tsv',sep="\t",header=TRUE)

data1 <- cbind(data$type,data$conservations)
#讨论四种突变类型与保守性之间的相关性

#Disrupting
a <- nrow(data[data$type=='Disrupting' & data$conservations==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Disrupting' & data$conservations==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Disrupting' & data$conservations==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Disrupting' & data$conservations==0,])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR <- result$estimate
p <- result$p.value

#Decreasing
a <- nrow(data[data$type=='Decreasing' & data$conservations==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Decreasing' & data$conservations==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Decreasing' & data$conservations==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Decreasing' & data$conservations==0,])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR1 <- result$estimate
p1 <- result$p.value

#Disrupting
a <- nrow(data[data$type=='Increasing' & data$conservations==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Increasing' & data$conservations==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Increasing' & data$conservations==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Increasing' & data$conservations==0,])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR2 <- result$estimate
p2 <- result$p.value

#Disrupting
a <- nrow(data[data$type=='No effect' & data$conservations==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='No effect' & data$conservations==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'No effect' & data$conservations==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'No effect' & data$conservations==0,])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR3 <- result$estimate
p3 <- result$p.value



##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#t检验
setwd("C:/Users/jiang/Desktop/chayitongji/conservations")
data1 <- read.table('Reviewed_provean_score.tsv',sep="\t",header=TRUE)
data=data1[,c("conservations","provean_score","type")]
a=data[data$type=='Disrupting',]
b=data[data$type=='Decreasing',]
c=data[data$type=='Increasing',]
d=data[data$type=='No effect',]
a1=a[a$conservations==1,]
a2=a[a$conservations==0,]
b1=b[b$conservations==1,]
b2=b[b$conservations==0,]
c1=c[c$conservations==1,]
c2=c[c$conservations==0,]
d1=d[d$conservations==1,]
d2=d[d$conservations==0,]

#以保守性得分为数据来进行计算
result1=t.test(a1$provean_score,a2$provean_score,var.equal = FALSE)
result2=t.test(b1$provean_score,b2$provean_score,var.equal = FALSE)
result3=t.test(c1$provean_score,c2$provean_score,var.equal = FALSE)
result4=t.test(d1$provean_score,d2$provean_score,var.equal = FALSE)


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#wilcoxon检验
setwd("C:/Users/jiang/Desktop/chayitongji/conservations")
data1 <- read.table('Unreviewed_provean_score.tsv',sep="\t",header=TRUE)
data=data1[,c("conservations","provean_score","type")]
a=data[data$type=='Disrupting',]
b=data[data$type=='Decreasing',]
c=data[data$type=='Increasing',]
d=data[data$type=='No effect',]
a1=a[a$conservations==1,]
a2=a[a$conservations==0,]
b1=b[b$conservations==1,]
b2=b[b$conservations==0,]
c1=c[c$conservations==1,]
c2=c[c$conservations==0,]
d1=d[d$conservations==1,]
d2=d[d$conservations==0,]


result1=wilcox.test(a1$provean_score,a2$provean_score)
result2=wilcox.test(b1$provean_score,b2$provean_score)
result3=wilcox.test(c1$provean_score,c2$provean_score)
result4=wilcox.test(d1$provean_score,d2$provean_score)

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#卡方检验
setwd("C:/Users/jiang/Desktop/chayitongji/conservations")
data <- read.table('Reviewed_provean_score.tsv',sep="\t",header=TRUE)

#讨论四种突变类型与保守性之间的相关性

#Disrupting
a <- nrow(data[data$type=='Disrupting' & data$conservations==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Disrupting' & data$conservations==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Disrupting' & data$conservations==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Disrupting' & data$conservations==0,])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
p <- result$p.value

#Decreasing
a <- nrow(data[data$type=='Decreasing' & data$conservations==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Decreasing' & data$conservations==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Decreasing' & data$conservations==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Decreasing' & data$conservations==0,])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
p1 <- result$p.value

#Disrupting
a <- nrow(data[data$type=='Increasing' & data$conservations==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Increasing' & data$conservations==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Increasing' & data$conservations==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Increasing' & data$conservations==0,])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
p2 <- result$p.value

#Disrupting
a <- nrow(data[data$type=='No effect' & data$conservations==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='No effect' & data$conservations==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'No effect' & data$conservations==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'No effect' & data$conservations==0,])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
p3 <- result$p.value










#######################################################################################################################################




##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#OR值
#首先对保守性数据，进行差异表达分析
setwd("C:/Users/jiang/Desktop/PNG1")
data <- read.table('Unreviewed_dssp2.tsv',sep="\t",header=TRUE)



#Disrupting
a <- nrow(data[data$type=='Disrupting' & data$residue_location=='Core',])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Disrupting' & data$residue_location=='Surface',])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Disrupting' & data$residue_location=='Core',])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Disrupting' & data$residue_location=='Surface',])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR <- result$estimate
p <- result$p.value

#Decreasing
a <- nrow(data[data$type=='Decreasing' & data$residue_location=='Core',])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Decreasing' & data$residue_location=='Surface',])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Decreasing' & data$residue_location=='Core',])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Decreasing' & data$residue_location=='Surface',])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR1 <- result$estimate
p1 <- result$p.value

#Disrupting
a <- nrow(data[data$type=='Increasing' & data$residue_location=='Core',])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Increasing' & data$residue_location=='Surface',])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Increasing' & data$residue_location=='Core',])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Increasing' & data$residue_location=='Surface',])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR2 <- result$estimate
p2 <- result$p.value

#Disrupting
a <- nrow(data[data$type=='No effect' & data$residue_location=='Core',])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='No effect' & data$residue_location=='Surface',])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'No effect' & data$residue_location=='Core',])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'No effect' & data$residue_location=='Surface',])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR3 <- result$estimate
p3 <- result$p.value



##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#t检验
setwd("C:/Users/jiang/Desktop/PNG1")
data <- read.table('Unreviewed_dssp2.tsv',sep="\t",header=TRUE)
a=data[data$type=='Disrupting',]
b=data[data$type=='Decreasing',]
c=data[data$type=='Increasing',]
d=data[data$type=='No effect',]
a1=a[a$residue_location=='Core',]
a2=a[a$residue_location=='Surface',]
b1=b[b$residue_location=='Core',]
b2=b[b$residue_location=='Surface',]
c1=c[c$residue_location=='Core',]
c2=c[c$residue_location=='Surface',]
d1=d[d$residue_location=='Core',]
d2=d[d$residue_location=='Surface',]

#以保守性得分为数据来进行计算
result1=t.test(a1$dssp_acc,a2$dssp_acc,var.equal = FALSE)
result2=t.test(b1$dssp_acc,b2$dssp_acc,var.equal = FALSE)
result3=t.test(c1$dssp_acc,c2$dssp_acc,var.equal = FALSE)
result4=t.test(d1$dssp_acc,d2$dssp_acc,var.equal = FALSE)
result1$p.value
result2$p.value
result3$p.value
result4$p.value

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#wilcoxon检验
setwd("C:/Users/jiang/Desktop/PNG1")
data <- read.table('Reviewed_dssp2.tsv',sep="\t",header=TRUE)
a=data[data$type=='Disrupting',]
b=data[data$type=='Decreasing',]
c=data[data$type=='Increasing',]
d=data[data$type=='No effect',]
a1=a[a$residue_location=='Core',]
a2=a[a$residue_location=='Surface',]
b1=b[b$residue_location=='Core',]
b2=b[b$residue_location=='Surface',]
c1=c[c$residue_location=='Core',]
c2=c[c$residue_location=='Surface',]
d1=d[d$residue_location=='Core',]
d2=d[d$residue_location=='Surface',]


result1=wilcox.test(a1$dssp_acc,a2$dssp_acc,var.equal = FALSE)
result2=wilcox.test(b1$dssp_acc,b2$dssp_acc,var.equal = FALSE)
result3=wilcox.test(c1$dssp_acc,c2$dssp_acc,var.equal = FALSE)
result4=wilcox.test(d1$dssp_acc,d2$dssp_acc,var.equal = FALSE)
result1$p.value
result2$p.value
result3$p.value
result4$p.value
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#卡方检验
setwd("C:/Users/jiang/Desktop/PNG1")
data <- read.table('Unreviewed_dssp2.tsv',sep="\t",header=TRUE)



#Disrupting
a <- nrow(data[data$type=='Disrupting' & data$residue_location=='Core',])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Disrupting' & data$residue_location=='Surface',])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Disrupting' & data$residue_location=='Core',])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Disrupting' & data$residue_location=='Surface',])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
result$p.value

#Decreasing
a <- nrow(data[data$type=='Decreasing' & data$residue_location=='Core',])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Decreasing' & data$residue_location=='Surface',])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Decreasing' & data$residue_location=='Core',])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Decreasing' & data$residue_location=='Surface',])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
result$p.value

#Disrupting
a <- nrow(data[data$type=='Increasing' & data$residue_location=='Core',])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Increasing' & data$residue_location=='Surface',])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Increasing' & data$residue_location=='Core',])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Increasing' & data$residue_location=='Surface',])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
result$p.value

#Disrupting
a <- nrow(data[data$type=='No effect' & data$residue_location=='Core',])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='No effect' & data$residue_location=='Surface',])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'No effect' & data$residue_location=='Core',])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'No effect' & data$residue_location=='Surface',])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
result$p.value









#######################################################################################################################################




##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#OR值
#首先对保守性数据，进行差异表达分析
setwd("C:/Users/jiang/Desktop/PNG1")
data <- read.table('Unreviewed_scannet2.tsv',sep="\t",header=TRUE)



#Disrupting
a <- nrow(data[data$type=='Disrupting' & data[,11]==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Disrupting' & data[,11]==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Disrupting' & data[,11]==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Disrupting' & data[,11]==0,])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR <- result$estimate
p <- result$p.value

#Decreasing
a <- nrow(data[data$type=='Decreasing' & data[,11]==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Decreasing' & data[,11]==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Decreasing' & data[,11]==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Decreasing' & data[,11]==0,])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR1 <- result$estimate
p1 <- result$p.value

#Disrupting
a <- nrow(data[data$type=='Increasing' & data[,11]==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Increasing' & data[,11]==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Increasing' & data[,11]==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Increasing' & data[,11]==0,])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR2 <- result$estimate
p2 <- result$p.value

#Disrupting
a <- nrow(data[data$type=='No effect' & data[,11]==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='No effect' & data[,11]==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'No effect' & data[,11]==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'No effect' & data[,11]==0,])    #既不是Disrupting，又不是保守性残基
result <- fisher.test(matrix(c(a,b,c,d),nrow=2))
OR3 <- result$estimate
p3 <- result$p.value



##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#t检验
setwd("C:/Users/jiang/Desktop/PNG1")
data <- read.table('Unreviewed_scannet2.tsv',sep="\t",header=TRUE)
a=data[data$type=='Disrupting',]
b=data[data$type=='Decreasing',]
c=data[data$type=='Increasing',]
d=data[data$type=='No effect',]
a1=a[a[,11]==1,]
a2=a[a[,11]==0,]
b1=b[a[,11]==1,]
b2=b[a[,11]==0,]
c1=c[a[,11]==1,]
c2=c[a[,11]==0,]
d1=d[a[,11]==1,]
d2=d[a[,11]==0,]

#以保守性得分为数据来进行计算
result1=t.test(a1$probability,a2$probability,var.equal = FALSE)
result2=t.test(b1$probability,b2$probability,var.equal = FALSE)
result3=t.test(c1$probability,c2$probability,var.equal = FALSE)
result4=t.test(d1$probability,d2$probability,var.equal = FALSE)
result1$p.value
result2$p.value
result3$p.value
result4$p.value

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#wilcoxon检验
setwd("C:/Users/jiang/Desktop/PNG1")
data <- read.table('Unreviewed_scannet2.tsv',sep="\t",header=TRUE)
a=data[data$type=='Disrupting',]
b=data[data$type=='Decreasing',]
c=data[data$type=='Increasing',]
d=data[data$type=='No effect',]
a1=a[a[,11]==1,]
a2=a[a[,11]==0,]
b1=b[a[,11]==1,]
b2=b[a[,11]==0,]
c1=c[a[,11]==1,]
c2=c[a[,11]==0,]
d1=d[a[,11]==1,]
d2=d[a[,11]==0,]


result1=wilcox.test(a1$probability,a2$probability,var.equal = FALSE)
result2=wilcox.test(b1$probability,b2$probability,var.equal = FALSE)
result3=wilcox.test(c1$probability,c2$probability,var.equal = FALSE)
result4=wilcox.test(d1$probability,d2$probability,var.equal = FALSE)
result1$p.value
result2$p.value
result3$p.value
result4$p.value
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
#卡方检验
setwd("C:/Users/jiang/Desktop/PNG1")
data <- read.table('Unreviewed_scannet2.tsv',sep="\t",header=TRUE)



#Disrupting
a <- nrow(data[data$type=='Disrupting' & data[,11]==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Disrupting' & data[,11]==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Disrupting' & data[,11]==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Disrupting' & data[,11]==0,])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
result$p.value

#Decreasing
a <- nrow(data[data$type=='Decreasing' & data[,11]==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Decreasing' & data[,11]==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Decreasing' & data[,11]==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Decreasing' & data[,11]==0,])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
result$p.value

#Disrupting
a <- nrow(data[data$type=='Increasing' & data[,11]==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='Increasing' & data[,11]==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'Increasing' & data[,11]==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'Increasing' & data[,11]==0,])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
result$p.value

#Disrupting
a <- nrow(data[data$type=='No effect' & data[,11]==1,])    #既是Disrupting，又是保守性残基
b <- nrow(data[data$type=='No effect' & data[,11]==0,])    #是Disrupting，不是保守性残基
c <- nrow(data[data$type != 'No effect' & data[,11]==1,])    #不是Disrupting，但是保守性残基
d <- nrow(data[data$type != 'No effect' & data[,11]==0,])    #既不是Disrupting，又不是保守性残基
result <- chisq.test(matrix(c(a,b,c,d),nrow=2))
result$p.value


find /data/liuyang/bishe_data/conservations/Reviewed/bin/Unreviewed -type f -newermt "2023-05-06 20:00:00" ! -newermt "2023-05-8 8:00:00"
find . -type f -newermt "2022-05-06 21:00:00" ! -newermt "2023-05-37 8:40:00" > result1.txt


#################################################################################
#PremPS结果方差分析，检验四个之间的差异
setwd('C://Users//jiang//Desktop//PS')
a <-read.table('Reviewed_PremPS_result.tsv',sep="\t",header=TRUE)
b <-read.table('Unreviewed_PremPS_result.tsv',sep="\t",header=TRUE)

a1 <-a[,c(3,20)]
b1 <-b[,c(3,20)]


data<-rbind(a1,b1)
colnames(data) <- c("type", "PremPS")
m=data[data$type=='Disrupting',]
n=data[data$type=='Decreasing',]
p=data[data$type=='Increasing',]
q=data[data$type=='No effect',]

data1=c(m$PremPS,n$PremPS,p$PremPS,q$PremPS)
PremPS <-data.frame(X=data1,A=factor(rep(1:4,c(2856,2160,265,4042))))
PremPS.aov <- aov(X~A,data=PremPS)
summary(PremPS.aov)



#########################################################################################

#R绘制ROC曲线
library(QuantPsyc)
library(car)
library('gplots')
library('ROCR')
library(caTools)
library(PRROC)
library('mccr')
library('pROC')

setwd('C://Users//jiang//Desktop//PNG')
a <-read.table('PremPS_1.tsv',sep="\t",header=TRUE)
#### Fig S3. ROC curve for predicting deleterious mutations by applying PremPDI on the training set of “Prempdi” ####
roc <- paste('C://Users//jiang//Desktop//PNG//', 'R_PremPS_ROC1.png', sep='')
pred <- prediction(a$PremPS, a$true)
perf <- performance(pred,"tpr","fpr")
auc1 <- performance(pred, "auc")@y.values[[1]]
png(file=roc,width=15,height=15,units="cm",res=300) 
par(cex.axis = 1.5)
plot(perf,col='black',lwd=4,xlab="", ylab="",cex.main=1.5)
mtext('True positive rate',side=2,line=2.5,cex=2.3,col="black")
mtext('False positive rate',side=1,line=3,cex=2.3,col="black")
auc1_1=paste('AUC = ' , as.character(round(auc1,2)))
text(0.7,0.2,labels=auc1_1,,cex=1.5)
# ROC optimal threshold 
cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], tpr=perf@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
sqrtMin1 <- cutoffs[which.min(sqrt((1-cutoffs$tpr)^2+(cutoffs$fpr)^2)),]  
# print('cutoff:')
# print(sqrtMin)
abline(h=sqrtMin1$tpr, lty=2, lwd=1)   # h，添加水平线
abline(v=sqrtMin1$fpr, lty=2, lwd=1)   # v，添加垂直线
points(sqrtMin1$fpr,sqrtMin1$tpr,pch=21, col='red', bg='red',cex = 1.6)
abline(a = 0, b = 1, lty = 2, col = "gray",lwd=3)
dev.off()




b <-read.table('PremPS_2.tsv',sep="\t",header=TRUE)
roc <- paste('C://Users//jiang//Desktop//PNG//', 'R_PremPS_ROC2.png', sep='')
pred <- prediction(b$PremPS, b$true)
perf <- performance(pred,"tpr","fpr")
png(file=roc,width=15,height=15,units="cm",res=300) 
par(cex.axis = 1.5)
plot(perf,col='black',lwd=4,xlab="", ylab="",cex.main=1.5)
mtext('True positive rate',side=2,line=2.5,cex=2.3,col="black")
mtext('False positive rate',side=1,line=3,cex=2.3,col="black")
auc2 <- performance(pred, "auc")@y.values[[1]]
auc2_1=paste('AUC = ' , as.character(round(auc2,2)))
text(0.7,0.2,labels=auc2_1,cex=1.5)
# ROC optimal threshold 
cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], tpr=perf@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
sqrtMin2 <- cutoffs[which.min(sqrt((1-cutoffs$tpr)^2+(cutoffs$fpr)^2)),]  
# print('cutoff:')
# print(sqrtMin)
abline(h=sqrtMin2$tpr, lty=2, lwd=1)   # h，添加水平线
abline(v=sqrtMin2$fpr, lty=2, lwd=1)   # v，添加垂直线
points(sqrtMin2$fpr,sqrtMin2$tpr,pch=21, col='red', bg='red',cex = 1.6)
abline(a = 0, b = 1, lty = 2, col = "gray",lwd=3)
dev.off()



c <-read.table('Mutabind2_1.tsv',sep="\t",header=TRUE)
roc <- paste('C://Users//jiang//Desktop//PNG//', 'R_Mutabind2_ROC1.png', sep='')
pred <- prediction(c$Mutabind2, c$true)
perf <- performance(pred,"tpr","fpr")
png(file=roc,width=15,height=15,units="cm",res=300) 
par(cex.axis = 1.5)
plot(perf,col='black',lwd=4,xlab="", ylab="",cex.main=1.5)
mtext('True positive rate',side=2,line=2.5,cex=2.3,col="black")
mtext('False positive rate',side=1,line=3,cex=2.3,col="black")
auc3 <- performance(pred, "auc")@y.values[[1]]
auc3_1=paste('AUC = ' , as.character(round(auc3,2)))
text(0.7,0.2,labels=auc3_1,cex=1.5)
# ROC optimal threshold 
cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], tpr=perf@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
sqrtMin3 <- cutoffs[which.min(sqrt((1-cutoffs$tpr)^2+(cutoffs$fpr)^2)),]  
# print('cutoff:')
# print(sqrtMin)
abline(h=sqrtMin3$tpr, lty=2, lwd=1)   # h，添加水平线
abline(v=sqrtMin3$fpr, lty=2, lwd=1)   # v，添加垂直线
points(sqrtMin3$fpr,sqrtMin3$tpr,pch=21, col='red', bg='red',cex = 1.6)
abline(a = 0, b = 1, lty = 2, col = "gray",lwd=3)
dev.off()




d <-read.table('Mutabind2.tsv',sep="\t",header=TRUE)
roc <- paste('C://Users//jiang//Desktop//PNG//', 'R_Mutabind2_ROC2.png', sep='')
pred <- prediction(d$Mutabind2, d$true)
perf <- performance(pred,"tpr","fpr")
png(file=roc,width=15,height=15,units="cm",res=300) 
par(cex.axis = 1.5)
plot(perf,col='black',lwd=4,xlab="", ylab="",cex.main=1.5)
mtext('True positive rate',side=2,line=2.5,cex=2.3,col="black")
mtext('False positive rate',side=1,line=3,cex=2.3,col="black")
auc4 <- performance(pred, "auc")@y.values[[1]]
auc4_1=paste('AUC = ' , as.character(round(auc4,2)))
text(0.7,0.2,labels=auc4_1,cex=1.5)
# ROC optimal threshold 
cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], tpr=perf@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
sqrtMin4 <- cutoffs[which.min(sqrt((1-cutoffs$tpr)^2+(cutoffs$fpr)^2)),]  
# print('cutoff:')
# print(sqrtMin)
abline(h=sqrtMin4$tpr, lty=2, lwd=1)   # h，添加水平线
abline(v=sqrtMin4$fpr, lty=2, lwd=1)   # v，添加垂直线
points(sqrtMin4$fpr,sqrtMin4$tpr,pch=21, col='red', bg='red',cex = 1.6)
abline(a = 0, b = 1, lty = 2, col = "gray",lwd=3)
dev.off()






