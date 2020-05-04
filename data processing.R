
##########
##data processing asthma project
##case-control study in Langone Hospital, 2019
##Jiaqian Xing
##########




setwd("/Users/scarlett/Desktop/asthma")
source('Untitled.R')
library(readxl)
library(ggplot2)
asthma<-read_excel("filtered.xlsx")
library(dplyr)
asthma<-asthma%>%mutate(age=AGE,sex=GENDER,educ=Education,race=New.Race,bmi=BMI,
                        cc=Group_2,pre_FVC_L=as.numeric(pre_FVC_L),post_FVC_L=as.numeric(post_FVC_L),
                        pre_FEV1_L=as.numeric(pre_FEV1_L),post_FEV1_L=as.numeric(post_FEV1_L),pre_FVC_pred_pct=as.numeric(pre_FVC_pred_pct),
                        post_FVC_pred_pct=as.numeric(post_FVC_pred_pct),pre_FEV1_pred_pct=as.numeric(pre_FEV1_pred_pct),post_FEV1_pred_pct=as.numeric(post_FEV1_pred_pct),
                        pre_FEV1.FVC=as.numeric(pre_FEV1.FVC),post_FVC_pred_pct=as.numeric(post_FVC_pred_pct),pre_FEV1_pred_pct=as.numeric(pre_FEV1_pred_pct),post_FEV1_pred_pct=as.numeric(post_FEV1_pred_pct),
                        pre_FEV1.FVC=as.numeric(pre_FEV1.FVC),post_FEV1.FVC=as.numeric(post_FEV1.FVC),Eos_pct=as.numeric(Eos_pct))%>%
  mutate(cc=case_when(cc=='Asthma'~1,
                      cc=='Control'~0),
         Group=Group_2)%>%
  select(-AGE,-BMI,-GENDER,-Education,-New.Race,-Eos_pct)

asthma<-asthma%>%mutate(sex=as.factor(case_when(sex=='F'~1,sex=='M'~0)),
                              educ=as.factor(case_when(educ=="High school (12th grade)"~2,
                                             educ=="More than high school" ~3,
                                             educ== "Grade school (up to 6th grade)" ~1)),
                              race=as.factor(case_when(race=="Black"~1,
                                             race=="Latino"~2,
                                             race=='White'~3)),
                        cc=as.factor(cc))


#dummy coding
asthma<-asthma%>%mutate(edu1=case_when(educ=='1'~1,
                                       educ=='2'~0,
                                       educ=='3'~0),
                        edu2=case_when(educ=='1'~0,
                                       educ=='2'~1,
                                       educ=='3'~0),
                        edu3=case_when(educ=='1'~0,
                                       educ=='2'~0,
                                       educ=='3'~1))


#decide which comes in final model
names(asthma)
asthma
asthma.cov<-asthma[,c(25:30,32)]
glm.cov<-glm(cc~., data=asthma.cov,family='binomial')
stepglm<-step(glm.cov)
stepglm
glm.cov
summary(stepglm)




#final dataset:  asthma   204*29   13 protein,  5 demo, 1cc, 



#exclude protein pairs with high correlation in both case and control group
ind1<-which(asthma$cc==1)
ind0<-which(asthma$cc==0)
cor_matrix<-matrix(nrow=0,ncol=4)
for (i in 1:13){
  for (j in 1:13){
    if (i<j){
      cortest1<-cor.test(as.vector(as.matrix(asthma[ind1,i])),as.vector(as.matrix(asthma[ind1,j])),method='kendall')
      cortest0<-cor.test(as.vector(as.matrix(asthma[ind0,i])),as.vector(as.matrix(asthma[ind0,j])),method='kendall')
      cor_res<-c(i,j,cortest0$p.value,cortest1$p.value)
      cor_matrix<-rbind(cor_matrix,cor_res)
      }
  }
}
edge_index<-which(cor_matrix[,3]<0.05|cor_matrix[,4]<0.05)
EDGE_list<-cor_matrix[edge_index,1:2]

#EDGE would be the final protein pair list that would be put in the model.


#############further exclude protein pairs with high computational cost
#DATA preparation

protein<-data.matrix(asthma[,1:13])
cov<-data.matrix(asthma[,c(25,32:34)])
cl<-data.matrix(asthma$cc)
set.seed(7)
#run the denpre2d model to exclude high computational cost protein pairs
for (i in 1:59){
  denPre <- denPre2D(edge=edge_list[i,],classLabel=cLabel,DataFit=DataFit[splitid,],
                     newData=rbind(DataFit[-splitid,],DataPre),method=c("integers","bandwidth")[2])
  denX<-cbind(denX,denPre)
  print(i)
}
edge_list<-EDGE_list[-c(4,11,34:37,45,53),]
edge_list<-edge_list[-c(22,51),]


#run step by step JDINAC.z model
preY <- NULL
vset <- NULL
classLabel<-cl
DataFit<-protein
DataPre<-protein
zFit<-cov
zPre<-cov
nfolds<-5
EDGE<-edge_list
set.seed(7)
for(j in 1:20){
size0 <- sum(classLabel==0)
size1 <- sum(classLabel==1)
sn0 <- round(size0/2)
sn1 <- round(size1/2)
splitid <- c(sample(1:size0,sn0),sample((1:size1)+size0,sn1-1))

cLabel <- classLabel[splitid]
denX<-NULL
y <- classLabel[-splitid]
zfit <- zFit[-splitid, ,drop=F]
cv.fit <- cv.glmnet(x=cbind(zfit,denX[1:length(y),]), y=y, family = "binomial", nfolds=nfolds) 

yp <- predict(cv.fit,newx=cbind(zPre,denX[-(1:length(y)),]), s="lambda.min",type="response")
preY <- cbind(preY,yp) 
coefs <- which(coef(cv.fit,s="lambda.min")[-(1:(1+ncol(zFit)))] !=0)
vset <- c(vset,coefs)
print(coefs)


cLabel <- classLabel[-splitid]  
denX<-NULL
y <- classLabel[splitid]
zfit <- zFit[splitid, ,drop=F]
cv.fit <- cv.glmnet(x=cbind(zfit,denX[1:length(y),]), y=y, family = "binomial", nfolds=nfolds) 
yp <- predict(cv.fit,newx=cbind(zPre,denX[-(1:length(y)),]), s="lambda.min",type="response")
preY <- cbind(preY,yp) 
coefs <- which(coef(cv.fit,s="lambda.min")[-(1:(1+ncol(zFit)))] !=0)
vset <- c(vset,coefs) 
print(coefs)
print(j*100)
}

yPre <- rowMeans(preY) 
numb <- table(vset)
Vid <- as.numeric(rownames(numb))  
Eset <- cbind(EDGE[Vid,],numb)
Eset <- Eset[order(Eset[,3],decreasing=T),]
colnames(Eset) <- c("row","col","numb")

list(yPre=yPre,Eset=Eset) 




###############univariable analysis#################

summary(asthma$age)
sqrt(var(asthma$age))
summary(asthma$bmi)
sqrt(var(asthma$bmi))
#
table(asthma$sex)
prop.table(table(asthma$sex))
table(asthma$race)
prop.table(table(asthma$race))
table(asthma$educ)
prop.table(table(asthma$educ))
table(asthma$bmire)
prop.table(table(asthma$bmire))
table(asthma$cc)
prop.table(table(asthma$cc))
a<-asthma[,14]



###############bivariate analysis####################

index_case<-which(asthma$cc==T)
index_control<-which(asthma$cc!=T)


p<-ggplot(asthma,aes(x=MMP.1))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=CCL17.TARC))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(Angiogenin))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(Angiopoietin.2))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=CRP))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=ST2))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=CHI3L1.YKL.40))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=IL.1ra))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=CCL20.MIP.alpha))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=RAGE))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=Periostin))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=VCAM.1))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=TREM.1))+
  geom_density(aes(color=Group))
p


p<-ggplot(asthma,aes(x=pre_FVC_L))+
  geom_density(aes(color=Group))+
  geom_density(aes(x=post_FVC_L))+
  geom_density(aes(color=Group))

p
p<-ggplot(asthma,aes(x=post_FVC_L))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=pre_FEV1_L))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=post_FEV1_L))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=pre_FVC_pred_pct))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=post_FVC_pred_pct))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=pre_FEV1_pred_pct))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=post_FEV1_pred_pct))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=pre_FEV1.FVC))+
  geom_density(aes(color=Group))
p
p<-ggplot(asthma,aes(x=post_FEV1.FVC))+
  geom_density(aes(color=Group))
p


table(asthma$sex,asthma$cc)
prop.table(table(asthma$sex,asthma$cc),margin=2)
table(asthma$educ,asthma$cc)
prop.table(table(asthma$educ,asthma$cc),margin=2)
table(asthma$race,asthma$cc)
prop.table(table(asthma$race,asthma$cc),margin=2)


mean(asthma$bmi[index_case])
mean(asthma$bmi[index_control])
mean(asthma$age[index_case])
mean(asthma$age[index_control])

summary(asthma$bmi[index_case])
summary(asthma$bmi[index_control])
summary(asthma$age[index_case])
summary(asthma$age[index_control])

sqrt(var(asthma$age[index_control]))
sqrt(var(asthma$age[index_case]))
sqrt(var(asthma$bmi[index_control]))
sqrt(var(asthma$bmi[index_case]))

asthma<-asthma%>%mutate(bmirec=case_when(bmi<18.5~1,
                                         bmi>=18.5&bmi<=24.9~2,
                                         bmi>=25&bmi<=29.9~3,
                                         bmi>29.9~4))



tablesex<-matrix(c(20,35,60,91),nrow=2)
chisq.test(tablesex)
tablerace<-matrix(c(9,36,10,37,84,30),nrow=3)
chisq.test(tablerace)
tableedu<-matrix(c(2,18,35,17,68,66),nrow=2)
chisq.test(tableedu)
tablebmirec<-matrix(c(0,17,17,21,1,29,47,74),nrow=4)
chisq.test(tablebmirec)
wilcox.test(asthma$age[asthma$cc==1], asthma$age[asthma$cc==0], alternative = "two.sided")
wilcox.test(asthma$bmi[asthma$cc==1], asthma$bmi[asthma$cc==0], alternative = "two.sided")




