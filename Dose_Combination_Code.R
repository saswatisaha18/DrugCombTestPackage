#Source Code
library(mvtnorm)
library(multcomp)
library(prodlim)
library(parallel)
library(e1071)
library(plyr)
library(limSolve)
library(doSNOW)
library(bifactorial)
#package can be found here
library(nloptr)
#'\code{eval_f} This function is for computing the SSE while fitting the model
#'@param Y obsservation for response
#'@param len Variable contain the number of patients in each dose combination group
#'@param x mean estimates under each dose group
#'@return the sum
eval_f <- function(Y,len,x) {
  for(i in 1:length(x))
    x[i]<-x[i]
  crossprod(Y-rep(x,len))
}

#'\code{eval_f_eq} This function evaluates the estimate under the null constraint (later used in function Est_eqrest)
#'@param x mean estimates under each dose group
#'The contraint c=MAX(a,b) is evaluated as (a+b-2c+|a-b|)=0
#'@param C1 constraint matrix which evaluates a+b-2c
#'@param C2 constraint matrix which evaluates a-b
#'@return a vector
eval_f_eq<-function(x,C1,C2)
{
  for(i in 1:length(x))
    x[i]<-x[i]
  h<-rep(NA,ncol(C1))
  h<-sapply(1:ncol(C1),function(i){crossprod(C1[,i],x)+abs(crossprod(C2[,i],x))})
  return(h)
}
#'\code{gr} function evaluating the gradient of the results from eval_f function with respect to the mean estimates
#'gr(x)
gr <- function(x) nl.grad(x, eval_f)

#'\code{Est_eqrest} function estimating the mean under the null constraints for the bootstrap method
#'@param Y obsservation for response
#'@param model Annova model fitted to the data
#'@param ind all the drug dose combination groups (index including monotherapies and placebo control) in our study
#'@param indcomb all the drug dose combination groups (index excluding monotherapies and placebo control) in our study
#'@param len Variable contain the number of patients in each dose combination group
#'@param levelsDrugA no of dose groups for Drug A (r)
#'@param levelsDrugB no of dose groups for Drug B (s)
#'@return a vector
#'@references \url{http://stackoverflow.com/a/5577647}
#'@export
Est_eqrest<-function(Y,model,ind,indcomb,len,levelsDrugA,levelsDrugB)
{
  #C1: Evaluates the contrast matrix corresponding to the difference with the first monotherapy
  #C2: Evaluates the contrast matrix corresponding to the difference with the second monotherapy
  
  C1<-C2<-matrix(0,ncol=nrow(ind),nrow=nrow(indcomb))
  for(i in 1:nrow(indcomb))
  { index<- row.match(indcomb[i,], ind) 
  C1[i,index]=C2[i,index]=1
  vec1<-row.match(c(0,indcomb[i,2]),ind)
  vec2<-row.match(c(indcomb[i,1],0),ind)
  C1[i,vec1]=-1
  C2[i,vec2]=-1
  }
  rm(index)
  rm(vec1)
  rm(vec2)
  
  #All possible drug combination indexed by a K digit binary variable
  comb<-bincombinations(levelsDrugA*levelsDrugB)
  #For a particular row (i) in comb, if a particular column (j) entry is:
  #0- it indicates that the combined mean (i,j combination) is equal to the first monotherapy and 
  #1- it indicates that the combined mean (i,j combination)  is equal to the second monotherapy
  
  Contrast<-list()
  Contcomb<-function(c,k)
  { Contmat1<-Contmat2<-matrix(0,nrow(C1),ncol(C1))
  Contmat1=sapply(1:length(c),function(i){if(c[i]==0)return(C1[i,])else return(C2[i,])})
  Contmat2=sapply(1:length(c),function(i){if(c[i]==0)return(C2[i,])else return(C1[i,])})
  if(k==1)
    return(Contmat1)
  else
    return(Contmat2)
  }  
  Contrast1<- lapply(alply(comb,1),function(i)Contcomb(i,1))
  Contrast2<- lapply(alply(comb,1),function(i)Contcomb(i,2))
  #'LSEI FUNCTION  gets the best estimates under the constrain Ex=F and Gx<=F
  #'comb will lead to 2^nrow(comb cases)
  #'powell1 and powell2 evlaluates the mean estimate under the different comb scenarios which gives the minimum SSE.
  #'powell1 evaluates the mean under the null boundary
  #'powell2 evaluates the mean under the global null
  x0=coef(model)
  powell<-matrix(0,nrow=nrow(comb),ncol=length(coef(model)))
  powell1=sapply(1:nrow(comb),function(i)lsei(A=sqrt(diag(len)),B=sqrt(diag(len))%*%x0,E=t(-Contrast1[[i]]),F=rep(0,nrow(indcomb)))$X)
   rm(comb);rm(Contrast1);rm(Contrast2);rm(x0)
  return(powell1[,which.min(apply(powell1,2,function(i)eval_f(Y,len,i)))])
}

#'\code{DrugCombtest} drug combination tests are conducted here.
#'@param Y obsservation for response
#'@param ind all the drug dose combination groups (index including monotherapies and placebo control) in our study
#'@param indcomb all the drug dose combination groups (index excluding monotherapies and placebo control) in our study
#'@param len Variable contain the number of patients in each dose combination group
#'@param levelsDrugA no of dose groups for Drug A (r)
#'@param levelsDrugB no of dose groups for Drug B (s)
#'@return The object save in .
#'@references \url{http://stackoverflow.com/a/5577647}
#'@export
DrugCombtest<-function(Y,len,levelsDrugB,levelsDrugA,ind,indcomb)
{
  data<-data.frame(Treatment=rep(1:((levelsDrugB+1)*(levelsDrugA+1)),len),response=Y)
  C<-list()
  #'Contrast Matrix and Observed Tstat calculation(with individual pvalue)
  model = lm(response ~ factor(Treatment)-1, data = data)
  for(i in 1:nrow(indcomb))
  {
    C[[i]]<-matrix(0,nrow=2,ncol=length(coef(model)))
    index<- row.match(indcomb[i,], ind) 
    C[[i]][,index]=1
    index1<-row.match(c(0,indcomb[i,2]), ind)
    C[[i]][2,index1]=-1
    index2<-row.match(c(indcomb[i,1],0), ind)
    C[[i]][1,index2]=-1
  }
  v<-vcov(model)
  tStat<-sapply(1:nrow(indcomb),function(i)min(crossprod(coef(model),t(C[[i]]))/sqrt(diag(C[[i]]%*%v%*%t(C[[i]])))))
  #'unadjust p values generated from the minimum test
  pStat<-sapply(1:nrow(indcomb),function(i)1-pt(tStat[i], ncp = 0,df=model$df))
  #'Bootstrap method
  #'Mean Estimation 1:
  C1<-C2<-matrix(0,ncol=nrow(ind),nrow=nrow(indcomb))
  for(i in 1:nrow(indcomb))
  { index<- row.match(indcomb[i,], ind) 
  C1[i,index]=-2
  vec1<-row.match(c(0,indcomb[i,2]),ind)
  vec2<-row.match(c(indcomb[i,1],0),ind)
  C1[i,vec1]=C1[i,vec2]=1
  C2[i,vec1]=1
  C2[i,vec2]=-1
  }
  C1<-t(C1)
  C2<-t(C2)
  
  #'the auglag function estimates the mean under the null boundary
  
  #'the Esr_eqrest is performing a similar action as auglag but it is written by the user.
  r<-Est_eqrest(Y,model,ind,indcomb,len,levelsDrugA,levelsDrugB)
  rm(index2)
  rm(vec1)
  rm(vec2)
  rm(index1)
  
  nSim=5000
  tStat_sim<-pStat_sim<-list()
  model1<-list()
  simVec11<-rep(r,len)
  sig0<-crossprod(Y-rep(coef(model),len))/(sum(len)-((levelsDrugA+1)*(levelsDrugB+1)))
  
  #This function evaluates the quantile under the bootstrap method
  boot_combtest<- function(simVec,levelsDrugB,levelsDrugA,len)
  {
    
    Y_sim<- Data_Sim<-list()
    Y_sim<-simVec+rnorm(length(simVec),0,sqrt(sig0))
    Data_Sim<-data.frame(Treatment=rep(1:((levelsDrugB+1)*(levelsDrugA+1)),len),response=Y_sim)
    est<-tapply(Data_Sim$response,Data_Sim$Treatment,mean)
    vc<-as.numeric(crossprod(Y_sim-rep(est,len))/(length(Y_sim)-length(len)))*diag(1/len)
    return(sapply(1:nrow(indcomb),function(i)min(crossprod(est,t(C[[i]]))/sqrt(diag(C[[i]]%*%vc%*%t(C[[i]]))))))
  }
  

  sim1<-sapply(1:nSim,function(i)boot_combtest(simVec11,levelsDrugB,levelsDrugA,len))

  if(nrow(indcomb)==1)
  tStat_sim1<-sim1
  else
  tStat_sim1<-sapply(1:nSim,function(i)max(sim1[,i]))

  m<-list()
  m$Tstat<-tStat

  m$pValBoot1adj<-sapply(1:(levelsDrugB*levelsDrugA),function(l)mean(tStat_sim1>tStat[l]))
  m$pValunadj<-pStat
  m$pValBonfadj<-round(p.adjust(pStat,method="bonferroni"),3)
  return (m)
}
#'\code{Dist_restnull2} This function evaluates the p value under the least favourable null approach 
#'@param ind all the drug dose combination groups (index including monotherapies and placebo control) in our study
#'@param indcomb all the drug dose combination groups (index excluding monotherapies and placebo control) in our study
#'@param len Variable contain the number of patients in each dose combination group
#'@param Crit value of obs test statistics where the p value should be computed
#'@return returns the p-value vector
#'@references \url{http://stackoverflow.com/a/5577647}
#'@export
Dist_restnull2<-function(levelsDrugA,levelsDrugB,len,ind,indcomb,Crit)
{
  C1<-C2<-matrix(0,ncol=nrow(ind),nrow=nrow(indcomb))
  for(i in 1:nrow(indcomb))
  { index<- row.match(indcomb[i,], ind) 
  C1[i,index]=C2[i,index]=1
  vec1<-row.match(c(0,indcomb[i,2]),ind)
  vec2<-row.match(c(indcomb[i,1],0),ind)
  C1[i,vec1]=-1
  C2[i,vec2]=-1
  }
  Sigma<-diag(1/len)
  comb<-bincombinations(levelsDrugA*levelsDrugB)
  #compute the test statistics T11 and T12 using the following function
  prod<-function(a,b){a%*%Sigma%*%b/{sqrt(a%*%Sigma%*%a)%*%sqrt(b%*%Sigma%*%b)}}
  
  Contrast<-list()
  Contcomb<-function(c)
  { 
    A<-matrix(0,nrow=length(c)-1,ncol=length(c)-1)
    
    for(i in 1:(length(c)-1))
    {k=1; 
    while(i+k<=length(c))
    {if(c[i]==0)
    {A[i,k]=ifelse(c[i+k]==0,prod(C1[i,],C1[i+k,]),prod(C1[i,],C2[i+k,]))}
      else 
      {A[i,k]=ifelse(c[i+k]==0,prod(C1[i+k,],C2[i,]),prod(C2[i+k,],C2[i,]))};
      k=k+1
    }
    }
    CovMat<-diag(length(c)) 
    for(k in 1:(nrow(CovMat)-1))
      CovMat[(k+1):nrow(CovMat),k]=CovMat[k,(k+1):nrow(CovMat)]=A[k,1:(nrow(CovMat)-k)]
    
    return(CovMat)
  }  
  #This evaluates the covariance matrix under different combination of the test statistic
  ContMat<- suppressWarnings(lapply(alply(comb,1),function(i)Contcomb(i)))
  #fn1<-function(x)sapply(ContMat,function(y)qmvt(x,tail="upper.tail",sigma=y,df=sum(len)-nrow(ind)))
  fn2<-function(y){(max(sapply(ContMat,function(i){(1-pmvt(upper=rep(y,ncol(comb)),lower=-Inf,delta=rep(0,ncol(comb)),sigma=i,df=sum(len)-nrow(ind)))})))}
  return(fn2(Crit))
}
#'\code{ApplyDrugComb} This function evaluates the p value under the different methos
#'@param levelsDrugA no of dose groups for Drug A (r)
#'@param levelsDrugB no of dose groups for Drug B (s)
#'@param len Variable contain the number of patients in each dose combination group
#'@param Y raw data
#'@param m dose-response means
#'@references \url{http://stackoverflow.com/a/5577647}
#'@export
ApplyDrugComb<-function(levelsDrugA,levelsDrugB,Y,len)
{
  ind<-expand.grid(c(0:levelsDrugA),c(0:levelsDrugB))
  indcomb<-expand.grid(c(1:levelsDrugA),c(1:levelsDrugB))
  nlev=(levelsDrugA+1)*(levelsDrugB+1)
  m<-DrugCombtest(Y,len=ln,levelsDrugB,levelsDrugA,ind,indcomb)
  m$LF_Null<-sapply(m$Tstat,function(i)Dist_restnull2(levelsDrugA,levelsDrugB,len,ind,indcomb,i))
  res<-do.call("cbind",m)
  
  #res<-res[,c(1,4,5,2,6)]
  colnames(res)<-c("Tstat","Boot_adj","Unadj","Bonf","LFC")
  
  #Final result
  DoseInd<-do.call("paste",c(indcomb,sep=","))
  DoseInd<- DoseInd[order(res[,1])]
  res<-res[order(res[,1]),]
  return(res)
  # print(as.data.frame(cbind(Dose=DoseInd,formatC( round( res, 3 ),format='f', digits=3 ))))
  
  
}

#RealData
levelsDrugB=2
levelsDrugA=3

ind<-expand.grid(c(0:levelsDrugA),c(0:levelsDrugB))
indcomb<-expand.grid(c(1:levelsDrugA),c(1:levelsDrugB))

Y<-list()
len<-n<-c(75,75,74,48,74,75,74,49,48,50,48,48)
m<-c(0,1.4,2.7,4.6,1.8,2.8,5.7,8.2,2.8,4.5,7.2,10.9)
s<-rep(7.07,12)
x<-list(12)
for(i in 1:12){
  x[[i]]<-rnorm(n[i],mean=0,sd=1)
  x[[i]]<-((x[[i]]-mean(x[[i]]))*(s[i]/sd(x[[i]])))+m[i]
}

hung<-carpet(x,D=c(2,3))
avetest(hung,test="ttest",nboot=20000)
maxtest(hung,test="ttest",nboot=20000)
Y<-unlist(x)
m<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len)

res<-do.call("cbind",m)
res<-res[order(res[,1]),]
round(res,3)



P <- detectCores(logical =TRUE) # physical cores
cl <- makeCluster(P)
clusterEvalQ(cl,library(DoseFinding))
clusterEvalQ(cl,library(multcomp))
clusterEvalQ(cl,library(mvtnorm))
clusterEvalQ(cl,library(prodlim))
clusterEvalQ(cl,library(parallel))
clusterEvalQ(cl,library(doSNOW))
clusterEvalQ(cl,library(e1071))
clusterEvalQ(cl,library(plyr))
clusterEvalQ(cl,library(limSolve))
clusterEvalQ(cl,library(bifactorial))
registerDoSNOW(cl)
levelsDrugB=1
levelsDrugA=2
nsim=10
res<-list()
ind<-expand.grid(c(0:levelsDrugA),c(0:levelsDrugB))
indcomb<-expand.grid(c(1:levelsDrugA),c(1:levelsDrugB))
ssize=c(10,25,50,75,100)
iteration=rep(ssize,each=nsim)
clusterExport(cl,c("eval_f","DrugCombtest","Est_eqrest","ApplyDrugComb"))  #make required functions known to cluster
clusterExport(cl, c("ssize","ind","indcomb","levelsDrugA","levelsDrugB"))
#pvalue_ Scenario-1
res2<-list()
del=2
ind<-expand.grid(c(0:levelsDrugA),c(0:levelsDrugB))
indcomb<-expand.grid(c(1:levelsDrugA),c(1:levelsDrugB))
comb<-bincombinations(nrow(indcomb))
source_data('https://github.com/saswatisaha18/DrugCombTestPackage/blob/master/Seed_Seq.Rdata?raw=TRUE')

#seed_seq<-readRDS("C:/Users/saha/Desktop/Reading Material/New Papers Implementation/Confidence Interval/Secondment/Novartis/Seed_Seq.RDS")
#Scenario 1 and Scenario 2
Al=c(9999,0.7)
system.time(res<-foreach(al = Al)%:%
foreach(i= 1:length(iteration),.combine="rbind")%dopar%
{ 
  set.seed(seed_seq[i])
  rand<-sample(1:nrow(comb),1)
  vec=switch(rand,"1"=c(del,del+al,del+al,del,del+al,del+al),
             "2"=c(del,del,del+2*al,del+al,del+al,del+2*al),
             "3"=c(del,del+2*al,del,del+al,del+2*al,del+al),
             "4"=c(del,del,del,del+al,del+al,del+al))
  ln=rep(iteration[i],nrow(ind))
  set.seed(seed_seq[i])
  Y<-rep(vec,each=iteration[i])+rnorm(sum(ln),0,1)
   Z<-list()
   Z[[1]]<-Y[1:ln[1]]
   for(j in 2:length(vec))
     Z[[j]]<-Y[(sum(ln[1:(j-1)])+1):(sum(ln[1:j]))]
   hung<-carpet(Z,D=c(levelsDrugB,levelsDrugA))
   p_max2<-maxtest(hung,method="hung",test="ttest")@p 
   m1<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len=ln)
  return(c(Tstat=max(m1[,1]),apply(m1[,-1],2,min),Hung=p_max2))
 })

res3<-apply(res[[1]][,-1],2,function(i)ifelse(i<=0.05,1,0))
t(cbind(apply(res3[1:5000,],2,sum)/nsim,apply(res3[5001:10000,],2,sum)/nsim,apply(res3[10001:15000,],2,sum)/nsim,apply(res3[15001:20000,],2,sum)/nsim,apply(res3[20001:25000,],2,sum)/nsim))
#Scenario 3
del=2
al=0.5
vec=c(del,del,del,del,del+al,del+al)
system.time(res1<-foreach(i= 1:length(iteration),.combine="rbind")%dopar%
{ 
  ln=rep(iteration[i],nrow(ind))
  set.seed(seed_seq[i])
  Y<-rep(vec,each=iteration[i])+rnorm(sum(ln),0,1)
  Z<-list()
  Z[[1]]<-Y[1:ln[1]]
  for(j in 2:length(vec))
    Z[[j]]<-Y[(sum(ln[1:(j-1)])+1):(sum(ln[1:j]))]
  hung<-carpet(Z,D=c(levelsDrugB,levelsDrugA))
  p_max2<-maxtest(hung,method="hung",test="ttest")@p 
  m1<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len=ln)
  return(c(Tstat=max(m1[,1]),apply(m1[,-1],2,min),Hung=p_max2))
})
res3<-apply(res1[,-1],2,function(i)ifelse(i<=0.05,1,0))
t(cbind(apply(res3[1:5000,],2,sum)/nsim,apply(res3[5001:10000,],2,sum)/nsim,apply(res3[10001:15000,],2,sum)/nsim,apply(res3[15001:20000,],2,sum)/nsim,apply(res3[20001:25000,],2,sum)/nsim))

#Scenario 4
levelsDrugA=3
levelsDrugB=2
nsim=5000
indcomb<-expand.grid(c(1:levelsDrugA),c(1:levelsDrugB))
ind<-expand.grid(c(0:levelsDrugA),c(0:levelsDrugB))

vec=c(0,0.1,0.3,0.6,0.2,0.5,0.6,0.9,0.5,0.8,0.8,0.9)
#vec=c(0,0.1,0.3,0.6,0.2,0.25,0.7,1.0,0.5,0.65,0.90,1.0)
ln=rep(50,length(vec))
system.time(res2<-foreach(i= 1:nsim)%dopar%
{
  set.seed(seed_seq[i])
  Y<-rep(vec,ln)+rnorm(sum(ln),0,1)
  Z<-list()
  Z[[1]]<-Y[1:ln[1]]
  for(j in 2:length(vec))
    Z[[j]]<-Y[(sum(ln[1:(j-1)])+1):(sum(ln[1:j]))]
  hung<-carpet(Z,D=c(levelsDrugB,levelsDrugA))
  p_max2<-maxtest(hung,method="hung",test="ttest")@p 
  m1<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len=ln)
  return(c(Tstat=max(m1[,1]),apply(m1[,-1],2,min),Hung=p_max2))
})
Res<-do.call("rbind",res2)
res3<-apply(Res[,-1],2,function(i)ifelse(i<=0.05,1,0))
apply(res3,2,sum)/nsim
#scenario 5
ln=c(50,35,90,35,90,35,90,35,35,35,35,35)
system.time(res3<-foreach(i= 1:nsim)%dopar%
{
  set.seed(seed_seq[i])
  Y<-rep(vec,ln)+rnorm(sum(ln),0,1)
  Z<-list()
  Z[[1]]<-Y[1:ln[1]]
  for(j in 2:length(vec))
    Z[[j]]<-Y[(sum(ln[1:(j-1)])+1):(sum(ln[1:j]))]
  hung<-carpet(Z,D=c(levelsDrugB,levelsDrugA))
  p_max2<-maxtest(hung,method="hung",test="ttest")@p 
  m1<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len=ln)
  return(c(Tstat=max(m1[,1]),apply(m1[,-1],2,min),Hung=p_max2))
})
Res<-do.call("rbind",res3)
re3<-apply(Res[,-1],2,function(i)ifelse(i<=0.05,1,0))
apply(re3,2,sum)/nsim


#Scenario 6
ln=c(50,70,70,70,20,50,50,50,20,50,50,50)
system.time(res4<-foreach(i= 1:nsim)%dopar%
{
  set.seed(seed_seq[i])
  Y<-rep(vec,ln)+rnorm(sum(ln),0,1)
  Z<-list()
  Z[[1]]<-Y[1:ln[1]]
  for(j in 2:length(vec))
    Z[[j]]<-Y[(sum(ln[1:(j-1)])+1):(sum(ln[1:j]))]
  hung<-carpet(Z,D=c(levelsDrugB,levelsDrugA))
  p_max2<-maxtest(hung,method="hung",test="ttest")@p 
  m1<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len=ln)
  return(c(Tstat=max(m1[,1]),apply(m1[,-1],2,min),Hung=p_max2))
})

Res<-do.call("rbind",res4)
re<-apply(Res[,-1],2,function(i)ifelse(i<=0.05,1,0))
apply(re,2,sum)/nsim

#Scenario 7
ln=c(50,30,58,58,56,30,58,58,56,30,58,58)
system.time(res5<-foreach(i= 1:nsim)%dopar%
{
  set.seed(seed[i])
  Y<-rep(vec,ln)+rnorm(sum(ln),0,1)
  Z<-list()
  Z[[1]]<-Y[1:ln[1]]
  for(j in 2:length(vec))
    Z[[j]]<-Y[(sum(ln[1:(j-1)])+1):(sum(ln[1:j]))]
  hung<-carpet(Z,D=c(levelsDrugB,levelsDrugA))
  p_max2<-maxtest(hung,method="hung",test="ttest")@p 
  m1<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len=ln)
  return(c(Tstat=max(m1[,1]),apply(m1[,-1],2,min),Hung=p_max2))
})

Res<-do.call("rbind",res5)
re<-apply(Res[,-1],2,function(i)ifelse(i<=0.05,1,0))
apply(re,2,sum)/nsim

#Scenario 8
vec=c(0,0.1,0.3,0.6,0.2,0.25,0.7,1.0,0.5,0.65,0.90,1.0)
ln=rep(50,length(vec))
system.time(res6<-foreach(i= 1:nsim)%dopar%
{
  set.seed(seed[i])
  Y<-rep(vec,ln)+rnorm(sum(ln),0,1)
  Z<-list()
  Z[[1]]<-Y[1:ln[1]]
  for(j in 2:length(vec))
    Z[[j]]<-Y[(sum(ln[1:(j-1)])+1):(sum(ln[1:j]))]
  hung<-carpet(Z,D=c(levelsDrugB,levelsDrugA))
  p_max2<-maxtest(hung,method="hung",test="ttest")@p 
  m1<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len=ln)
  return(c(Tstat=max(m1[,1]),apply(m1[,-1],2,min),Hung=p-max2))
})

Res<-do.call("rbind",res6)
re<-apply(Res[,-1],2,function(i)ifelse(i<=0.05,1,0))
apply(re,2,sum)/nsim
#scenario 9
ln=c(50,35,90,35,90,35,90,35,35,35,35,35)
system.time(res7<-foreach(i= 1:nsim)%dopar%
{
  set.seed(seed[i])
  Y<-rep(vec,ln)+rnorm(sum(ln),0,1)
  Z<-list()
  Z[[1]]<-Y[1:ln[1]]
  for(j in 2:length(vec))
    Z[[j]]<-Y[(sum(ln[1:(j-1)])+1):(sum(ln[1:j]))]
  hung<-carpet(Z,D=c(levelsDrugB,levelsDrugA))
  p_max2<-maxtest(hung,method="hung",test="ttest")@p
  m1<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len=ln)
  return(c(Tstat=max(m1[,1]),apply(m1[,-1],2,min),Hung=p_max2))
})

Res<-do.call("rbind",res7)
re<-apply(Res[,-1],2,function(i)ifelse(i<=0.05,1,0))
apply(re,2,sum)/nsim
#Scenario 10
ln=c(50,70,70,70,20,50,50,50,20,50,50,50)
system.time(res8<-foreach(i= 1:nsim)%dopar%
{
  set.seed(seed[i])
  Y<-rep(vec,ln)+rnorm(sum(ln),0,1)
  Z<-list()
  Z[[1]]<-Y[1:ln[1]]
  for(j in 2:length(vec))
    Z[[j]]<-Y[(sum(ln[1:(j-1)])+1):(sum(ln[1:j]))]
  hung<-carpet(Z,D=c(levelsDrugB,levelsDrugA))
  p_max2<-maxtest(hung,method="hung",test="ttest")@p
  m1<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len=ln)
  return(c(Tstat=max(m1[,1]),apply(m1[,-1],2,min),Hung=p_max2))
})

Res<-do.call("rbind",res8)
re<-apply(Res[,-1],2,function(i)ifelse(i<=0.05,1,0))
apply(re,2,sum)/nsim
#Scenario 11
ln=c(50,30,58,58,56,30,58,58,56,30,58,58)
system.time(res9<-foreach(i= 1:nsim)%dopar%
{
  set.seed(seed[i])
  Y<-rep(vec,ln)+rnorm(sum(ln),0,1)
  Z<-list()
  Z[[1]]<-Y[1:ln[1]]
  for(j in 2:length(vec))
    Z[[j]]<-Y[(sum(ln[1:(j-1)])+1):(sum(ln[1:j]))]
  hung<-carpet(Z,D=c(levelsDrugB,levelsDrugA))
  m1<-ApplyDrugComb(levelsDrugA,levelsDrugB,Y,len=ln)
  return(c(Tstat=max(m1[,1]),apply(m1[,-1],2,min),Hung=p_max2))
})
Res<-do.call("rbind",res9)
re<-apply(Res[,-1],2,function(i)ifelse(i<=0.05,1,0))
apply(re,2,sum)/nsim
