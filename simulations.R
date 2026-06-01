#########################################################################################
#
#    Simulations for 95% CI on genotypes and traits
#
#    24 April 2025
#
#
#########################################################################################

#set up workspace
setwd("")
library(readxl)
library(dplyr)
library(tidyr,purrr)

#bring in data files
summary_byrep<-read.csv("./data_summary_clean.csv")
traits<-read.csv("traits.csv")
summary_edges<-read.csv("./data_summaryedges_clean.csv")

# simulating genotypes
n<-c(5,8:20) #all possible sample sizes
genotype_list<-c("LJ01","LJ02","LJ03","LJ04","LM01",
                 "LM02","LM03","LM04","LM05","LM06")
CI<-matrix(NA,nrow=14,ncol=3)
freq<-list()
colnames(CI)<-c("sample","lower","upper")
CI[,1]<-c(5,8:20)

for (j in 1:14){
  bootstrapped<-matrix(NA,nrow=1000,ncol=n[j])
  for (i in 1:1000){
    bootstrapped[i,]<- sample(x=genotype_list,
                              size=n[j],
                              replace=TRUE)
  }
  freq[[j]]<-t(apply(bootstrapped,1,function(row_data){
    table(factor(row_data,level=genotype_list))/n[j]
  }))
  bounds<-matrix(NA,nrow=2,ncol=ncol(freq[[j]]))
  for (i in 1:10){
    bounds[,i]<-quantile(freq[[j]][,i],probs=c(0.0025,0.9975))
  }
  CI[j,2]<-mean(bounds[1,])
  CI[j,3]<-mean(bounds[2,])
}

CI<-as.data.frame(CI)
summary_byrep$CI_freq<-CI$upper[match(summary_byrep$n,CI$sample)]

# bootstrapping traits
traitsims<-list()
traitsim<-data.frame(
  raft.z = rep.int(NA,1000),
  SLA.z = rep.int(NA,1000),
  root.z = rep.int(NA,1000),
  growth.z = rep.int(NA,1000)
)
CI_traits<-data.frame(
  sample=c(5,8:20),
  lower.rafts=rep.int(NA,14),
  upper.rafts=rep.int(NA,14),
  lower.SLA=rep.int(NA,14),
  upper.SLA=rep.int(NA,14),
  lower.root=rep.int(NA,14),
  upper.root=rep.int(NA,14),
  lower.growth=rep.int(NA,14),
  upper.growth=rep.int(NA,14)
)

traits<-traits[order(traits$genotype.ID),]

for (j in 1:length(freq)){
  for (i in 1:1000){
    traitsim$raft.z[i]<-
      sum(freq[[j]][i,]*traits$rafts.z)
    traitsim$SLA.z[i]<-
      sum(freq[[j]][i,]*traits$SLA.z)
    traitsim$root.z[i]<-
      sum(freq[[j]][i,]*traits$root.z)
    traitsim$growth.z[i]<-
      sum(freq[[j]][i,]*traits$growth.z)
  }
  CI_traits[j,2:3]<-quantile(traitsim$raft.z,probs=c(0.0125,0.9875))
  CI_traits[j,4:5]<-quantile(traitsim$SLA.z,probs=c(0.0125,0.9875))
  CI_traits[j,6:7]<-quantile(traitsim$root.z,probs=c(0.0125,0.9875))
  CI_traits[j,8:9]<-quantile(traitsim$growth.z,probs=c(0.0125,0.9875))
  traitsims[[j]]<-traitsim 
}

## get 95% CI for edge-core differences

## make 1000 pair-wise differences for each sample size

CI_diffs<-matrix(NA,nrow=14,ncol=3)
CI_diffs[,1]<-n
colnames(CI_diffs)<-c("sample","lower","upper")
freq_pairs<-matrix(NA,nrow=1000,ncol=14)
colnames(freq_pairs)<-n
for (j in 1:14){
  for (i in 1:1000){
    x<-sample(freq[[j]][,1],2)
    freq_pairs[i,j]<-x[1]-x[2]
  }
  bounds<-matrix(NA,nrow=2,ncol=14)
  bounds[,j]<-quantile(freq_pairs[,j],probs=c(0.0025,0.9975))
  CI_diffs[j,2]<-bounds[1,j]
  CI_diffs[j,3]<-bounds[2,j]
}


CI_diff<-as.data.frame(CI_diffs)
summary_byrep$CI_diff
summary_byrep$CIup_diff<-CI_diff$upper[match(summary_byrep$n,CI_diff$sample)]
summary_byrep$CIlow_diff<-CI_diff$lower[match(summary_byrep$n,CI_diff$sample)]

# core - edge diffs for traits

CI_td<-data.frame(
  sample=c(5,8:20),
  lower.rafts=rep.int(NA,14),
  upper.rafts=rep.int(NA,14),
  lower.SLA=rep.int(NA,14),
  upper.SLA=rep.int(NA,14),
  lower.root=rep.int(NA,14),
  upper.root=rep.int(NA,14),
  lower.growth=rep.int(NA,14),
  upper.growth=rep.int(NA,14)
)

for (j in 1:length(freq)){
  trait_pm<-matrix(NA,nrow=1000,ncol=4)
  for (i in 1:1000){
    x<-sample(traitsims[[j]][,1],2)
    y<-sample(traitsims[[j]][,2],2)
    z<-sample(traitsims[[j]][,3],2)
    v<-sample(traitsims[[j]][,4],2)
    trait_pm[i,]<-c(
      x[1]-x[2],y[1]-y[2],z[1]-z[2],v[1]-v[2]
    )
  }
  CI_td[j,2:3]<-quantile(trait_pm[,1],probs=c(0.0125,0.9875))
  CI_td[j,4:5]<-quantile(trait_pm[,2],probs=c(0.0125,0.9875))
  CI_td[j,6:7]<-quantile(trait_pm[,3],probs=c(0.0125,0.9875))
  CI_td[j,8:9]<-quantile(trait_pm[,4],probs=c(0.0125,0.9875))
}

#add sig to summary files
summary_byrep$sig_LJ01<-summary_byrep$prop_LJ01 > summary_byrep$CI_freq
summary_byrep$sig_LJ02<-summary_byrep$prop_LJ02 > summary_byrep$CI_freq
summary_byrep$sig_LJ03<-summary_byrep$prop_LJ03 > summary_byrep$CI_freq
summary_byrep$sig_LJ04<-summary_byrep$prop_LJ04 > summary_byrep$CI_freq
summary_byrep$sig_LM01<-summary_byrep$prop_LM01 > summary_byrep$CI_freq
summary_byrep$sig_LM02<-summary_byrep$prop_LM02 > summary_byrep$CI_freq
summary_byrep$sig_LM03<-summary_byrep$prop_LM03 > summary_byrep$CI_freq
summary_byrep$sig_LM04<-summary_byrep$prop_LM04 > summary_byrep$CI_freq
summary_byrep$sig_LM05<-summary_byrep$prop_LM05 > summary_byrep$CI_freq
summary_byrep$sig_LM06<-summary_byrep$prop_LM06 > summary_byrep$CI_freq

summary_byrep$CI_raftUP<-CI_traits$upper.rafts[match(summary_byrep$n,CI_traits$sample)]
summary_byrep$CI_raftLOW<-CI_traits$lower.rafts[match(summary_byrep$n,CI_traits$sample)]
summary_byrep$CI_SLAUP<-CI_traits$upper.SLA[match(summary_byrep$n,CI_traits$sample)]
summary_byrep$CI_SLALOW<-CI_traits$lower.SLA[match(summary_byrep$n,CI_traits$sample)]
summary_byrep$CI_rootUP<-CI_traits$upper.root[match(summary_byrep$n,CI_traits$sample)]
summary_byrep$CI_rootLOW<-CI_traits$lower.root[match(summary_byrep$n,CI_traits$sample)]
summary_byrep$CI_growthUP<-CI_traits$upper.growth[match(summary_byrep$n,CI_traits$sample)]
summary_byrep$CI_growthLOW<-CI_traits$lower.growth[match(summary_byrep$n,CI_traits$sample)]

summary_byrep$CIdiff_raftUP<-CI_td$upper.rafts[match(summary_byrep$n,CI_td$sample)]
summary_byrep$CIdiff_raftLOW<-CI_td$lower.rafts[match(summary_byrep$n,CI_td$sample)]
summary_byrep$CIdiff_SLAUP<-CI_td$upper.SLA[match(summary_byrep$n,CI_td$sample)]
summary_byrep$CIdiff_SLALOW<-CI_td$lower.SLA[match(summary_byrep$n,CI_td$sample)]
summary_byrep$CIdiff_rootUP<-CI_td$upper.root[match(summary_byrep$n,CI_td$sample)]
summary_byrep$CIdiff_rootLOW<-CI_td$lower.root[match(summary_byrep$n,CI_td$sample)]
summary_byrep$CIdiff_growthUP<-CI_td$upper.growth[match(summary_byrep$n,CI_td$sample)]
summary_byrep$CIdiff_growthLOW<-CI_td$lower.growth[match(summary_byrep$n,CI_td$sample)]


#change over space
summary_edges$LJ01_diff<-c(summary_edges$prop_LJ01-summary_cores$prop_LJ01)
summary_edges$LJ02_diff<-c(summary_edges$prop_LJ02-summary_cores$prop_LJ02)
summary_edges$LJ03_diff<-c(summary_edges$prop_LJ03-summary_cores$prop_LJ03)
summary_edges$LJ04_diff<-c(summary_edges$prop_LJ04-summary_cores$prop_LJ04)
summary_edges$LM01_diff<-c(summary_edges$prop_LM01-summary_cores$prop_LM01)
summary_edges$LM02_diff<-c(summary_edges$prop_LM02-summary_cores$prop_LM02)
summary_edges$LM03_diff<-c(summary_edges$prop_LM03-summary_cores$prop_LM03)
summary_edges$LM04_diff<-c(summary_edges$prop_LM04-summary_cores$prop_LM04)
summary_edges$LM05_diff<-c(summary_edges$prop_LM05-summary_cores$prop_LM05)
summary_edges$LM06_diff<-c(summary_edges$prop_LM06-summary_cores$prop_LM06)

summary_edges$LJ01_sigdiff<-
  summary_edges$LJ01_diff > summary_edges$CIup_diff | summary_edges$LJ01_diff < summary_edges$CIlow_diff
summary_edges$LJ02_sigdiff<-
  summary_edges$LJ02_diff > summary_edges$CIup_diff | summary_edges$LJ02_diff < summary_edges$CIlow_diff
summary_edges$LJ03_sigdiff<-
  summary_edges$LJ03_diff > summary_edges$CIup_diff | summary_edges$LJ03_diff < summary_edges$CIlow_diff
summary_edges$LJ04_sigdiff<-
  summary_edges$LJ04_diff > summary_edges$CIup_diff | summary_edges$LJ04_diff < summary_edges$CIlow_diff
summary_edges$LM01_sigdiff<-
  summary_edges$LM01_diff > summary_edges$CIup_diff | summary_edges$LM01_diff < summary_edges$CIlow_diff
summary_edges$LM02_sigdiff<-
  summary_edges$LM02_diff > summary_edges$CIup_diff | summary_edges$LM02_diff < summary_edges$CIlow_diff
summary_edges$LM03_sigdiff<-
  summary_edges$LM03_diff > summary_edges$CIup_diff | summary_edges$LM03_diff < summary_edges$CIlow_diff
summary_edges$LM04_sigdiff<-
  summary_edges$LM04_diff > summary_edges$CIup_diff | summary_edges$LM04_diff < summary_edges$CIlow_diff
summary_edges$LM05_sigdiff<-
  summary_edges$LM05_diff > summary_edges$CIup_diff | summary_edges$LM05_diff < summary_edges$CIlow_diff
summary_edges$LM06_sigdiff<-
  summary_edges$LM06_diff > summary_edges$CIup_diff | summary_edges$LM06_diff < summary_edges$CIlow_diff

summary_edges$SLA_diff<-c(summary_edges$SLA.z-summary_cores$SLA.z)
summary_edges$root_diff<-c(summary_edges$root.z-summary_cores$root.z)
summary_edges$raft_diff<-c(summary_edges$raft.z-summary_cores$raft.z)
summary_edges$growth_diff<-c(summary_edges$growth.z-summary_cores$growth.z)

