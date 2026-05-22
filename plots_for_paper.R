##############################################################################
#
#  PLOTS  --  Use clean data (see data_cleanup for orignal files + processing)
#  24 April 2025
#
##############################################################################

# Set up workspace

setwd("~/Documents/PhD/THESIS/DUCKWEED")
library(dplyr)
library(tidyr,purrr)


###############################################################################

# bring in data

data_ext<-read.csv("./extent-over-time_clean.csv")
densities<-read.csv("./density_data_clean.csv")
peak_edge_dist<-read.csv("./steepness_clean.csv")
peak_edge_pop<-read.csv("./peak_edge_pop_clean.csv")
genotypes<-read.csv("./genotypes_clean.csv")
summary_byrep<-read.csv("./data_summary_clean.csv")
summary_edges<-read.csv("./data_summaryedges_clean.csv")
summary_cores<-read.csv("./data_summarycores_clean.csv")
traits<-read.csv("traits.csv")

###############################################################################

# Colour scheme for plots

colorvect<-c("#fb9f9f", #LJ01
             "#FF8811", #LJ02
             "#fb5858", #LJ03
             "#ff4500", #LJ04
             
             "#b2d8d8", #LM01
             "#66b2b2", #LM02
             "#004c4c", #LM03
             "#89ecda", #LM04
             "#008080", #LM05
             "#006666" #LM06
)
#ind colors 
LJ01<-"#fb9f9f" #LJ01
LJ02<-"#FF8811" #LJ02
LJ03<-"#fb5858" #LJ03
LJ04<-"#ff4500" #LJ04

LM01<-"#b2d8d8" #LM01
LM02<-"#66b2b2" #LM02
LM03<-"#004c4c" #LM03
LM04<-"#89ecda" #LM04
LM05<-"#008080" #LM05
LM06<-"#006666" #LM06


ctrlcol<-"#62C2F9"
spirocol<-"#A4DA9E"
cmeancol<-"#0571B0"
smeancol<-"#4EAF43"
cfill<-"#3387d6"
sfill<-"#5e9961"

cols16<-c("#89CFF0","#7393B3","#0096FF","#6495ED","#00FFFF","#6F8FAF",
          "#7DF9FF","#6082B6","#5D3FD3","#CCCCFF","#87CEEB","#4682B4",
          "#4169E1","#3F00FF","#0047AB","#A7C7E7")
cols16trans<-paste0(cols16,"88")
cols19<-c("#7FFFD4","#454B1B","#AAFF00","#AFE1AF","#50C878","#5F8575",
          "#4F7942","#228B22","#7CFC00","#008000","#2AAA8A","#4CBB17",
          "#90EE90","#32CD32","#0BDA51","#98FB98","#93C572","#009E60",
          "#00FF7F")
cols19trans<-paste0(cols19,"88")

################################################################################

# Extent over time (Fig 1)

################################################################################

#split by treatment
extentC<-subset(data_ext,Treatment == "C")
extentS<-subset(data_ext,Treatment == "S")

# separate them by replicate
extC_byrep<- list(
  subset(data_ext, Treatment =="C" & Rep ==1),
  subset(data_ext, Treatment =="C" & Rep ==2),
  subset(data_ext, Treatment =="C" & Rep ==4),
  subset(data_ext, Treatment =="C" & Rep ==5),
  subset(data_ext, Treatment =="C" & Rep ==8),
  subset(data_ext, Treatment =="C" & Rep ==11),
  subset(data_ext, Treatment =="C" & Rep ==12),
  subset(data_ext, Treatment =="C" & Rep ==13),
  subset(data_ext, Treatment =="C" & Rep ==14),
  subset(data_ext, Treatment =="C" & Rep ==15),
  subset(data_ext, Treatment =="C" & Rep ==16),
  subset(data_ext, Treatment =="C" & Rep ==18),
  subset(data_ext, Treatment =="C" & Rep ==22),
  subset(data_ext, Treatment =="C" & Rep ==23),
  subset(data_ext, Treatment =="C" & Rep ==24),
  subset(data_ext, Treatment =="C" & Rep ==25)
)

extS_byrep<-list(
  subset(data_ext, Treatment =="S" & Rep ==1),
  subset(data_ext, Treatment =="S" & Rep ==2),
  subset(data_ext, Treatment =="S" & Rep ==3),
  subset(data_ext, Treatment =="S" & Rep ==5),
  subset(data_ext, Treatment =="S" & Rep ==7),
  subset(data_ext, Treatment =="S" & Rep ==9),
  subset(data_ext, Treatment =="S" & Rep ==10),
  subset(data_ext, Treatment =="S" & Rep ==11),
  subset(data_ext, Treatment =="S" & Rep ==12),
  subset(data_ext, Treatment =="S" & Rep ==13),
  subset(data_ext, Treatment =="S" & Rep ==14),
  subset(data_ext, Treatment =="S" & Rep ==16),
  subset(data_ext, Treatment =="S" & Rep ==17),
  subset(data_ext, Treatment =="S" & Rep ==18),
  subset(data_ext, Treatment =="S" & Rep ==19),
  subset(data_ext, Treatment =="S" & Rep ==20),
  subset(data_ext, Treatment =="S" & Rep ==21),
  subset(data_ext, Treatment =="S" & Rep ==24),
  subset(data_ext, Treatment =="S" & Rep ==25)
)

#sensitivity -- without S outlier
extentS_nO<-subset(extentS,extentS$Rep != 20)

#model extent over time
lm_ext_C<-lm(extent~day,extentC)
lm_ext_C<-lm(I(extent-8)~0+day,extentC) # force intercept through 8
summary(lm_ext_C)
lm_ext_S<-lm(extent~day,extentS)
lm_ext_S<-lm(I(extent-8)~0+day,extentS) # force intercept through 8
summary(lm_ext_S)

#get data at day 30 for boxplots
endextent<-subset(data_ext,day==30)

#sensitivity -- without S outlier
endextent_nO<-endextent %>% 
  filter(!(Rep == 20 & endextent$Treatment == "S"))

#draw plot
par(fig=c(0.01,0.8,0,1))
plot(NULL,ylim=c(5,273),xlim=c(1,29),cex.axis=1.5,ylab="Extent (cm)",xlab="Day",cex.lab=1.5)
for (i in 1:length(extC_byrep)){
  data<-extC_byrep[[i]]
  lines(x=data$day,y=data$extent,col=ctrlcol,lwd=2.5)
}
for (i in 1:length(extS_byrep)){
  data<-extS_byrep[[i]]
  lines(x=data$day,y=data$extent,col=spirocol,lwd=2.5)
}
#abline(a=coef(lm_ext_C[1]),coef(lm_ext_C[2]),col=cmeancol,lwd=6)
#abline(a=coef(lm_ext_S[1]),b=coef(lm_ext_S[2]),col=smeancol,lwd=6)

#forced ones
abline(a=8,b=coef(lm_ext_C),col=cmeancol,lwd=6)
abline(a=8,b=coef(lm_ext_S),col=smeancol,lwd=6)

#boxplots
par(fig=c(0.7,0.9,0,1),new=TRUE)
boxplot(endextent$extent[endextent$Treatment=="C"],
        axes=FALSE,
        col=ctrlcol,lwd=2,ylim=c(5,273))

par(fig=c(0.75,0.95,0,1),new=TRUE)
boxplot(endextent$extent[endextent$Treatment=="S"],
        axes=FALSE,
        col=spirocol,ylim=c(5,273),lwd=2)

#stats
sd.ctrl<-sd(endextent$extent[endextent$Treatment=="C"])
sd.spiro<-sd(endextent$extent[endextent$Treatment=="S"])
cv.ctrl<-sd(endextent$extent[endextent$Treatment=="C"])/
  mean(endextent$extent[endextent$Treatment=="C"])

cv.spiro<-sd(endextent$extent[endextent$Treatment=="S"])/
  mean(endextent$extent[endextent$Treatment=="S"])

lnVR_ext<-log(sd.spiro/sd.ctrl)+1/(2*length(endextent$extent[endextent$Treatment=="S"])-1)-1/(2*length(endextent$extent[endextent$Treatment=="C"])-1)


lnCVR_ext<-log(cv.spiro/cv.ctrl)+1/(2*length(endextent$extent[endextent$Treatment=="S"])-1)-1/(2*length(endextent$extent[endextent$Treatment=="C"])-1)

# without outlier
sd.ctrl<-sd(endextent_nO$extent[endextent_nO$Treatment=="C"])
sd.spiro<-sd(endextent_nO$extent[endextent_nO$Treatment=="S"])
cv.ctrl<-sd(endextent_nO$extent[endextent_nO$Treatment=="C"])/
  mean(endextent_nO$extent[endextent_nO$Treatment=="C"])

cv.spiro<-sd(endextent_nO$extent[endextent_nO$Treatment=="S"])/
  mean(endextent_nO$extent[endextent_nO$Treatment=="S"])

lnVR_extNO<-log(sd.spiro/sd.ctrl)+1/(2*length(endextent_nO$extent[endextent_nO$Treatment=="S"])-1)-1/(2*length(endextent_nO$extent[endextent_nO$Treatment=="C"])-1)

lnCVR_extNO<-log(cv.spiro/cv.ctrl)+1/(2*length(endextent_nO$extent[endextent_nO$Treatment=="S"])-1)-1/(2*length(endextent_nO$extent[endextent_nO$Treatment=="C"])-1)


#diff in extent
summary(aov(data_ext$extent[data_ext$day==6]~data_ext$Treatment[data_ext$day==6]))#***
summary(aov(data_ext$extent[data_ext$day==30]~data_ext$Treatment[data_ext$day==30])) 

#########################################################################

# LE Density (Fig 2)

#########################################################################

densC_byrep<- list(
  subset(densities, treatment =="C" & Rep ==1),
  subset(densities, treatment =="C" & Rep ==2),
  subset(densities, treatment =="C" & Rep ==4),
  subset(densities, treatment =="C" & Rep ==5),
  subset(densities, treatment =="C" & Rep ==8),
  subset(densities, treatment =="C" & Rep ==11),
  subset(densities, treatment =="C" & Rep ==12),
  subset(densities, treatment =="C" & Rep ==13),
  subset(densities, treatment =="C" & Rep ==14),
  subset(densities, treatment =="C" & Rep ==15),
  subset(densities, treatment =="C" & Rep ==16),
  subset(densities, treatment =="C" & Rep ==18),
  subset(densities, treatment =="C" & Rep ==22),
  subset(densities, treatment =="C" & Rep ==23),
  subset(densities, treatment =="C" & Rep ==24),
  subset(densities, treatment =="C" & Rep ==25)
)

densS_byrep<-list(
  subset(densities, treatment =="S" & Rep ==1),
  subset(densities, treatment =="S" & Rep ==2),
  subset(densities, treatment =="S" & Rep ==3),
  subset(densities, treatment =="S" & Rep ==5),
  subset(densities, treatment =="S" & Rep ==7),
  subset(densities, treatment =="S" & Rep ==9),
  subset(densities, treatment =="S" & Rep ==10),
  subset(densities, treatment =="S" & Rep ==11),
  subset(densities, treatment =="S" & Rep ==12),
  subset(densities, treatment =="S" & Rep ==13),
  subset(densities, treatment =="S" & Rep ==14),
  subset(densities, treatment =="S" & Rep ==16),
  subset(densities, treatment =="S" & Rep ==17),
  subset(densities, treatment =="S" & Rep ==18),
  subset(densities, treatment =="S" & Rep ==19),
  subset(densities, treatment =="S" & Rep ==20),
  subset(densities, treatment =="S" & Rep ==21),
  subset(densities, treatment =="S" & Rep ==24),
  subset(densities, treatment =="S" & Rep ==25)
)


# density reps + steepness
layout(matrix(c(1,2,3,3,4,4),2,3),widths =c(5,2,2))
par(mar=c(2.6, 4.1, 4.1, 2.1))
plot(NULL, xlim = c(3, 200), ylim = c(0, 80), 
     xlab = "", ylab = "", type = "n",xaxt="n",cex.axis=1.5) 
mtext("a)",side=3,line=0.5,at=1,cex=1.5)
mtext("Frond density",side=2,line=2.5,cex=1.3)
axis(1,at=seq(from=0,to=200,by=10),cex.axis=1.5)
for (i in 1:length(densC_byrep)){
  data<-densC_byrep[[i]]
  lines(x=data$cm_plus,y=data$frond_no,col=cols16trans[i],lwd=2,type="b")
}
par(mar=c(5.1, 4.1, 1.6, 2.1))
# boxplot(peak_edge_pop$peak[peak_edge_pop$treatment=="S"],
#         peak_edge_pop$edge[peak_edge_pop$treatment=="S"],at=c(1.2,74),xlim=c(0,75),
#         boxwex=3.4,xlab="",ylab="", xaxt="n",ylim=c(0,80),col=smeancol)
plot(NULL, xlim = c(3, 70), ylim = c(0, 80), 
     xlab = "", ylab = "", type = "n",xaxt="n",cex.axis=1.5) 
for (i in 1:length(densS_byrep)){
  data<-densS_byrep[[i]]
  lines(x=data$cm_plus,y=data$frond_no,col=cols19trans[i],lwd=2,type="b")
}
axis(1,at=seq(from=5,to=70,by=5),cex.axis=1.5)
mtext("b)",side=3,line=0.5,cex=1.5,at=1)
mtext("Frond density",side=2,line=2.5,cex=1.3)
mtext("Distance (cm)",side=1,line=3,cex=1.3)
par(mar=c(5.1, 3.1, 4.1, 2.1))
boxplot(peak_edge_pop$peak[peak_edge_pop$treatment=="C"],
        peak_edge_pop$peak[peak_edge_pop$treatment=="S"],
        peak_edge_pop$edge[peak_edge_pop$treatment=="C"],
        peak_edge_pop$edge[peak_edge_pop$treatment=="S"],
        at=c(1,2,4,5),xaxt="n",
        col=c(cmeancol,smeancol),cex.axis=1.5)
mtext("c)",side=3,line=0.5,at=0.5,cex=1.5)
mtext("Population density",side=2,line=2.5,cex=1.3)
mtext("Cores",side=1,line=2,cex=1.3,at=1.5)
mtext("Edges",side=1,line=2,cex=1.3,at=4.5)
par(mar=c(5.1, 3.1, 4.1, 2.1))
boxplot(1/peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="C"],
        c(1,1/peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="S"]), # add 1 for one rep with 0 dist, comes out NA otherwise 
        col=c(cmeancol,smeancol),xaxt="n",boxwex=0.8,cex.axis=1.5)
mtext("Wave steepness",side=2,line=2.5,cex=1.5)
mtext("d)",side=3,line=0.5,at=0.25,cex=1.5)

#stats
sd.cntrl<-sd(1/peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="C"])
sd.spiro<-sd(1/peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="S"])
cv.ctrl<-sd(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="C"])/
  mean(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="C"])

cv.spiro<-sd(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="S"])/
  mean(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="S"])

lnVR_steep<-log(sd.spiro/sd.cntrl)+1/(2*length(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="S"])-1)-1/(2*length(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="C"])-1)



lnCVR_steep<-log(cv.spiro/cv.ctrl)+
  1/(2*length(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="S"])-1)-1/(2*length(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="C"])-1)


cv.peak<-sd(peak_edge_dist$peak[peak_edge_dist$treatment=="C"])/
  mean(peak_edge_dist$peak[peak_edge_dist$treatment=="C"])

cv.edgeC<-sd(peak_edge_dist$edge[peak_edge_dist$treatment=="C"])/
  mean(peak_edge_dist$edge[peak_edge_dist$treatment=="C"])

lnVR_steep<-log(sd.spiro/sd.cntrl)+1/(2*length(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="S"])-1)-1/(2*length(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="C"])-1)



lnCVR_steep<-log(cv.spiro/cv.ctrl)+
  1/(2*length(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="S"])-1)-1/(2*length(peak_edge_dist$fifty_pct_edge_dist[peak_edge_dist$treatment=="C"])-1)



cv.peak<-sd(peak_edge_dist$peak[peak_edge_dist$treatment=="S"])/
  mean(peak_edge_dist$peak[peak_edge_dist$treatment=="S"])

cv.edgeS<-sd(peak_edge_dist$edge[peak_edge_dist$treatment=="S"])/
  mean(peak_edge_dist$edge[peak_edge_dist$treatment=="S"])

lnCVR_peS<-log(cv.peak/cv.edge) 
lnCVR_edges<-log(cv.edgeS/cv.edgeC)

summary(aov(peak_edge_dist$fifty_pct_edge_dist~peak_edge_dist$treatment)) #***

##############################################################################################

# LE genotypes and traits (Fig 3)

##############################################################################################

#pooled genotypes
ctrl_gen<-subset(genotypes,genotypes$Treatment == "CTRL")
spiro_gen<-subset(genotypes,genotypes$Treatment == "SPIRO")

c_core_gen<-subset(ctrl_gen,ctrl_gen$Position == "CORE")
s_core_gen<-subset(spiro_gen,spiro_gen$Position == "CORE")

c_edge_gen<-subset(ctrl_gen,ctrl_gen$Position == "EDGE")
s_edge_gen<-subset(spiro_gen,spiro_gen$Position == "EDGE")

#split for plots

summ_ctrl_edge<-subset(summary_edges, summary_edges$Treatment == "C")
summ_ctrl_core<-subset(summary_cores, summary_edges$Treatment == "C")
summ_spiro_edge<-subset(summary_edges, summary_edges$Treatment == "S")
summ_spiro_core<-subset(summary_cores, summary_edges$Treatment == "S")

# bootstrapped SL

freqsigsCE<-summary_byrep %>%
  filter(Treatment=="C" & position=="EDGE") %>%
  select(Treatment,Rep,position,final_ext,prop_LJ01,prop_LJ02, prop_LJ03,prop_LJ04,
         prop_LM01,prop_LM02,prop_LM03,prop_LM04,prop_LM05,prop_LM06,
         sig_LJ01,sig_LJ02,sig_LJ03,sig_LJ04,sig_LM01,sig_LM02,
         sig_LM03,sig_LM04,sig_LM05,sig_LM06) %>%
  pivot_longer(cols=prop_LJ01:sig_LM06,
               names_to = c(".value","genotype"),
               names_sep = "_")


freqsigsSE<-summary_byrep %>%
  filter(Treatment=="S" & position=="EDGE") %>%
  select(Treatment,Rep,position,final_ext,prop_LJ01,prop_LJ02, prop_LJ03,prop_LJ04,
         prop_LM01,prop_LM02,prop_LM03,prop_LM04,prop_LM05,prop_LM06,
         sig_LJ01,sig_LJ02,sig_LJ03,sig_LJ04,sig_LM01,sig_LM02,
         sig_LM03,sig_LM04,sig_LM05,sig_LM06) %>%
  pivot_longer(cols=prop_LJ01:sig_LM06,
               names_to = c(".value","genotype"),
               names_sep = "_")

summary_traitfreq<-summary_byrep %>%
  select(Treatment,Rep,position,final_ext,n,CI_raftUP,CI_raftLOW,
         CI_SLAUP,CI_SLALOW,CI_rootUP,CI_rootLOW,CI_growthLOW,CI_growthUP,CI_growthLOW,raft.z,
         root.z,SLA.z,growth.z)

summary_traitfreq$raft.sig<-
  summary_traitfreq$raft.z > summary_traitfreq$CI_raftUP | summary_traitfreq$raft.z < summary_traitfreq$CI_raftLOW
summary_traitfreq$SLA.sig<-
  summary_traitfreq$SLA.z > summary_traitfreq$CI_SLAUP | summary_traitfreq$SLA.z < summary_traitfreq$CI_SLALOW
summary_traitfreq$root.sig<-
  summary_traitfreq$root.z > summary_traitfreq$CI_rootUP | summary_traitfreq$root.z < summary_traitfreq$CI_rootLOW
summary_traitfreq$growth.sig<-
  summary_traitfreq$growth.z > summary_traitfreq$CI_growthUP | summary_traitfreq$growth.z < summary_traitfreq$CI_growthLOW

traitsigsCE<-summary_traitfreq %>%
  filter(Treatment=="C" & position=="EDGE") %>%
  select(Treatment,Rep,position,final_ext,raft.z,SLA.z,root.z,growth.z,
         raft.sig,SLA.sig,root.sig,growth.sig) %>%
  pivot_longer(cols=raft.z:growth.sig,
               names_to = c("trait",".value"),
               names_sep = "\\.")


traitsigsSE<-summary_traitfreq %>%
  filter(Treatment=="S" & position=="EDGE") %>%
  select(Treatment,Rep,position,final_ext,raft.z,SLA.z,root.z,growth.z,
         raft.sig,SLA.sig,root.sig,growth.sig) %>%
  pivot_longer(cols=raft.z:growth.sig,
               names_to = c("trait",".value"),
               names_sep = "\\.")


freqsigsCE$genotype<-as.factor(freqsigsCE$genotype)
freqsigsSE$genotype<-as.factor(freqsigsSE$genotype)
traitsigsCE$trait<-as.factor(traitsigsCE$trait)
traitsigsSE$trait<-as.factor(traitsigsSE$trait)

#grouped pies
#png("fig3_pies_ctrl.png",width=750,height=850)
#layout(matrix(c(1,1,2,2,3,3,4,4,5,5,5,5,5,6,6,6,7,7,7,7,7,8,8,8),8,3,byrow=FALSE))
png("fig3_pies_ctrl-core.png",width=250,height=200) #stacked is 370 height
#par(mfrow=c(2,1))
par(mar=c(0,0,0,0))
pie(table(c_core_gen$GenotypeID),col=colorvect,cex=1.2)
dev.off()

png("fig3_pies_ctrl-edge.png",width=250,height=200)
par(mar=c(0,0,0,0))
pie(table(c_edge_gen$GenotypeID),col=colorvect,cex=1.2)
dev.off()

png("fig3_pies_spir-core.png",width=250,height=200)
#par(mfrow=c(2,1))
par(mar=c(0,0,0,0))
pie(table(s_core_gen$GenotypeID),col=colorvect,cex=1.2)
dev.off()

png("fig3_pies_spir-edge.png",width=250,height=200)
par(mar=c(0,0,0,0))
pie(table(s_edge_gen$GenotypeID),col=colorvect,cex=1.2)
dev.off()

#Le genfreq (ctrl)

# numeric positions for each level
pos <- c(1:10)

png("fig3_genfreq.png",width=500,height=620)
par(mfrow=c(1,2))
par(mar=c(3.1,5.1,0.5,0))
stripchart(freqsigsCE$prop[freqsigsCE$sig==FALSE]~freqsigsCE$genotype[freqsigsCE$sig==FALSE],
           method="jitter",pch=9,col="darkgray",xlim=c(0,1),xlab="",cex=2,las=1,cex.axis=1.5)
stripchart(freqsigsCE$prop[freqsigsCE$sig==TRUE]~freqsigsCE$genotype[freqsigsCE$sig==TRUE],
           method="jitter",pch=18,col=cmeancol,add=TRUE,cex=3.5,at=pos)
abline(v=0.1,lty=2)

#genfreq (spiro)
par(mar=c(3.1,0.6,0.5,3.1))
stripchart(freqsigsSE$prop[freqsigsSE$sig==FALSE]~freqsigsSE$genotype[freqsigsSE$sig==FALSE],
           method="jitter",pch=10,col="darkgray",xlim=c(0,1),xlab="",cex=2,cex.axis=1.5,yaxt="n")
stripchart(freqsigsSE$prop[freqsigsSE$sig==TRUE]~freqsigsSE$genotype[freqsigsSE$sig==TRUE],
           method="jitter",pch=16,col=smeancol,add=TRUE,cex=2.5,at=pos)
abline(v=0.1,lty=2)
dev.off()

pos<-c(1:4)
#LE trait mean (ctrl)
png("fig3_traitfreq.png",width=540,height=270)
par(mfrow=c(1,2))
par(mar=c(3.1,5.1,1.1,0))
stripchart(traitsigsCE$z[traitsigsCE$sig==FALSE]~traitsigsCE$trait[traitsigsCE$sig==FALSE],
           method="jitter",pch=9,col="darkgray",xlim=c(-2,2),xlab="",cex=2,las=1,cex.axis=1.5,
           group.names=rev(c("SLA","roots","raft sz","fec.")))
stripchart(traitsigsCE$z[traitsigsCE$sig==TRUE]~traitsigsCE$trait[traitsigsCE$sig==TRUE],
           method="jitter",pch=18,col=cmeancol,add=TRUE,cex=3.5,at=pos)
abline(v=0,lty=2)

#LE trait mean (spiro)
par(mar=c(3.1,0.6,1.1,3.1))
stripchart(traitsigsSE$z[traitsigsSE$sig==FALSE]~traitsigsSE$trait[traitsigsSE$sig==FALSE],
           method="jitter",pch=10,col="darkgray",xlim=c(-2,2),xlab="",cex=2,cex.axis=1.5,yaxt="n")
stripchart(traitsigsSE$z[traitsigsSE$sig==TRUE]~traitsigsSE$trait[traitsigsSE$sig==TRUE],
           method="jitter",pch=16,col=smeancol,add=TRUE,cex=2.5,at=pos)
abline(v=0,lty=2)
dev.off()

#####################################################################################################

# Edge - core differences (geno, traits, diversity)

#####################################################################################################


# bootstrapped SL

diffsigsCE<-summary_edges %>%
  filter(Treatment=="C") %>%
  select(Treatment,Rep,position,final_ext,LJ01_diff,LJ02_diff,LJ03_diff,LJ04_diff,
         LM01_diff,LM02_diff,LM03_diff,LM04_diff,LM05_diff,LM06_diff,
         LJ01_sigdiff,LJ02_sigdiff,LJ03_sigdiff,LJ04_sigdiff,
         LM01_sigdiff,LM02_sigdiff,LM03_sigdiff,LM04_sigdiff,LM05_sigdiff,LM06_sigdiff) %>%
  pivot_longer(cols=LJ01_diff:LM06_sigdiff,
               names_to = c("genotype",".value"),
               names_sep = "_")


diffsigsSE<-summary_edges %>%
  filter(Treatment=="S") %>%
  select(Treatment,Rep,position,final_ext,LJ01_diff,LJ02_diff,LJ03_diff,LJ04_diff,
         LM01_diff,LM02_diff,LM03_diff,LM04_diff,LM05_diff,LM06_diff,
         LJ01_sigdiff,LJ02_sigdiff,LJ03_sigdiff,LJ04_sigdiff,
         LM01_sigdiff,LM02_sigdiff,LM03_sigdiff,LM04_sigdiff,LM05_sigdiff,LM06_sigdiff) %>%
  pivot_longer(cols=LJ01_diff:LM06_sigdiff,
               names_to = c("genotype",".value"),
               names_sep = "_")

summary_trait_diff<-summary_edges %>%
  select(Treatment,Rep,position,final_ext,n,CIdiff_raftUP,CIdiff_raftLOW,
         CIdiff_SLAUP,CIdiff_SLALOW,CIdiff_rootUP,CIdiff_rootLOW,CIdiff_growthLOW,CIdiff_growthUP,CIdiff_growthLOW,raft_diff,
         root_diff,SLA_diff,growth_diff)

summary_trait_diff$raft_sigdiff<-
  summary_trait_diff$raft_diff > summary_trait_diff$CIdiff_raftUP | summary_trait_diff$raft_diff < summary_trait_diff$CIdiff_raftLOW
summary_trait_diff$SLA_sigdiff<-
  summary_trait_diff$SLA_diff > summary_trait_diff$CIdiff_SLAUP | summary_trait_diff$SLA_diff < summary_trait_diff$CIdiff_SLALOW
summary_trait_diff$root_sigdiff<-
  summary_trait_diff$root_diff > summary_trait_diff$CIdiff_rootUP | summary_trait_diff$root_diff < summary_trait_diff$CIdiff_rootLOW
summary_trait_diff$growth_sigdiff<-
  summary_trait_diff$growth_diff > summary_trait_diff$CIdiff_growthUP | summary_trait_diff$growth_diff < summary_trait_diff$CIdiff_growthLOW

traitsigdiffsCE<-summary_trait_diff %>%
  filter(Treatment=="C") %>%
  select(Treatment,Rep,position,final_ext,raft_diff,SLA_diff,root_diff,growth_diff,
         raft_sigdiff,SLA_sigdiff,root_sigdiff,growth_sigdiff) %>%
  pivot_longer(cols=raft_diff:growth_sigdiff,
               names_to = c("trait",".value"),
               names_sep = "_")


traitsigdiffsSE<-summary_trait_diff %>%
  filter(Treatment=="S") %>%
  select(Treatment,Rep,position,final_ext,raft_diff,SLA_diff,root_diff,growth_diff,
         raft_sigdiff,SLA_sigdiff,root_sigdiff,growth_sigdiff) %>%
  pivot_longer(cols=raft_diff:growth_sigdiff,
               names_to = c("trait",".value"),
               names_sep = "_")

diffsigsCE$genotype<-as.factor(diffsigsCE$genotype)
diffsigsSE$genotype<-as.factor(diffsigsSE$genotype)
traitsigdiffsCE$trait<-as.factor(traitsigdiffsCE$trait)
traitsigdiffsSE$trait<-as.factor(traitsigdiffsSE$trait)


#png("fig4_gendiff.png",width=700,height=900)
#par(mfrow=c(1,2))
#layout.matrix<-matrix(c(rep.int(1,16),rep.int(2,16),rep.int(0,8),rep.int(3,5),0,5,5,rep.int(3,5),0,5,5,
                       # rep.int(4,5),0,5,5,rep.int(4,5),0,5,5),nrow=8,ncol=9,byrow=FALSE)
layout(matrix(c(1,1,1,2,2,2,0,0,0,3,0,5,4,0,5),nrow=3,ncol=5),widths=c(10,10,0.2,10,10),
       heights=c(9,0.7,4))
par(mar=c(3.1,5.1,2.6,0))
par(oma=c(1,2,1,1))
pos<-c(1:10)
stripchart(diffsigsCE$diff[diffsigsCE$sigdiff==FALSE]~diffsigsCE$genotype[diffsigsCE$sigdiff==FALSE],
           method="jitter",pch=9,col="darkgray",xlim=c(-1.1,1.1),xlab="",cex=2.5,las=1,cex.axis=2)
stripchart(diffsigsCE$diff[diffsigsCE$sigdiff==TRUE]~diffsigsCE$genotype[diffsigsCE$sigdiff==TRUE],
           method="jitter",pch=18,col=cmeancol,add=TRUE,cex=3.5,at=pos)
abline(v=0,lty=2)
mtext("A. Genotype frequencies",side=3,line=0.5,at=1,cex=1.5)

par(mar=c(3.1,0.8,2.6,4.1))


stripchart(diffsigsSE$diff[diffsigsSE$sigdiff==TRUE]~diffsigsSE$genotype[diffsigsSE$sigdiff==TRUE],
           method="jitter",pch=16,col=smeancol,xlim=c(-1.1,1.1),xlab="",cex.axis=2,yaxt="n",cex=2.5)
stripchart(diffsigsSE$diff[diffsigsSE$sigdiff==FALSE]~diffsigsSE$genotype[diffsigsSE$sigdiff==FALSE],
           method="jitter",pch=10,col="darkgray",cex=2,add=TRUE,at=pos)
abline(v=0,lty=2)
#mtext("B.",side=3,line=0.5,at=-1,cex=1.5)
#dev.off()

#traits
#png("fig4_traitdiff.png",width=700,height=500)
#par(mfrow=c(1,2))
pos<-c(1:4)
par(mar=c(3.1,3.1,2.6,1))
stripchart(traitsigdiffsCE$diff[traitsigdiffsCE$sigdiff==FALSE]~traitsigdiffsCE$trait[traitsigdiffsCE$sigdiff==FALSE],
           method="jitter",pch=9,col="darkgray",xlim=c(-2.2,2.2),xlab="",cex=2.5,las=1,cex.axis=2,
           group.names=rev(c("SLA","roots","raft sz","fec."))
           )
stripchart(traitsigdiffsCE$diff[traitsigdiffsCE$sigdiff==TRUE]~traitsigdiffsCE$trait[traitsigdiffsCE$sigdiff==TRUE],
           method="jitter",pch=18,col=cmeancol,add=TRUE,cex=3.5,at=pos)
abline(v=0,lty=2)
mtext("B. Genotype-weighted trait means",side=3,line=0.5,at=2,cex=1.5)

#LE trait mean (spiro)
par(mar=c(3.1,0.1,2.6,3.6))
stripchart(traitsigdiffsSE$diff[traitsigdiffsSE$sigdiff==FALSE]~traitsigdiffsSE$trait[traitsigdiffsSE$sigdiff==FALSE],
           method="jitter",pch=10,col="darkgray",xlim=c(-2.2,2.2),xlab="",cex=2.5,cex.axis=2,yaxt="n")
stripchart(traitsigdiffsSE$diff[traitsigdiffsSE$sigdiff==TRUE]~traitsigdiffsSE$trait[traitsigdiffsSE$sigdiff==TRUE],
           method="jitter",pch=16,col=smeancol,add=TRUE,cex=2.5,at=pos)
abline(v=0,lty=2)
#mtext("d)",side=3,line=0.5,at=-2.5,cex=1.5)
#dev.off()

#png("fig4_divloss.png",width=700,height=200)
par(mar=c(3.1,0,0.1,3.6))
boxplot(summary_edges$div_diff[summary_edges$Treatment=="C"],summary_edges$div_diff[summary_edges$Treatment=="S"],
        horizontal = TRUE,col=c(cmeancol,smeancol),yaxt="n",cex.axis=1.8,outcex=2)
mtext("C. Genotype diversity",side=3,line=0.5,at=0,cex=1.5)
#dev.off()

#PCA

##########################################################################

## Fig 5 Trait correlations / speed

##########################################################################

# steepness (50%)
summary_edges$steepness<-peak_edge_distance$steepness_50

# population density (5 cm behind the edge)
summary_edges$pop_size<-peak_edge_pop$edge


# 4 by 2 panels -- trait dist by genotype + extent by trait


LM.roots<-lm(summary_edges$final_ext~
               summary_edges$Treatment*summary_edges$root.z)
LM.roots<-lm(summary_edges$final_ext~
               summary_edges$Treatment*summary_edges$root.z)
LM.rafts<-lm(summary_edges$final_ext~
               summary_edges$Treatment*summary_edges$raft.z)
LM.SLA<-lm(summary_edges$final_ext~
               summary_edges$Treatment*summary_edges$SLA.z)
LM.growth<-lm(summary_edges$final_ext~
             summary_edges$Treatment*summary_edges$growth.z)

#layout(matrix(c(1,2,3,4,5,6,7,8),nrow=2,ncol=4),heights = c(2,1.5),widths=c(2,2,2,2))
par(mfrow=c(1,4))
par(mar=c(2.1,2.1,1.1,1.1)) 
plot(summary_edges$final_ext[summary_edges$Treatment == "C"]~
       summary_edges$SLA.z[summary_edges$Treatment == "C"],col=ctrlcol,
     pch=16,cex=3,xlab="",ylab="Final extent",ylim=c(0,280),xlim=c(-2,2),cex.axis=1.5)
points(x=summary_edges$SLA.z[summary_edges$Treatment == "S"],
       y=summary_edges$final_ext[summary_edges$Treatment == "S"],
       pch=18,cex=3.5,col=spirocol)
abline(h=275,col="grey")
points(x=traits$SLA.z,y=rep.int(275,times=10),col=colorvect,pch=4,lwd=2,cex=1.5)
lines(range(summary_edges$SLA.z),coef(LM.SLA)[1]+coef(LM.SLA)[3]*range(summary_edges$SLA.z),
                                                                      col="darkgrey",lwd=2,lty=2)
lines(range(summary_edges$SLA.z),coef(LM.SLA)[1]+coef(LM.SLA)[2]
      +(coef(LM.SLA)[3]+coef(LM.SLA)[4])*range(summary_edges$SLA.z),
      col="darkgrey",lwd=2,lty=2)

plot(summary_edges$final_ext[summary_edges$Treatment == "C"]~summary_edges$root.z[summary_edges$Treatment == "C"],
     col=ctrlcol,pch=16,cex=3,xlim=c(-2,2),ylim=c(0,280),xlab="",
     ylab="",cex.axis=1.5)
points(x=summary_edges$root.z[summary_edges$Treatment == "S"],
       y=summary_edges$final_ext[summary_edges$Treatment == "S"],cex=3.5,pch=18,col=spirocol)
abline(h=275,col="grey")
points(x=traits$root.z,y=rep.int(275,times=10),col=colorvect,pch=4,lwd=2,cex=1.5)
lines(range(summary_edges$root.z),coef(LM.roots)[1]+coef(LM.roots)[3]*range(summary_edges$root.z),
      col=cmeancol,lwd=2)
lines(range(summary_edges$root.z),coef(LM.roots)[1]+coef(LM.roots)[2]
      +(coef(LM.roots)[3]+coef(LM.roots)[4])*range(summary_edges$root.z),
      col=smeancol,lwd=2)

plot(summary_edges$final_ext[summary_edges$Treatment == "C"]~summary_edges$raft.z[summary_edges$Treatment == "C"],
     col=ctrlcol,pch=16,cex=3,xlim=c(-2,2),ylim=c(0,280),xlab="",
     ylab="",cex.axis=1.5)
points(x=summary_edges$raft.z[summary_edges$Treatment == "S"],
       y=summary_edges$final_ext[summary_edges$Treatment == "S"],cex=3.5,pch=18,col=spirocol)
abline(h=275,col="grey")
points(x=traits$rafts.z,y=rep.int(275,times=10),col=colorvect,pch=4,lwd=2,cex=1.5)
lines(range(summary_edges$raft.z),coef(LM.rafts)[1]+coef(LM.rafts)[3]*range(summary_edges$raft.z),
      col="darkgrey",lwd=2,lty=2)
lines(range(summary_edges$raft.z),coef(LM.rafts)[1]+coef(LM.rafts)[2]
      +(coef(LM.rafts)[3]+coef(LM.rafts)[4])*range(summary_edges$raft.z),
      col="darkgrey",lwd=2,lty=2)

plot(summary_edges$final_ext[summary_edges$Treatment == "C"]~summary_edges$growth.z[summary_edges$Treatment == "C"],
     col=ctrlcol,pch=16,cex=3,xlim=c(-2.5,2.5),ylim=c(0,280),xlab="",
     ylab="",cex.axis=1.5)
points(x=summary_edges$growth.z[summary_edges$Treatment == "S"],
       y=summary_edges$final_ext[summary_edges$Treatment == "S"],cex=3.5,pch=18,col=spirocol)
abline(h=275,col="grey")
points(x=traits$growth.z,y=c(rep.int(275,times=8),269,281),col=colorvect,pch=4,lwd=2,cex=1.5)
lines(range(summary_edges$growth.z),coef(LM.growth)[1]+coef(LM.growth)[3]*range(summary_edges$growth.z),
      col="darkgrey",lwd=2,lty=2)
lines(range(summary_edges$growth.z),coef(LM.growth)[1]+coef(LM.growth)[2]
      +(coef(LM.growth)[3]+coef(LM.growth)[4])*range(summary_edges$growth.z),
      col="darkgrey",lwd=2,lty=2)
