#This file contains the code to explore how sensitive is the propagation to changes in the values of dependence on 
#the mutualisms (R) assigned to species.

########### 1.  Estimation of the propagation for each treatment after reduced and increased the values of R of each species

#-----------------------------------------------------------------------------------------------------------------
###### 1. Estimation of the propagation for each treatment --

library(tidyverse)
library(dplyr)
source('/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Extra_Functions/Extrafunctions.R')

#Upload matrix
NI_disp <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_disp.csv", sep=",",row.names=1)
NI_pol  <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_pol.csv", sep=",",row.names=1)
I_disp<- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_disp.csv", sep=",",row.names=1)
I_pol <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_pol.csv", sep=";",row.names=1)

# Load Ri values 
Rdisp_NI <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rdisp_NI.csv")
Rdisp_I <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rdisp_I.csv")
Rpol_NI <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rpol_NI.csv")
Rpol_I <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rpol_I.csv")
Rplant_pol_NI <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rplant_pol_NI.csv")
Rplant_pol_I <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rplant_pol_I.csv")
Rplant_disp_NI <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rplant_disp_NI.csv")
Rplant_disp_I <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/Rplant_disp_I.csv")

#Nodes ID
Nodes_NI<- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_nodelist.csv", sep=",",row.names=1)
Nodes_I<- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_nodelist.csv", sep=",",row.names=1)

#roles
NI_sprole <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Modularity/NI_super_list.csv", sep=",",row.names=1)
I_sprole <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Modularity/I_super_list.csv", sep=",",row.names=1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      NI TREATMENT                            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Ndisp=NI_disp
  Npol=NI_pol
  Rdisp=Rdisp_NI
  Rpol=Rpol_NI
  rplants_pol=Rplant_pol_NI
  rplants_disp=Rplant_disp_NI
  sprole=NI_sprole
  Nodes=Nodes_NI

Npol  <- Npol [ order(row.names(Npol)), ] #order pollinators

#Ri values for each species and trophic group
rdisp_base=Rdisp$x
rpol_base=Rpol$x
rplants_pol_base=rplants_pol$x
rplants_disp_base=rplants_disp$x

#Create scenarios where we reduced and increased the values of R of each species used in the simulations by 0.1 until all species had either no or full dependence on the mutualisms
k<-seq(-0.9,0.9,by=0.1) 

res3=data.frame()

for(t in 1:length(k)){
  rdisp<-rdisp_base+k[t]
  rpol<-rpol_base+k[t]
  rplants_pol<-rplants_pol_base+k[t]
  rplants_disp<-rplants_disp_base+k[t]
  rdisp<-replace(rdisp,rdisp<0,0)
  rpol<-replace(rpol,rpol<0,0)
  rplants_pol<-replace(rplants_pol,rplants_pol<0,0)
  rplants_disp<-replace(rplants_disp,rplants_disp<0,0)
  rdisp<-replace(rdisp,rdisp>1,1)
  rpol<-replace(rpol,rpol>1,1)
  rplants_pol<-replace(rplants_pol,rplants_pol>1,1)
  rplants_disp<-replace(rplants_disp,rplants_disp>1,1)
  ndisp=dim(Ndisp)[1]
  npol=dim(Npol)[1]
  nplants=dim(Ndisp)[2]
  nanim=npol+ndisp
  ntot=nplants+npol+ndisp 
  
  # create IM 
  pol_indx=1:npol
  disp_indx=(npol+1):(ndisp+npol)
  
  imatrix=rbind(Npol,Ndisp) 
  IM= as.matrix(imatrix)
  
  ####################
  ####  Scenarios ####
  ####################
  
  # Removing species 
  
  # a. random plants 
  nsims_p = nplants*100
  
  PL.degree=numeric(nsims_p)
  PL.dead.plants=numeric(nsims_p)
  PL.dead.pols=numeric(nsims_p)
  PL.dead.disp=numeric(nsims_p)
  
  set.seed(999)
  bye.plants=rep(1:nplants,100) #
  
  for (i in 1:nsims_p){
    NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='plant',target=bye.plants[i],return.matrix=F)
    PL.degree[i]=max(NC$cascade_data$degree)
    if(is.data.frame(NC$plant_species_data)) {
      PL.dead.plants[i]=length(NC$plant_species_data$lost_plants)
    }
    if (is.data.frame(NC$pol_species_data)) {
      PL.dead.pols[i]=length(NC$pol_species_data$lost_pols)
    }
    if (is.data.frame(NC$disp_species_data)){
      PL.dead.disp[i]=length(NC$disp_species_data$lost_disps)
    } 
  }
  
  #b. seed dispersers
  nsims_d = ndisp*100
  
  DP.degree=numeric(nsims_d)
  DP.dead.plants=numeric(nsims_d)
  DP.dead.pols=numeric(nsims_d)
  DP.dead.disp=numeric(nsims_d)
  
  set.seed(9999)
  bye.disp=rep(disp_indx,100)
  
  
  for (i in 1:nsims_d)
  {
    NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='disperser',target=bye.disp[i],return.matrix=F)
    DP.degree[i]=max(NC$cascade_data$degree)
    if(is.data.frame(NC$plant_species_data)) {
      DP.dead.plants[i]=length(NC$plant_species_data$lost_plants)
    }
    if (is.data.frame(NC$pol_species_data)) {
      DP.dead.pols[i]=length(NC$pol_species_data$lost_pols)
    }
    if (is.data.frame(NC$disp_species_data)){
      DP.dead.disp[i]=length(NC$disp_species_data$lost_disps)
    } 
  }
  
  #c. pollinators
  nsims_pol = npol*100
  
  PO.degree=numeric(nsims_pol)
  PO.dead.plants=numeric(nsims_pol)
  PO.dead.pols=numeric(nsims_pol)
  PO.dead.disp=numeric(nsims_pol)
  
  set.seed(99999)
  bye.pol= rep(pol_indx,100) 
  
  for (i in 1:nsims_pol)
  {
    NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='pollinator',target=bye.pol[i],return.matrix=F)
    PO.degree[i]=max(NC$cascade_data$degree)
    if(is.data.frame(NC$plant_species_data)) {
      PO.dead.plants[i]=length(NC$plant_species_data$lost_plants)
    }
    if (is.data.frame(NC$pol_species_data)) {
      PO.dead.pols[i]=length(NC$pol_species_data$lost_pols)
    }
    if (is.data.frame(NC$disp_species_data)){
      PO.dead.disp[i]=length(NC$disp_species_data$lost_disps)
    } 
  }
  
  #######################################################
  ####### Create data frame and summarize results #######
  #######################################################

  dead.plants=c(PL.dead.plants,DP.dead.plants,PO.dead.plants)
  dead.disps=c(PL.dead.disp,DP.dead.disp,PO.dead.disp)
  dead.pols=c(PL.dead.pols,DP.dead.pols,PO.dead.pols)
  degrees=c(PL.degree,DP.degree,PO.degree)
  initial.random=c(rep('Plant',nsims_p),rep('Disp',nsims_d),rep('Pol',nsims_pol))
  Species_ID = c(rep(1:nplants,100), rep(disp_indx,100), rep(pol_indx,100))
  Results=data.frame(Species_ID,initial.random,dead.plants,dead.disps,dead.pols,degrees)
  Results=Results %>% mutate(porcDeadTot= ((dead.plants+dead.disps+dead.pols)/ntot)*100, porcDeadPlants=(dead.plants/nplants)*100,porcDeadDisp=(dead.disps/ndisp)*100,porcDeadPols=(dead.pols/npol)*100)
  

  Results_perspeciesID<- Results %>% 
    group_by(Species_ID,initial.random) %>% 
    summarise_at(vars(c("dead.plants","dead.disps","dead.pols","degrees","porcDeadPlants","porcDeadTot",
                        "porcDeadDisp","porcDeadPols")),funs(mean(.))) %>% 
    arrange(initial.random)
  Treat = rep("NI",ntot)
  
  Results_final<- cbind(Treat = Treat,Results_perspeciesID)
  Results_final2<- Results_final %>% 
    arrange(factor(initial.random, levels = c("Plant","Pol","Disp"))) %>% 
    mutate(Species_ID = 1:ntot) #we replace the Species ID of pollinator and seed disperser number to match their name in the node list
  
  #change species ID for species name and add species role
  Results_final2$Species_ID<-as.factor(Results_final2$Species_ID)
  levels(Results_final2$Species_ID) <- Nodes$name_node
  factor(Results_final2$Species_ID)
  
  #assign the species role to node id 
  roles<- sprole %>%
    group_by(node_id)%>% 
    slice(1) %>% 
    select(role)
  
  final1<-cbind(Results_final2,roles = roles$role)
  
  res<-final1 %>% 
    group_by(Treat,roles) %>% 
    summarise (dg_mean = mean(degrees),
               dg_sd = sd(degrees),
               porc_dead_tot_mean = mean(porcDeadTot),
               porc_dead_tot_ = sd(porcDeadTot),
               porc_dead_plants_mean = mean(porcDeadPlants),
               porc_dead_disp_mean = mean(porcDeadDisp),
               porc_dead_pol_mean = mean(porcDeadPols))
  
  R = rep(k[t],length(unique(sprole$role)))
  res2<-rbind(data.frame(R = rep(k[t],length(unique(sprole$role))), res))
  res3 = rbind(res3,res2)
}

res_NI = res3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      I TREATMENT                            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ndisp=I_disp
Npol=I_pol
Rdisp=Rdisp_I
Rpol=Rpol_I
rplants_pol=Rplant_pol_I
rplants_disp=Rplant_disp_I
sprole=I_sprole
Nodes=Nodes_I

Npol  <- Npol [ order(row.names(Npol)), ] #order pollinators

#Ri values for each species and trophic group
rdisp_base=Rdisp$x
rpol_base=Rpol$x
rplants_pol_base=rplants_pol$x
rplants_disp_base=rplants_disp$x

#Create scenarios where we reduced and increased the values of R of each species used in the simulations by 0.1 until all species had either no or full dependence on the mutualisms
k<-seq(-0.9,0.9,by=0.1) 

res3=data.frame()

for(t in 1:length(k)){
  rdisp<-rdisp_base+k[t]
  rpol<-rpol_base+k[t]
  rplants_pol<-rplants_pol_base+k[t]
  rplants_disp<-rplants_disp_base+k[t]
  rdisp<-replace(rdisp,rdisp<0,0)
  rpol<-replace(rpol,rpol<0,0)
  rplants_pol<-replace(rplants_pol,rplants_pol<0,0)
  rplants_disp<-replace(rplants_disp,rplants_disp<0,0)
  rdisp<-replace(rdisp,rdisp>1,1)
  rpol<-replace(rpol,rpol>1,1)
  rplants_pol<-replace(rplants_pol,rplants_pol>1,1)
  rplants_disp<-replace(rplants_disp,rplants_disp>1,1)
  ndisp=dim(Ndisp)[1]
  npol=dim(Npol)[1]
  nplants=dim(Ndisp)[2]
  nanim=npol+ndisp
  ntot=nplants+npol+ndisp 
  
  # create IM 
  pol_indx=1:npol
  disp_indx=(npol+1):(ndisp+npol)
  
  imatrix=rbind(Npol,Ndisp) 
  IM= as.matrix(imatrix)
  
  ####################
  ####  Scenarios ####
  ####################
  
  # Removing species 
  
  # a. random plants 
  nsims_p = nplants*100
  
  PL.degree=numeric(nsims_p)
  PL.dead.plants=numeric(nsims_p)
  PL.dead.pols=numeric(nsims_p)
  PL.dead.disp=numeric(nsims_p)
  
  set.seed(999)
  bye.plants=rep(1:nplants,100) #
  
  for (i in 1:nsims_p){
    NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='plant',target=bye.plants[i],return.matrix=F)
    PL.degree[i]=max(NC$cascade_data$degree)
    if(is.data.frame(NC$plant_species_data)) {
      PL.dead.plants[i]=length(NC$plant_species_data$lost_plants)
    }
    if (is.data.frame(NC$pol_species_data)) {
      PL.dead.pols[i]=length(NC$pol_species_data$lost_pols)
    }
    if (is.data.frame(NC$disp_species_data)){
      PL.dead.disp[i]=length(NC$disp_species_data$lost_disps)
    } 
  }
  
  #b. seed dispersers
  nsims_d = ndisp*100
  
  DP.degree=numeric(nsims_d)
  DP.dead.plants=numeric(nsims_d)
  DP.dead.pols=numeric(nsims_d)
  DP.dead.disp=numeric(nsims_d)
  
  set.seed(9999)
  bye.disp=rep(disp_indx,100)
  
  
  for (i in 1:nsims_d)
  {
    NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='disperser',target=bye.disp[i],return.matrix=F)
    DP.degree[i]=max(NC$cascade_data$degree)
    if(is.data.frame(NC$plant_species_data)) {
      DP.dead.plants[i]=length(NC$plant_species_data$lost_plants)
    }
    if (is.data.frame(NC$pol_species_data)) {
      DP.dead.pols[i]=length(NC$pol_species_data$lost_pols)
    }
    if (is.data.frame(NC$disp_species_data)){
      DP.dead.disp[i]=length(NC$disp_species_data$lost_disps)
    } 
  }
  
  #c. pollinators
  nsims_pol = npol*100
  
  PO.degree=numeric(nsims_pol)
  PO.dead.plants=numeric(nsims_pol)
  PO.dead.pols=numeric(nsims_pol)
  PO.dead.disp=numeric(nsims_pol)
  
  set.seed(99999)
  bye.pol= rep(pol_indx,100) 
  
  for (i in 1:nsims_pol)
  {
    NC=netcascade_multi(as.matrix(IM),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='pollinator',target=bye.pol[i],return.matrix=F)
    PO.degree[i]=max(NC$cascade_data$degree)
    if(is.data.frame(NC$plant_species_data)) {
      PO.dead.plants[i]=length(NC$plant_species_data$lost_plants)
    }
    if (is.data.frame(NC$pol_species_data)) {
      PO.dead.pols[i]=length(NC$pol_species_data$lost_pols)
    }
    if (is.data.frame(NC$disp_species_data)){
      PO.dead.disp[i]=length(NC$disp_species_data$lost_disps)
    } 
  }
  
  #######################################################
  ####### Create data frame and summarize results #######
  #######################################################
  
  dead.plants=c(PL.dead.plants,DP.dead.plants,PO.dead.plants)
  dead.disps=c(PL.dead.disp,DP.dead.disp,PO.dead.disp)
  dead.pols=c(PL.dead.pols,DP.dead.pols,PO.dead.pols)
  degrees=c(PL.degree,DP.degree,PO.degree)
  initial.random=c(rep('Plant',nsims_p),rep('Disp',nsims_d),rep('Pol',nsims_pol))
  Species_ID = c(rep(1:nplants,100), rep(disp_indx,100), rep(pol_indx,100))
  Results=data.frame(Species_ID,initial.random,dead.plants,dead.disps,dead.pols,degrees)
  Results=Results %>% mutate(porcDeadTot= ((dead.plants+dead.disps+dead.pols)/ntot)*100, porcDeadPlants=(dead.plants/nplants)*100,porcDeadDisp=(dead.disps/ndisp)*100,porcDeadPols=(dead.pols/npol)*100)
  
  
  Results_perspeciesID<- Results %>% 
    group_by(Species_ID,initial.random) %>% 
    summarise_at(vars(c("dead.plants","dead.disps","dead.pols","degrees","porcDeadPlants","porcDeadTot",
                        "porcDeadDisp","porcDeadPols")),funs(mean(.))) %>% 
    arrange(initial.random)
  Treat = rep("I",ntot)
  
  Results_final<- cbind(Treat = Treat,Results_perspeciesID)
  Results_final2<- Results_final %>% 
    arrange(factor(initial.random, levels = c("Plant","Pol","Disp"))) %>% 
    mutate(Species_ID = 1:ntot) #we replace the Species ID of pollinators and seed disperser number to match their name in the node list
  
  #change species ID for species name and add species role
  Results_final2$Species_ID<-as.factor(Results_final2$Species_ID)
  levels(Results_final2$Species_ID) <- Nodes$name_node
  factor(Results_final2$Species_ID)
  
  #assign the species role to node id 
  roles<- sprole %>%
    group_by(node_id)%>% 
    slice(1) %>% 
    select(role)
  
  final1<-cbind(Results_final2,roles = roles$role)
  
  res<-final1 %>% 
    group_by(Treat,roles) %>% 
    summarise (dg_mean = mean(degrees),
               dg_sd = sd(degrees),
               porc_dead_tot_mean = mean(porcDeadTot),
               porc_dead_tot_ = sd(porcDeadTot),
               porc_dead_plants_mean = mean(porcDeadPlants),
               porc_dead_disp_mean = mean(porcDeadDisp),
               porc_dead_pol_mean = mean(porcDeadPols))
  
  R = rep(k[t],length(unique(sprole$role)))
  res2<-rbind(data.frame(R = rep(k[t],length(unique(sprole$role))), res))
  res3 = rbind(res3,res2)
}

res_I = res3

#######################################################
####### Create FINAL data frame  #######
#######################################################

Prop_sens<- rbind(res_NI,res_I)

#Plot (Appendix S3)--
library(ggplot2)

Prop_sens$Treat<-factor(Prop_sens$Treat, levels=c("NI","I"))

ggplot(Prop_sens) + geom_point(aes(R,porc_dead_tot_mean, color = roles))+ 
  geom_line(aes(R,porc_dead_tot_mean, color = roles))+
  geom_ribbon(aes(x= R, y = porc_dead_tot_mean, ymin = porc_dead_tot_mean - porc_dead_tot_ ,ymax = porc_dead_tot_mean + porc_dead_tot_, fill = roles), alpha = 0.3)+
  geom_vline(xintercept = 0.0, linetype="dotted", 
             color = "red", size=1)+
  scale_color_manual(values = c("#e94c8b","#7f59b0","#fcb424","#00b0f0"))+
  scale_fill_manual(values =  c("#e94c8b","#7f59b0","#fcb424","#00b0f0"))+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black",linetype = "solid", size = 1.4 ),
        axis.title = element_text(size= 15),
        axis.ticks = element_line(colour = "black", size = 1.4),
        axis.text = element_text (colour = "black", size = 14 ))+
  ylab("Extinct species")+
  facet_grid(~Treat)+theme_classic()

