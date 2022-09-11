########### 1. Estimation of the propagation of disturbances in networks controlled by size (NI treatment) 

########### 2. Comparison with matrix of NI without controlling by size

#-----------------------------------------------------------------------------------------------------------------

########### 1. Estimation of the propagation of disturbances in networks controlled by size (NI treatment) 

library(tidyverse)
library(dplyr)
source('/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Extra_Functions/Extrafunctions.R')
source("Simulation_control_size.R") #simulated networks controlled by size

# Load Ri values to assign of each species in each network
tdisp <- read.csv("traits_frug.csv", sep = ";")
tpol<- read.csv("traits_pol.csv", sep = ",")
tplants_pol<- read.csv("traits_plants_pol_2.csv")
tplants_disp<- read.csv("traits_plants_disp_2.csv")
tpol %>% arrange(Species_Ref) #order the pollinators in the R list according alphabetic order

# List of roles to assign of each species
sprole <- read.csv("/Users/agustin/Desktop/Papers/Paper_multicapa/Matrices_consecos/ISRAEL/Analisis_pooled/Modularity/Modularity_sproles/Change_sprole/super_list_treat.csv", sep=",",row.names=1)

#loops
Final=data.frame()
for(m in res){
  
###we assign the R for each species in each simulated matrix --
  
#plants_pollinators
  plant.common=colnames(m)[colnames(m)%in%tplants_pol[,1]]
  Rplants_pol<- tplants_pol %>% 
    filter(Plant_name%in%plant.common)
  rplants_pol<- Rplants_pol[,5]
  
#plants_dispersers
  plant.common=colnames(m)[colnames(m)%in%tplants_disp[,1]]
  Rplants_disp<- tplants_disp %>% 
    filter(Plant_name%in%plant.common)
  rplants_disp<- Rplants_disp[,4]
  
#pollinators
  pol.common = rownames(m)[rownames(m)%in%tpol[,2]]
  Rpol<-tpol %>% 
    filter(Species_Ref%in%pol.common) %>% 
    arrange(Species_Ref)
  rpol<-Rpol[,3]
  
#dispersers
  disp.common=row.names(m)[row.names(m)%in%tdisp[,1]]
  Rdisp<- tdisp %>% 
    filter(Seeddisp_name%in%disp.common)
  rdisp<-Rdisp[,7]
  
  ndisp= length(rdisp)
  npol= length(rpol)
  nplants=ncol(m)
  nanim= nrow(m)
  ntot=nplants+npol+ndisp #Total number of nodes
  pol_indx=1:npol
  disp_indx=(npol+1):(ndisp+npol)
  
####################
####  Scenarios ####
####################
  
  # Removing species 
  
  # a. random plants 
  nsims_p = nplants
  
  PL.degree=numeric(nsims_p)
  PL.dead.plants=numeric(nsims_p)
  PL.dead.pols=numeric(nsims_p)
  PL.dead.disp=numeric(nsims_p)
  
  set.seed(99999)
  bye.plants=1:nplants#
  
  for (i in 1:nsims_p){#Call the function Stochastic coextinction model
    NC=netcascade_multi(as.matrix(m),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='plant',target=bye.plants[i],return.matrix=F)
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
  nsims_d = length(disp_indx)
  
  DP.degree=numeric(nsims_d)
  DP.dead.plants=numeric(nsims_d)
  DP.dead.pols=numeric(nsims_d)
  DP.dead.disp=numeric(nsims_d)
  
  set.seed(9999)
  bye.disp=disp_indx
  
  
  for (i in 1:nsims_d){
    NC=netcascade_multi(as.matrix(m),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='disperser',target=bye.disp[i],return.matrix=F)
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
  nsims_pol = npol
  
  PO.degree=numeric(nsims_pol)
  PO.dead.plants=numeric(nsims_pol)
  PO.dead.pols=numeric(nsims_pol)
  PO.dead.disp=numeric(nsims_pol)
  
  set.seed(99999)
  bye.pol= pol_indx 
  
  for (i in 1:nsims_pol)
  {
    NC=netcascade_multi(as.matrix(m),rplants_disp,rplants_pol, rpol, rdisp,deadPlants=NULL, deadPols=NULL, deadDisps=NULL,targetGuild='pollinator',target=bye.pol[i],return.matrix=F)
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
  
  dead.plants=c(PL.dead.plants,PO.dead.plants,DP.dead.plants) 
  dead.disps=c(PL.dead.disp,PO.dead.disp,DP.dead.disp)
  dead.pols=c(PL.dead.pols,PO.dead.pols,DP.dead.pols)
  degrees=c(PL.degree,PO.degree,DP.degree)
  initial.random=c(rep('Plant',nsims_p),rep('Pol',nsims_pol),rep('Disp',nsims_d))
  Species_ID= c(colnames(m),rownames(m))
  Results=data.frame(Species_ID,initial.random,dead.plants,dead.disps,dead.pols,degrees)
  Results=Results %>% 
    mutate(porcDeadTot= ((dead.plants+dead.disps+dead.pols)/ntot)*100, porcDeadPlants=(dead.plants/nplants)*100,
           porcDeadDisp=(dead.disps/ndisp)*100,porcDeadPols=(dead.pols/npol)*100) %>% 
    arrange(Species_ID)
  
  
  #assign the species role to node id 
  roles<- sprole %>%
    filter(Treat == "NI") %>% 
    group_by(node_id)%>% 
    slice(1) %>% 
    select(node_id,role)
  
  roles_m<- roles %>% 
    filter(node_id%in%Results$Species_ID) %>% 
    select(role)
  
Results_2<-cbind(Results,role = roles_m$role)

Final= rbind(Final,Results_2)
}

Prop_NI_size = Final #Propagation of 300 networks controlled by size

#######################################################
####### Create data frame and summarize results #######
#######################################################
NI<-Prop_NI_size%>% 
  group_by(Species_ID, initial.random,role) %>% 
  summarise(
        porc_dead_tot_mean = mean(porcDeadTot),
        porc_dead_plants_mean = mean(porcDeadPlants),
        porc_dead_disp_mean = mean(porcDeadDisp),
        porc_dead_pol_mean = mean(porcDeadPols))
NI_c<-cbind(Treat = rep("NI", nrow(NI)),NI)

names(NI_c)<- c("Treat","Species_ID","initial.random","role","porcDeadTot","porcDeadPlants","porcDeadDisp","porcDeadPols")

#-------------------------------------------------------------------------------
  ########### 2. Comparison with matrix of NI without controlling by size

#upload matrix of NI without controlling for size

NI_nc<-read.csv()#should upload the file containing the propagation in NI, generated in the code titled "Extinction_Propagation"

NI_sprole <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Modularity/NI_super_list.csv", sep=",",row.names=1)

#assign the species role to node id 
NI_roles<- NI_sprole %>%
  group_by(node_id)%>% 
  dplyr::slice(1) %>% 
  dplyr::select(role)

NI_ncr<-cbind(NI_nc,role = NI_roles$role)

NI_ncr2<-NI_ncr %>% 
  select(-c(dead.plants,dead.disps,dead.pols,degrees))

NI_ncr2<-NI_ncr2[,c(1,2,3,8,5,4,6,7)]#reorder columns

#create final dataframe with controlled and not controlled data
Final<- rbind(NI_c,NI_ncr2)
Final<-Final [,-1]
Final2<-cbind(Treat = rep (c("Controlled", "Uncontrolled"),each = 136), Final)

Final2$Species_ID = as.factor(Final2$Species_ID)
Final2$Treat<- factor(Final2$Treat, levels=c("Controlled", "Uncontrolled"))
Final2$role <- factor(Final2$role, levels=c("network hub","module hub","connector", "peripheral"))
Final2<-Final2[,-9]

#Model to compare percentage of extinct species according to uncontrolled vs controlled data 
library(lme4)
library("stats4")
library("bbmle")
library(car)
library(lsmeans)
library(multcomp)

m1<- glmer (porcDeadTot~ Treat +( 1| Species_ID), family = Gamma, data = Final2)
summary(m1)
anova(m1)

#Homogeneity
E_lme<-resid(m1, type= "deviance") 
F_lme<-fitted(m1) 
plot(x=F_lme, y=E_lme, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3)

##Independence 
boxplot(E_lme~Final2$Treat, main="Tratamiento")
boxplot(E_lme~Final2$role, main="roles") 

#Plot (Appendix S3)-
ext_order=ggplot(Final2)+geom_boxplot(aes(y=porcDeadTot,x=role,fill=Treat))+
  scale_fill_manual(values=c("#a69880", "cornflowerblue"))+
  scale_x_discrete(name="Topological role")+
  scale_y_continuous(name = "Percentage of extinct species")+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line (colour = "black",linetype = "solid", size = 1.4 ),
        axis.title = element_text(size= 15),
        axis.ticks = element_line(colour = "black", size = 1.4),
        axis.text = element_text (colour = "black", size = 14 ),
        legend.background = element_blank(),
        legend.key = element_rect(colour = "white", fill = NA))
ext_order

