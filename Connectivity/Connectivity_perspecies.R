###########Calculate the presence of plant species connecting layers and their number and proportion
# of paths (1)
########### Compare the presence and ability of plants to connect both layers between treatments (2)

library(dplyr)
library(vegan)
source('/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Extra_Functions/Extrafunctions.R') 

################## (1) ###############

# NI TREATMENT --
NI_disp <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_disp.csv", sep=",",row.names=1)
NI_pol<- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_pol.csv", sep=",",row.names=1)

Npol<-NI_pol
Ndisp<-NI_disp

#check 
Ndisp=Ndisp[apply(Ndisp[,], 1, function(x) !all(x==0)),] 
Npol=Npol[apply(Npol[,], 1, function(x) !all(x==0)),]
NpolE=empty(Npol)
NdispE=empty(Ndisp)

#Call function to change the network' edges to proportions
PropNetPol=ToProp(Npol) 
PropNetDisp=ToProp(Ndisp)

#Call function to create the network edge list
files=CreateFiles(PropNetDisp,PropNetPol) 
Edges=files$Edges.info
Nodes=files$Nodes.info
Layers=files$Layers.info

#Call function to calculate the number of path per plant species
plantesName= colnames(Ndisp)
wp=which(plantesName %in% Nodes$name_node)
Edges2=Edges %>% select(node_from,layer_from,node_to,layer_to,weight)
chau=which(Edges2$layer_from!=Edges2$layer_to)
Edges.new=Edges2[-chau,]#Remove the interlayer edges from the edge lost because
#the function "Connectivityindex" already do that

con_sp=Connectivityindex(wp,edges.list=Edges.new) #call Connectivityindex function
conn_sp_NI<-cbind(Treat = rep("NI",length(plantesName)), plantesName,con_sp) 
conn_sp_NI<-as.data.frame(conn_sp_NI)
conn_sp_NI[is.na(conn_sp_NI)]<-0
conn_sp_NI$N_pol<- as.integer(conn_sp_NI$N_pol)
conn_sp_NI$N_disp<- as.integer(conn_sp_NI$N_disp)
conn_sp_NI$P<- as.integer(conn_sp_NI$P)
conn_sp_NI<- conn_sp_NI %>% 
  mutate (Connect= ifelse(P == 0, 0,1), #Calculate the presence of plant connecting layer
          Prop_path = P / (nrow(NI_pol) * nrow(NI_disp) )) %>% #calculate theproportion of paths per plant sp
  select(Treat,plantesName,P,Connect,Prop_path)


# I TREATMENT --
I_disp <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_disp.csv", sep=",",row.names=1)
I_pol <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_pol.csv", sep=";",row.names=1)

Npol<-I_pol
Ndisp<-I_disp

#check
Ndisp=Ndisp[apply(Ndisp[,], 1, function(x) !all(x==0)),] 
Npol=Npol[apply(Npol[,], 1, function(x) !all(x==0)),]
NpolE=empty(Npol) 
NdispE=empty(Ndisp)

#Call function to change the network' edges to proportions
PropNetPol=ToProp(Npol) 
PropNetDisp=ToProp(Ndisp)

#Call function to create the network edge list
files=CreateFiles(PropNetDisp,PropNetPol) 
Edges=files$Edges.info
Nodes=files$Nodes.info
Layers=files$Layers.info

#Call function to the number of path per plant species
plantesName= colnames(Ndisp)
wp=which(plantesName %in% Nodes$name_node)
Edges2=Edges %>% select(node_from,layer_from,node_to,layer_to,weight)
chau=which(Edges2$layer_from!=Edges2$layer_to)
Edges.new=Edges2[-chau,]#Remove the interlayer edges from the edge lost because
#the function "Connectivityindex" already do that

con_sp=Connectivityindex(wp,edges.list=Edges.new) #call Connectivityindex function
conn_sp_I<-cbind(Treat = rep("I",length(plantesName)), plantesName,con_sp) 
conn_sp_I<-as.data.frame(conn_sp_I)
conn_sp_I[is.na(conn_sp_I)]<-0
conn_sp_I$N_pol<- as.integer(conn_sp_I$N_pol)
conn_sp_I$N_disp<- as.integer(conn_sp_I$N_disp)
conn_sp_I$P<- as.integer(conn_sp_I$P)
conn_sp_I<- conn_sp_I %>% 
  mutate (Connect= ifelse(P == 0, 0,1), #Calculate the presence of plants connecting layers
          Prop_path = P / (nrow(I_pol) * nrow(I_disp) )) %>% #calculate the proportion of paths per plant sp
  select(Treat,plantesName,P,Connect,Prop_path)

#final matrix
conn_sp<-rbind(conn_sp_NI,conn_sp_I)


################## (2) ###############
library(tidyverse)
library(lme4)
library("stats4")
library("bbmle")

# test if presence of plants connecting layers change between treatments ---
m1<-glmer (Connect ~ Treat + ( 1| plantesName), family = binomial(link="logit"), data = conn_sp)
summary(m1)

#Homogeneity
EM<-resid(m1, type= "deviance") 
FM<-fitted(m1) 
plot(x=FM, y=EM, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3) 

#independence 
E1_lme<-resid(m1, type= "deviance") 
boxplot(E1_lme~conn_sp$Treat, main="Tratamiento")


# test if Number and prop of path per species connecting layer change between treatment---
conn_sp <-conn_sp%>% 
  filter(P>0)#filter only plants connecting layers

#Number of path per species
m2<- glmer.nb (P ~ Treat + ( 1| plantesName), data = conn_sp)
summary(m2)

#overdispersion
overdisp_fun<- function(model) {
  vpars<- function(m) {
    nrow(m)*(nrow(m)+1)/2 
  }
  model.df<- sum(sapply(VarCorr(model),vpars))+length(fixef(model)) 
  (rdf<- nrow(model@frame)-model.df)
  rp<- residuals(model)
  Pearson.chisq<- sum(rp^2)
  prat<- Pearson.chisq/rdf
  pval<- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE,log.p=TRUE) 
  c(chisq=Pearson.chisq,ratio=prat,p=exp(pval))
}
overdisp_fun(m2)

#Homogeneity 
EM<-resid(m2, type= "deviance") 
FM<-fitted(m2) 
plot(x=FM, y=EM, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3) 

#Independence 
boxplot(EM~conn_sp$Treat, main="Tratamiento") 


#Prop of paths per species 
library(glmmTMB)

m3<-glmmTMB(Prop_path~ Treat + ( 1| plantesName),beta_family(link = "logit"), data = conn_sp)
summary(m3) 

#Homogeneity 
EM<-resid(m3, type= "pearson") 
FM<-fitted(m3)
plot(x=FM, y=EM, xlab = "Ajustados", ylab = "Residuales normalizados")
abline(0,0, col="red", lwd= 3) 

#Independence
boxplot(EM~conn_sp$Treat, main="Tratamiento")

