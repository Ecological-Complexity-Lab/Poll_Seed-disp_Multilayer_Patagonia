#This file contains the code necesarry to simulate networks in NI (non-invaded treatment) with the same size as networks in I (invaded treatment)

##########Simulation of networks in non invaded treatment with the same size as networks in invaded treatment
#We removed species at random in the non-invaded site to reach the number of pollinators and plants (we have the same number of seed disperser)

library(tidyverse)

#Upload matrix
NI_disp <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_disp.csv", sep=",",row.names=1)
NI_pol  <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_pol.csv", sep=",",row.names=1)
I_disp<- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_disp.csv", sep=",",row.names=1)
I_pol <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_pol.csv", sep=";",row.names=1)


#Create the matrix of pollinators,plant, seed dispersers
Ndisp=NI_disp
Npol=NI_pol
Npol  <- Npol [ order(row.names(Npol)), ] #order pollinators
imatrix=rbind(Npol,Ndisp) 
full_matrix= as.matrix(imatrix)


#Remove species of each trophic group at random until get the number of species in the other treatment (we have the same number of seed disperser)

nsims = 450
row.names_matrix=c()
col.names_matrix=c()
res=vector(mode = "list", length = 300)
NN=1 

set.seed(15000)
for (s in 1:nsims){
  if (NN<=300) 
  {
  full_matrix_plants<-full_matrix[, sample(ncol(full_matrix), 24)] #we sampled random plants 

  out=full_matrix_plants[apply(full_matrix_plants[,], 1, function(x) all(x==0)),] #detect pollinators without interactions after removing plants
  pol.removed = nrow(out)#pollinators removed after sampling plants because don't have interactions
  to.continue.removing=full_matrix_plants[apply(full_matrix_plants[,], 1, function(x) !all(x==0)),] #matrix to  continue removing pollinators
  
  matrix_pol<-to.continue.removing[1:(nrow(NI_pol)-pol.removed),]#subset pollinators and plants
  matrix_pol<-as.data.frame(matrix_pol)
  full_matrix_pol<-matrix_pol[sample(nrow(matrix_pol),67),]  #we sampled pollinators
    
    full_matrix_pol  <- full_matrix_pol[order(row.names(full_matrix_pol)), ] #order pollinators
  
  full_matrix_disp <- NI_disp[,colnames(full_matrix_pol)]
  reduced_matrix <- rbind(full_matrix_pol,full_matrix_disp)
  reduced_matrix2<-reduced_matrix%>% 
    select(sort(colnames(.)))
  
  row_with_0= reduced_matrix2[apply(reduced_matrix2 [,], 1, function(x) all(x==0)),]#count rows having sum = 0

  if (nrow(row_with_0) == 0) 
  {
    res[[NN]]= reduced_matrix2
    NN=NN+1
  }
  
  }
  
  else
  {break}
  
}

