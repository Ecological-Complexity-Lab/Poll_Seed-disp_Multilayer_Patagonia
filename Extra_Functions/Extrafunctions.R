#This file contain

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      1. CONNECTIVITY INDEX                            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function calculates the number of path per plant species connecting layers and the number of pollinators and seed dispersers 
#interacting with the plant 

Connectivityindex=function(wp,edges.list)
{
  Nplants= length(wp) #vector containing plant species
  P=numeric(Nplants) #Number of paths per plant
  N_pol= numeric(Nplants) #Number of pollinators interacting with the plant
  N_disp= numeric(Nplants)#Number of seed dispersers interacting with the plant
  
  for (j in 1:Nplants) 
  {
    nodeP=wp[j]
    Nperlayer=edges.list %>% filter(node_from==nodeP) %>% group_by(layer_to) %>% dplyr::summarise(N=n())
    Nperlayer=Nperlayer$N 
    if (length(Nperlayer)>1) #if the plant interacts with more than one layer
    {
      P[j]=Nperlayer[1]*Nperlayer[2] #calculate its number of paths between layers
    }
    N_pol[j]= Nperlayer[1]
    N_disp[j]= Nperlayer[2]
  }
  conectivity_sp<-cbind(P,N_pol,N_disp)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #                      2. TO PROPORTION FUNCTION                      
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function change the weight of links between species in the network to proportions

ToProp=function(IntNetwork) 
{
  NN=sum(IntNetwork)
  out=IntNetwork/NN
  return(out)
}  



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #                      3. CREATE FILES FUNCTION                      
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function create an edge list containing intralayer and interlayer interaction between species
  
CreateFiles=function(networkD,networkP)

{ 
  require(tidyverse)
  Plantas=names(networkP)
  Polinizadores=row.names(networkP)
  Dispersores=row.names(networkD)
  
  networkD=networkD %>% select(Plantas)
  
  n.plantas=length(Plantas)
  n.polinizadores=length(Polinizadores)
  n.dispersores=length(Dispersores)
  n.nodos=n.plantas+n.polinizadores+n.dispersores
  
  #1. Create edge list ---------------------------
  
  #Create edges list of plant-pollinator layer ----
  nodeFrom=c()
  nodeTo=c()
  weight=c()
  
  for (i in 1:n.plantas)
  { #Save pollinators that have at least one interaction
    edges.tpm=networkP[,i]>0 
    quienes.tmp=which(edges.tpm)
    
    for (j in quienes.tmp) 
    {#Save the name of the plant interacting with each pollinator and the weight of the interaction
      nodeFrom=c(nodeFrom,i)
      nodeTo=c(nodeTo,j)
      weight=c(weight,as.numeric(networkP[,i][j]))
    }
    
  }
  
  nodeTo=nodeTo+n.plantas #Correct for the pollinator id
  
  Edges.List1=as.data.frame(cbind(nodeFrom,nodeTo,weight))#data frame containing intralayer edges of plant-pollination layer
  
#Create edge list of plant-seed disperser layer ----
  nodeFrom=c()
  nodeTo=c()
  weight=c()
  
  for (i in 1:n.plantas)
  {#Save seed dispersers that have at least one interaction
    edges.tpm=networkD[,i]>0
    quienes.tmp=which(edges.tpm)
    
    for (j in quienes.tmp)
    {#Save the name of the plant interacting with each seed disperser and the weight of the interaction
      nodeFrom=c(nodeFrom,i)
      nodeTo=c(nodeTo,j)
      weight=c(weight,as.numeric(networkD[,i][j]))
    }
    
  }
  
  nodeTo=nodeTo+n.plantas+n.polinizadores #Correct for the seed disperser id
  
  Edges.List2=as.data.frame(cbind(nodeFrom,nodeTo,weight)) #data frame containing intralayer edges of plant-seed disperser layer
  
  
  #we create vectors with the ID of each layer (pollination layer = 1, seed disperser layer = 2)
  Nlayer1=rep(1,nrow(Edges.List1)) 
  Nlayer2=rep(2,nrow(Edges.List2))
  
  #add the vectors to the edge lists created before
  Extend1=data.frame('layer_from'=as.integer(Nlayer1),
                     "node_from"=as.integer(Edges.List1[,1]),
                     'layer_to'=as.integer(Nlayer1),
                     'node_to'=as.integer(Edges.List1[,2]),
                     'weight'=Edges.List1[,3])

  Extend2=data.frame('layer_from'=as.integer(Nlayer2),
                     "node_from"=as.integer(Edges.List2[,1]),
                     'layer_to'=as.integer(Nlayer2),
                     'node_to'=as.integer(Edges.List2[,2]),
                     'weight'=Edges.List2[,3])
  
 
  #Create interlayer edge list ----
  
  #select those plants  interacting within each layer
  Qplan1=which(vapply(networkD,sum,1)>0)
  Qplan2=which(vapply(networkP,sum,1)>0)
  
  #identify those plants common to both layers
  quienes.inter=Qplan2[Qplan2%in%Qplan1] 
  P = vector(mode="numeric", length((quienes.inter)))
  D = vector(mode="numeric", length((quienes.inter)))
  
  for (j in seq_along(quienes.inter)){#calculate the number of pol and seed disperser interacting with plants connecting layers
    nodeP=quienes.inter[j] 
    P[j]=as.numeric(Extend1 %>% filter(node_from ==nodeP)  %>% summarise(n = n())) 
    D[j]=as.numeric(Extend2 %>% filter(node_from ==nodeP)  %>% summarise(n = n())) 
  } 
  
  W = P*D #number of path per plant species connecting both layers
  
  cuan.inter=length(quienes.inter)
  for (i in seq_along(W)){
    cuan.inter[i] = W[i]/ (n.polinizadores*n.dispersores)#calculate the proportion of path per plant species 
  }
  
  #dataframe containing interlayer edges
  extendInter=data.frame("layer_from"=as.integer(rep(1,length(cuan.inter))),'node_from'=as.integer(quienes.inter),
                         'layer_to'=as.integer(rep(2,length(cuan.inter))),'node_to'=as.integer(quienes.inter),
                         "weight"=as.numeric(cuan.inter))
  
  #join all dataframe together
  AllExtend=rbind(Extend1,Extend2,extendInter)
  
  
  #2. layer info ---------------------------
  layerID=c(1,2)
  layerLabel=c("PlantaPol","PlantaDisp")
  LayerInfo=data.frame('layer_id'=as.integer(layerID),'name_layer'=layerLabel)
  
  #3. layer info ---------------------------
  nodeID=seq(1,n.nodos) 
  nodeLabel=c(Plantas,Polinizadores,Dispersores)
  NodesInfo=data.frame('node_id'=nodeID,'name_node'=nodeLabel)
  
  return(list(Edges.info=AllExtend,Layers.info=LayerInfo,Nodes.info=NodesInfo))




}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #                      4. SHUFFLING FUNCTION                      
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#This function shuffle the interactions of the empirical network to create 1000 random networks according to the "r2dtable" algorithm (see main manuscript for more details). 
#For each generated network, we calculated the number of modules and L value.

Shuffle <- function(Npol,Ndisp)
  
{ require(tidyverse)
  require(infomapecology)
  
  nsim=1000
  
  nmpol <- vegan::nullmodel(Npol, method = "r2dtable") #create the objects to calculate the parameters
  nmdisp <- vegan::nullmodel(Ndisp, method = "r2dtable")
  
  set.seed(11)
  
  #Simulation
  simupol <- simulate(nmpol, nsim=nsim+100) 
  simudisp <- simulate(nmdisp, nsim=nsim+100)
  
  #Create vectors to store the output of each simulation
  L=numeric(nsim)
  M=numeric(nsim)
  x_modules = list()
  NN=1 
  
  for (i in 1:(nsim+100)) 
  {
    if (NN<=nsim) #if the number of outputs is smaller than 1000, continue calculating modularity
    {
      IntNetPol=simupol[,,i] 
      IntNetDisp=simudisp[,,i] 
      
      IntNetPol=as.data.frame(IntNetPol)
      IntNetPol=filter_all(IntNetPol, any_vars(. != 0)) #check of row with 0 
      
      IntNetDisp=as.data.frame(IntNetDisp)
      IntNetDisp=filter_all(IntNetDisp, any_vars(. != 0))  
      
      #Call function to change the network' edges to proportions
      PNetPol=ToProp(IntNetPol)
      PNetDisp=ToProp(IntNetDisp)
      
      #Call function to create the network edge list
      files=CreateFiles(PNetDisp,PNetPol)
      Edges=files$Edges.info
      Nodes=files$Nodes.info
      Layers=files$Layers.info
      
      Nodes<-cbind(Nodes, Type = rep(c("Plant","Pol","Disp"),times =c(ncol(IntNetPol), nrow(IntNetPol),nrow(IntNetDisp)))) #we added the trophic group of species to the dataframe
      
      #Calculate number of modules and L
      Network= create_multilayer_object(extended = Edges,nodes=Nodes,
                                        intra_output_extended = TRUE,layers=Layers) #create a multilayer object using the edge list generated before
      
      salida=run_infomap_multilayer(Network, relax = F, silent = T, trials = 100, 
                                    temporal_network = F,flow_model='undirected') #output containing list of modules with the species integrating each module and L
      
     
      if (!is.na(salida$m)) #Save the value of number of modules and L when the number of modules is not Na
      {
        L[NN]=salida$L 
        M[NN]=salida$m
        x_modules = salida$modules
        NN=NN+1
      }
      
    }
    
    else
    {break} 
  }
  return(list(Sim_Results=data.frame("M"=as.integer(M),"L"=L),
           Role = data.frame (x_modules)))     
  }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      5.  Stochastic coextinction model                       
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#This function contain the stochastic coextinction model (SCM), modified from (Vieira & Almeida-Neto 2015) to include three different trophic groups,
#to simulate species extinction (see main manuscript for more detail)


netcascade_multi <- function(imatrix,rplants_disp,rplants_pol, rpol, rdisp, deadPlants=NULL, deadPols=NULL, 
                             
                             deadDisps=NULL, targetGuild,target,return.matrix=F){
  
  
  #---------ARGUMENT CHECKS-----------------------------------
  
  if(class(imatrix)!="matrix" || (nrow(imatrix)+ncol(imatrix))<3){stop("'imatrix' must be an object of class 'matrix', with animal species on rows, plant species on columns and at least three species overall")}
  if(class(rplants_disp)!="numeric" || class(rplants_pol) != "numeric" || class(rpol)!="numeric" || class(rdisp) != "numeric" ||  max(c(max(rplants_disp),max(rplants_pol),max(rdisp),max(rpol)))>1 || min(c(min(rplants_disp),min(rplants_pol),min(rdisp),min(rpol)))<0) {stop("'r' must be numeric vectors with values ranging between 0 and 1")}
  if(targetGuild%in%c("pollinator","disperser","plant")==F){stop('Invalid target guild for primary extinction. Valid targets guilds are "pol", "disp" and "plant"')}
  if(is.numeric(target)==F){stop('Invalid value for the "target" argument. You may specify a single species by entering its row or column number or you may use a vector of relative probabilites for all species in the target guild.')}
  if(is.null(deadPols)==F && class(deadPols)!= "integer"){stop("deadPols must be either NULL or an integer vector specifying the row numbers of pollinators considered to be extinct on the original matrix")}
  if(is.null(deadPlants)==F && class(deadPlants)!= "integer"){stop("deadPlants must be either NULL an integer vector specifying the column numbers of plants considered to be extinct on the original matrix")}
  if(is.null(deadDisps)==F && class(deadDisps)!= "integer"){stop("deadDisp must be either NULL an integer vector specifying the seed disperser considered to be extinct on the original matrix")}
  if(length(rpol)!= length(pol_indx)){stop("The length of vector'rpol' must be equal to the number of pollinator species")}
  if(length(rdisp)!= length(disp_indx)){stop("The length of vector'rdisp' must be equal to the number of seed disperser species")}   
  if(length(rplants_pol)!= ncol(imatrix)){stop("The length of vector'rplants_pol' must be equal to number of columns (i.e. plant species) in 'imatrix'")}
  if(length(rplants_disp)!= ncol(imatrix)){stop("The length of vector'rplants_disp' must be equal to number of columns (i.e. plant species) in 'imatrix'")}
  
  #---------DEFINING SOME VARIABLES---------------------------
  poll = pol_indx
  disp = disp_indx 
  plants = 1:nplants
  
  pollNA = pol_indx
  dispNA = disp_indx
  plantsNA = 1:nplants
  
  pollNA[deadPols] <- NA
  dispNA[deadDisps] <- NA
  plantsNA[deadPlants] <- NA
  
  degree_when_lost_pols <- c()
  degree_when_lost_disps<- c()
  degree_when_lost_plants<- c()
  
  #----------CALCULATING DEPENDENCE MATRICES------------------- 
  M_pp<- matrix(0, npol,nplants)#create a matrix of dependence between plants and pol
  M_pd<- matrix(0, ndisp,nplants)#create a matrix of dependence between plant- seed disperser
  M_anim<- matrix(0, nanim,nplants)#create a matrix of dependence between animals and plants
  
  for(i in 1:nplants){
    M_pp[,i]<-imatrix[pol_indx,i]/sum(imatrix[pol_indx,i])
  }#matrix of plant dependence on each pollinator
  for(i in 1:nplants){
    M_pd[,i]<-imatrix[disp_indx,i]/sum(imatrix[disp_indx,i])
  }#matrix of plant dependence on seed disperser
  
  M_pp[is.nan(M_pp)] <- 0 
  M_pd[is.nan(M_pd)] <- 0 
  M_plants<-rbind(M_pp,M_pd) #matrix of plant dependence on animal species
  
  for(i in 1:nanim){
    M_anim[i,]<-imatrix[i,]/sum(imatrix[i,])
  }
  
  #-----------CHOOSING TARGET SPECIES FOR PRIMARY EXTINCTION---
  coext_pols <- c()
  coext_disps <- c()
  coext_plants <- c()
  
  if(length(target)==1){
    if(targetGuild=="pollinator"){
      if(target %in% deadPols){stop('Specified target species for the primary extinction is already extinct')}
      coext_pols <- target
      degree_when_lost_pols <- 1 #stores the degree of the extinction event of every pollinator species lost during the coextinction cascade. 
    }
    if(targetGuild=="disperser"){
      if(target %in% deadDisps){stop('Specified target species for the primary extinction is already extinct')}
      coext_disps <- target
      degree_when_lost_disps <- 1 #stores the degree of the extinction event of every seed disperser species lost during the coextinction cascade. 
    }
    if(targetGuild=="plant"){
      if(target %in% deadPlants){stop('Specified target species for the primary extinction is already extinct')}
      coext_plants <- target
      degree_when_lost_plants <- 1#stores the degree of the extinction event of every seed plant species lost during the coextinction cascade. 
    }
  }else{
    nspecies <- switch(targetGuild,pollinator = npol, disperser = ndisp, plant = nplants)
    if(length(target)==nspecies){
      if(targetGuild =="pollinator"){
        alive <- poll[is.na(pollNA)==F]
        coext_pols <- sample(c(alive,0),1,prob = c(target[is.na(pollNA)==F],0))
        degree_when_lost_pols <- 1
      }
      if(targetGuild =="disperser"){
        alive <- disp[is.na(dispNA)==F]
        coext_disps <- sample(c(alive,0),1,prob = c(target[is.na(dispNA)==F],0))
        degree_when_lost_disps <- 1
      }
      
      if(targetGuild =="plant"){
        alive <- plants[is.na(plantsNA)==F]
        coext_plants <- sample(c(alive,0),1,prob = c(target[is.na(plantsNA)==F],0))
        degree_when_lost_plants <- 1
      }
    }else{
      stop('Length of "target" must be 1 (specifying a single species within the target guild) or else be equal to the number of species in the target guild (specifying probabilities of primary extinction for each species in the target guild)')
    }
  }
  
  
  imatrix[coext_pols,] <- 0 
  imatrix[coext_disps,]<- 0 
  imatrix[,coext_plants] <- 0
  
  #final list of species (pollinators, seed dispersers and plants) which were "alive" in the original community but became extinct during this primary extinction + extinction cascade
  lostpols <- coext_pols #pollinators
  lostdisps<- coext_disps#seed dispersers
  lostplants <- coext_plants#plants
  
  
  #-------------------CASCADE LOOP--------------------------- 
  equilibrium <- FALSE
  degree <- 1
  degree_table <- data.frame(degree,guild=factor(targetGuild,levels=c("pollinator","disperser","plant")),n_extinctions=1)
  while(equilibrium == FALSE){
    ext_pols<- coext_pols
    ext_disps <- coext_disps
    ext_plants <- coext_plants
    
    pollNA[ext_pols] <- NA
    dispNA[ext_disps] <- NA
    plantsNA[ext_plants] <- NA
    
    polleft <- poll[is.na(pollNA)==F]
    dleft <- disp[is.na(dispNA)==F]
    pleft <- plants[is.na(plantsNA)==F]
    
    #hacer esto para un polinzadores y calcular las plantas que mueren (leer mejor estos dos loops) 
    
    #If a pollinator get extinct, then next species in going extinct is a plant
    if(length(ext_pols)>0){
      for(i in 1:length(ext_pols)){
        unlucky <- rplants_pol[pleft]*M_pp[ext_pols[i],pleft] > runif(length(pleft))
        coext_plants = c(coext_plants,pleft[unlucky])
      }
      coext_pols <- c()
      coext_plants <- unique(coext_plants)
      plantsNA[coext_plants] <- NA
      lostplants <- c(lostplants,coext_plants)
      imatrix[,coext_plants] <- 0
      
      for(i in 1:nplants){
        if(sum(imatrix[pol_indx,i])==0){
          M_plants[,i] <- 0#if the plant went extinct, assign a 0
        }else{
          M_plants[pol_indx,i] <- imatrix[pol_indx,i]/sum(imatrix[pol_indx,i])#otherwise, calculate the dependence again
        }
      } 
      if(length(coext_plants)>0){
        degree <- degree + 1
        degree_when_lost_plants <- c(degree_when_lost_plants, rep(degree,length(coext_plants)))
        degree_table[degree,] <- data.frame(degree,"plant",length(coext_plants))
      }
    }
    
    #If a seed disperser get extinct, then next species in going extinct is a plant
    if (length(ext_disps)>0){
      for(i in 1:length(ext_disps)){
        unlucky <- rplants_disp[pleft]*M_pd[ext_disps[i]-npol,pleft] > runif(length(pleft))
        coext_plants = c(coext_plants,pleft[unlucky])
      }
      coext_disps <- c()
      coext_plants <- unique(coext_plants)
      plantsNA[coext_plants] <- NA
      lostplants <- c(lostplants,coext_plants)
      imatrix[,coext_plants] <- 0
      
      for(i in 1:nplants){
        if(sum(imatrix[disp_indx,i])==0){
          M_plants[,i] <- 0#if the plant went extinct, assign a 0
        }else{
          M_plants[disp_indx,i] <- imatrix[disp_indx,i]/sum(imatrix[disp_indx,i])#otherwise, calculate the dependence again
        }
      }  
      if(length(coext_plants)>0){
        degree <- degree + 1
        degree_when_lost_plants <- c(degree_when_lost_plants, rep(degree,length(coext_plants)))
        degree_table[degree,] <- data.frame(degree,"plant",length(coext_plants))
      }
    }
    
    #If a plant get extinct, then next species in going extinct is a pollinator and/or a seed disperser
    else{ 
      for(i in 1:length(ext_plants)){
        unlucky_pol <- rpol[polleft]*M_anim[polleft,ext_plants[i]] > runif(length(polleft))
        unlucky_disp <- rdisp[dleft-npol]*M_anim[dleft,ext_plants[i]] > runif(length(dleft-npol))
        coext_pols <- c(coext_pols, polleft[unlucky_pol])
        coext_disps <- c(coext_disps, dleft[unlucky_disp])
      }

      coext_plants = c()
      coext_pols <- unique(coext_pols) 
      coext_disps <- unique(coext_disps)
      lostpols <- c(lostpols,coext_pols)
      lostdisps <- c(lostdisps,coext_disps)
      pollNA[coext_pols] <- NA
      dispNA[coext_disps] <- NA
      imatrix[coext_pols,] <- 0 
      imatrix[coext_disps,] <- 0 
      
      for(i in 1:nanim){
        if(sum(imatrix[i,])==0){
          M_anim[i,] <- 0 #if the animal went extinct, assign a 0
        }else{
          M_anim[i,] <- imatrix[i,]/sum(imatrix[i,])#otherwise, calculate the dependence again
        }
      }
      if(length(coext_pols)>0){
        degree <- degree + 1
        degree_when_lost_pols <- c(degree_when_lost_pols, rep(degree,length(coext_pols)))
        degree_table[degree,] <- data.frame(degree,"pollinator",length(coext_pols))  
      }
      if(length(coext_disps)>0){
        degree <- degree + 1
        degree_when_lost_disps <- c(degree_when_lost_disps, rep(degree,length(coext_disps)))
        degree_table[degree,] <- data.frame(degree,"disperser",length(coext_disps))  
      }
    }
    equilibrium <- equilibrium + (length(coext_plants)+length(coext_pols)+length(coext_disps))==0
  }
  
  #-------------------OUTPUT---------------------------
  
  if(return.matrix==T){
    return(list(interaction_matrix = imatrix, lost_pols = lostpols, lost_disps = lostdisps,lost_plants = lostplants))  
  }else{
    if(length(lostpols)>0){
      spp_data_pols<- data.frame(lost_pols = lostpols,degree_of_extinction=degree_when_lost_pols)
    }else{
      spp_data_pols <- "No pollinator species were lost"
    }
    if(length(lostdisps)>0){
      spp_data_disps <- data.frame(lost_disps = lostdisps,degree_of_extinction=degree_when_lost_disps)
    }else{
      spp_data_disps <- "No seed disperser species were lost"
    }
    
    if(length(lostplants)>0){
      spp_data_plants <- data.frame(lost_plants = lostplants,degree_of_extinction=degree_when_lost_plants)
    }else{
      spp_data_plants <- "No plant species were lost"
    }
    return(list(cascade_data=degree_table,pol_species_data=spp_data_pols,disp_species_data=spp_data_disps,
                plant_species_data=spp_data_plants))
  }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      6. AUC FUNCTION                           
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#This function calculate the area under the curve (AUC, proxy the robustness) of the species lost against species present in the network
auc=function(rem,pre){
  y <- pre
  x <- rem
  ext.curve <- splinefun(x, y)
  ext.area <- integrate(ext.curve, 0, 1)
  return(as.numeric(ext.area[[1]]))
}