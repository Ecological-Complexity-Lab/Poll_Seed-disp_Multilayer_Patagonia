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



############################# AUC FUNCTION ###################################

auc=function(rem,pre){
  y <- pre
  x <- rem
  ext.curve <- splinefun(x, y)
  ext.area <- integrate(ext.curve, 0, 1)
  return(as.numeric(ext.area[[1]]))
}