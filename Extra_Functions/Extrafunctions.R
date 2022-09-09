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
  
  #mismo orden para la plantas
  networkD=networkD %>% select(Plantas)
  
  n.plantas=length(Plantas)
  n.polinizadores=length(Polinizadores)
  n.dispersores=length(Dispersores)
  n.nodos=n.plantas+n.polinizadores+n.dispersores
  
  #Create edges list of plant-pollinator layer --
  
  nodeFrom=c()
  nodeTo=c()
  weight=c()
  
  for (i in 1:n.plantas)#para cada palnta for 1 a n, fijate los pesos que son positivos #y guardate quienes son esos polinizadores "quienes"
  {
    edges.tpm=networkP[,i]>0 #para la fila uno, ver todos los valores de las columnas que son amyores a 1
    quienes.tmp=which(edges.tpm)
    
    for (j in quienes.tmp) #para cada polinizadores que esta en "quienes",  #guardo el nombre dle nodo y el peso
    {
      nodeFrom=c(nodeFrom,i)
      nodeTo=c(nodeTo,j)
      weight=c(weight,as.numeric(networkP[,i][j]))
    }
    
  }
  
  nodeTo=nodeTo+n.plantas #lo que le digo ac? es dejarle un espacio para las plantas y
  #empezar del x (polinizadores van del x al tanto), y lo copiamos en un dataframe:
  
  Edges.List1=as.data.frame(cbind(nodeFrom,nodeTo,weight))
  
#Create edge list of plant-seed disperser layer --
  
  nodeFrom=c()
  nodeTo=c()
  weight=c()
  
  for (i in 1:n.plantas)
  {
    edges.tpm=networkD[,i]>0
    quienes.tmp=which(edges.tpm)
    
    for (j in quienes.tmp)
    {
      nodeFrom=c(nodeFrom,i)
      nodeTo=c(nodeTo,j)
      weight=c(weight,as.numeric(networkD[,i][j]))
    }
    
  }
  
  nodeTo=nodeTo+n.plantas+n.polinizadores
  
  Edges.List2=as.data.frame(cbind(nodeFrom,nodeTo,weight))
  
  
  ###agregamos el ID de las capas a cada lista intra
  
  Nlayer1=rep(1,nrow(Edges.List1)) #creamos el vector queva a corresponder a la columna layer 1 (polinizacion) (rep significa que se repite el valor 1 -capa 1- en todas las filas; nrow es la cantidad filas que va a tener el vector ( = cant. de conexiones ne la primer capa
  Nlayer2=rep(2,nrow(Edges.List2))
  
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
  
 
  #Create weigthed interlayer edges (Prop path = Npath/Total path) --
  
  #select the plant present in both layers 
  Qplan1=which(vapply(networkD,sum,1)>0)
  Qplan2=which(vapply(networkP,sum,1)>0)
  
  #creamos el vector quienes inter con la identidad de las plantas y cuantificamos
  #las plantas que conectan ambos procesos
  quienes.inter=Qplan2[Qplan2%in%Qplan1] 
  P = vector(mode="numeric", length((quienes.inter)))# aca creamos el vector para guardar el numero de pol de cada planta que une ambas capas
  D = vector(mode="numeric", length((quienes.inter)))#vector donde guardamos el Ndisp con los que interactua cada planta que unen ambas capas
  
   #Loop para calcular el peso entre capa de cada planta compartidas

  for (j in seq_along(quienes.inter)){
    nodeP=quienes.inter[j]
    P[j]=as.numeric(Extend1 %>% filter(node_from ==nodeP)  %>% summarise(n = n())) #calculamos el num de pol que interactua con cada planta que conecta ambas capas
    D[j]=as.numeric(Extend2 %>% filter(node_from ==nodeP)  %>% summarise(n = n())) #calculamos el num de disp que interactua con cada planta que conecta ambas capas
  } 
  
  W = P*D #Calculate the number of path per plant species connecting both layers
  
  cuan.inter=length(quienes.inter)# create the vector containing the Num path/Total path
  for (i in seq_along(W)){
    cuan.inter[i] = W[i]/ (n.polinizadores*n.dispersores) 
  }
  
  #creamos el archivo de interacciones inter
  extendInter=data.frame("layer_from"=as.integer(rep(1,length(cuan.inter))),'node_from'=as.integer(quienes.inter),
                         'layer_to'=as.integer(rep(2,length(cuan.inter))),'node_to'=as.integer(quienes.inter),
                         "weight"=as.numeric(cuan.inter))
  
  #unimos enlaces intra y inter en un archivo
  AllExtend=rbind(Extend1,Extend2,extendInter)
  
  ### Layers info 
  layerID=c(1,2)
  layerLabel=c("PlantaPol","PlantaDisp")
  LayerInfo=data.frame('layer_id'=as.integer(layerID),'name_layer'=layerLabel)
  
  
  ### Nodes info 
  
  nodeID=seq(1,n.nodos) #  Primeros plantas , luego polinizadores , luego dispersores 
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