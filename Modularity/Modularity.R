########### 1. Estimation of the number of modules in the network for each treatment 

########### 2. Comparison of the observed number of modules to values generated by shuffling 
# the network with a null model 

library(infomapecology)
library(dplyr)
library(vegan)
source('./Extra_Functions/Extrafunctions.R')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      NI TREATMENT                            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###### 1. Estimation of the number of modules in the network for each treatment--
NI_disp <- read.csv("./Data/NI_disp.csv", sep=",",row.names=1)
NI_pol  <- read.csv("./Data/NI_pol.csv", sep=",",row.names=1)

Npol<-NI_pol
Ndisp<-NI_disp

#Call function to change the network' edges to proportions
PropNetPol=ToProp(Npol) 
PropNetDisp=ToProp(Ndisp)

#Call function to create the network edge list
files=CreateFiles(PropNetDisp,PropNetPol)
Edges=files$Edges.info
Nodes=files$Nodes.info
Layers=files$Layers.info

#Calculate the real value of map equation (number of modules and L) using infomapecology package.
Multi_object= create_multilayer_object(extended = Edges,nodes=Nodes,
                                      intra_output_extended = TRUE,layers=Layers) #create a multilayer object using the edge list generated before

Salida_mod=run_infomap_multilayer(Multi_object, relax = F, silent = T, trials = 100,seed=1111, 
                                  temporal_network = F,flow_model='undirected') #output containing list of modules with the species integrating each module and L.
NI_Num_Modules=Salida_mod$m #number of modules
NI_L = Salida_mod$L #L

#Plot--
#this plot contains information about the occurrence of a species module in a layer and the proportion 
#of species in the network that were assigned to the module. 

#Preparing dataframe-

#Calculate the number and proportion of species in each module
N_species_mod<-Salida_mod$modules %>% 
  group_by(module) %>% 
  summarize(module_size=n_distinct(node_id))
N_species_mod<-N_species_mod %>% 
  mutate(Prop_sp = module_size / sum (N_species_mod$module_size))

#Select just the plants and add the identity of layers
Species_modules<- Salida_mod$modules %>%
  filter(between (node_id,1,ncol(NI_pol))) %>% 
  group_by(module) %>% 
  mutate(Layer = if_else (layer_id==1, "Pol","Disp"))%>% 
  group_by(module) %>% 
  mutate(Shared = n_distinct(Layer)) %>% 
  mutate(Layer = ifelse(Shared == 2, "Shared",Layer))
Species_modules$module<- as.numeric(Species_modules$module)

#add proportions of species in each module to the data frame
Species_modules2<-Species_modules %>% 
  mutate(Prop = case_when(
    module == 1 ~ as.numeric(N_species_mod[1,3]),
    module == 2 ~ as.numeric(N_species_mod[2,3]),
    module == 3 ~ as.numeric(N_species_mod[3,3]),
    module == 4 ~ as.numeric(N_species_mod[4,3]),
    module == 5 ~ as.numeric(N_species_mod[5,3]),
    module == 6~ as.numeric(N_species_mod[6,3]),
    module == 7 ~ as.numeric(N_species_mod[7,3]),
    module == 8 ~ as.numeric(N_species_mod[8,3]),  
    module == 9 ~ as.numeric(N_species_mod[9,3]),
    module == 10~ as.numeric(N_species_mod[10,3]),
    module == 11 ~ as.numeric(N_species_mod[11,3]),
    module == 12 ~ as.numeric(N_species_mod[12,3]),
    module == 13 ~ as.numeric(N_species_mod[13,3]),
    module == 14 ~ as.numeric(N_species_mod[14,3]),
    module == 15 ~ as.numeric(N_species_mod[15,3]),
    module == 16 ~ as.numeric(N_species_mod[16,3]), 
    module == 17 ~ as.numeric(N_species_mod[17,3]),
    module == 18~ as.numeric(N_species_mod[18,3]), 
  ))

#change the order of x axis according to module size
Species_modules3<-Species_modules2 %>% 
  select(module,layer_id,Prop) %>% 
  arrange(desc(Prop))
Species_modules3$module <- factor(Species_modules3$module, level = c("3","1","5","7","12","4","14","17","15","13","8","16",
                                                                     "11","10","2","9","6","18"))

#Plot-
NI<-ggplot(Species_modules3,aes(x = module, y = layer_id, fill = Prop))+
  geom_point(size=9, shape = 22)+
  scico::scale_fill_scico(limits = c(0, 0.31),palette = "lajolla")+
  scale_y_continuous(limits = c(0.8,2.2),breaks=c(1, 2))+
  labs(x='Module ID', y='Layer')+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line (colour = "black",linetype = "solid", size = 1 ),
        axis.text = element_text(size=18, color='black'),
        axis.title = element_text(size=25, color='black'), 
        title = element_text(size = 10, color = 'black'))+
  coord_fixed(ratio = 5)
NI

###### 2. Comparison of the observed number of modules to values generated by shuffling the network with a null model-- 

#Call the function "Shuffle" to generated shuffled networks and calculate modularity for each simulation
Sim_Results<-Shuffle(Npol,Ndisp)

shapiro.test(Sim_Results$Sim_Results$M) #testing for normality

t.test(Sim_Results$Sim_Results$M, mu = NI_Num_Modules) #two tailed - t test

#Plot histogram (Appendix S1)--
Mod.null_NI<- ggplot(Sim_Results$Sim_Results)+geom_bar(aes(x=M), fill = "cornflowerblue")+
  theme_bw() +geom_vline(xintercept =NI_Num_Modules,linetype = "dashed")+
  theme(axis.line = element_line (colour = "black",linetype = "solid", size = 1 ),
        axis.title = element_text(size= 25),
        axis.text= element_text(size= 18, colour = "black"))
Mod.null_NI


# Plot NI Multilayer networks ----                         
#~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot using igraph a quapartide network for both biprtite layers
# polls -> plants1 -- plants2 -> destrib

# The groups:
# plants: 1-37
# pollinators: 38-132
# seed distributes: 133-136

pols <- Nodes %>% filter(node_id > 37 & node_id < 133)  %>% 
        add_column(colour="orange") %>%
        mutate(label = node_id, node_id = paste("1_", node_id, sep = ""))
plnt1 <- Nodes %>% filter(node_id <= 37) %>% 
        add_column(colour="palegreen") %>%
        mutate(label = node_id, node_id = paste("1_", node_id, sep = ""))
plnt2 <- Nodes %>% filter(node_id <= 37) %>% 
        add_column(colour="palegreen") %>%
        mutate(label = node_id, node_id = paste("2_", node_id, sep = ""))
dist <- Nodes %>% filter(node_id >= 133)  %>% 
        add_column(colour="plum1")%>%
        mutate(label = node_id, node_id = paste("2_", node_id, sep = ""))

pols$x <- 0
plnt1$x <- 0.5
plnt2$x <- 1
dist$x <- 1.5

pols$y <- (nrow(pols):1)*20/nrow(pols) 
plnt1$y <- (nrow(plnt1):1)*20/nrow(plnt1)
plnt2$y <- (nrow(plnt2):1)*20/nrow(plnt2)
dist$y <- (nrow(dist):1)*5/nrow(dist) + 7.5

flat_nodes <- rbind(plnt1, plnt2, pols, dist)

flat_edges <- Edges %>% 
        mutate(colour=case_when(layer_from==layer_to ~ 'black',
                                layer_from!=layer_to ~ 'steelblue1')) %>%
        mutate(from = paste(layer_from, "_", node_from, sep = "")) %>%
        mutate(to = paste(layer_to, "_", node_to, sep = "")) %>% 
        select(from, to, weight, colour)

net <- igraph::graph_from_data_frame(flat_edges, directed = FALSE, vertices = flat_nodes)
pdf("./Modularity/ni_mln.pdf")
plot.igraph(net,  axes = FALSE, vertex.frame.color = "gray20",
            vertex.label = NA,
            vertex.size = 5,
            vertex.color = flat_nodes$colour,
            edge.width = flat_edges$weight*30,
            edge.color = flat_edges$colour,
            main = "Not Invaded Network")#,
#            margin = c(-0.3,-0.15,-0.2,-0.15))
dev.off()

# ----

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      I TREATMENT                            
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###### 1. Estimation of the number of modules in the network for each treatment--
I_disp<- read.csv("./Data/I_disp.csv", sep=",",row.names=1)
I_pol <- read.csv("./Data/I_pol.csv", sep=";",row.names=1)

Npol<-I_pol
Ndisp<-I_disp

#Call function to change the network' edges to proportions
PropNetPol=ToProp(Npol) 
PropNetDisp=ToProp(Ndisp)

#Call function to create the network edge list
files=CreateFiles(PropNetDisp,PropNetPol) 
Edges=files$Edges.info
Nodes=files$Nodes.info
Layers=files$Layers.info

#Calculate the real value of map equation (number of modules and L) using infomapecology package.
Multi_object= create_multilayer_object(extended = Edges,nodes=Nodes,
                                       intra_output_extended = TRUE,layers=Layers)#create a multilayer object using the edge list generated before

Salida_mod=run_infomap_multilayer(Multi_object, relax = F, silent = T, trials = 100,seed=1111, 
                                  temporal_network = F,flow_model='undirected')#output containing list of modules with the species integrating each module and L.

I_Num_Modules=Salida_mod$m #Number of modules
I_L = Salida_mod$L #L

#Plot--
#this plot contains information about the occurrence of a species module in a layer and the proportion 
#of species in the network that were assigned to the module. 

#Preparing dataframe-

#Calculate the number and proportion of species in each module
N_species_mod<-Salida_mod$modules %>% 
  group_by(module) %>% 
  summarize(module_size=n_distinct(node_id))
N_species_mod<-N_species_mod %>% 
  mutate(Prop_sp = module_size / sum (N_species_mod$module_size))

#Select just the plants and add the identity of layers
Species_modules_I<- Salida_mod$modules %>%
  filter(between (node_id,1,ncol(I_pol))) %>% 
  group_by(module) %>% 
  mutate(Layer = if_else (layer_id==1, "Pol","Disp"))%>% 
  group_by(module) %>% 
  mutate(Shared = n_distinct(Layer)) %>% 
  mutate(Layer = ifelse(Shared == 2, "Shared",Layer))
Species_modules_I$module<- as.numeric(Species_modules_I$module)

#add proportions of species in each module to the data frame
Species_modules2_I<-Species_modules_I %>% 
  mutate(Prop = case_when(
    module == 1 ~ as.numeric(N_species_mod[1,3]),
    module == 2 ~ as.numeric(N_species_mod[2,3]),
    module == 3 ~ as.numeric(N_species_mod[3,3]),
    module == 4 ~ as.numeric(N_species_mod[4,3]),
    module == 5 ~ as.numeric(N_species_mod[5,3]),
    module == 6~ as.numeric(N_species_mod[6,3]),
    module == 7 ~ as.numeric(N_species_mod[7,3]),
    module == 8 ~ as.numeric(N_species_mod[8,3]),  
    module == 9 ~ as.numeric(N_species_mod[9,3]),
    module == 10~ as.numeric(N_species_mod[10,3]),
    module == 11 ~ as.numeric(N_species_mod[11,3]),
    module == 12 ~ as.numeric(N_species_mod[12,3]),
    module == 13 ~ as.numeric(N_species_mod[13,3]),
    module == 14 ~ as.numeric(N_species_mod[14,3]),
    module == 15 ~ as.numeric(N_species_mod[15,3]),
    module == 16 ~ as.numeric(N_species_mod[16,3]), 
    module == 17 ~ as.numeric(N_species_mod[17,3]),
    module == 18~ as.numeric(N_species_mod[18,3]), 
  ))



#Change the order of x axis according to module size
Species_modules3_I<-Species_modules2_I %>% 
  select(module,layer_id,Prop) %>% 
  arrange(desc(Prop))

mod_10<- setNames(data.frame(t(c(10,1,0.01052632))),colnames(Species_modules3_I))
Species_modules4_I<- rbind(Species_modules3_I, mod_10)
Species_modules4_I$module <- factor(Species_modules4_I$module, level = c("1","5","4","3","7","2","8","6","11","9","10"))

#Plot-
I<-ggplot(Species_modules4_I,aes(x = module, y = layer_id, fill = Prop))+
  geom_point(size=9, shape = 22)+
  xlim(c("1","5","3","7","2","8","6","11","9","10","12","13","14","15","16","17","18"))+
  scico::scale_fill_scico(limits = c(0, 0.31),palette = "lajolla")+
  scale_y_continuous(limits = c(0.8,2.2),breaks=c(1, 2))+
  labs(x='Module ID', y='Layer')+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line (colour = "black",linetype = "solid", size = 1 ),
        axis.text = element_text(size=18, color='black'),
        axis.title = element_text(size=25, color='black'), 
        title = element_text(size = 10, color = 'black'))+
  coord_fixed(ratio = 5)
I

###### 2. Comparison of the observed number of modules to values generated by shuffling the network with a null model--

#Call the function "Shuffle" to generated shuffled networks and calculate modularity for each simulation
Sim_Results<-Shuffle(Npol,Ndisp)

shapiro.test(Sim_Results$Sim_Results$M) #testing for normality

t.test(Sim_Results$Sim_Results$M, mu = I_Num_Modules) #two tailed - t test

#Plot histogram (Appendix S1)--
Mod.null_I<- ggplot(Sim_Results$Sim_Results)+geom_bar(aes(x=M), fill = "orange2")+theme_bw() +
  geom_vline(xintercept =I_Num_Modules,linetype = "dashed")+
  theme(axis.line = element_line (colour = "black",linetype = "solid", size = 1 ),
        axis.title = element_text(size= 25),
        axis.text= element_text(size= 18, colour = "black"))
Mod.null_I


# Plot I Multilayer networks ----                         
#~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot using igraph a quapartide network for both biprtite layers
# polls -> plants1 -- plants2 -> destrib

# Prepare the groups:
# plants: 1-24
# pollinators: 25-91
# seed distributes: 92-94

pols <- Nodes %>% filter(node_id > 24 & node_id < 92)  %>% 
  add_column(colour="orange") %>%
  mutate(label = node_id, node_id = paste("1_", node_id, sep = ""))
plnt1 <- Nodes %>% filter(node_id <= 24) %>% 
  add_column(colour="palegreen") %>%
  mutate(label = node_id, node_id = paste("1_", node_id, sep = ""))
plnt2 <- Nodes %>% filter(node_id <= 24) %>% 
  add_column(colour="palegreen") %>%
  mutate(label = node_id, node_id = paste("2_", node_id, sep = ""))
dist <- Nodes %>% filter(node_id >= 92)  %>% 
  add_column(colour="plum1")%>%
  mutate(label = node_id, node_id = paste("2_", node_id, sep = ""))

pols$x <- 0
plnt1$x <- 0.5
plnt2$x <- 1
dist$x <- 1.5

pols$y <- (nrow(pols):1)*20/nrow(pols) 
plnt1$y <- (nrow(plnt1):1)*20/nrow(plnt1)
plnt2$y <- (nrow(plnt2):1)*20/nrow(plnt2)
dist$y <- (nrow(dist):1)*5/nrow(dist) + 7.5

flat_nodes <- rbind(plnt1, plnt2, pols, dist)

flat_edges <- Edges %>% 
  mutate(colour=case_when(layer_from==layer_to ~ 'black',
                          layer_from!=layer_to ~ 'steelblue1')) %>%
  mutate(from = paste(layer_from, "_", node_from, sep = "")) %>%
  mutate(to = paste(layer_to, "_", node_to, sep = "")) %>% 
  select(from, to, weight, colour)

net <- igraph::graph_from_data_frame(flat_edges, directed = FALSE, vertices = flat_nodes)
pdf("./Modularity/i_mln.pdf")
plot.igraph(net,  axes = FALSE, vertex.frame.color = "gray20",
            vertex.label = NA,# V(net)$label, #vertex.label.color = "gray20",
            vertex.size = 5, #vertex.size2 = 30,
            vertex.color = flat_nodes$colour,
            edge.width = flat_edges$weight*30,
            edge.color = flat_edges$colour,
            #edge.curved = T,
            main = "Invaded Network")#,
#            margin = c(-0.3,-0.15,-0.2,-0.15))
dev.off()
