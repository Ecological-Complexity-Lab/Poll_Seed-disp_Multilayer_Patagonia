#In this file we calculated the structural role of species in each treatment

library(infomapecology)
library(dplyr)
library(vegan)
source('/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Extra_Functions/Extrafunctions.R')

#Upload matrix 
NI_disp<- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_disp.csv", sep=",",row.names=1)
NI_pol <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/NI_pol.csv", sep=",",row.names=1)

I_disp<- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_disp.csv", sep=",",row.names=1)
I_pol <- read.csv("/Users/agustin/Documents/GitHub/Multilayer_Ecology-letters/Data/I_pol.csv", sep=";",row.names=1)

# set treatment (we run the model one time per treatment)
Tre='NI'

if (Tre=="NI")
{
  Ndisp=NI_disp
  Npol=NI_pol
}

if (Tre=="I")
{
  Ndisp=I_disp
  Npol=I_pol
}

#Call function to change the network' edges to proportions
PropNetPol=ToProp(Npol) 
PropNetDisp=ToProp(Ndisp)

#Call function to create the network edge list
files=CreateFiles(PropNetDisp,PropNetPol) 
Edges=files$Edges.info
Nodes=files$Nodes.info
Layers=files$Layers.info

Nodes<-cbind(Nodes, Type = rep(c("Plant","Pol","Disp"),times =c(ncol(Npol), nrow(Npol),nrow(Ndisp)))) 
                             
intra_edges<- Edges %>% #filter intra_edges
  filter(layer_from == layer_to) %>% 
  select(layer = layer_from, node_from, node_to, weight)
inter_edges <- filter(Edges, layer_from != layer_to)# filter inter_edges

#Create functions to calculate z and c score to assign a role for each species (see manuscript for more details)
z_score_multi <- function(i, modules) {
  m_i <- modules$module[modules$state_node == i]
  k_m_i <- modules$k_m[modules$state_node == i]
  avg_k <- k_m_avg$k_m_avg[k_m_avg$module == m_i]
  sd_k <- k_m_avg$k_m_sd[k_m_avg$module == m_i]
  return ((k_m_i - avg_k)/sd_k)
}

c_score_multi <- function(i, modules, edges) {
  c_sum_i <- 0
  c_sum_i <- c()
  if (i %in% edges$sn_from) {
    i_edges <- filter(edges, sn_from == i)
    for (t in unique(i_edges$m_to)) {
      k_i_t <- sum(i_edges$m_to == t)
      k_i <- as.numeric(modules$k_total[modules$state_node == i])
      c_sum_i <- c(c_sum_i, (k_i_t / k_i)**2)
    }
  }
  if (i %in% edges$sn_to) {
    i_edges <- filter(edges, sn_to == i)
    for (t in unique(i_edges$m_from)) {
      k_i_t <- sum(i_edges$m_from == t)
      k_i <- as.numeric(modules$k_total[modules$state_node == i])
      c_sum_i <- c(c_sum_i, (k_i_t / k_i)**2)
    }
  }
  return(c_score_i <- as.numeric(1-sum(c_sum_i)))
}

#Calculate a module list 
Multi_object <- create_multilayer_object(extended = Edges,
                                     nodes = Nodes, layers = Layers)#create a multilayer object using the edge list generated before
Salida_mod <- run_infomap_multilayer(Multi_object, relax = F, silent = T, trials = 100,seed=1111, 
                                     temporal_network = F,flow_model='undirected') #output containing list of modules with the species integrating each module and L.
module_list <- Salida_mod$modules 
module_list  <- module_list %>% 
  mutate(state_node = c(1:length(module_list$node_id)), 
         type = Nodes$Type[match(module_list$node_id, Nodes$node_id)]) %>% 
  dplyr::select(state_node, node_id, layer_id, module, type) %>%
  arrange(module)


# add unique IDs for state nodes
expanded_list <- Edges %>% 
  left_join(module_list, by = c('node_from' = 'node_id', 'layer_from' = 'layer_id')) %>% #merge a new dataframe searching on edges when node from == node id
  left_join(module_list, by = c('node_to' = 'node_id', 'layer_to' = 'layer_id')) %>%   #same with node_to
  dplyr::select(sn_from = state_node.x, sn_to = state_node.y, node_from, node_to, 
                layer_from, layer_to, m_from = module.x,
                m_to = module.y, weight) %>% 
  filter(is.na(sn_from) == F & is.na(sn_to) == F) # removes interlayer "edges" with zero weight
intra_only <- filter(expanded_list, layer_from == layer_to)

# total degree of state nodes *not including inter edges*
main_graph <- graph.data.frame(intra_only, directed = F, vertices = module_list)
module_list %<>% mutate(k_total = igraph::degree(main_graph)) # degree

# degree within modules
k_m <- data.frame(matrix(nrow = 0, ncol = 2))
names(k_m) <- c('state_node', 'k_m')
for (m in unique(module_list$module)) {
  m_state_nodes <- filter(module_list, module == m) 
# edges where both nodes are in m 
  m_edges <- filter(intra_only, m_from == m & m_to == m) 
  m_graph <- graph.data.frame(m_edges, directed = F, vertices = m_state_nodes)
  m_degrees <- data.frame(state_node = m_state_nodes$state_node, k_m = igraph::degree(m_graph)) # degree
  k_m <- rbind(k_m, m_degrees)
}
module_list %<>% inner_join(k_m, by = "state_node") 

# average degree in each module
k_m_avg <- data.frame(matrix(nrow = 0, ncol = 3))
for (m in unique(module_list$module)) {
  k_m_avg <- rbind(k_m_avg,
                   c(m, mean(module_list$k_m[module_list$module == m]),
                     sd(module_list$k_m[module_list$module == m])))  
}
names(k_m_avg) <- c('module', 'k_m_avg', 'k_m_sd')

#we calculate the Z and C scores of each species with the function
module_list %<>% arrange(state_node)
Z <- lapply(module_list$state_node, module_list, FUN = z_score_multi) #chequear esto!
module_list %<>% mutate(z_score = unlist(Z))
C <- lapply(sort(module_list$state_node), module_list, intra_only, FUN = c_score_multi)
module_list %<>% mutate(c_score = unlist(C)) 

super_list <- mutate(module_list, 
                     role = case_when(module_list$z_score <= 2.5 & module_list$c_score <= 0.62 ~ 'peripheral',
                                      module_list$z_score <= 2.5 & module_list$c_score > 0.62 ~ 'connector',
                                      module_list$z_score > 2.5 & module_list$c_score <= 0.62 ~ 'module hub',
                                      module_list$z_score > 2.5 & module_list$c_score > 0.62 ~ 'network hub',
                                      is.nan(module_list$z_score) == T & module_list$c_score <= 0.62 ~ 'peripheral',
                                      is.nan(module_list$z_score) == T & module_list$c_score > 0.62 ~ 'connector',
                                      is.na(module_list$z_score) == T & module_list$c_score <= 0.62 ~ 'peripheral',
                                      is.na(module_list$z_score) ==T & module_list$c_score > 0.62 ~ 'connector'))
#NI_sprole_list <- write.csv(super_list, 'NI_super_list.csv') 


#Plot indicating the distribution of species according to their structural role (Appendix S3)--
ggplot(super_list, aes(x=c_score, y=z_score, color= type)) +
  geom_point(size = 1.8) +geom_hline(yintercept = 2.5)+geom_vline(xintercept = 0.62)+
  scale_color_manual(values= c("blue","darkgreen","red"))+
  xlim(c(0,1)) + ylim(c(-1.65,4))+
  xlab("Among-module connectivity (c)")+ ylab('Within-module degree (z)')+
  theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = "white", size = 4, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "bottom",
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank (),
    panel.grid.minor = element_blank()
  )



