
setwd("SET WORKING DIRECTORY")

## Load needed libraries

library(readr)
library(dplyr)
library(igraph)
library(ggplot2)
#spatial analysis
library(geosphere)
library(sf)


### FUNCTIONS

########################################
## Add attributes from files to vertices

incorporate_att <- function(g,att_df){  # Incorporate new node attribute to network
  
  num <- length(V(g))
  
  for(pos in 1:num){
    idx <- which(att_df$site == V(g)[pos]$name)
    print(pos)
    g <- set_vertex_attr(g, "lat", pos, att_df$lat[idx])
    g <- set_vertex_attr(g, "long", pos, att_df$long[idx])
  }
  return(g)
}

#######################################
## Get the giant component of a network

giant_component <- function(g) {
  components <- components(g)
  which_max <- which.max(components$csize)
  induced_subgraph(g, which(components$membership == which_max))
}


###########################################
## Perform a basic network characterisation

metrics <- function(g){
  metrics <- list(Num_nodes =length(V(g)), Density=graph.density(g), CC=transitivity(g, type="average"), APL=average.path.length(g), DegCentr=centr_degree(g)$centralization, BetCentr=centr_betw(g, directed = FALSE)$centralization)
  return(metrics)
}





###############################################################################
## Spatial analysis to compare the geographical dispersion of two sets of nodes 
## (used to determine whether discarded nodes might be relevant)

check_nodes <- function(file_A,file_B){ # file_A and file_B are names of files in .CSV format
  
  # Import both sets of sites
  set_A <- read_delim(file_A, delim = ";", escape_double = FALSE, trim_ws = TRUE) # Sites used as nodes
  set_A <- set_A[,-1] # used nodes

  set_B <- read_delim(file_B, delim = ";", escape_double = FALSE, trim_ws = TRUE) # Sites to be checked
  set_B$lat <- as.numeric(set_B$lat)
  set_B$long <- as.numeric(set_B$long)
  
  set_A <- setdiff(set_A, set_B)
  
  
  # Convert to sf object type for statistical analysis
  pts_sf_A <- st_as_sf(set_A, coords = c("lat", "long"), crs = 4326)
  pts_sf_B <- st_as_sf(set_B, coords = c("lat", "long"), crs = 4326)
  

  # Plotting

  plot(pts_sf_A)
  plot(pts_sf_B)

  col_A <- adjustcolor("blue", alpha.f = 0.3)  # alpha: 0 = fully transparent, 1 = fully opaque
  col_B <- adjustcolor("red",  alpha.f = 0.3)

  plot(st_geometry(pts_sf_A), col = col_A, pch = 16, cex = 0.8)
  plot(st_geometry(pts_sf_B), col = col_B, pch = 16, cex = 0.8, add = TRUE)
  legend("topright", legend = c("Used", "Not used"), col = c(col_A, col_B), pch = 16)

}




#########################################################################################
## Network robustness test 
## implementation of the procedure for archaeological networks proposed by Peeples and Brughmans)

net_robustness_test <- function(g,metric=c("deg","bet"), type=c("node","edge"), num_iterations=500, jump, end){
  
  metric <- match.arg(metric)
  type <- match.arg(type)
  jump <- as.integer(length(V(g))*jump)
  end <- length(V(g))*end
  
  result <- data.frame()
  set.seed(1)
  
  if(metric=="deg"){g<- set_vertex_attr(g, "cent", V(g), degree(g))}   # Centrality measure: Degree
  else{g <- set_vertex_attr(g, "cent", V(g), betweenness(g))}   # Centrality measure: Betweenness
  
  for (i in 1:num_iterations){
    for (n in seq(jump,end,jump)){
      if(type=="node"){  # Removing nodes
        out <- sample(1:length(V(g)),n)
        test_net <- delete_vertices(g, out)}
      else{  # Removing edges
        out <- sample(1:length(E(g)),n)
        test_net <- delete_edges(g, out)}
      
      if(metric=="deg") {x <- degree(test_net)}
      else{x <- betweenness(test_net)}
      y <- V(test_net)$cent
      correlation <- cor(x, y, method = 'spearman')
      
      if(type=="node"){max=length(V(g))}
      else{max=length(E(g))}
      
      result <- rbind(result, data.frame(erased=as.character(round(n*100/max)), metric=correlation))    }  
  }
  result <- result[order(result$erased),]
  
  # Plotting test results
  if(type=="node"){label="Removed nodes (%)"}
  else{label="Removed edges (%)"}
    ggplot(data=result, aes(x=erased, y=metric)) + geom_boxplot() + 
    xlab(label) + 
    ylab("Correlation") + 
    ylim(c(0.0, 1.0))
}




################################################################
## Network model 1 (Random) + small-world-ness index calculation 

random_nm <- function(test_net, all_net, num_rep){
  
  path=c()  # Average path length
  cc <- c()  # Average Clustering Coefficient 
  cent_bet <- c() # Betweenness centralization
  cent_deg <- c() # Degree centralization

  set.seed(1)


  # Extract list of non-overlapping edges (so, those ones over land)
  extra_edges <- as_data_frame(difference(all_net, test_net), what = "edges")
  edge_ids <- get.edge.ids(all_net, vp = t(as.matrix(extra_edges)))

  # Obtain geodistance matrix. This will be used to check distances covered by edges
  xy <- data.frame(long=as.numeric(V(all_net)$long),lat=as.numeric(V(all_net)$lat))
  dist_mat <- distm(xy)

  # Get the maximum distance in the set of edges to be rewired
  edges_df <- as.data.frame(ends(all_net,edge_ids,names=FALSE))

  max_dist = 0
  for (i in 1:length(edge_ids)){
    dist_i <- dist_mat[edges_df$V1[i],edges_df$V2[i]]
    max_dist <- max(c(max_dist,dist_i))
  }

  # Random networks' generation

  for (num in 1:num_rep){ 
    print(num)
    working_net <- all_net
  
    for (edge_pos in edge_ids){
      pos <- ends(working_net,edge_pos,names=FALSE)[1,1] # First node at the edge to be rewired. This one is kept.
      neigh <- sample.int(length(V(working_net)),1) # Second node at the edge to be rewired. This one is changed.
      #while((are_adjacent(working_net,pos,neigh)) || (pos==neigh)) {neigh <- sample.int(length(V(working_net)),1)} # WITHOUT geographical restriction
      while((are_adjacent(working_net,pos,neigh)) || (pos==neigh) || (dist_mat[pos,neigh]>max_dist)) {neigh <- sample.int(length(V(working_net)),1)} # WITH geographical restriction
    
      # Once the alternative node is identified, we proceed to rewiring 
      edge_list <- as_data_frame(working_net, what = "edges")
      edge_list[edge_pos, "to"] <- V(working_net)$name[neigh]
      working_net <- graph_from_data_frame(edge_list, directed = is_directed(working_net), vertices = as_data_frame(working_net, what = "vertices"))
    }
    working_net <- giant_component(working_net)
    cent_bet<- c(cent_bet,centr_betw(working_net,directed=FALSE)$centralization)
    cent_deg<- c(cent_deg, centr_degree(working_net)$centralization)
    path <- c(path,mean_distance(working_net, unconnected=TRUE))
    cc <- c(cc,transitivity(working_net, type="average"))
  }

  # Calculation of the small-world-ness index a la (Humphries & Gurney, 2008)

  norm_path <- mean_distance(all_net, directed=FALSE, unconnected=TRUE) / mean(path)
  norm_cc <- transitivity(all_net, type="average")/mean(cc)

  # Generate results
  results <- list("cc"=cc, "apl"=path, "cent_bet"=cent_bet, "cent_deg"=cent_deg, "S"=norm_cc/norm_path)
  return(results)
}


######################################
## Network model 2 (Nearest neighbors)

nearneigh_nm <- function(test_net, all_net, num_rep){
  
  path=c()  # Average path length
  cc <- c()  # Average Clustering Coefficient 
  cent_bet <- c() # Betweenness centralization
  cent_deg <- c() # Degree centralization
  
  set.seed(1)
  


  # Extract list of non-overlapping edges (so, those one over land)
  extra_edges <- as_data_frame(difference(all_net, test_net), what = "edges")
  edge_ids <- get.edge.ids(all_net, vp = t(as.matrix(extra_edges))) # Rewiring ONLY extra edges


  # Obtain geodistance matrix. This will be used to identify links to closest neighbors
  xy <- data.frame(long=as.numeric(V(all_net)$long),lat=as.numeric(V(all_net)$lat))

  dist_mat <- distm(xy)
  diag(dist_mat) <- 9999999

  for (num in 1:num_rep){
    print(num)
  
    working_net <- all_net
  
    for (edge_pos in edge_ids){
      pos = sample.int(2,1)
      pos <- ends(working_net,edge_pos,names=FALSE)[1,pos] # First node at the edge to be rewired. This one is kept.
      pot_neigh <- dist_mat[pos,]
      neigh <- which.min(pot_neigh) # Second node is re-defined by choosing the closest one.
      while((are_adjacent(working_net,pos,neigh))) {
        pot_neigh[neigh] <- 9999999
        neigh <- which.min(pot_neigh)
      }
      # Once the alternative node is identified, we proceed to rewiring 
      edge_list <- as_data_frame(working_net, what = "edges")
      edge_list[edge_pos, "to"] <- V(working_net)$name[neigh]
      working_net <- graph_from_data_frame(edge_list, directed = is_directed(working_net), vertices = as_data_frame(working_net, what = "vertices"))
    
  }
    working_net <- giant_component(working_net)
    cent_bet<- c(cent_bet,centr_betw(working_net,directed=FALSE)$centralization)
    cent_deg<- c(cent_deg, centr_degree(working_net)$centralization)
    path <- c(path,mean_distance(working_net, unconnected=TRUE))
    cc <- c(cc,transitivity(working_net, type="average"))
  
  }
  # Generate results
  results <- list("cc"=cc, "apl"=path, "cent_bet"=cent_bet, "cent_deg"=cent_deg)
  return(results)
}


########################################################################
## Plot comparison between empirical and synthetic results (requires ggplot2)

plot_compare <- function(empirical, random, nn){

  # Build a long-format dataframe for the two null models
  df <- data.frame(values = c(random, nn),
    model = c(rep("Random", length(random)),
              rep("Nearest-neighbor", length(nn)))
  )

ggplot(data=df, aes(x = model, y = values, fill = model)) +
  
  # Boxplot for null models
  geom_boxplot(alpha = 0.7, width = 0.4, outlier.shape = 21) +
  
  # Empirical value as a horizontal reference line
  geom_hline(yintercept = empirical, 
             linetype = "dashed", 
             color = "red", 
             linewidth = 0.8) +
  
  
  # Labels and theme
  labs(x = NULL, y = NULL, title = NULL) +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 13))
}





#######################################
#######################################
############     MAIN      ############
#######################################

### NETWORK CONSTRUCTION

## Load edgelist

LC3_water_edges <- read_delim("Data/LC3_water_v6.csv", delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
LC3_all_edges <- read_delim("Data/LC3_all_v6.csv", delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

LC4_water_edges <- read_delim("Data/LC4_water_v6.csv", delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
LC4_all_edges <- read_delim("Data/LC4_all_v6.csv", delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

LC5_water_edges <- read_delim("Data/LC5_water_v6.csv", delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
LC5_all_edges <- read_delim("Data/LC5_all_v6.csv", delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


## Create graph object

LC3_rivers<-graph_from_data_frame(LC3_water_edges,directed=FALSE)
LC3_all<-graph_from_data_frame(LC3_all_edges,directed=FALSE)

LC4_rivers<-graph_from_data_frame(LC4_water_edges,directed=FALSE)
LC4_all<-graph_from_data_frame(LC4_all_edges,directed=FALSE)

LC5_rivers<-graph_from_data_frame(LC5_water_edges,directed=FALSE)
LC5_all<-graph_from_data_frame(LC5_all_edges,directed=FALSE)


## Add geographic position as a vertex attribute

LC_geo <- read_delim("Data/LC_lat-long_v4.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE) # Read hierarchies

LC3_all<- incorporate_att(LC3_all,LC_geo)
LC4_all<- incorporate_att(LC4_all,LC_geo)
LC5_all<- incorporate_att(LC5_all,LC_geo)


## Generate .GML files for visualisation (as inputs for Gephi, for instance)

write_graph(LC3_all,"LC3_all_v6.gml",format="gml")
write_graph(LC4_all,"LC4_all_v6.gml",format="gml")
write_graph(LC5_all,"LC5_all_v6.gml",format="gml")



### TESTS


## Spatial analysis to assess the potential impact of including discarded nodes

check_nodes("Data/LC_lat-long_v4.csv","Data/discarded.csv") # Sites disregarded because of low chronological precision
check_nodes("Data/LC_lat-long_v4.csv","Data/Non-Uruk.csv")  # Sites disregarded because of non-Uruk assemblage

## Network reliability/robustness tests (results presented as plots)


# LC3
net_robustness_test(LC3_all,metric="deg", type="node", num_iterations=500, jump = 0.01, end=0.09) #Removing Nodes/Impact on Degree
net_robustness_test(LC3_all,metric="bet", type="node", num_iterations=500, jump = 0.01, end=0.09) #Removing Nodes/Impact on Betweenness
net_robustness_test(LC3_all,metric="deg", type="edge", num_iterations=500, jump = 0.01, end=0.09) #Removing Edges/Impact on Degree
net_robustness_test(LC3_all,metric="bet", type="edge", num_iterations=500, jump = 0.01, end=0.09) #Removing Edges/Impact on Betweenness

#LC4
net_robustness_test(LC4_all,metric="deg", type="node", num_iterations=500, jump = 0.01, end=0.09) #Removing Nodes/Impact on Degree
net_robustness_test(LC4_all,metric="bet", type="node", num_iterations=500, jump = 0.01, end=0.09) #Removing Nodes/Impact on Betweenness
net_robustness_test(LC4_all,metric="deg", type="edge", num_iterations=500, jump = 0.01, end=0.09) #Removing Edges/Impact on Degree
net_robustness_test(LC4_all,metric="bet", type="edge", num_iterations=500, jump = 0.01, end=0.09) #Removing Edges/Impact on Betweenness

#LC5
net_robustness_test(LC5_all,metric="deg", type="node", num_iterations=500, jump = 0.01, end=0.09) #Removing Nodes/Impact on Degree
net_robustness_test(LC5_all,metric="bet", type="node", num_iterations=500, jump = 0.01, end=0.09) #Removing Nodes/Impact on Betweenness
net_robustness_test(LC5_all,metric="deg", type="edge", num_iterations=500, jump = 0.01, end=0.09) #Removing Edges/Impact on Degree
net_robustness_test(LC5_all,metric="bet", type="edge", num_iterations=500, jump = 0.01, end=0.09) #Removing Edges/Impact on Betweenness




### NETWORK ANALYSIS

## Basic network characterisation (only reported metrics)

metrics(LC3_all)
metrics(LC4_all)
metrics(LC5_all)



## Computational experiments 


# Generate values corresponding to network model 1 (Random)

random_LC3 <- random_nm(LC3_rivers,LC3_all, 100)
random_LC4 <- random_nm(LC4_rivers,LC4_all, 100)
random_LC5 <- random_nm(LC5_rivers,LC5_all, 100)


# Generate values corresponding to network model 2 (nearest neighbors)

nn_LC3 <- nearneigh_nm(LC3_rivers,LC3_all, 100)
nn_LC4 <- nearneigh_nm(LC4_rivers,LC4_all, 100)
nn_LC5 <- nearneigh_nm(LC5_rivers,LC5_all, 100) 



# Plotting the experiments' results (to feed Figs. 3 and S6)

# LC3
plot_compare(transitivity(LC3_all, type = "average"),random_LC3$cc, nn_LC3$cc) # Clustering coefficient
plot_compare(mean_distance(LC3_all, directed=FALSE, unconnected=TRUE),random_LC3$apl, nn_LC3$apl) # Average Path Length
plot_compare(centr_betw(LC3_all,directed=FALSE)$centralization,random_LC3$cent_bet, nn_LC3$cent_bet) # Betweenness Centralization
plot_compare(centr_degree(LC3_all)$centralization,random_LC3$cent_deg, nn_LC3$cent_deg) # Degree Centralization

# LC4
plot_compare(transitivity(LC4_all, type = "average"),random_LC4$cc, nn_LC4$cc) # Clustering coefficient
plot_compare(mean_distance(LC4_all, directed=FALSE, unconnected=TRUE),random_LC4$apl, nn_LC4$apl) # Average Path Length
plot_compare(centr_betw(LC4_all,directed=FALSE)$centralization,random_LC4$cent_bet, nn_LC4$cent_bet) # Betweenness Centralization
plot_compare(centr_degree(LC4_all)$centralization,random_LC4$cent_deg, nn_LC4$cent_deg) # Degree Centralization

# LC5
plot_compare(transitivity(LC5_all, type = "average"),random_LC5$cc, nn_LC5$cc) # Clustering coefficient
plot_compare(mean_distance(LC5_all, directed=FALSE, unconnected=TRUE),random_LC5$apl, nn_LC5$apl) # Average Path Length
plot_compare(centr_betw(LC5_all,directed=FALSE)$centralization,random_LC5$cent_bet, nn_LC5$cent_bet) # Betweenness Centralization
plot_compare(centr_degree(LC5_all)$centralization,random_LC5$cent_deg, nn_LC5$cent_deg) # Degree Centralization

