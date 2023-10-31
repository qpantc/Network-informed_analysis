library(dplyr)
library(tidygraph)
library(igraph)
library(raster)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(factoextra)
library(ggpmisc) #加载ggpmisc包
library(mgcv)


# 0 Affiliated function ------------------------------------------------------------
## Correlation to edges 
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  return(data.frame( row = rownames(cormat)[row(cormat)[ut]],
                     column = rownames(cormat)[col(cormat)[ut]],
                     cor  =(cormat)[ut],
                     p = pmat[ut]))
}

data_to_edges <- function(data,min.p=0.05,min.r=0.2){
  temp_corr_data <- Hmisc::rcorr(as.matrix(data,type="pearson"))
  temp_corr_data <-flattenCorrMatrix(temp_corr_data$r, temp_corr_data$P)
  temp_corr_data <-temp_corr_data[(temp_corr_data$p<=min.p),]
  colnames(temp_corr_data) <- c("to","from","estimate","p.value")
  all_edges <- temp_corr_data[,c("to","from","estimate")]
  all_edges$estimate <- abs(all_edges$estimate)
  all_edges <-all_edges[(all_edges$estimate>=min.r),]
  colnames(all_edges)<- c("from","to","weight")
  return(all_edges)
}
#··································· -------------------------------------------

## Number --> trait combination

# Functions ---------------------------------------------------------------
PowerSetsBinary <- function(number,all_node_list=c("LA","Lthick","Llen","Lwid","SLA","LWC","LFM","L15N","LDMC","LCC.m","LNC.m","LPC.m","LCNR","SSD","SCDen","SCDia","SCElen","WDL","RD","SRL","SdM","SdNum","SdL","SdGE","DispL","Height","SD")){
  ## 根据一个10进制的数和一个性状名称列表，算出该数字转化为二进制对应的性状组合。
  one_network_nodes = c()
  for (j in c(0:(length(all_node_list)-1))) {
    # print(bitwShiftR(number,j))
    if (bitwShiftR(number,j) %% 2 == 1){
      # print(all_node_list[j+1])
      one_network_nodes <- c(one_network_nodes,all_node_list[j+1])
      # print(one_network_nodes)
    }
  }
  if (length(one_network_nodes)>1){
    return(one_network_nodes)
  } else{
    return(c())
  }
}

## trait combination --> Number 
trait_list_to_number <- function(trait_list,full_node_list=c("LA","Lthick","Llen","Lwid","SLA","LWC","LFM","L15N","LDMC","LCC.m","LNC.m","LPC.m","LCNR","SSD","SCDen","SCDia","SCElen","WDL","RD","SRL","SdM","SdNum","SdL","SdGE","DispL","Height","SD")){
  x <- full_node_list %in% trait_list
  x[x==TRUE] = 1
  x[x==FALSE] = 0
  trait_number_str <- paste0(rev(x),collapse = '')
  reduced_network_number <- strtoi(trait_number_str, base = 2)
  return(reduced_network_number)
}

## get_different_traits in two given trait list
get_different_traits <- function(all_traits,concise_traits)
{
  different_traits <- c()
  for (trait in all_traits) {
    if (!trait %in% concise_traits){
      different_traits <- c(different_traits,trait)
    }
  }
  return(different_traits)
}
#··································· -------------------------------------------

## Nomalized_tif_save
Nomalized_tif_save <- function(raw_matrix,vairable_name){
  raw_matrix[raw_matrix == 99999999] <- NA
  normailze_matrix <-  BBmisc::normalize(raw_matrix, method = "standardize")
  save_file = paste0('./_Results/Data/Metrics_tif/standardize_',vairable_name,'.tif')
  writeRaster(raster(normailze_matrix),save_file,overwrite=TRUE)
}

calculate_distance_and_save <- function(raw_matrix,max_metrics_index,vairable_name,method='standardize'){
  raw_matrix[raw_matrix == 99999999] <- NA
  if (method == 'raw'){
    normailze_matrix <-  raw_matrix
    save_file = paste0('./_Results/Data/Distance_tif/raw_',vairable_name,'_distance.tif')
  }
  else if (method == 'standardize'){
    normailze_matrix <-  BBmisc::normalize(raw_matrix, method = "standardize")#scale_in_1
    save_file = paste0('./_Results/Data/Distance_tif/standardize_',vairable_name,'_distance.tif')
  }else if (method == '0-1'){
    normailze_matrix <-  scale_in_1(raw_matrix,b = 0)#
    save_file = paste0('./_Results/Data/Distance_tif/0-1_',vairable_name,'_distance.tif')
  }else if (method == 'scale'){
    normailze_matrix <-  BBmisc::normalize(raw_matrix, method = "scale")#
    save_file = paste0('./_Results/Data/Distance_tif/scale_',vairable_name,'_distance.tif')
  }
  max_metrics_value <- normailze_matrix[max_metrics_index]
  distance_matricx <- (normailze_matrix-max_metrics_value)^2
  writeRaster(raster(distance_matricx),save_file,overwrite=TRUE)
  return(max_metrics_value^2)
}
#··································· -------------------------------------------
## Module box: calculate_nested_module 
calculate_nested_box <- function(df){
  #对于一个包含模块和性状百分比的数据
  # df <- Module.composition.walktrap
  module_df <- df %>% group_by(module) %>% summarise(variance = sum(variance))
  # module_df <- df %>% group_by(module) %>% summarise(variance = mean(variance))
  module_df$Nodes <- '.module'
  module_box_df <- calculate_box(module_df)
  for (each_module in module_box_df$module) {
    module_trait_df <- subset(df,module==each_module)[,c('module','variance','Nodes')]
    module_bottom <- module_box_df[which(module_box_df$module==each_module),]$low_line[[1]]
    module_up <- module_box_df[which(module_box_df$module==each_module),]$top_line[[1]]
    one_module_trait_df <- calculate_box(module_trait_df,inter=0,low_est=module_bottom,top_est=module_up)
    module_box_df <- rbind(module_box_df,one_module_trait_df)
  }
  return(module_box_df)
}
## calculate one module 
calculate_box <- function(df,inter=0.04,low_est=0,top_est=1){
  # 要求传入的df只有两列，一列名称，一列variance,必须要有一列是 variance
  df$percentage <- df$variance/sum(df$variance)
  df <- df[order(df$percentage),]#  
  # df <- df[order(df$Nodes),]
  # 对于每个box计算其上顶，下底，中间值，和高度
  space = ifelse(nrow(df)>1,inter/(nrow(df)-1),0)
  interval = c()
  for (i in seq(0,nrow(df))) {
    interval = c(interval,sum(df$percentage[0:i]))
  }
  box_line = c()
  for (i in seq(1,length(interval)-1)) {
    box_length <- interval[i+1] - interval[i]
    low_line = ifelse(i== 1,low_est,interval[i]*(top_est-low_est) + space/2 +low_est)
    top_line = ifelse(i+1==length(interval),top_est,interval[i+1]*(top_est-low_est)- space/2+low_est)
    mid_line = (low_line+top_line)/2
    box_line = append(box_line,c(low_line,mid_line,top_line)) 
    # print(interval)
    # print(c(low_line,mid_line,top_line,top_line-low_line))
  } 
  box_line_df = as.data.frame(matrix(box_line,ncol = 3,byrow = TRUE))
  colnames(box_line_df) = c('low_line','mid_line','top_line')
  # print(box_line_df)
  df <- cbind(df,box_line_df)
  # print(df)
  return(df)
}

stratified_network_sample <- function(distance.tif,size.tif,number.tif,
                                      interval=0.01,sample_size=100,
                                      strat_size=4,end_size=26){
  result_df <- data.frame()
  for(size in c(strat_size:end_size)){
    # size filter
    size_filter <- size.tif
    size_filter[which(size_filter!=size)] = NA
    size_filter[which(size_filter==size)] = 1 
    distance_with_size_filter <- size_filter * distance.tif
    percentage_min = 0
    for (percentage_max in seq(0,1,interval)[-1]) {
      distance_min = quantile(distance_with_size_filter,percentage_min,na.rm=TRUE)
      distance_max = quantile(distance_with_size_filter,percentage_max,na.rm=TRUE)
      print(paste0(size,'  ',percentage_max,': ',distance_min,'-',distance_max))
      distance_filter <- size_filter
      distance_filter[which(distance_with_size_filter>=distance_min & distance_with_size_filter<distance_max)] = TRUE
      distance_filter[which(distance_with_size_filter<distance_min | distance_with_size_filter>=distance_max)] = FALSE
      
      if(sum(distance_filter,na.rm=TRUE) >= sample_size){
        network_number_df <- data.frame(number = sample(number.tif[which(distance_filter==TRUE)],sample_size))
      }
      else{
        network_number_df <- data.frame(number = number.tif[which(distance_filter==TRUE)])
      }
      if (nrow(network_number_df)>0){
        network_number_df$size = size
        network_number_df$value = (distance_min+distance_max)/2
        network_number_df$percentage = percentage_max
        result_df <- rbind(result_df,network_number_df)
      }
      percentage_min = percentage_max
    }
  }
  return(result_df)
}


get_network_cost <- function(number){
  trait_list <- PowerSetsBinary(number)
  cost_df <- subset(Vertices,Nodes %in% trait_list)
  network_cost <- sum(cost_df$cost)
  return(network_cost)
}

get_network_cost2 <- function(number){
  trait_list <- PowerSetsBinary(number,all_node_list = new_node_list)
  cost_df <- subset(Vertices_applicaiton,node_new %in% trait_list)
  network_cost <- sum(cost_df$cost)
  return(network_cost)
}

#··································· -------------------------------------------

# 1 Basical -- Network calculating -------------------------
#  data + trait combination --> network 
## 1.1 Calculate everything --------
## number of network 
data_to_network <- function(data,number = NA,trait_combination=NA,r_min = 0.2) {
  # 根据给定的性状组合trait_combination，从数据data中筛选对应的性状列进行网络构建
  if (!is.na(number)){
    node_list = PowerSetsBinary(number)
  } else if (!is.na(trait_combination)){
    node_list <- strsplit(trait_combination,'_')[[1]]
  } else{
    print("no network number, neither trait combination")
    return(NA)
  }
  
  data_subset <- data %>% dplyr::select(node_list)

  ### node mean sd calculate
  node_value_sd <-   data_subset %>% 
    summarise_all(list(mean=mean,variance = sd),na.rm = TRUE) %>% as.data.frame()  %>% 
    pivot_longer(everything(),names_sep ="_",names_to = c("Nodes", ".value"))
  
  ### calculate the edges
  temp_corr <- Hmisc::rcorr(as.matrix(data_subset,type="pearson"))
  temp_corr<-flattenCorrMatrix(temp_corr$r, temp_corr$P)
  corr_list <-temp_corr[(temp_corr$p<=0.05),]
  colnames(corr_list) <- c("to","from","estimate","p.value")
  edges<-corr_list[,c("to","from","estimate")]
  edges$estimate <- abs(edges$estimate)
  edges <-edges[(edges$estimate>=r_min),]
  colnames(edges)<- c("from","to","weight")
  
  ### build the network
  one_network <- tbl_graph(nodes = node_value_sd, edges = edges,directed = FALSE) 
  Net <- graph_from_data_frame(d=edges,vertices=node_value_sd, directed = FALSE)
  
  output <- list (node_value_sd,
                  Net,
                  one_network)
  names(output) <- c("One.Node.metrics",
                     "Net",
                     "tidynet")
  return (output)
}


## data,trait_combination --> network metrics; node parameters; module composition
Calculate_network <- function(data,trait_combination,r_min = 0.2) {
  # 根据给定的性状组合trait_combination，从数据data中筛选对应的性状列进行网络构建
  # print(trait_combination)
  node_list <- strsplit(trait_combination,'_')[[1]]
  data_subset <- data %>% dplyr::select(node_list)
  
  ### node mean sd calculate
  node_value_sd <- data.frame(Nodes = node_list)
  
  ### calculate the edges
  temp_corr <- Hmisc::rcorr(as.matrix(data_subset,type="pearson"))
  temp_corr<-flattenCorrMatrix(temp_corr$r, temp_corr$P)
  corr_list <-temp_corr[(temp_corr$p<=0.05),]
  colnames(corr_list) <- c("to","from","estimate","p.value")
  edges<-corr_list[,c("to","from","estimate")]
  edges$estimate <- abs(edges$estimate)
  edges <-edges[(edges$estimate>=r_min),]
  colnames(edges)<- c("from","to","weight")

  ### build the network
  one_network <- tbl_graph(nodes = node_value_sd, edges = edges,directed = FALSE) 
  Net <- graph_from_data_frame(d=edges,vertices=node_value_sd, directed = FALSE)
  
  #--------------------------------------------------         [network metrics] 
  #1. connectance/ weighted density
  density <- 2*nrow(edges)/(length(node_list)*(length(node_list)-1))
  density_weighted <- 2*sum(edges$weight)/(length(node_list)*(length(node_list)-1))
  #2. diameter # 所有点之间最短路径的最大值,由于权重是相关系数，距离不好衡量，不能使用权重
  diameter_ <- diameter(Net, weights=NA) # 与graph_diameter() 同
  #3. average path length (mean of shortest path between all nodes)
  # 距离应该默认为1
  average_path_length <- mean_distance(Net, weights=NA)
  
  #1. Modulatity
  cluster_walktrap_ <- NA  #默认有加权 cluster_edge_betweenness(Net,weights=abs(E(Net)$weight))
  cluster_walktrap_ <- modularity(cluster_walktrap(Net))
  modularity_optimal <- modularity(cluster_optimal(Net))
  
  #2. average clustering coefficient (AC)
  node_ac <- transitivity(Net, type="local")
  avera_cc <- sum(node_ac,na.rm = TRUE)/length(node_ac)
  
  one_combination_network <- c(length(node_list),round(density,4),round(density_weighted,12),round(diameter_,4),
                               round(average_path_length,12),round(cluster_walktrap_,12),
                               round(modularity_optimal,12),round(avera_cc,12),trait_combination)
  names(one_combination_network) <- c('size','density','density_weighted','diameter',
                                      'average_path_length','cluster_walktrap_',
                                      'modularity_optimal','avera_cc','name')
  One_network_metrics <- as.data.frame(t(one_combination_network))
  
  #--------------------------------------------------            [node metrics] 
  #1. degree ###
  node_value_sd$node_degree <- degree(Net) # 与centrality_degree同
  #1. weighted degree ###
  node_value_sd$weighted_degree <- strength(Net) # 与centrality_degree同，加weighted
  #2. closenness: reciprocal of the mean shortest path between a focal node trait and all other nodes 距离应该默认为1
  node_value_sd$node_closenness <- closeness(Net)
  #3. Node betweenness (B): the number of shortest paths between the focal node trait 距离应该默认为1
  node_value_sd$node_betweenness <- betweenness(Net) 
  #4.1 Clustering coefficient (CC) #邻居节点之间的连接程度(有连接的数量/总可能数量), 不要weighted，另一种含义
  node_value_sd$local_transitivity <- transitivity(Net, type="local") # 与clustering coefficient有不同
  #5 hub value
  node_value_sd$hub_score <- hub_score(Net)$vector #默认weighted 与centrality_hub相同，
  #4.2 Clustering ability (CA)
  node_value_sd$clustering_ability = (node_value_sd$node_degree-1)*(node_value_sd$node_degree-0)*node_value_sd$local_transitivity/2
  
  node_value_sd$variance = 1#node_value_sd$variance*node_value_sd$variance
  
  #--------------------------------------------------[Module composition betweenness]
  Module.composition.louvain0 <- node_value_sd
  Module.composition.louvain0$module <- cluster_louvain(Net)$membership
  Module.composition.louvain <- calculate_nested_box(Module.composition.louvain0)
  
  #-----------------------------------------------  [Module composition optimal]
  trait_module_optimal <- node_value_sd
  trait_module_optimal$module <- cluster_optimal(Net)$membership
  Module.composition.optimal <- calculate_nested_box(trait_module_optimal)
  
  output <- list (One_network_metrics,
                  node_value_sd,
                  Module.composition.louvain,
                  Module.composition.optimal,
                  Net)
  names(output) <- c("One_network_metrics",
                     "One.Node.metrics",
                     'Module.composition.louvain',
                     'Module.composition.optimal',
                     "Net")
  return (output)
}
## 1.2 Calculate network metrics only and but also export igraph  -------
Calculate_network_metrics <- function(data,trait_combination,r_min = 0.2) {
  # 根据给定的性状组合trait_combination，从数据data中筛选对应的性状列进行网络构建
  # print(trait_combination)
  node_list <- strsplit(trait_combination,'_')[[1]]
  data_subset <- data %>% dplyr::select(node_list)
  
  ### node mean sd calculate
  node_value_sd <- data_subset %>% 
    summarise(across(everything(), list(mean=mean,variance = sd)),na.rm = TRUE) %>% as.data.frame()  %>% 
    dplyr::select(-'na.rm') %>% pivot_longer(everything(),names_sep ="_",names_to = c("Nodes", ".value"))
  
  ### calculate the edges
  temp_corr <- Hmisc::rcorr(as.matrix(data_subset,type="pearson"))
  temp_corr<-flattenCorrMatrix(temp_corr$r, temp_corr$P)
  corr_list <-temp_corr[(temp_corr$p<=0.05),]
  colnames(corr_list) <- c("to","from","estimate","p.value")
  edges<-corr_list[,c("to","from","estimate")]
  edges$estimate <- abs(edges$estimate)
  edges <-edges[(edges$estimate>=r_min),]
  
  colnames(edges)<- c("from","to","weight")
  # print(node_value_sd)
  ### build the network
  Net <- graph_from_data_frame(d=edges,vertices=node_value_sd, directed = FALSE)
  
  #--------------------------------------------------     [network metrics] 
  #1. connectance/ weighted density
  density <- 2*nrow(edges)/(length(node_list)*(length(node_list)-1))
  density_weighted <- 2*sum(edges$weight)/(length(node_list)*(length(node_list)-1))
  #2. diameter # 所有点之间最短路径的最大值,由于权重是相关系数，距离不好衡量，不能使用权重
  diameter_ <- diameter(Net, weights=NA) # 与graph_diameter() 同
  #3. average path length (mean of shortest path between all nodes)
  # 距离应该默认为1
  average_path_length <- mean_distance(Net, weights=NA)
  
  #1. Modulatity
  modularity_optimal <- modularity(cluster_optimal(Net))
  modularity_louvain <- modularity(cluster_louvain(Net))
  modularity_walktrap <- modularity(cluster_walktrap(Net))
  
  
  #2. average clustering coefficient (AC)
  node_ac <- transitivity(Net, type="local")
  avera_cc <- sum(node_ac,na.rm = TRUE)/length(node_ac)
  
  one_combination_network <- c(length(node_list),round(density,4),round(density_weighted,6),round(diameter_,4),
                               round(average_path_length,6),round(avera_cc,6),
                               round(modularity_optimal,6),round(modularity_louvain,6),round(modularity_walktrap,6),
                               trait_combination)
  names(one_combination_network) <- c('size','density','density_weighted','diameter',
                                      'average_path_length','avera_cc',
                                      'modularity_optimal','modularity_louvain','modularity_walktrap',
                                      'name')
  One_network_metrics <- as.data.frame(t(one_combination_network))

  output <- list (One_network_metrics,
                  Net)
  names(output) <- c("One_network_metrics",
                     "Net")
  return (output)
}
Calculate_network_metrics_inner <- function(data,trait_combination,r_min = 0.2,network_number) {
  # 根据给定的性状组合trait_combination，从数据data中筛选对应的性状列进行网络构建
  # print(trait_combination)
  node_list <- strsplit(trait_combination,'_')[[1]]
  data_subset <- data %>% dplyr::select(node_list)
  
  ### node mean sd calculate
  node_value_sd <- as.data.frame(node_list)
  
  ### calculate the edges
  temp_corr <- Hmisc::rcorr(as.matrix(data_subset,type="pearson"))
  temp_corr<-flattenCorrMatrix(temp_corr$r, temp_corr$P)
  corr_list <-temp_corr[(temp_corr$p<=0.05),]
  colnames(corr_list) <- c("to","from","estimate","p.value")
  edges<-corr_list[,c("to","from","estimate")]
  edges$estimate <- abs(edges$estimate)
  edges <-edges[(edges$estimate>=r_min),]
  
  colnames(edges)<- c("from","to","weight")

  Net <- graph_from_data_frame(d=edges,vertices=node_value_sd, directed = FALSE)
  
  #--------------------------------------------------     [network metrics] 
  density <- 2*nrow(edges)/(length(node_list)*(length(node_list)-1))
  density_weighted <- 2*sum(edges$weight)/(length(node_list)*(length(node_list)-1))    
  #1. Modulatity
  modularity_optimal <- modularity(cluster_optimal(Net))
  modularity_louvain <- modularity(cluster_louvain(Net))
  modularity_walktrap <- modularity(cluster_walktrap(Net))

  one_combination_network <- c(length(node_list),round(density,4),round(density_weighted,6),
                               round(modularity_optimal,6),round(modularity_louvain,6),round(modularity_walktrap,6),
                               network_number)
  names(one_combination_network) <- c('size','density','density_weighted',
                                      'modularity_optimal','modularity_louvain','modularity_walktrap',
                                      'number')
  One_network_metrics <- as.data.frame(t(one_combination_network))
  return(One_network_metrics)
}

## 1.3 Calculate network metrics only -------
network_metrics_with_combination <- function(data,trait_combination,r_min = 0.2){
  # 根据给定的性状组合trait_combination，从数据data中筛选对应的性状列进行网络构建
  # print(trait_combination)
  node_list <- strsplit(trait_combination,'_')[[1]]
  data_subset <- data %>% dplyr::select(node_list)
  
  ### node mean sd calculate
  node_value_sd <- data_subset %>% 
    summarise(across(everything(), list(mean=mean,variance = sd)),na.rm = TRUE) %>% as.data.frame()  %>% 
    dplyr::select(-'na.rm') %>% pivot_longer(everything(),names_sep ="_",names_to = c("Nodes", ".value"))
  
  ### calculate the edges
  temp_corr <- Hmisc::rcorr(as.matrix(data_subset,type="pearson"))
  temp_corr<-flattenCorrMatrix(temp_corr$r, temp_corr$P)
  corr_list <-temp_corr[(temp_corr$p<=0.05),]
  colnames(corr_list) <- c("to","from","estimate","p.value")
  edges<-corr_list[,c("to","from","estimate")]
  edges$estimate <- abs(edges$estimate)
  all_edges <-all_edges[(all_edges$estimate>=r_min),]
  
  colnames(edges)<- c("from","to","weight")
  ### build the network
  Net <- graph_from_data_frame(d=edges,vertices=node_value_sd, directed = FALSE)
  
  ### network metrics
  #1. connectance/ weighted density
  density <- 2*nrow(edges)/(length(node_list)*(length(node_list)-1))
  density_weighted <- 2*sum(edges$weight)/(length(node_list)*(length(node_list)-1))
  #2. diameter # 所有点之间最短路径的最大值,由于权重是相关系数，距离不好衡量，不能使用权重
  diameter_ <- diameter(Net, weights=NA) # 与graph_diameter() 同
  #3. average path length (mean of shortest path between all nodes)
  # 距离应该默认为1
  average_path_length <- mean_distance(Net, weights=NA)
  #4. evenness
  H <- vegan::diversity(abs(edges [,"weight"]))# with absolute values
  S <- dim (edges)[1]
  evenness <- H/log(S)
  
  # Complexity ###
  #1. Modulatity
  modularity_edge_betweenness <- modularity(cluster_edge_betweenness(Net))
  modularity_optimal <- modularity(cluster_optimal(Net))
  
  #2. average clustering coefficient (AC)
  node_ac <- transitivity(Net, type="local")
  avera_cc <- sum(node_ac,na.rm = TRUE)/length(node_ac)
  
  one_combination_network <- c(length(node_list),round(density,4),round(density_weighted,4),round(diameter_,4),
                               round(average_path_length,4),round(evenness,4),round(modularity_edge_betweenness,4),
                               round(modularity_optimal,4),round(avera_cc,4),trait_combination)
  names(one_combination_network) <- c('size','density','density_weighted','diameter',
                                      'average_path_length','evenness','modularity_edge_betweenness',
                                      'modularity_optimal','avera_cc','name')
  one_result <- as.data.frame(t(one_combination_network))
  return(one_result)
}
#··································· -------------------------------------------

# 2 Gradient -- Network calculating -------------------------
## 2.1 Combination list: Calculate everything -------------------------
## network metrics; nodes parameters; module composition
Network_metrics_along_combination <- function(df,combination_list){
  
  Net.metrics.final <- data.frame()
  Node.metrics.final <-data.frame()
  Module.composition.walktrap.final <-data.frame()
  Module.composition.optimal.final <-data.frame()
  Net.list <- list()
  
  for (combination in combination_list) {
    node_list <- strsplit(combination,'_')[[1]]
    one_subset_data <- data %>% dplyr::select(node_list)
    # calculate the full network
    Calculated_network <- Calculate_network(data=one_subset_data %>% dplyr::select(node_list),
                                            trait_combination=combination)
    
    One_network_metrics = Calculated_network$One_network_metrics
    One.Node.metrics = Calculated_network$One.Node.metrics
    Module.composition.walktrap = Calculated_network$Module.composition.walktrap
    Module.composition.optimal = Calculated_network$Module.composition.optimal
    
    One_network_metrics['size'] = length(node_list)
    One.Node.metrics['size'] =  length(node_list)
    Module.composition.walktrap['size'] =  length(node_list)
    Module.composition.optimal['size'] =  length(node_list)
    
    One_network_metrics['network'] = combination
    One.Node.metrics['network'] =  combination
    Module.composition.walktrap['network'] = combination
    Module.composition.optimal['network'] =  combination
    
    Net.metrics.final <- rbind(Net.metrics.final,One_network_metrics)
    Node.metrics.final <-rbind(Node.metrics.final,One.Node.metrics)
    Module.composition.walktrap.final <-rbind(Module.composition.walktrap.final,Module.composition.walktrap)
    Module.composition.optimal.final <-rbind(Module.composition.optimal.final,Module.composition.optimal)
    
    Net.list[[length(Net.list)+1]] <- Calculated_network$Net
  }
  output <- list (Net.metrics.final,
                  Node.metrics.final,
                  Module.composition.walktrap.final,
                  Module.composition.optimal.final,
                  Net.list)
  names(output) <- c('Net.metrics.final',
                     'Node.metrics.final',
                     'Module.composition.walktrap.final',
                     'Module.composition.optimal.final',
                     'Net.list')
  return (output)
}
Network_metrics_along_number <- function(df,number_list){
  
  Net.metrics.final <- data.frame()
  Node.metrics.final <-data.frame()
  Module.composition.walktrap.final <-data.frame()
  Module.composition.optimal.final <-data.frame()
  Net.list <- list()
  
  for (number in number_list) {
    node_list = PowerSetsBinary(number)
    combination = paste0(node_list,collapse = '_')
    one_subset_data <- df %>% dplyr::select(node_list)
    # calculate the full network
    Calculated_network <- Calculate_network(data=one_subset_data,
                                            trait_combination=combination)
    
    One_network_metrics = Calculated_network$One_network_metrics
    One.Node.metrics = Calculated_network$One.Node.metrics
    Module.composition.walktrap = Calculated_network$Module.composition.walktrap
    Module.composition.optimal = Calculated_network$Module.composition.optimal
    
    One_network_metrics['size'] = length(node_list)
    One.Node.metrics['size'] =  length(node_list)
    Module.composition.walktrap['size'] =  length(node_list)
    Module.composition.optimal['size'] =  length(node_list)
    
    One_network_metrics['network'] = combination
    One.Node.metrics['network'] =  combination
    Module.composition.walktrap['network'] = combination
    Module.composition.optimal['network'] =  combination
    
    Net.metrics.final <- rbind(Net.metrics.final,One_network_metrics)
    Node.metrics.final <-rbind(Node.metrics.final,One.Node.metrics)
    Module.composition.walktrap.final <-rbind(Module.composition.walktrap.final,Module.composition.walktrap)
    Module.composition.optimal.final <-rbind(Module.composition.optimal.final,Module.composition.optimal)
    
    Net.list[[length(Net.list)+1]] <- Calculated_network$Net
  }
  output <- list (Net.metrics.final,
                  Node.metrics.final,
                  Module.composition.walktrap.final,
                  Module.composition.optimal.final,
                  Net.list)
  names(output) <- c('Net.metrics.final',
                     'Node.metrics.final',
                     'Module.composition.walktrap.final',
                     'Module.composition.optimal.final',
                     'Net.list')
  return (output)
}

## 2.2 Gradient: Calculate network metrics and export igraph -------------------------
# network metrics; igraph objective
Network_metrics_gradient <- function(df,trait_combination,gradient='gradient',
                                     min_species = 20,overlab=0,interval = 1,
                                     interval.n=0){
  node_list <- strsplit(trait_combination,'_')[[1]]
  Net.metrics.final <- data.frame()
  Net.list <- list()
  
  gradient.min = min(df[gradient]);gradient.max = max(df[gradient])
  if (interval.n != 0){ # 按照数量划分
    gradient_interval <- (gradient.max-gradient.min)/(interval.n+1)
  }
  else{#按照间隔为1的划分
    gradient_interval <- interval
  }
  gradient_list <- seq(gradient.min,gradient.max,by = gradient_interval)
  for (gradient.start in gradient_list) {
    gradient.end = gradient.start + gradient_interval*(overlab+1)
    if (gradient.end <= max(gradient_list)){
      one_subset_data <- df[which(df[gradient]>=gradient.start & df[gradient]<gradient.end),]
      species_number <- sum(!duplicated(one_subset_data$Species))
      
      if (species_number>=min_species){
        print(c(trait_combination,gradient.start,gradient.end,nrow(one_subset_data)))
        # calculate the full network
        Calculated_network <- Calculate_network_metrics(data=one_subset_data %>% dplyr::select(node_list),
                                                        trait_combination=trait_combination)
        
        One_network_metrics = Calculated_network$One_network_metrics
        One_network_metrics[gradient] = gradient.start
        One_network_metrics['species_number'] = species_number
        Net.metrics.final <- rbind(Net.metrics.final,One_network_metrics)
        
        Net.list[[gradient.start]] <- Calculated_network$Net
      }
    }
  }
  output <- list (Net.metrics.final,
                  Net.list)
  names(output) <- c('Net.metrics.final',
                     'Net.list')
  return (output)
}
#··································· -------------------------------------------

# 3 Node parameters --------------------------------------------------------
## 3.1 Network sensitivity to each node in the network ---------------------
## with data trait combination
Net_metrics_sensitivity_to_node <- function(trait_df,trait_combination){
  Node_sensitive.final <-data.frame()
  full_network_metrics.final <- data.frame()
  
  node_list <- strsplit(trait_combination,'_')[[1]]
  # calculate the full network
  full_network_metrics <- 
    network_metrics_with_combination(data=trait_df %>% dplyr::select(node_list),
                                     trait_combination=trait_combination)
  full_network_metrics.final <- rbind(full_network_metrics.final,full_network_metrics)
  
  for (omit_trait in node_list) {
    omint_node_list <- node_list[node_list != omit_trait]
    omint_network_metrics <- 
      network_metrics_with_combination(data=trait_df %>% dplyr::select(omint_node_list),
                                       trait_combination=paste0(omint_node_list,collapse = '_'))
    delta_df <- rbind(as.numeric(full_network_metrics[,2:9]),as.numeric(omint_network_metrics[,2:9] ))
    colnames(delta_df) <- colnames(full_network_metrics)[2:9]
    # calculte total distance with density weighted_density diameter modularity
    Dis_density_2 = (delta_df[2,1]-delta_df[1,1])^2
    Dis_density_weighted_2 = (delta_df[2,2]-delta_df[1,2])^2
    Dis_diameter_2 = (delta_df[2,3]-delta_df[1,3])^2
    Dis_modularity_optimal_2 = (delta_df[2,7]-delta_df[1,7])^2
    distance = (0.5*(Dis_density_2+Dis_density_weighted_2) + 
                  Dis_modularity_optimal_2)^0.5
    
    
    sensitivity_df = as.data.frame(t(delta_df[2,]-delta_df[1,]))
    sensitivity_df$distance =distance
    sensitivity_df$name = omit_trait
    sensitivity_df$type = 'value'
    
    
    sensitivity_precentage = as.data.frame(t((delta_df[2,]-delta_df[1,])/delta_df[1,]))
    sensitivity_precentage$distance = (distance)/
      ((delta_df[1,1]^2 + delta_df[1,2]^2)/2  + delta_df[1,7]^2)^0.5
    sensitivity_precentage$name = omit_trait
    sensitivity_precentage$type = 'precentage'
    
    Node_sensitive.final <- rbind(Node_sensitive.final,sensitivity_df,sensitivity_precentage)
  }
  full_network_metrics.final$network_name <- trait_combination
  Node_sensitive.final$network_name <- trait_combination
  # print(Node_sensitive.final)
  output <- list (full_network_metrics.final,
                  Node_sensitive.final)
  names(output) <- c("full_network_metrics",
                     "Node_sensitive")
  return (output)
}
## with edges/combination number
Net_metrics_sensitivity_to_nodes <- function(all_edges,network_number){
  all_node_list <- c("LA","Lthick","Llen","Lwid","SLA","LWC","LFM","L15N","LDMC",
                     "LCC.m","LNC.m","LPC.m","LCNR","SSD","SCDen","SCDia","SCElen",
                     "WDL","RD","SRL","SdM","SdNum","SdL","SdGE","DispL","Height","SD")
  node_list <- PowerSetsBinary(network_number,all_node_list)
  
  sensitivity_df <- data.frame()
  
  # define function and calculate the trait list with network_number
  print(node_list)
  network_metrics_with_combination_inner <- function(node_list,all_edges = all_edges){
    ## 根据给定的性状列表，所有边的关系，计算该网络的density，density_weighted，modularity_optimal
    # 根据数字确定需要计算的网络
    if (length(node_list)>1){
      # 根据给定的性状组合trait_combination，从数据data中筛选对应的性状列进行网络构建
      ### node mean sd calculate
      node_value_sd <- as.data.frame(node_list)
      colnames(node_value_sd) <- 'Nodes'
      ### select all edges
      edges <- subset(all_edges,from %in% node_list & to %in% node_list) 
      ### build the network
      Net <- graph_from_data_frame(d=edges,vertices=node_value_sd, directed = FALSE)
      
      ### network metrics
      #1. connectance/ weighted density
      len_node_list <- length(node_list)
      density <- 2*nrow(edges)/(len_node_list*(len_node_list-1))
      density_weighted <- 2*sum(edges$weight)/(len_node_list*(len_node_list-1))
      
      # Complexity ###
      #1. Modulatity
      modularity_optimal <- modularity(cluster_optimal(Net))
      return(c(round(density,4),round(density_weighted,4),round(modularity_optimal,4)))
    }
  }
  if (length(node_list)>3){
    # calculate the full network with full node list
    full_network_metrics <- network_metrics_with_combination_inner(node_list,all_edges = all_edges)
    
    full_network_density <- full_network_metrics[1]
    full_network_weighted_density <- full_network_metrics[2]
    full_network_modularity <- full_network_metrics[3]
    # print(c(full_network_density,full_network_weighted_density,full_network_modularity))
    
    value_of_full_network <- (0.5*full_network_metrics[1]^2 + 0.5*full_network_metrics[2]^2 +
                                full_network_metrics[3]^2)^0.5
    print(value_of_full_network)
    if (!is.na(value_of_full_network) & value_of_full_network>0){
      for (trait in node_list) {
        omint_node_list <- node_list[node_list != trait]
        reduced_network_metrics <- network_metrics_with_combination_inner(omint_node_list,all_edges = all_edges)
        
        reduced_network_density <- reduced_network_metrics[1]
        reduced_network_weighted_density <- reduced_network_metrics[2]
        reduced_network_modularity <- reduced_network_metrics[3]
        # print(c(reduced_network_density,reduced_network_weighted_density,reduced_network_modularity))
        
        distance_to_full_network <- (0.5*(full_network_density - reduced_network_density)^2 +
                                       0.5*(full_network_weighted_density - reduced_network_weighted_density)^2 +
                                       (full_network_modularity - reduced_network_modularity)^2)^0.5
        
        print(c(distance_to_full_network,value_of_full_network))
        one_result <- c(distance_to_full_network/value_of_full_network,
                        (full_network_density - reduced_network_density)/full_network_density,
                        (full_network_weighted_density - reduced_network_weighted_density)/full_network_weighted_density,
                        (full_network_modularity - reduced_network_modularity)/full_network_modularity,
                        full_network_density,full_network_weighted_density,full_network_modularity,length(node_list),network_number)
        one_result <- as.data.frame(t(one_result))
        one_result$Trait <- trait
        one_result$network_name <- paste0(node_list,collapse = '_')
        # print(one_result)
        sensitivity_df <- rbind(sensitivity_df,one_result)
      }
    }
    else{
      for (trait in node_list) {
        one_result <- c(0,0,0,0,
                        full_network_density,full_network_weighted_density,full_network_modularity,length(node_list),network_number)
        one_result <- as.data.frame(t(one_result))
        one_result$Trait <- trait
        one_result$network_name <- paste0(node_list,collapse = '_')
        # print(one_result)
        sensitivity_df <- rbind(sensitivity_df,one_result)
      }
    }
  }
  colnames(sensitivity_df) <- c('sensitivity_all','sensitivity_density','sensitivity_connectance','sensitivity_modularity',
                                'density','connectance','modularity','size','network_number','Trait','network_name')
  
  return(sensitivity_df)
}
Net_metrics_sensitivity_to_nodes2 <- function(all_edges,network_number){
  all_node_list <-  c("Aarea","L15N", "LA","LCC.m","LCNR", "LD","LDM","LDMC","LL",
                      "LNC.m","LPC.m","LSC", "LSLA","Lthick", "LWid","Max.H","Rarea",
                      "RD","Rdia","RMI","RNC", "RSRL","RTD", "SD","SdM","SHV","SSD")
  node_list <- PowerSetsBinary(network_number,all_node_list)
  
  sensitivity_df <- data.frame()
  
  # define function and calculate the trait list with network_number
  print(node_list)
  network_metrics_with_combination_inner <- function(node_list,all_edges = all_edges){
    ## 根据给定的性状列表，所有边的关系，计算该网络的density，density_weighted，modularity_optimal
    # 根据数字确定需要计算的网络
    if (length(node_list)>1){
      # 根据给定的性状组合trait_combination，从数据data中筛选对应的性状列进行网络构建
      ### node mean sd calculate
      node_value_sd <- as.data.frame(node_list)
      colnames(node_value_sd) <- 'Nodes'
      ### select all edges
      edges <- subset(all_edges,from %in% node_list & to %in% node_list) 
      ### build the network
      Net <- graph_from_data_frame(d=edges,vertices=node_value_sd, directed = FALSE)
      
      ### network metrics
      #1. connectance/ weighted density
      len_node_list <- length(node_list)
      density <- 2*nrow(edges)/(len_node_list*(len_node_list-1))
      density_weighted <- 2*sum(edges$weight)/(len_node_list*(len_node_list-1))
      
      # Complexity ###
      #1. Modulatity
      modularity_optimal <- modularity(cluster_optimal(Net))
      return(c(round(density,4),round(density_weighted,4),round(modularity_optimal,4)))
    }
  }
  if (length(node_list)>3){
    # calculate the full network with full node list
    full_network_metrics <- network_metrics_with_combination_inner(node_list,all_edges = all_edges)
    
    full_network_density <- full_network_metrics[1]
    full_network_weighted_density <- full_network_metrics[2]
    full_network_modularity <- full_network_metrics[3]
    # print(c(full_network_density,full_network_weighted_density,full_network_modularity))
    
    value_of_full_network <- (0.5*full_network_metrics[1]^2 + 0.5*full_network_metrics[2]^2 +
                                full_network_metrics[3]^2)^0.5
    print(value_of_full_network)
    if (!is.na(value_of_full_network) & value_of_full_network>0){
      for (trait in node_list) {
        omint_node_list <- node_list[node_list != trait]
        reduced_network_metrics <- network_metrics_with_combination_inner(omint_node_list,all_edges = all_edges)
        
        reduced_network_density <- reduced_network_metrics[1]
        reduced_network_weighted_density <- reduced_network_metrics[2]
        reduced_network_modularity <- reduced_network_metrics[3]
        # print(c(reduced_network_density,reduced_network_weighted_density,reduced_network_modularity))
        
        distance_to_full_network <- (0.5*(full_network_density - reduced_network_density)^2 +
                                       0.5*(full_network_weighted_density - reduced_network_weighted_density)^2 +
                                       (full_network_modularity - reduced_network_modularity)^2)^0.5
        
        print(c(distance_to_full_network,value_of_full_network))
        one_result <- c(distance_to_full_network/value_of_full_network,
                        (full_network_density - reduced_network_density)/full_network_density,
                        (full_network_weighted_density - reduced_network_weighted_density)/full_network_weighted_density,
                        (full_network_modularity - reduced_network_modularity)/full_network_modularity,
                        full_network_density,full_network_weighted_density,full_network_modularity,length(node_list),network_number)
        one_result <- as.data.frame(t(one_result))
        one_result$Trait <- trait
        one_result$network_name <- paste0(node_list,collapse = '_')
        # print(one_result)
        sensitivity_df <- rbind(sensitivity_df,one_result)
      }
    }
    else{
      for (trait in node_list) {
        one_result <- c(0,0,0,0,
                        full_network_density,full_network_weighted_density,full_network_modularity,length(node_list),network_number)
        one_result <- as.data.frame(t(one_result))
        one_result$Trait <- trait
        one_result$network_name <- paste0(node_list,collapse = '_')
        # print(one_result)
        sensitivity_df <- rbind(sensitivity_df,one_result)
      }
    }
  }
  colnames(sensitivity_df) <- c('sensitivity_all','sensitivity_density','sensitivity_connectance','sensitivity_modularity',
                                'density','connectance','modularity','size','network_number','Trait','network_name')
  
  return(sensitivity_df)
}

## 3.2 Network sensitivity to each node in the network ---------------------
Node_sensitivity_gradient <- function(df,node_list,gradient='gradient',
                                      min_species = 20,overlab=0,interval = 0,
                                      interval.n=0){
    
    full_network_metrics.final <- data.frame()
    Node_sensitive.final <-data.frame()
    
    gradient.min = min(df[gradient]);gradient.max = max(df[gradient])
    if (interval == 0 & interval.n != 0){
      gradient_interval <- (gradient.max-gradient.min)/(interval.n+1)
    }
    else if (interval != 0 & interval.n == 0) {
      gradient_interval <- interval
    }
    else{
      print('set value to interval or interval.n, only one!')
    }
    gradient_list <- seq(gradient.min,gradient.max,by = gradient_interval)
    for (gradient.start in gradient_list) {
      gradient.end = gradient.start + gradient_interval*(overlab+1)
      if (gradient.end <= max(gradient_list)){
        
        one_subset_data <- df[which(df[gradient]>=gradient.start & df[gradient]<gradient.end),]
        species_number <- sum(!duplicated(one_subset_data$Species))
        
        if (species_number>=min_species){
          print(c(gradient.start,gradient.end,nrow(one_subset_data)))
          # calculate the full network
          full_network_metrics <- 
            network_metrics_with_combination(data=one_subset_data %>% dplyr::select(node_list),
                                             trait_combination=paste0(node_list,collapse = '_'))
          full_network_metrics[gradient] = gradient.start
          full_network_metrics.final <- rbind(full_network_metrics.final,full_network_metrics)
          
          for (omit_trait in node_list) {
            omint_node_list <- node_list[node_list != omit_trait]
            omint_network_metrics <- 
              network_metrics_with_combination(data=one_subset_data %>% dplyr::select(omint_node_list),
                                               trait_combination=paste0(omint_node_list,collapse = '_'))
            det_df <- rbind(as.numeric(full_network_metrics[,2:9]),as.numeric(omint_network_metrics[,2:9] ))
            colnames(det_df) <- colnames(full_network_metrics)[2:9]
            sensitivity_df = as.data.frame(t((det_df[1,]-det_df[2,])/det_df[1,]))
            sensitivity_df$name = omit_trait
            sensitivity_df[gradient] = gradient.start
            Node_sensitive.final <- rbind(Node_sensitive.final,sensitivity_df)
            print(sensitivity_df)
          }
        }
      }
    }
    output <- list (full_network_metrics.final,
                    Node_sensitive.final)
    names(output) <- c("full_network_metrics.final",
                       "Node_sensitive.final")
    return (output)
  }

#··································· -------------------------------------------
# 4 Data analysis -----------------------------------------------------------
## 4.1 Value range: Calculate network parameters range in different size 4-27 -----------------------------------------------------------
mean_sd_on_size_from_tif <- function(size.tif,data.tif,min_size=4,variable='value'){
  size.tif[size.tif == 99999999] <- NA
  max_size = max(size.tif,na.rm = TRUE)
  result <- data.frame()
  for (size in c(min_size:max_size)) {
    
    size.tif_mask <- ifelse(size.tif==size,1,NA)
    mean_value = mean(size.tif_mask*data.tif,na.rm = TRUE)
    sd_value = sd(size.tif_mask*data.tif,na.rm = TRUE)
    one_result = c(size,mean_value,sd_value)
    result = rbind(result,t(as.data.frame(one_result)))
    print(one_result)
  }
  colnames(result) = c('size','mean','sd')
  rownames(result) = NULL
  result$variable = variable
  return(result)
}

## 4.2 RDA and plots -----------------------------------------------------------
Rda_plot <- function(trait_combination,climate_factor_list,Trait_ecoregion_merged){
  node_list <- strsplit(trait_combination,'_')[[1]]
  datatemp<-Trait_ecoregion_merged[,node_list]
  dataENV<-Trait_ecoregion_merged[,climate_factor_list]#c(29:59) c(29,30,31,37,40,41,42,44,45)
  
  p=vegan::rda(datatemp,dataENV)
  sm=summary(p,scaling=0)
  # plot(p, display=c("bp","si","sp"),type="points",scaling=3,xlim=c(-2,4),choice=c(1,3))
  bi=sm$species
  bs=sm$sites
  bp=sm$biplot
  p$tot.chi
  jpeg(filename = paste0('./_Results/Figures/RDA/',length(node_list),'rda_env_factors_',trait_combination,'.jpg'),width = 3000, height = 3000, units = "px", pointsize = 60)
  axis_length <- max(abs(min(bp[,2])),abs(max(bp[,2])),
                     abs(min(bp[,1])),abs(max(bp[,1])),
                     abs(min(bs[,2])),abs(max(bs[,2])),
                     abs(min(bs[,1])),abs(max(bs[,1])))
  plot(bs[,1],bs[,2],col="gray63",xlim=c(-1.5*axis_length,1.5*axis_length),
       ylim=c(-1.5*axis_length,1.5*axis_length),
       main=paste0(length(node_list),' - rda_env_factors'),
       xlab=paste0("RDA1 ",scales::percent(sm$cont$importance[3,1],0.01)),
       ylab=paste0("RDA2 ",scales::percent(sm$cont$importance[3,2]-sm$cont$importance[3,1],0.01)))
  
  abline(h=c(0), v=c(0), col="gray33", lty=3,lwd=1.5)
  axis(1, at=c(-10,10),lwd = 2)
  axis(2, at=c(-10,10),lab=FALSE,lwd = 2)
  
  for (i in c(1:nrow(bp))) {
    trait_name <- rownames(bp)[i]
    arrows(0,0,bp[i,1],bp[i,2],length=0.4,lwd =1.5,col="black")
    text_x <- ifelse(bp[i,1]>0,bp[i,1]+0.05*axis_length,bp[i,1]-0.05*axis_length)
    text_y <- ifelse(bp[i,2]>0,bp[i,2]+0.05*axis_length,bp[i,2]-0.05*axis_length)
    text(text_x,text_y,trait_name,col="black")
  }
  dev.off()
  
  jpeg(filename = paste0('./_Results/Figures/RDA/',length(node_list),'rda_trait_',trait_combination,'.jpg'),width = 3000, height = 3000, units = "px", pointsize = 60)

  plot(bs[,1],bs[,2],col="gray63",xlim=c(-1*axis_length,1*axis_length),
       ylim=c(-1*axis_length,1*axis_length),
       main=paste0(length(node_list),' - rda_trait '),
       xlab=paste0("RDA1 ",scales::percent(sm$cont$importance[3,1],0.01)),
       ylab=paste0("RDA2 ",scales::percent(sm$cont$importance[3,2]-sm$cont$importance[3,1],0.01)))
  
  abline(h=c(0), v=c(0), col="gray33", lty=3,lwd=1.5)
  

  # #Climate variables
  # axis(3, col = "blue",lab=FALSE,at=c(-round(axis_length,2)*2,-round(axis_length,2),round(axis_length,2),round(axis_length,2)*2),lwd = 2)
  # axis(3, col = "blue",lty = 1, lwd = 2,at=c(-round(axis_length,2),0,round(axis_length,2)),labels=c(-round(axis_length,2),0,round(axis_length,2)) ,col.axis = "blue")
  # axis(4, col = "blue",lab=FALSE,at=c(-round(axis_length,2)*2,-round(axis_length,2),round(axis_length,2),round(axis_length,2)*2),lwd = 2)
  # axis(4, col = "blue", lty = 1, lwd = 2,at=c(-round(axis_length,2),0,round(axis_length,2)),labels=c(-round(axis_length,2),0,round(axis_length,2)) ,col.axis = "blue")
  # 
  for (i in c(1:nrow(bi))) {
    trait_name <- rownames(bi)[i]
    arrows(0,0,bi[i,1],bi[i,2],length=0.4,lwd =1.5,col="blue")
    text_x <- ifelse(bi[i,1]>0,bi[i,1]+0.05*axis_length,bi[i,1]-0.05*axis_length)
    text_y <- ifelse(bi[i,2]>0,bi[i,2]+0.05*axis_length,bi[i,2]-0.05*axis_length)
    text(text_x,text_y,trait_name,col="blue")
  }
  
  dev.off()
}

Rda_plot_number <- function(network_number,climate_factor_list,Trait_ecoregion_merged){
  node_list <- PowerSetsBinary(network_number)
  trait_combination <- paste0(node_list,collapse = '_')
  datatemp<-Trait_ecoregion_merged[,node_list]
  dataENV<-Trait_ecoregion_merged[,climate_factor_list]#c(29:59) c(29,30,31,37,40,41,42,44,45)
  
  p=vegan::rda(datatemp,dataENV)
  sm=summary(p,scaling=0)
  # plot(p, display=c("bp","si","sp"),type="points",scaling=3,xlim=c(-2,4),choice=c(1,3))
  bi=sm$species
  bs=sm$sites
  bp=sm$biplot
  p$tot.chi
  jpeg(filename = paste0('./_Results/Figures/RDA/',length(node_list),'rda_env','.jpg'),width = 3000, height = 3000, units = "px", pointsize = 60)
  axis_length <- max(abs(min(bp[,2])),abs(max(bp[,2])),
                     abs(min(bp[,1])),abs(max(bp[,1])),
                     abs(min(bs[,2])),abs(max(bs[,2])),
                     abs(min(bs[,1])),abs(max(bs[,1])))
  plot(bs[,1],bs[,2],col="gray63",xlim=c(-1.5*axis_length,1.5*axis_length),
       ylim=c(-1.5*axis_length,1.5*axis_length),
       main=paste0(length(node_list),' - rda_env_factors'),
       xlab=paste0("RDA1 ",scales::percent(sm$cont$importance[3,1],0.01)),
       ylab=paste0("RDA2 ",scales::percent(sm$cont$importance[3,2]-sm$cont$importance[3,1],0.01)))
  
  abline(h=c(0), v=c(0), col="gray33", lty=3,lwd=1.5)
  axis(1, at=c(-10,10),lwd = 2)
  axis(2, at=c(-10,10),lab=FALSE,lwd = 2)
  
  for (i in c(1:nrow(bp))) {
    trait_name <- rownames(bp)[i]
    arrows(0,0,bp[i,1],bp[i,2],length=0.4,lwd =1.5,col="black")
    text_x <- ifelse(bp[i,1]>0,bp[i,1]+0.05*axis_length,bp[i,1]-0.05*axis_length)
    text_y <- ifelse(bp[i,2]>0,bp[i,2]+0.05*axis_length,bp[i,2]-0.05*axis_length)
    text(text_x,text_y,trait_name,col="black")
  }
  dev.off()
  
  jpeg(filename = paste0('./_Results/Figures/RDA/',length(node_list),'rda_trait','.jpg'),width = 3000, height = 3000, units = "px", pointsize = 60)
  
  plot(bs[,1],bs[,2],col="gray63",xlim=c(-1*axis_length,1*axis_length),
       ylim=c(-1*axis_length,1*axis_length),
       main=paste0(length(node_list),' - rda_trait '),
       xlab=paste0("RDA1 ",scales::percent(sm$cont$importance[3,1],0.01)),
       ylab=paste0("RDA2 ",scales::percent(sm$cont$importance[3,2]-sm$cont$importance[3,1],0.01)))
  
  abline(h=c(0), v=c(0), col="gray33", lty=3,lwd=1.5)

  for (i in c(1:nrow(bi))) {
    trait_name <- rownames(bi)[i]
    arrows(0,0,bi[i,1],bi[i,2],length=0.4,lwd =1.5,col="blue")
    text_x <- ifelse(bi[i,1]>0,bi[i,1]+0.05*axis_length,bi[i,1]-0.05*axis_length)
    text_y <- ifelse(bi[i,2]>0,bi[i,2]+0.05*axis_length,bi[i,2]-0.05*axis_length)
    text(text_x,text_y,trait_name,col="blue")
  }
  
  dev.off()
}

#··································· -------------------------------------------

# 5. Application of network --------------------------------------------------
network_metrics_with_edges <- function(node_list,all_edges = all_edges){
  ## 根据给定的性状列表，所有边的关系，计算该网络的density，density_weighted，modularity_optimal
  # 根据数字确定需要计算的网络
  if (length(node_list)>1){
    # 根据给定的性状组合trait_combination，从数据data中筛选对应的性状列进行网络构建
    data_subset <- data[,node_list]
    
    ### node mean sd calculate
    node_value_sd <- as.data.frame(node_list)
    colnames(node_value_sd) <- 'Nodes'
    ### select all edges
    edges <- subset(all_edges,from %in% node_list & to %in% node_list) 
    ### build the network
    Net <- graph_from_data_frame(d=edges,vertices=node_value_sd, directed = FALSE)
    
    ### network metrics
    #1. connectance/ weighted density
    len_node_list <- length(node_list)
    density <- 2*nrow(edges)/(len_node_list*(len_node_list-1))
    density_weighted <- 2*sum(edges$weight)/(len_node_list*(len_node_list-1))
    
    # Complexity ###
    #1. Modulatity
    modularity_optimal <- modularity(cluster_optimal(Net))
    return(c(round(density,4),round(density_weighted,4),round(modularity_optimal,4)))
  }
}
## 5.1 Identify redundant traits ----------------------------------------------
# provide the full network with all traits
# provide the concise network only with reduced traits

Redundant_traits <- function(all_edges,all_traits,concise_traits){
  # all_edges <- all_edges
  # all_traits <- c("LA","Lthick","Llen","Lwid","SLA","LWC","LFM",
  #                 "L15N","LDMC","LCC.m","LNC.m","LPC.m","LCNR",
  #                 "SSD","SCDen","SCDia","SCElen","WDL","RD","SRL",
  #                 "SdM","SdNum","SdL","SdGE","DispL","Height","SD")
  # concise_traits <- c("LA","Lthick",
  #                     "L15N","LDMC","LCC.m",
  #                     "SSD","SCDen","WDL","RD","SRL",
  #                     "SdM","SdNum","Height","SD")
  # 
  # Get the parameters of the network with all traits
  full_network_metrics <- network_metrics_with_edges(all_traits,all_edges = all_edges)
  full_network_density <- full_network_metrics[1]
  full_network_weighted_density <- full_network_metrics[2]
  full_network_modularity <- full_network_metrics[3]
  print(c('full_network: ',full_network_density,full_network_weighted_density,full_network_modularity))
  value_of_full_network <- (0.5*full_network_metrics[1]^2 + 0.5*full_network_metrics[2]^2 +
                              full_network_metrics[3]^2)^0.5
  
  # Get the parameters of the network with concise traits
  concise_network_metrics <- network_metrics_with_edges(concise_traits,all_edges = all_edges)
  concise_network_density <- concise_network_metrics[1]
  concise_network_weighted_density <- concise_network_metrics[2]
  concise_network_modularity <- concise_network_metrics[3]
  print(c('concise_network: ',concise_network_density,concise_network_weighted_density,concise_network_modularity))
  value_of_concise_network <- (0.5*concise_network_metrics[1]^2 + 0.5*concise_network_metrics[2]^2 +
                                 concise_network_metrics[3]^2)^0.5
  
  distance_of_concise_to_full <- (0.5*(full_network_density - concise_network_density)^2 +
                                    0.5*(full_network_weighted_density - concise_network_weighted_density)^2 +
                                    (full_network_modularity - concise_network_modularity)^2)^0.5
  ini_result <- c('-',
                  '-',
                  round(distance_of_concise_to_full/value_of_full_network,6),
                  0,
                  round(distance_of_concise_to_full,6),round(value_of_full_network,6),
                  paste0(all_traits,collapse = '_'))
  # the list redundant traits (n)
  redundant_traits <- get_different_traits(all_traits,concise_traits)
  result <- rbind(data.frame(),t(as.data.frame(ini_result)))
  
  # The example/objective must be good 
  if (!is.na(value_of_full_network) & value_of_full_network>0){ 
    un_redundant_traits = c()
    # a loop or iteration to find the most redundant (m) traits combination
    while (length(redundant_traits)>0) {
      un_redundant_trait = ''
      Distance = 1
      
      
      for (trait_i in redundant_traits) {
        
        # 比较每个(m-1)个性状组合的与全性状网络的距离，找到最小距离的组合就是最冗余的组合
        # 
        # get a relatively un-redundant trait
        sub_trait_list <- redundant_traits[redundant_traits != trait_i]
        sub_trait_list <- c(concise_traits,un_redundant_traits,trait_i)
        # relativly unredundant trait network
        
        reduced_network_metrics <- network_metrics_with_edges(sub_trait_list,all_edges = all_edges)
        reduced_network_density <- reduced_network_metrics[1]
        reduced_network_weighted_density <- reduced_network_metrics[2]
        reduced_network_modularity <- reduced_network_metrics[3]
        
        
        distance_to_full_network <- (0.5*(full_network_density - reduced_network_density)^2 +
                                       0.5*(full_network_weighted_density - reduced_network_weighted_density)^2 +
                                       (full_network_modularity - reduced_network_modularity)^2)^0.5
        
        distance_to_concise_network <- (0.5*(concise_network_density - reduced_network_density)^2 +
                                          0.5*(concise_network_weighted_density - reduced_network_weighted_density)^2 +
                                          (concise_network_modularity - reduced_network_modularity)^2)^0.5
        
        # 判断当前组合是否为距离变化最小的组合。
        if (distance_to_full_network < Distance){
          # 如果是，更新距离、冗余的性状列表
          Distance =  round(distance_to_full_network,6)
          un_redundant_trait = trait_i
          lost_information <- round(Distance/value_of_full_network,6)
          add_information <-  round(distance_to_concise_network/value_of_concise_network,6)
        }
        
        # 一次回合结束之后，生成一次结果记录 [最大冗余 = 距离变化最小]
        # 最大冗余的性状组合 | 去除冗余性状组合后的与原始网络的距离 | 
        # 本次保留的较不冗余性状 | 该性状的冗余度 |
        # 去除冗余性状的网络参数 | 目标网络的参数 |
      }
      # print(c(un_redundant_trait,'Sensitivity to left trait combination: ',redundancy))
      un_redundant_traits <- c(un_redundant_traits,un_redundant_trait)
      redundant_traits <- redundant_traits[redundant_traits != un_redundant_trait]
      one_row_result <- c(paste0(un_redundant_traits,collapse = '_'),
                          un_redundant_trait,
                          lost_information,
                          add_information,
                          Distance,round(value_of_full_network,6),
                          paste0(all_traits,collapse = '_'))
      
      result <- rbind(result,as.data.frame(t(one_row_result)))
      
    }
  }
  colnames(result)<- c('Kept_traits','informative_trait','lost_information',
                       'add_information','Full_network',
                       'Distance','Full_traits')
  rownames(result) <- NULL
  return(result)
}

## 5.2 identify informative traits ---------------------------------------------
# provide the full network with all traits (as a best target)
# provide a initial small network which need to add new traits
# determine how many traits need to be added
Informative_traits <- function(all_traits,initial_traits,new_traits_number=3,
                               network_number.tif,
                               density.tif,density_weighted.tif,modularity.tif){

  
  # Get the main parameters of the network with all traits
  
  all_traits_number <- trait_list_to_number(trait_list=all_traits,full_node_list=all_traits)
  
  full_network_index <- which(network_number.tif==all_traits_number)
  full_network_density <- density.tif[full_network_index]
  full_network_weighted_density <- density_weighted.tif[full_network_index]
  full_network_modularity <- modularity.tif[full_network_index]
  
  print(c('full_network: ',all_traits_number,full_network_density,full_network_weighted_density,full_network_modularity))
  
  value_of_full_network <- (0.5*full_network_density^2 + 0.5*full_network_weighted_density^2 +
                              full_network_modularity^2)^0.5
  
  # the list of traits (m) which are redundant
  unknown_traits <- get_different_traits(all_traits,initial_traits)
  # new_traits_number = 5
  result <- data.frame()
  if (!is.na(value_of_full_network) & value_of_full_network>0){
    #  a loop or iteration to find the most redundant (m-1) traits combination
    while (new_traits_number>=0) {
      Informative_trait_list = ''
      lost_information = ''
      Distance = 1
      
      # 给定性状数量，生成所有可能的性状组合，找到使得原来网络更加接近最终网络的组合
      # Generate any possible trait combination with this 
      # unknown_traits 为需要添加的性状组合 从new_traits_number开始，生成一个具有随机生成5的性状一行的表
      possible_trait_combination <- as.data.frame(t(combn(unknown_traits,new_traits_number)))
      
      # 测试每一种可添加的性状组合，找到信息损失最少的性状组合，以及最冗余的性状
      for (i in c(1:nrow(possible_trait_combination))) {
        # 将该行性状转化为列表
        new_trait_list <- as.vector(unlist(possible_trait_combination[i,]))
        # relativly unredundant trait network, 将原始性状组合和新性状组合一同计算网络参数
        
        reduced_traits_number <- trait_list_to_number(trait_list=c(initial_traits,new_trait_list),full_node_list=all_traits)
        
        print(c(initial_traits,new_trait_list,reduced_traits_number))
        reduced_network_index <- which(network_number.tif==reduced_traits_number)
        reduced_network_density <- density.tif[reduced_network_index]
        reduced_network_weighted_density <- density_weighted.tif[reduced_network_index]
        reduced_network_modularity <- modularity.tif[reduced_network_index]
        
        # 计算当前组合的网络与目标网络的距离（一般是全网络）
        distance_to_full_network <- (0.5*(full_network_density - reduced_network_density)^2 +
                                       0.5*(full_network_weighted_density - reduced_network_weighted_density)^2 +
                                       (full_network_modularity - reduced_network_modularity)^2)^0.5
        
        # 比较这种性状组合是否为最接近目标网络的组合
        if (distance_to_full_network < Distance){# 不断往下，找到信息损失最少的性状
          Distance =  round(distance_to_full_network,4)
          Informative_trait_list = new_trait_list
          lost_information <- round(Distance/value_of_full_network,6)
          # print(Informative_trait_list)
        }
      }
      if (new_traits_number ==0){
        Informative_trait_list = c('-')
        reduced_traits_number <- trait_list_to_number(trait_list=c(initial_traits),full_node_list=all_traits)
        
        print(c(initial_traits,new_trait_list,reduced_traits_number))
        reduced_network_index <- which(network_number.tif==reduced_traits_number)
        reduced_network_density <- density.tif[reduced_network_index]
        reduced_network_weighted_density <- density_weighted.tif[reduced_network_index]
        reduced_network_modularity <- modularity.tif[reduced_network_index]
        
        # 计算当前组合的网络与目标网络的距离（一般是全网络）
        distance_to_full_network <- (0.5*(full_network_density - reduced_network_density)^2 +
                                       0.5*(full_network_weighted_density - reduced_network_weighted_density)^2 +
                                       (full_network_modularity - reduced_network_modularity)^2)^0.5
        
        # 比较这种性状组合是否为最接近目标网络的组合
        Distance =  round(distance_to_full_network,4)
        lost_information <- round(Distance/value_of_full_network,6)

        }
      one_row_result <- c(paste0(Informative_trait_list,collapse = '_'),
                          lost_information,
                          Distance,
                          paste0(initial_traits,collapse = '_'))
      
      result <- rbind(result,as.data.frame(t(one_row_result)))
      # print(one_row_result)#c(paste0(Informative_trait_list,collapse = '_'),'Distance to the Objective: ',Distance))
      unknown_traits = Informative_trait_list
      new_traits_number = new_traits_number -1 
    }
  }
  colnames(result)<- c('Trait_combination',
                       'lost_information',
                       'Distance to full network',
                       'initial_traits')
  rownames(result) <- NULL
  
  return(result)
}

get_filted_size_tif <- function(threshold=0.05,size.tif,distance.tif,name_number.tif){
  size.tif[size.tif == 99999999] <- NA
  min_size = min(size.tif,na.rm = TRUE)
  max_size = max(size.tif,na.rm = TRUE)
  
  # let the name number of the network which is far from the full network to NA
  for (size_i in c(min_size:max_size)) {
    print(size_i)
    size.tif_mask <- ifelse(size.tif==size_i,1,NA)
    size_distance.tif <- size.tif_mask*distance.tif
    max_in_size <- max(size_distance.tif,na.rm=TRUE)
    min_in_size <- min(size_distance.tif,na.rm=TRUE)
    
    distance_threshold <- min_in_size + threshold*(max_in_size-min_in_size)
    
    name_number_mask <- ifelse(size_distance.tif <= distance_threshold,1,0)
    size.tif[which(name_number_mask==0)] <- NA
  }
  return(size.tif)
}

find_next_network <- function(current_df,next_df){
  next_df_full <- data.frame()
  current_size = max(current_df$size)
  
  # 获取每一个新的网络（大小+1 的网络）
  for (next_i in c(1:nrow(next_df))) {
    next_network <- next_df$network[next_i]
    print(paste0('Check next network ',next_df$network_number[next_i], 
                 ' --- Finished: ',round(next_i/nrow(next_df),3),'---', Sys.time()))
    # 看每一个原来的网络的最后一个,防止出现多次重复的网络（最后一个网络可能一样），浪费计算资源
    for (current_network_number in unique(current_df$network_number)) {
      current_network <- current_df$network[which(current_df$network_number==current_network_number)][1]
      # 如果匹配上
      if(all(is.element(current_network[[1]], next_network[[1]]))){
        next_number = next_df$network_number[next_i]
        # 更新当前所以以这个网络结束的网络
        one_new_row <- subset(current_df,network_number==current_network_number)
        one_new_row$network_number <- next_number
        one_new_row$network <- lapply(one_new_row$network_number,PowerSetsBinary)
        one_new_row$size <- current_size+1
        one_new_row$consistency_id <- paste0(one_new_row$consistency_id,'_',next_number)
        # print(c(next_number,next_network[[1]]))
        
        next_df_full = rbind(next_df_full,one_new_row)
      }
    }
  }
  return(next_df_full)
}


# 99 Functions for data preprocess ---------------------------------------------------------
## 99.1 data filter and scaling with pipe ----------------------
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}
delete.zero <- function(DF, n=0) {
  DF[rowSums(DF==0) <= n,]
}
## 99.2 Scaling --- data filter and scaling with pipe ---------
scale_in_1 <- function(x,b=0.05){
  x_star <- (x-quantile(x,b,na.rm = T))/(quantile(x,1-b,na.rm = T)-quantile(x,b,na.rm = T))
  return(x_star)
}

log_scale <- function(x){
  x_star <- scale(log(x - min(x,na.rm = T) + 1),center = FALSE)[,1]
  return(x_star)
}

log_normalize <- function(x){
  x_star <- BBmisc::normalize(log(x - min(x,na.rm = T) + 1),method = "standardize")
  return(x_star)
}

## 99.3 Other  ---------
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

time_tracker<- function(last_time = Sys.time(),label=''){
  
  current_time =  Sys.time()
  x= difftime(current_time, last_time, tz,
              units = c("auto", "secs", "mins", "hours",
                        "days", "weeks")) # 转换成毫秒、秒、分钟、小时、天、周
  print(paste0('Current time for ',label,': ',x))
  return(current_time)
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
plot_top_performance <- function(variable,top_network_df){
  carss<- summarySE(top_network_df[c("size",variable,"Variable")], 
                    measurevar=variable, groupvars=c("size","Variable"))
  colnames(carss) <- c("size","Variable","N","performance","sd","se","ci")
  carss$d_performance <- NA
  carss$mean_scale <- NA
  for (i in unique(carss$size)) {
    mean_performance <- mean(carss$performance[which(carss$size==i)])
    carss$d_performance[which(carss$size==i)] <- carss$performance[which(carss$size==i)]-mean_performance
    carss$mean[which(carss$size==i)] <- mean_performance
  }
  pd <- position_dodge(0.4) # move them .05 to the left and right
  
  p <- ggplot(carss, aes(x=size, y=d_performance, colour=Variable)) + 
    # geom_errorbar(aes(ymin=d_performance-se, ymax=d_performance+se),
    #               width=1,position=pd) +
    geom_line(position=pd)+
    geom_line(aes(x=size, y=mean/20-0.025),position=pd,
              color='#000000',size=1.2,alpha=0.3) +
    scale_x_continuous(breaks = seq(4, 27,1)) + 
    scale_y_continuous(
      ~.*1,name = "Δ correlation coefficient",
      breaks=seq(-0.05,0.05,by=0.01),
      sec.axis = sec_axis(~.*20,name = 'Average correlation coefficient',
                          breaks=seq(-0.5,0.5,by=0.2),
                          labels=seq(0,1,0.2))) + 
    geom_point(position=pd,size=0.5) +
    labs(title = paste0("Top matched network (10%) in ", variable),
         x = "Number of traits") + theme_bw() + 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(axis.text.y.right = element_text(color = '#777777')) +
    theme(axis.text.y.left = element_text(color = '#997755'))
  
  return(p)
}
