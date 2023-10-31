rm(list = ls())
library(latex2exp)
source('./0_Loading_variables.R')


## Fig. 3-1  Modules composition of networks --------
best_matched_network_series <- read.csv('./_Results/Data/Network_series/best_performing_network_series_0.05.csv') %>% 
  mutate(performance = sqrt(sqrt(density_R2)*sqrt(density_weighted_R2))*sqrt(modularity_optimal_R2)) 

all_edges<- read.csv('./tips/all_edges.csv')

metric_df <- data.frame()
i = 5
for (number in best_matched_network_series$number) {
  
  if (!is.na(number)){
    node_list = PowerSetsBinary(number)
  } else if (!is.na(trait_combination)){
    node_list <- strsplit(trait_combination,'_')[[1]]
  } else{
    print("no network number, neither trait combination")
    return(NA)
  }
  node_value_sd <- as.data.frame(node_list)
  
  edges <- subset(all_edges,from %in% node_list & to %in% node_list) 
  
  
  network_i <- graph_from_data_frame(d=edges,vertices=node_value_sd, directed = FALSE)
  ceb <- cluster_optimal(network_i) #  (cluster_optimal
  # dendPlot(ceb)
  
  # jpeg(filename = paste0('./_Results/Figures/Network_module/all_optimal_',i,'_netowrk_module.jpg'),
  #      width = 1000, height = 1000, units = "px", pointsize = 30)
  pdf(paste0('./_Results/Figures/Network_module/all_optimal_',i,'_netowrk_module.pdf'),
      width = 10, height = 10,  pointsize = 20)

  
  size <- length(PowerSetsBinary(number))
  
  density <- 2*nrow(edges)/(size*(size-1))
  connecitivty <- round(2*sum(edges$weight)/(size*(size-1)),3)
  modularity <- round(modularity(ceb),3)
  
  
  metric_df <- rbind(metric_df,data.frame(
    size = size,
    density = density,
    connecitivty = connecitivty,
    modularity = modularity))
  plot(ceb, network_i)
  text(-1,1,size)
  # 
  # plot(ceb, network_i,main=paste0('global_',size),
  #      sub = paste0("Connecitivty: ",connecitivty,'   ',"Modularity: ",modularity))
  
  dev.off()
  
  jpeg(filename = paste0('./_Results/Figures/Network_module/best_performing_all_optimal_',i,'_netowrk_module.jpg'),
       width = 1000, height = 1000, units = "px", pointsize = 30)
  # pdf(paste0('./_Results/Figures/Network_module/all_optimal_',i,'_netowrk_module.pdf'),
  #     width = 10, height = 10,  pointsize = 20)

  size <- length(PowerSetsBinary(number))
  
  density <- 2*nrow(edges)/(size*(size-1))
  connecitivty <- round(2*sum(edges$weight)/(size*(size-1)),3)
  modularity <- round(modularity(ceb),3)
  
  
  metric_df <- rbind(metric_df,data.frame(
    size = size,
    density = density,
    connecitivty = connecitivty,
    modularity = modularity))
  plot(ceb, network_i)
  text(-1,1,size)
  # 
  # plot(ceb, network_i,main=paste0('global_',size),
  #      sub = paste0("Connecitivty: ",connecitivty,'   ',"Modularity: ",modularity))
  i = i +1
  
  dev.off()
}

# write.csv(metric_df,'./_Results/Data/Network_series/best_performing_network_series_metrics.csv',row.names = FALSE)

## ------------------------------------------------------------------------------------------------------
## PCA: Requires gap-filled plant trait data 
## ------------------------------------------------------------------------------------------------------
###  --------

# ecoregion_trait_all <- ecoregion_trait_with_env[c(seq(2,26),32,33)]
# Netowrk_number_all <- best_matched_network_series$number
# 
# # network_number <- c(675984,109793424,127621276,127916447,134217727)
# 
# color_df <- data.frame()
# for (trait in PowerSetsBinary(Netowrk_number_all[27 - 4])) {
#   one_trait_color = '#c0f0fF'
#   if (trait %in% PowerSetsBinary(Netowrk_number_all[20 - 4] )) {
#     one_trait_color <- '#8ff1f9'
#   }
#   if (trait %in% PowerSetsBinary(Netowrk_number_all[15 - 4] )) {
#     one_trait_color <- '#6baac1'
#   }
#   if (trait %in% PowerSetsBinary(Netowrk_number_all[10 - 4] )) {
#     one_trait_color <- '#5c6fa0'
#   }
#   if (trait %in% PowerSetsBinary(Netowrk_number_all[6 - 4])) {
#     one_trait_color <- '#99146d'
#   }
#   one_trait_df <- data.frame(trait = trait,color_ = one_trait_color)
#   color_df <- rbind(color_df,one_trait_df)
# }
# 
# 
# number= 134217727
# for (number in Netowrk_number_all) {
#   traits <- PowerSetsBinary(number)
#   
#   network_data <- ecoregion_trait_all[traits]
#   
#   # PCA计算
#   pca_result <- prcomp(network_data, scale.= TRUE)
#   
#   biplot_data <- as.data.frame(pca_result$x)  # Extract PCA scores (individuals' coordinates)
#   biplot_vars <- as.data.frame(pca_result$rotation)  # Extract PCA loadings (variables' coordinates)
#   biplot_vars$trait <- rownames(biplot_vars)
#   biplot_vars <- left_join(biplot_vars,color_df)
#   rownames(biplot_vars) <- biplot_vars$trait
#   
#   scale_times <- as.integer(1.5*mean(abs(pca_result$x[,c(1,2,3)]))/mean(abs(pca_result$rotation[,c(1,2,3)])))
#   plot_size <- 1.4 * max(scale_times* abs(pca_result$rotation[,c(1,2,3)])) 
#   p1 <- ggplot() +
#     geom_point(data = biplot_data, aes(x = PC1, y = PC2), size = 2, color = "gray",alpha=0.04) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_segment(data = biplot_vars, aes(x = 0, y = 0, xend = PC1*scale_times, yend = PC2*scale_times, color = color_),
#                  arrow = arrow(length = unit(0.1, "inches"))) +
#     geom_text(data = biplot_vars, aes(x = PC1*(scale_times + 0.1* plot_size), 
#                                       y = PC2*(scale_times + 0.1* plot_size), 
#                                       label = rownames(biplot_vars),color = color_),
#               size = 4) +
#     labs(x = paste0("PC1 (Hydraulic safety strategy",formattable::percent(pca_result$sdev[1]/sum(pca_result$sdev)),')'), 
#          y = paste0("PC2 (Leaf economic strategy ",formattable::percent(pca_result$sdev[2]/sum(pca_result$sdev)),')')) +
#     scale_x_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_y_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_color_identity()+
#     theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#                        legend.position = 'None',axis.text.x= element_blank(),axis.text.y= element_blank(),
#                        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
#   
#   p2 <- 
#     ggplot() +
#     geom_point(data = biplot_data, aes(x = PC2, y = PC3), size = 2, color = "gray",alpha=0.04) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_segment(data = biplot_vars, aes(x = 0, y = 0, xend = PC2*scale_times, yend = PC3*scale_times, color = color_),
#                  arrow = arrow(length = unit(0.1, "inches"))) +
#     geom_text(data = biplot_vars, aes(x = PC2*(scale_times+ 0.1* plot_size), 
#                                       y = PC3*(scale_times+ 0.1* plot_size), label = rownames(biplot_vars), color = color_),
#               size = 4) +
#     labs(x = paste0("PC2 (Leaf economic strategy ",formattable::percent(pca_result$sdev[2]/sum(pca_result$sdev)),')'), 
#          y = paste0("PC3 (Plant reproduction and competition ",formattable::percent(pca_result$sdev[3]/sum(pca_result$sdev)),')')) +
#     scale_x_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_y_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_color_identity()+
#     theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#                        legend.position = 'None',axis.text.x= element_blank(),axis.text.y= element_blank(),
#                        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
#   ggarrange(p1,p2,ncol = 2)
#   
#   ggsave(paste0('./_Results/Figures/PCA/best_performing_all_',length(traits),'.jpg'),width = 22 , height = 10, units = "cm",dpi = 600)
#   ggsave(paste0('./_Results/Figures/PCA/best_performing_all_',length(traits),'.pdf'),width = 22 , height = 10, units = "cm",dpi = 600)
# }
# 
# #arrow
# for (number in Netowrk_number_all) { # 
#   traits <- PowerSetsBinary(number)
#   
#   network_data <- ecoregion_trait_all[traits]
#   
#   # PCA计算
#   pca_result <- prcomp(network_data, scale.= TRUE)
#   
#   biplot_data <- as.data.frame(pca_result$x)  # Extract PCA scores (individuals' coordinates)
#   biplot_vars <- as.data.frame(pca_result$rotation)  # Extract PCA loadings (variables' coordinates)
#   biplot_vars$trait <- rownames(biplot_vars)
#   biplot_vars <- left_join(biplot_vars,color_df)
#   rownames(biplot_vars) <- biplot_vars$trait
#   
#   scale_times <- as.integer(1.5*mean(abs(pca_result$x[,c(1,2,3)]))/mean(abs(pca_result$rotation[,c(1,2,3)])))
#   plot_size <- 1.4 * max(scale_times* abs(pca_result$rotation[,c(1,2,3)])) 
#   p1 <- ggplot() +
#     # geom_point(data = biplot_data, aes(x = PC1, y = PC2), size = 2, color = "gray",alpha=0.04) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_segment(data = biplot_vars, aes(x = 0, y = 0, xend = PC1*scale_times, yend = PC2*scale_times, color = color_),
#                  arrow = arrow(length = unit(0.1, "inches"))) +
#     geom_text(data = biplot_vars, aes(x = PC1*(scale_times + 0.1* plot_size), 
#                                       y = PC2*(scale_times + 0.1* plot_size), 
#                                       label = rownames(biplot_vars),color = color_),
#               size = 4) +
#     labs(x = paste0("PC1 (Hydraulic safety strategy",formattable::percent(pca_result$sdev[1]/sum(pca_result$sdev)),')'), 
#          y = paste0("PC2 (Leaf economic strategy ",formattable::percent(pca_result$sdev[2]/sum(pca_result$sdev)),')')) +
#     scale_x_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_y_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_color_identity()+
#     theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#                        legend.position = 'None',axis.text.x= element_blank(),axis.text.y= element_blank(),
#                        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
#   
#   p2 <- 
#     ggplot() +
#     # geom_point(data = biplot_data, aes(x = PC2, y = PC3), size = 2, color = "gray",alpha=0.04) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_segment(data = biplot_vars, aes(x = 0, y = 0, xend = PC2*scale_times, yend = PC3*scale_times, color = color_),
#                  arrow = arrow(length = unit(0.1, "inches"))) +
#     geom_text(data = biplot_vars, aes(x = PC2*(scale_times+ 0.1* plot_size), 
#                                       y = PC3*(scale_times+ 0.1* plot_size), label = rownames(biplot_vars), color = color_),
#               size = 4) +
#     labs(x = paste0("PC2 (Leaf economic strategy ",formattable::percent(pca_result$sdev[2]/sum(pca_result$sdev)),')'), 
#          y = paste0("PC3 (Plant reproduction and competition ",formattable::percent(pca_result$sdev[3]/sum(pca_result$sdev)),')')) +
#     scale_x_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_y_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_color_identity()+
#     theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#                        legend.position = 'None',axis.text.x= element_blank(),axis.text.y= element_blank(),
#                        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
#   ggarrange(p1,p2,ncol = 2)
#   
#   ggsave(paste0('./_Results/Figures/PCA/best_performing_arrow_',length(traits),'.jpg'),width = 22 , height = 10, units = "cm",dpi = 600)
#   ggsave(paste0('./_Results/Figures/PCA/best_performing_arrow_',length(traits),'.pdf'),width = 22 , height = 10, units = "cm",dpi = 600)
# }
# p2
# # point
# for (number in Netowrk_number_all) {
#   traits <- PowerSetsBinary(number)
#   
#   network_data <- ecoregion_trait_all[traits]
#   
#   # PCA计算
#   pca_result <- prcomp(network_data, scale.= TRUE)
#   
#   biplot_data <- as.data.frame(pca_result$x)  # Extract PCA scores (individuals' coordinates)
#   biplot_vars <- as.data.frame(pca_result$rotation)  # Extract PCA loadings (variables' coordinates)
#   biplot_vars$trait <- rownames(biplot_vars)
#   biplot_vars <- left_join(biplot_vars,color_df)
#   rownames(biplot_vars) <- biplot_vars$trait
#   
#   scale_times <- as.integer(1.5*mean(abs(pca_result$x[,c(1,2)]))/mean(abs(pca_result$rotation[,c(1,2)])))
#   plot_size <- 1.4 * max(scale_times* abs(pca_result$rotation[,c(1,2)])) 
#   p1 <- ggplot() +
#     geom_point(data = biplot_data, aes(x = PC1, y = PC2), size = 2, color = "gray",alpha=0.04) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#     # geom_segment(data = biplot_vars, aes(x = 0, y = 0, xend = PC1*scale_times, yend = PC2*scale_times, color = color_),
#     #              arrow = arrow(length = unit(0.1, "inches"))) +
#     # geom_text(data = biplot_vars, aes(x = PC1*(scale_times + 0.1* plot_size), 
#     #                                   y = PC2*(scale_times + 0.1* plot_size), 
#     #                                   label = rownames(biplot_vars),color = color_),
#     #           size = 4) +
#     labs(x = paste0("PC1 (Hydraulic safety strategy",formattable::percent(pca_result$sdev[1]/sum(pca_result$sdev)),')'), 
#          y = paste0("PC2 (Leaf economic strategy ",formattable::percent(pca_result$sdev[2]/sum(pca_result$sdev)),')')) +
#     scale_x_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_y_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_color_identity()+
#     theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#                        legend.position = 'None',axis.text.x= element_blank(),axis.text.y= element_blank(),
#                        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
#   
#   p2 <- 
#     ggplot() +
#     geom_point(data = biplot_data, aes(x = PC2, y = PC3), size = 2, color = "gray",alpha=0.04) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#     # geom_segment(data = biplot_vars, aes(x = 0, y = 0, xend = PC2*scale_times, yend = PC3*scale_times, color = color_),
#     #              arrow = arrow(length = unit(0.1, "inches"))) +
#     # geom_text(data = biplot_vars, aes(x = PC2*(scale_times+ 0.1* plot_size), 
#     #                                   y = PC3*(scale_times+ 0.1* plot_size), label = rownames(biplot_vars), color = color_),
#     #           size = 4) +
#     labs(x = paste0("PC2 (Leaf economic strategy ",formattable::percent(pca_result$sdev[2]/sum(pca_result$sdev)),')'), 
#          y = paste0("PC3 (Plant reproduction and competition ",formattable::percent(pca_result$sdev[3]/sum(pca_result$sdev)),')')) +
#     scale_x_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_y_continuous(limits = c(0 - plot_size,plot_size))+
#     scale_color_identity()+
#     theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#                        legend.position = 'None',axis.text.x= element_blank(),axis.text.y= element_blank(),
#                        axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
#   ggarrange(p1,p2,ncol = 2)
#   
#   ggsave(paste0('./_Results/Figures/PCA/best_performing_points_',length(traits),'.jpg'),width = 22 , height = 10, units = "cm",dpi = 600)
#   ggsave(paste0('./_Results/Figures/PCA/best_performing_points_',length(traits),'.pdf'),width = 22 , height = 10, units = "cm",dpi = 600)
# }


