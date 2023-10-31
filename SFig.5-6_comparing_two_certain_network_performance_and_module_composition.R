rm(list = ls())
library(latex2exp)
source('./0_Loading_variables.R')

## ------------------------------------------------------------------------------------------------------
## PCA: Requires gap-filled plant trait data 
## ------------------------------------------------------------------------------------------------------

# 5 Performance and module composition in two certain networks----
## 5.1 commonly used 6 traits and best performing 6 traits -------
c("LA","Lthick","Llen","Lwid","SLA","LWC","LFM","L15N","LDMC","LCC.m","LNC.m",
  "LPC.m","LCNR","SSD","SCDen","SCDia","SCElen","WDL","RD","SRL","SdM","SdNum",
  "SdL","SdGE","DispL","Height","SD")

trait_common_number = 33563921
Diaz_trait_number <- trait_list_to_number(c('SSD','Height','SLA','SdM','LNC.m','LA'))
trait_best_performing_number = 675984

trait_best_performing <- paste0(PowerSetsBinary(trait_best_performing_number),collapse = '_') 
trait_common <- paste0(PowerSetsBinary(trait_common_number),collapse = '_')
trait_diaz <- paste0(PowerSetsBinary(Diaz_trait_number),collapse = '_')

trait_data <- ecoregion_trait_with_env[c(seq(2,26),32,33)]
###5.1.1 Common used 6 traits ------------------------------------------------------
network_info <- data_to_network(data = trait_data,trait_combination = trait_common)
network_i <- network_info$Net
ceb <- cluster_optimal(network_i) #  (
# dendPlot(ceb)

jpeg(filename = paste0('./_Results/Figures/6traitsNetowrk/',trait_common_number,'Network_module_trait_common.jpg'),
     width = 1000, height = 1000, units = "px", pointsize = 30)
plot(ceb, network_i)
text(0,-1.5,"Common used trait network")
dev.off()

pdf(paste0('./_Results/Figures/6traitsNetowrk/',trait_common_number,'Network_module_trait_common.pdf'),
    width = 10, height = 10, pointsize = 20)
plot(ceb, network_i)
text(0,-1.5,"Common used trait network")
dev.off()



###5.1.2 Diaz 6 traits ---------------------------------------------------
network_info <- data_to_network(data = trait_data,trait_combination = trait_diaz)
network_i <- network_info$Net
ceb <- cluster_optimal(network_i) #  (
# dendPlot(ceb)

jpeg(filename = paste0('./_Results/Figures/6traitsNetowrk/',Diaz_trait_number,'Network_module_trait_plant_form_spectrum.jpg'),
     width = 1000, height = 1000, units = "px", pointsize = 20)
plot(ceb, network_i)
text(0,-1.5,"Plant form spectrum trait network")
dev.off()

pdf(paste0('./_Results/Figures/6traitsNetowrk/',Diaz_trait_number,'Network_module_trait_plant_form_spectrum.pdf'),
    width = 10, height = 10, pointsize = 20)
plot(ceb, network_i)
text(0,-1.5,"Plant form spectrum trait network")
dev.off()

###5.1.3 best_performing 6 traits ------------------------------------------------------
network_info <- data_to_network(data = trait_data,trait_combination = trait_best_performing)
network_i <- network_info$Net
ceb <- cluster_optimal(network_i) #  (
# dendPlot(ceb)

jpeg(filename = paste0('./_Results/Figures/6traitsNetowrk/',trait_best_performing_number,'Network_module_trait_best_performing.jpg'),
     width = 1000, height = 1000, units = "px", pointsize = 20)
plot(ceb, network_i)
text(0,-1.5,"Best performing trait network")
dev.off()

pdf(paste0('./_Results/Figures/6traitsNetowrk/',trait_best_performing_number,'Network_module_trait_best_performing.pdf'),
     width = 10, height = 10, pointsize = 20)

plot(ceb, network_i)
text(0,-1.5,"Best performing trait network")
dev.off()

### Network performance --------------------------------------------------

full_network = read.csv('./tips/Network_metrics_all_ecoregion_27.csv')
colnames(full_network) <- c(paste0(colnames(full_network)[1:7],27),'gradient')
df = ecoregion_trait_with_env
gradient='gradient'
min_species = 20
overlab=0
interval = 1

for (number in c(trait_common_number,trait_best_performing_number,Diaz_trait_number)) {
  node_list = PowerSetsBinary(number)
  
  if (length(node_list)>2){
    trait_combination = paste0(node_list,collapse = '_')
    Net.metrics.final <- data.frame()
    
    gradient.min = min(df[gradient]);gradient.max = max(df[gradient])
    gradient_interval <- interval
    gradient_list <- seq(gradient.min,gradient.max,by = gradient_interval)
    for (gradient.start in gradient_list) {
      gradient.end = gradient.start + gradient_interval*(overlab+1)
      if (gradient.end <= max(gradient_list)){
        one_subset_data <- df[which(df[gradient]>=gradient.start & df[gradient]<gradient.end),]
        species_number <- sum(!duplicated(one_subset_data$Species))
        
        if (species_number>=min_species){
          # print(c(trait_combination,gradient.start,gradient.end,nrow(one_subset_data)))
          # calculate the full network
          network_metrics_One <- Calculate_network_metrics_inner(data=one_subset_data %>% dplyr::select(node_list),
                                                                 trait_combination=trait_combination,
                                                                 network_number = number)
          network_metrics_One[gradient] = gradient.start
          Net.metrics.final <- rbind(Net.metrics.final,network_metrics_One)
        }
      }
    }
    write.csv(Net.metrics.final,paste0('./_Results/Data/Network_series/',number,'_network_metrics.csv'),row.names = FALSE)
  }
  print(number)
}

for (number in c(trait_common_number,trait_best_performing_number,Diaz_trait_number)) {
  reduced_network <- read.csv(paste0('./_Results/Data/Network_series/',number,'_network_metrics.csv'))
  one_reduced_subset <- dplyr::inner_join(reduced_network,full_network)
  p_c <- ggplot(one_reduced_subset,aes(x=density_weighted,y=density_weighted27)) + geom_point(alpha=0.2,color = '#CC79A7') +
    geom_smooth(method = 'lm',se=FALSE,color = '#CC79A7',size=0.4) +
    stat_poly_eq(
      aes(label = paste(..rr.label.., sep = '~~~~')),
      formula = y ~ x,  parse = TRUE,
      size = 5, #公式字体大小
      label.x = 0.1,  #位置 ，0-1之间的比例
      label.y = 0.95,color = '#CC79A7')+
    labs(x = "Connectivity of redcued network",
         y = "Connectivity of full network") + theme_classic()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
          axis.title.x=element_blank(),axis.title.y=element_blank())
  p_m <- ggplot(one_reduced_subset,aes(x=modularity_optimal,y=modularity_optimal27)) + geom_point(alpha=0.2,color='#0072B2') +
    geom_smooth(method = 'lm',se=FALSE,color='#0072B2',size=0.4)+
    stat_poly_eq(
      aes(label = paste(..rr.label.., sep = '~~~~')),
      formula = y ~ x,  parse = TRUE,
      size = 5, #公式字体大小
      label.x = 0.1,  #位置 ，0-1之间的比例
      label.y = 0.95,color='#0072B2')+
    labs(x = "Modularity of redcued network",
         y = "Modularity of full network") + theme_classic()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
          axis.title.x=element_blank(),axis.title.y=element_blank())
  
  P_CM <- ggarrange(p_c,p_m,ncol = 1)
  
  ggsave(paste0('./_Results/Figures/6traitsNetowrk/',number,'_performance.jpg'),P_CM,device = "jpeg",
         width = 6,height = 6,units ="cm",dpi = 300)
  ggsave(paste0('./_Results/Figures/6traitsNetowrk/',number,'_performance.pdf'),P_CM,device = "pdf",
         width = 6,height = 6,units ="cm",dpi = 300)
}


### Network PCA --------------------------------------------------
ecoregion_trait_all <- ecoregion_trait_with_env[c(seq(2,26),32,33)]

for (number in c(trait_common_number,trait_best_performing_number,Diaz_trait_number)) {
  traits <- PowerSetsBinary(number)
  
  network_data <- ecoregion_trait_all[traits]
  
  # PCA计算
  pca_result <- prcomp(network_data, scale.= TRUE)
  
  biplot_data <- as.data.frame(pca_result$x)  # Extract PCA scores (individuals' coordinates)
  biplot_vars <- as.data.frame(pca_result$rotation)  # Extract PCA loadings (variables' coordinates)
  biplot_vars$trait <- rownames(biplot_vars)
  biplot_vars$color_ <- "#565656"
  rownames(biplot_vars) <- biplot_vars$trait
  
  scale_times <- as.integer(1.5*mean(abs(pca_result$x[,c(1,2)]))/mean(abs(pca_result$rotation[,c(1,2)])))
  plot_size <- 1.4 * max(scale_times* abs(pca_result$rotation[,c(1,2)])) 
  p1 <- ggplot() +
    geom_point(data = biplot_data, aes(x = PC1, y = PC2), size = 2, color = "gray",alpha=0.04) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(data = biplot_vars, aes(x = 0, y = 0, xend = PC1*scale_times, yend = PC2*scale_times, color = color_),
                 arrow = arrow(length = unit(0.1, "inches"))) +
    geom_text(data = biplot_vars, aes(x = PC1*(scale_times + 0.1* plot_size), 
                                      y = PC2*(scale_times + 0.1* plot_size), 
                                      label = rownames(biplot_vars),color = color_),
              size = 4) +
    labs(x = paste0("PC1  (",formattable::percent(pca_result$sdev[1]/sum(pca_result$sdev)),')'), 
         y = paste0("PC2  (",formattable::percent(pca_result$sdev[2]/sum(pca_result$sdev)),')')) +
    scale_x_continuous(limits = c(0 - plot_size,plot_size))+
    scale_y_continuous(limits = c(0 - plot_size,plot_size))+
    scale_color_identity()+
    theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       legend.position = 'None',axis.text.x= element_blank(),axis.text.y= element_blank(),
                       axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
  
  p2 <- 
    ggplot() +
    geom_point(data = biplot_data, aes(x = PC2, y = PC3), size = 2, color = "gray",alpha=0.04) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(data = biplot_vars, aes(x = 0, y = 0, xend = PC2*scale_times, yend = PC3*scale_times, color = color_),
                 arrow = arrow(length = unit(0.1, "inches"))) +
    geom_text(data = biplot_vars, aes(x = PC2*(scale_times+ 0.1* plot_size), 
                                      y = PC3*(scale_times+ 0.05* plot_size), label = rownames(biplot_vars), color = color_),
              size = 4) +
    labs(x = paste0("PC2  (",formattable::percent(pca_result$sdev[2]/sum(pca_result$sdev)),')'), 
         y = paste0("PC3  (",formattable::percent(pca_result$sdev[3]/sum(pca_result$sdev)),')')) +
    scale_x_continuous(limits = c(0 - plot_size,plot_size))+
    scale_y_continuous(limits = c(0 - plot_size,plot_size))+
    scale_color_identity()+
    theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       legend.position = 'None',axis.text.x= element_blank(),axis.text.y= element_blank(),
                       axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
  ggarrange(p1,p2,ncol = 2)
  
  ggsave(paste0('./_Results/Figures/6traitsNetowrk/',number,'.jpg'),width = 22 , height = 10, units = "cm",dpi = 600)
  ggsave(paste0('./_Results/Figures/6traitsNetowrk/',number,'.pdf'),width = 22 , height = 10, units = "cm",dpi = 600)
}


## 5.2 best worthy 10 traits and best performing 10 traits -------
c("LA","Lthick","Llen","Lwid","SLA","LWC","LFM","L15N","LDMC","LCC.m","LNC.m",
  "LPC.m","LCNR","SSD","SCDen","SCDia","SCElen","WDL","RD","SRL","SdM","SdNum",
  "SdL","SdGE","DispL","Height","SD")


trait_data <- ecoregion_trait_with_env[c(seq(2,26),32,33)]
###5.2.1 worthy used traits ------------------------------------------------------

worthy_number = 110355717# 110265605# 110355717 

for (worthy_number in best_worthy_networks_series_marginal_benefit$number) {
  

  trait_worthy <- paste0(PowerSetsBinary(worthy_number),collapse = '_')
  
  network_info <- data_to_network(data = trait_data,trait_combination = trait_worthy)
  network_i <- network_info$Net
  ceb <- cluster_optimal(network_i) #  cluster_louvain cluster_walktrap
  # dendPlot(ceb)
  
  jpeg(filename = paste0('./_Results/Figures/Best_worthy_networks/worthy',length(PowerSetsBinary(worthy_number)),'.jpg'),
       width = 1000, height = 1000, units = "px", pointsize = 20)
  plot(ceb, network_i)
  text(0,-1.5,"Cost-efficient trait network")
  dev.off()
  
  pdf(paste0('./_Results/Figures/Best_worthy_networks/worthy',length(PowerSetsBinary(worthy_number)),'.pdf'),
       width = 10, height = 10, pointsize = 20)
  plot(ceb, network_i)
  text(0,-1.5,"Cost-efficient trait network")
  dev.off()

}
###5.2.2 best_performing used traits ------------------------------------------------------
best_number = 126572696# 109793424# 126572696

trait_best_performing <- paste0(PowerSetsBinary(best_number),collapse = '_') 

network_info <- data_to_network(data = trait_data,trait_combination = trait_best_performing)
network_i <- network_info$Net
ceb <- cluster_optimal(network_i) #  (
# dendPlot(ceb)

jpeg(filename = paste0('./_Results/Figures/16_traits_cost_efficient_comparing/Network_module_trait_best_performing',best_number,'.jpg'),
     width = 1000, height = 1000, units = "px", pointsize = 20)
plot(ceb, network_i)
text(0,-1.5,"Best performing trait network")
dev.off()

pdf(paste0('./_Results/Figures/16_traits_cost_efficient_comparing/Network_module_trait_best_performing',best_number,'.pdf'),
     width = 10, height = 10, pointsize = 20)
plot(ceb, network_i)
text(0,-1.5,"Best performing trait network")
dev.off()

###5.2.3 Network performance -------------------------------------------------------
full_network = read.csv('./tips/Network_metrics_all_ecoregion_27.csv')
colnames(full_network) <- c(paste0(colnames(full_network)[1:7],27),'gradient')
df = ecoregion_trait_with_env
gradient='gradient'
min_species = 20
overlab=0
interval = 1

for (number in c(best_number,worthy_number)) {
  node_list = PowerSetsBinary(number)
  
  if (length(node_list)>2){
    trait_combination = paste0(node_list,collapse = '_')
    Net.metrics.final <- data.frame()
    
    gradient.min = min(df[gradient]);gradient.max = max(df[gradient])
    gradient_interval <- interval
    gradient_list <- seq(gradient.min,gradient.max,by = gradient_interval)
    for (gradient.start in gradient_list) {
      gradient.end = gradient.start + gradient_interval*(overlab+1)
      if (gradient.end <= max(gradient_list)){
        one_subset_data <- df[which(df[gradient]>=gradient.start & df[gradient]<gradient.end),]
        species_number <- sum(!duplicated(one_subset_data$Species))
        
        if (species_number>=min_species){
          # print(c(trait_combination,gradient.start,gradient.end,nrow(one_subset_data)))
          # calculate the full network
          network_metrics_One <- Calculate_network_metrics_inner(data=one_subset_data %>% dplyr::select(node_list),
                                                                 trait_combination=trait_combination,
                                                                 network_number = number)
          network_metrics_One[gradient] = gradient.start
          Net.metrics.final <- rbind(Net.metrics.final,network_metrics_One)
        }
      }
    }
    write.csv(Net.metrics.final,paste0('./_Results/Data/Network_series/',number,'_network_metrics.csv'))
  }
  
  print(number)
}

for (number in c(best_number,worthy_number)) {
  reduced_network <- read.csv(paste0('./_Results/Data/Network_series/',number,'_network_metrics.csv'))
  one_reduced_subset <- dplyr::inner_join(reduced_network,full_network)
  p_c <- ggplot(one_reduced_subset,aes(x=density_weighted,y=density_weighted27)) + geom_point(alpha=0.2,color = '#CC79A7') +
    geom_smooth(method = 'lm',se=FALSE,color = '#CC79A7',size=0.4) +
    stat_poly_eq(
      aes(label = paste(..rr.label.., sep = '~~~~')),
      formula = y ~ x,  parse = TRUE,
      size = 5, #公式字体大小
      label.x = 0.1,  #位置 ，0-1之间的比例
      label.y = 0.95,color = '#CC79A7')+
    labs(x = "Connectivity of redcued network",
         y = "Connectivity of full network") + theme_classic()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
          axis.title.x=element_blank(),axis.title.y=element_blank())
  p_m <- ggplot(one_reduced_subset,aes(x=modularity_optimal,y=modularity_optimal27)) + geom_point(alpha=0.2,color='#0072B2') +
    geom_smooth(method = 'lm',se=FALSE,color='#0072B2',size=0.4)+
    stat_poly_eq(
      aes(label = paste(..rr.label.., sep = '~~~~')),
      formula = y ~ x,  parse = TRUE,
      size = 5, #公式字体大小
      label.x = 0.1,  #位置 ，0-1之间的比例
      label.y = 0.95,color='#0072B2')+
    labs(x = "Modularity of redcued network",
         y = "Modularity of full network") + theme_classic()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
          axis.title.x=element_blank(),axis.title.y=element_blank())
  
  P_CM <- ggarrange(p_c,p_m,ncol = 1)
  
  ggsave(paste0('./_Results/Figures/16_traits_cost_efficient_comparing/',number,'_performance.jpg'),P_CM,device = "jpeg",
         width = 6,height = 6,units ="cm",dpi = 300)
  ggsave(paste0('./_Results/Figures/16_traits_cost_efficient_comparing/',number,'_performance.pdf'),P_CM,device = "pdf",
         width = 6,height = 6,units ="cm",dpi = 300)
}


### Network PCA --------------------------------------------------
ecoregion_trait_all <- ecoregion_trait_with_env[c(seq(2,26),32,33)]

for (number in c(best_number,worthy_number)) {
  traits <- PowerSetsBinary(number)
  
  network_data <- ecoregion_trait_all[traits]
  
  # PCA计算
  pca_result <- prcomp(network_data, scale.= TRUE)
  
  biplot_data <- as.data.frame(pca_result$x)  # Extract PCA scores (individuals' coordinates)
  biplot_vars <- as.data.frame(pca_result$rotation)  # Extract PCA loadings (variables' coordinates)
  biplot_vars$trait <- rownames(biplot_vars)
  biplot_vars$color_ <- "#565656"
  rownames(biplot_vars) <- biplot_vars$trait
  
  scale_times <- as.integer(1.5*mean(abs(pca_result$x[,c(1,2)]))/mean(abs(pca_result$rotation[,c(1,2)])))
  plot_size <- 1.4 * max(scale_times* abs(pca_result$rotation[,c(1,2)])) 
  p1 <- ggplot() +
    geom_point(data = biplot_data, aes(x = PC1, y = PC2), size = 2, color = "gray",alpha=0.04) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(data = biplot_vars, aes(x = 0, y = 0, xend = PC1*scale_times, yend = PC2*scale_times, color = color_),
                 arrow = arrow(length = unit(0.1, "inches"))) +
    geom_text(data = biplot_vars, aes(x = PC1*(scale_times + 0.1* plot_size), 
                                      y = PC2*(scale_times + 0.1* plot_size), 
                                      label = rownames(biplot_vars),color = color_),
              size = 4) +
    labs(x = paste0("PC1  (",formattable::percent(pca_result$sdev[1]/sum(pca_result$sdev)),')'), 
         y = paste0("PC2  (",formattable::percent(pca_result$sdev[2]/sum(pca_result$sdev)),')')) +
    scale_x_continuous(limits = c(0 - plot_size,plot_size))+
    scale_y_continuous(limits = c(0 - plot_size,plot_size))+
    scale_color_identity()+
    theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       legend.position = 'None',axis.text.x= element_blank(),axis.text.y= element_blank(),
                       axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
  
  p2 <- 
    ggplot() +
    geom_point(data = biplot_data, aes(x = PC2, y = PC3), size = 2, color = "gray",alpha=0.04) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(data = biplot_vars, aes(x = 0, y = 0, xend = PC2*scale_times, yend = PC3*scale_times, color = color_),
                 arrow = arrow(length = unit(0.1, "inches"))) +
    geom_text(data = biplot_vars, aes(x = PC2*(scale_times+ 0.1* plot_size), 
                                      y = PC3*(scale_times+ 0.1* plot_size), label = rownames(biplot_vars), color = color_),
              size = 4) +
    labs(x = paste0("PC2  (",formattable::percent(pca_result$sdev[2]/sum(pca_result$sdev)),')'), 
         y = paste0("PC3  (",formattable::percent(pca_result$sdev[3]/sum(pca_result$sdev)),')')) +
    scale_x_continuous(limits = c(0 - plot_size,plot_size))+
    scale_y_continuous(limits = c(0 - plot_size,plot_size))+
    scale_color_identity()+
    theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       legend.position = 'None',axis.text.x= element_blank(),axis.text.y= element_blank(),
                       axis.ticks.x = element_blank(),axis.ticks.y = element_blank())
  ggarrange(p1,p2,ncol = 2)
  
  ggsave(paste0('./_Results/Figures/16_traits_cost_efficient_comparing/',number,'.jpg'),width = 22 , height = 10, units = "cm",dpi = 600)
  ggsave(paste0('./_Results/Figures/16_traits_cost_efficient_comparing/',number,'.pdf'),width = 22 , height = 10, units = "cm",dpi = 600)
}


