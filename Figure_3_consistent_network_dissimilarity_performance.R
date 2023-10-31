rm(list = ls())
library(latex2exp)
source('./0_Loading_variables.R')


# 2 The consistent networks ------
# best_performing_network_series_0.05.csv
# best_matched_network_series_0.05.csv
best_matched_network <- read.csv('./_Results/Data/Network_series/best_matched_network_series_0.05.csv') %>% 
  mutate(performance = sqrt(sqrt(density_R2)*sqrt(density_weighted_R2))*sqrt(modularity_optimal_R2)) 
best_matched_network$type = 'Best matched network series'
best_matched_network$rgb = '#a5a5a5'
best_performing_network <- read.csv('./_Results/Data/Network_series/best_performing_network_series_0.05.csv')
best_performing_network$type = 'Best performing network series'
best_performing_network$rgb = '#c9cba3'

best_matched_and_performing_network_series <- rbind(best_matched_network[c('size','wd','performance','type','rgb')],
                                                    best_performing_network[c('size','wd','performance','type','rgb')])

best_matched_and_performing_network_series <- best_performing_network[c('size','wd','performance','type','rgb')]
## Fig. 3.a The Dissimilarity  ----
top_matched_network <- read.csv('./_Results/Data/best_matching_networks/best_matched_network_metrics_0.05.csv',header = TRUE) %>% 
  group_by(size) %>% summarise(max_distance = max(wd,na.rm = TRUE),
                               min_distance = min(wd,na.rm = TRUE),
                               mean_distance= mean(wd,na.rm = TRUE),
                               sd_distance= sd(wd,na.rm = TRUE))
top_matched_network[which(top_matched_network$size ==27),2:5] =0

consistent_network <- read.csv('./_Results/Data/Network_series/consistent_network_0.05.csv') %>% 
  group_by(size) %>% summarise(max_distance = max(wd,na.rm = TRUE),
                               min_distance = min(wd,na.rm = TRUE),
                               mean_distance= mean(wd,na.rm = TRUE),
                               sd_distance= sd(wd,na.rm = TRUE))


top_matched_network$type <- 'top'
top_matched_network$alpha <- 0.4
top_matched_network$color <- '#684e94'

consistent_network$type <- 'zconsistent'
consistent_network$alpha <- 0.7
consistent_network$color <- '#392D69'

distance_df <- rbind(top_matched_network,consistent_network)
distance_df$diff <- distance_df$max_distance - distance_df$min_distance
distance_df %>% group_by(type) %>% summarise(dii_all = sum(diff))

p0 <- ggplot(distance_df)+
  geom_line(aes(x=size, y=min_distance,color = color)) +
  geom_line(aes(x=size, y=max_distance, color = color)) +
  geom_ribbon(aes(x=size,ymin=min_distance,ymax=max_distance,alpha = alpha,fill = type))+
  
  geom_line(data=best_matched_and_performing_network_series,aes(x=size, y=wd,color=rgb),size = 0.6,linetype=1) +
  geom_point(data=best_matched_and_performing_network_series,aes(x=size, y=wd,color=rgb),size = 0.8) +
  
  
  scale_color_identity()+#scale_fill_identity()+#scale_y_log10()+
  scale_x_continuous(breaks = seq(5, 27,2)) + # scale_y_log10()+
  scale_y_continuous(breaks = seq(0, 0.6,0.02)) +
  scale_alpha_identity()+
  scale_fill_manual(values  = c('#684e94','#392D69'))+
  labs(# title = "The distance to full trait network (percentage)",
    x = "Number of traits",
    y = "Network dissimilarity") + theme_classic() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")
p0

## Fig. 3.b The Performance ----
performance_top_matched <- read.csv('./_Results/Data/best_matching_networks/best_matched_network_performance_0.05.csv') %>% 
  mutate(performance = sqrt(sqrt(density_R2)*sqrt(density_weighted_R2))*sqrt(modularity_optimal_R2)) %>% 
  group_by(size) %>% summarise(mean_performance = mean(performance,na.rm = TRUE),
                               max_performance = max(performance,na.rm = TRUE),
                               min_performance =min(performance,na.rm = TRUE))

performance_top_matched$type = 'top'
performance_top_matched$color = '#684e94'
performance_top_matched$alpha = 0.4

consistent_network_performance <- read.csv('./_Results/Data/Network_series/consistent_network_0.05.csv')  %>% 
  mutate(performance = sqrt(sqrt(density_R2)*sqrt(density_weighted_R2))*sqrt(modularity_optimal_R2)) %>% 
  group_by(size) %>% summarise(mean_performance = mean(performance,na.rm = TRUE),
                               max_performance = max(performance,na.rm = TRUE),
                               min_performance =min(performance,na.rm = TRUE))


consistent_network_performance$type = 'zconsistent'
consistent_network_performance$color = '#392D69'
consistent_network_performance$alpha = 0.7

mean(consistent_network_performance$mean_performance)
sd(consistent_network_performance$mean_performance)
mean(consistent_network_performance$min_performance)
sd(consistent_network_performance$min_performance)


sum(consistent_network_performance$mean_performance)/sum(performance_top_matched$mean_performance)

performance_df <- rbind(performance_top_matched,consistent_network_performance)
# best_matched_network <- read.csv('./_Results/Data/Network_series/best_matched_network_series_0.05.csv')

p1<-ggplot(performance_df)+
  geom_line(aes(x=size, y=min_performance,color = color)) +
  geom_line(aes(x=size, y=max_performance, color = color)) +
  geom_line(aes(x=size, y=mean_performance, color = color),size=0.8) +
  geom_ribbon(aes(x=size,ymin=min_performance,ymax=max_performance,
                  alpha = alpha, fill = type))+
  
  geom_line(data=best_matched_and_performing_network_series,aes(x=size, y=performance,color=rgb),size = 0.6,linetype=1) +
  geom_point(data=best_matched_and_performing_network_series,aes(x=size, y=performance,color=rgb),size = 0.8) +
  
  scale_color_identity()+#scale_fill_identity()+#scale_y_log10()+
  scale_x_continuous(breaks = seq(5, 27,2)) + # scale_y_log10()+
  scale_y_continuous(breaks = seq(0, 1,0.2)) +
  scale_alpha_identity()+
  scale_fill_manual(values  = c('#684e94','#392D69'))+
  labs(# title = "Modularity",
    x = "Number of traits",
    y = "Network performance") + theme_classic() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

ggarrange(p0,p1,ncol = 2,nrow = 1,labels =c('a','b'))
ggsave('./_Results/Figures/Fig. 3 The dissimilarity and network performance of the consistent reduced networks and the best performing network series.jpg',
       width = 22, height = 9, units = "cm",dpi = 600)
# ggsave('./_Results/Figures/Fig. 3 The dissimilarity and network performance of the consistent reduced networks and the best performing network series.pdf',
#        width = 22, height = 9, units = "cm",dpi = 600)


### Supplement: Performance detail ----
performance_top_matched <- read.csv('./_Results/Data/best_matching_networks/best_matched_network_performance_0.05.csv') %>% 
  group_by(size) %>% summarise(mean_density_R2 = mean(density_R2,na.rm = TRUE),max_density_R2 = max(density_R2,na.rm = TRUE),min_density_R2 =min(density_R2,na.rm = TRUE),
                               mean_density_weighted_R2 = mean(density_weighted_R2,na.rm = TRUE),max_density_weighted_R2 = max(density_weighted_R2,na.rm = TRUE),min_density_weighted_R2 =min(density_weighted_R2,na.rm = TRUE),
                               mean_modularity_optimal_R2 = mean(modularity_optimal_R2,na.rm = TRUE),max_modularity_optimal_R2 = max(modularity_optimal_R2,na.rm = TRUE),min_modularity_optimal_R2 =min(modularity_optimal_R2,na.rm = TRUE))
performance_top_matched$type = 'top'
performance_top_matched$color = '#684e94'
performance_top_matched$alpha = 0.5

consistent_network_performance <- read.csv('./_Results/Data/Network_series/consistent_network_0.05.csv')  %>% 
  group_by(size) %>% summarise(mean_density_R2 = mean(density_R2,na.rm = TRUE),max_density_R2 = max(density_R2,na.rm = TRUE),min_density_R2 =min(density_R2,na.rm = TRUE),
                               mean_density_weighted_R2 = mean(density_weighted_R2,na.rm = TRUE),max_density_weighted_R2 = max(density_weighted_R2,na.rm = TRUE),min_density_weighted_R2 =min(density_weighted_R2,na.rm = TRUE),
                               mean_modularity_optimal_R2 = mean(modularity_optimal_R2,na.rm = TRUE),max_modularity_optimal_R2 = max(modularity_optimal_R2,na.rm = TRUE),min_modularity_optimal_R2 =min(modularity_optimal_R2,na.rm = TRUE))
consistent_network_performance$type = 'zconsistent'
consistent_network_performance$color = '#392D69'
consistent_network_performance$alpha = 0.7

# mean_density_R2
mean(consistent_network_performance$mean_density_R2)
sd(consistent_network_performance$mean_density_R2)
mean(consistent_network_performance$min_density_R2)
sd(consistent_network_performance$min_density_R2)


# density_weighted_R2
mean(consistent_network_performance$mean_density_weighted_R2)
sd(consistent_network_performance$mean_density_weighted_R2)
mean(consistent_network_performance$min_density_weighted_R2)
sd(consistent_network_performance$min_density_weighted_R2)


# modularity_optimal_R2
mean(consistent_network_performance$mean_modularity_optimal_R2)
sd(consistent_network_performance$mean_modularity_optimal_R2)
mean(consistent_network_performance$min_modularity_optimal_R2)
sd(consistent_network_performance$min_modularity_optimal_R2)


performance_df <- rbind(performance_top_matched,consistent_network_performance)
p1<-ggplot(performance_df)+
  geom_line(aes(x=size, y=min_density_R2,color = color)) +
  geom_line(aes(x=size, y=max_density_R2, color = color)) +
  geom_line(aes(x=size, y=mean_density_R2, color = color),size=0.8) +
  geom_ribbon(aes(x=size,ymin=min_density_R2,ymax=max_density_R2,
                  alpha = alpha, fill = type))+
  geom_line(data=best_performing_network,aes(x=size, y=density_R2),color='#c9cba3',size = 0.8)+
  
  scale_color_identity()+#scale_fill_identity()+#scale_y_log10()+
  scale_x_continuous(breaks = seq(5, 27,2)) + # scale_y_log10()+
  scale_y_continuous(breaks = seq(0, 1,0.1)) +
  scale_alpha_identity()+
  scale_fill_manual(values  = c('#684e94','#392D69'))+
  labs(# title = "Density",
    x = "",
    y = TeX("$R^2$ on Edge Density")) + theme_classic() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")# ,axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()

p2<-ggplot(performance_df)+
  geom_line(aes(x=size, y=min_density_weighted_R2,color = color)) +
  geom_line(aes(x=size, y=max_density_weighted_R2, color = color)) +
  geom_line(aes(x=size, y=mean_density_weighted_R2, color = color),size=0.8) +
  geom_ribbon(aes(x=size,ymin=min_density_weighted_R2,ymax=max_density_weighted_R2,
                  alpha = alpha, fill = type))+
  geom_line(data=best_performing_network,aes(x=size, y=density_weighted_R2),color='#c9cba3',size = 0.8)+
  
  scale_color_identity()+#scale_fill_identity()+#scale_y_log10()+
  scale_x_continuous(breaks = seq(5, 27,2)) + # scale_y_log10()+
  scale_y_continuous(breaks = seq(0, 1,0.1)) +
  scale_alpha_identity()+
  scale_fill_manual(values  = c('#684e94','#392D69'))+
  labs(# title = "Connectivity",
    x = "Number of traits",
    y = TeX("$R^2$ on Connectivity")) + theme_classic() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")


p3<-ggplot(performance_df)+
  geom_line(aes(x=size, y=min_modularity_optimal_R2,color = color)) +
  geom_line(aes(x=size, y=max_modularity_optimal_R2, color = color)) +
  geom_line(aes(x=size, y=mean_modularity_optimal_R2, color = color),size=0.8) +
  geom_ribbon(aes(x=size,ymin=min_modularity_optimal_R2,ymax=max_modularity_optimal_R2,
                  alpha = alpha, fill = type))+
  geom_line(data=best_performing_network,aes(x=size, y=modularity_optimal_R2),color='#c9cba3',size = 0.8)+
  
  scale_color_identity()+#scale_fill_identity()+#scale_y_log10()+
  scale_x_continuous(breaks = seq(5, 27,2)) + # scale_y_log10()+
  scale_y_continuous(breaks = seq(0, 1,0.1)) +
  scale_alpha_identity()+
  scale_fill_manual(values  = c('#684e94','#392D69'))+
  labs(# title = "Modularity",
    x = "Number of traits",
    y = TeX("$R^2$ on Modularity")) + theme_classic() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

ggarrange(p1,p2,p3,ncol = 3,nrow = 1,labels =c('a','b','c'))

ggsave('./_Results/Figures/Supplement Fig. 3 Network performance of the consistent reduced networks and the optimal  network series.jpg',
       width = 20, height = 6, units = "cm",dpi = 600)
# ggsave('./_Results/Figures/Supplement Fig. 3 Network performance of the consistent reduced networks and the optimal  network series.pdf',
#        width = 20, height = 6, units = "cm",dpi = 600)
#*********************************************------------------
