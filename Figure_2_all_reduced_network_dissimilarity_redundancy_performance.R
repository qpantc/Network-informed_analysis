rm(list = ls())
library(latex2exp)
source('./0_Loading_variables.R')

# 1 Reduced network and the redundancy of traits in theory ------
## Fig. 2.a The Dissimilarity ------

# Data: all potential network | random sampled network  | top matched network 
all_potential_network <- read.csv('./_Results/Tables/BMN_dw.csv',header = TRUE) %>% 
  dplyr::select(c('size','max_distance','min_distance','mean_distance','sd_distance'))
all_potential_network$type <- 'All reduced networks'
all_potential_network$alpha <- 0.1
all_potential_network$color <- '#000000'

ramdon_sampled_network <- read.csv('./_Results/Data/Random_sampled_network/Random_network_sample.csv',header = TRUE) %>% 
  group_by(size) %>% summarise(max_distance = max(wd,na.rm = TRUE),
                               min_distance = min(wd,na.rm = TRUE),
                               mean_distance= mean(wd,na.rm = TRUE),
                               sd_distance= sd(wd,na.rm = TRUE)) 
ramdon_sampled_network[which(ramdon_sampled_network$size ==27),2:5] =0

ramdon_sampled_network$type <- 'Randomly sampled networks'
ramdon_sampled_network$alpha <- 0.8
ramdon_sampled_network$color <- '#666666'

top_matched_network <- read.csv('./_Results/Data/best_matching_networks/best_matching_networks_number_0.05.csv',header = TRUE) %>% 
  group_by(size) %>% summarise(max_distance = max(wd,na.rm = TRUE),
                               min_distance = min(wd,na.rm = TRUE),
                               mean_distance= mean(wd,na.rm = TRUE),
                               sd_distance= sd(wd,na.rm = TRUE))
top_matched_network[which(top_matched_network$size ==27),2:5] =0

top_matched_network$type <- 'Optimal reduced networks'
top_matched_network$alpha <- 0.8
top_matched_network$color <- '#684e94'

distance_df <- rbind(all_potential_network,ramdon_sampled_network,top_matched_network)


p1 <- ggplot(distance_df)+
  geom_ribbon(aes(x=size,ymin=min_distance,ymax=max_distance,
                  alpha = alpha*0.4, fill = type))+
  geom_line(aes(x=size, y=min_distance,color = color,alpha = 1-alpha)) +
  geom_line(aes(x=size, y=max_distance, color = color)) +
  geom_line(aes(x=size, y=mean_distance, alpha=ifelse(type == 'All reduced networks',0,1),color = color),size=0.6) +
  
  geom_point(aes(x=size, y=mean_distance, color = color,shape = type,
                 alpha=ifelse(type == 'All reduced networks',0,1)),size=1.5) +
  scale_color_identity()+#scale_fill_identity()+#scale_y_log10()+
  scale_x_continuous(breaks = seq(5, 27,1)) + # scale_y_log10()+
  scale_y_continuous(breaks = seq(0, 0.6,0.1)) +
  scale_alpha_identity()+
  scale_fill_manual(values  = c('#ffffff','#aaaaaa','#684e94'))+
  scale_shape_manual(values = c(0,15,16))+
  labs(# title = "The distance to full trait network (percentage)",
    x = "Number of traits",
    y = "Network dissimilarity") + theme_classic() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = c(.8, .85))#
p1

## Fig. 2.b Trait redundancy------
#### redundancy for sampled network traits
random_metrics <- read.csv('./_Results/Data/Random_sampled_network/metrics_of_random_1781855.csv')
random_metrics[c(7:33)] <- random_metrics[c(7:33)]/10000

random_metrics_summary <- random_metrics[c(3,7:33)] %>% 
  pivot_longer(cols = LA:SD,names_to = 'Trait',values_to = 'Value') %>% 
  mutate(Value =1- Value) %>%   mutate(Value =ifelse(Value<0,0,Value)) %>% 
  # mutate(Value = (1+exp(-1))/(exp(6*Value -1)+1)) %>% #Value = 1/(27*Value+1) (1+exp(-1))/(exp(27*Value -1)+1)
  group_by(Trait,size) %>% summarise(value_mean = mean(Value,na.rm = TRUE),
                                     value_sd = sd(Value,na.rm = TRUE))
random_metrics_summary$type = 'random'
random_metrics_summary$color = '#999999'

colnames(random_metrics_summary)  <- c("Nodes","size","value_mean", "value_sd","type","color")
random_metrics_summary <- inner_join(random_metrics_summary,Vertices)
random_metrics_summary$alpha = 1-random_metrics_summary$value_mean

#### redundancy for best matched network traits
best_matched_metrics <- read.csv('./_Results/Data/best_matching_networks/best_matched_network_metrics_0.05.csv')
best_matched_metrics[c(7:33)] <- best_matched_metrics[c(7:33)]/10000
best_matched_metrics_summary <- best_matched_metrics[c(3,7:33)] %>% 
  pivot_longer(cols = LA:SD,names_to = 'Trait',values_to = 'Value') %>% 
  mutate(Value =1- Value) %>%   mutate(Value =ifelse(Value<0,0,Value)) %>% 
  group_by(Trait,size) %>% summarise(value_mean = mean(Value,na.rm = TRUE),
                                     value_sd = sd(Value,na.rm = TRUE))
best_matched_metrics_summary$type = 'matched'
best_matched_metrics_summary$color = '#009E73'

colnames(best_matched_metrics_summary)  <- c("Nodes","size","value_mean", "value_sd","type","color")
best_matched_metrics_summary <- inner_join(best_matched_metrics_summary,Vertices)
best_matched_metrics_summary$alpha = 1-best_matched_metrics_summary$value_mean

trait_redundancy_df <- rbind(best_matched_metrics_summary,random_metrics_summary)
trait_redundancy_summary <- trait_redundancy_df %>% group_by(size,type)%>% 
  summarise(value_max = max(value_mean,na.rm = TRUE),
            value_m = mean(value_mean,na.rm = TRUE),
            value_min = min(value_mean,na.rm = TRUE))

trait_redundancy_summary$color <- 'a'
trait_redundancy_summary$color[which(trait_redundancy_summary$type == 'matched')] <- 'b'

trait_redundancy_summary$breakpoint[which(trait_redundancy_summary$size >10 )] <- 'large'  
trait_redundancy_summary$breakpoint[which(trait_redundancy_summary$size <=10 )] <- 'small'  

sum(trait_redundancy_summary$value_m[which(trait_redundancy_summary$size <=10 & trait_redundancy_summary$type == 'matched')])/
  sum(trait_redundancy_summary$value_m[which(trait_redundancy_summary$size <=10 & trait_redundancy_summary$type == 'random')])


p2 <- ggplot(trait_redundancy_summary)+
  geom_ribbon(data = trait_redundancy_summary %>% filter(breakpoint =='small'),aes(x=size,ymin=value_min,ymax=value_max,
                                                                                   fill = color),alpha = 0.4)+
  geom_line(data = trait_redundancy_summary %>% filter(breakpoint =='small'),aes(x=size, y=value_min, color = color)) +
  geom_line(data = trait_redundancy_summary %>% filter(breakpoint =='small'),aes(x=size, y=value_max, color = color)) +
  
  
  geom_ribbon(data = trait_redundancy_summary %>% filter(breakpoint =='large'),aes(x=size,ymin=value_min,ymax=value_max,
                                                                                   fill = color),alpha = 0.4)+
  geom_line(data = trait_redundancy_summary %>% filter(breakpoint =='large'),aes(x=size, y=value_min, color = color)) +
  geom_line(data = trait_redundancy_summary %>% filter(breakpoint =='large'),aes(x=size, y=value_max, color = color)) +
  
  
  geom_point(data = trait_redundancy_summary %>% filter(breakpoint =='small'),
             aes(x=size,y=value_m, color = color,shape = type),size=1.5) +
  geom_smooth(data = trait_redundancy_summary %>% filter(breakpoint =='small'),
              aes(x=size,y=value_m, color = color,fill = color),method = 'loess',se = FALSE,size=0.6) +#,method = 'lm'
  
  
  geom_point(data = trait_redundancy_summary %>% filter(breakpoint =='large'),
             aes(x=size,y=value_m, color = color,shape = type),size=1.5)+
  geom_smooth(data = trait_redundancy_summary %>% filter(breakpoint =='large'),
              aes(x=size,y=value_m, color = color,fill = color),size=0.6) +
  
  scale_x_continuous(breaks = seq(5, 27,2)) + # scale_y_log10()+
  scale_y_continuous( limits = c(0.4,1),breaks = seq(0, 1,0.1)) +
  scale_color_manual(values =c('#666666','#684e94'))+
  scale_fill_manual(values  = c('#aaaaaa','#684e94'))+
  scale_alpha_identity()+
  scale_shape_manual(values = c(16,15))+
  labs(
    x = "Number of traits",y = "Trait redundancy") +
  theme_classic()+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = 'None')
p2

## Fig. 2.3 The performance ----
performance_random <- read.csv('./_Results/Data/Random_sampled_network/Performance_of_random_1781855.csv') %>% 
  mutate(performance = sqrt(sqrt(density_R2)*sqrt(density_weighted_R2))*sqrt(modularity_optimal_R2)) %>% 
  group_by(size) %>% summarise(mean_performance = mean(performance,na.rm = TRUE),
                               max_performance = max(performance,na.rm = TRUE),
                               min_performance =min(performance,na.rm = TRUE))
performance_random$type = 'random'
performance_random$color = '#666666'
performance_random$alpha = 0.8

performance_top_matched <- read.csv('./_Results/Data/best_matching_networks/best_matched_network_performance_0.05.csv')
performance_top_matched <- performance_top_matched  %>% 
  mutate(performance = sqrt(sqrt(density_R2)*sqrt(density_weighted_R2))*sqrt(modularity_optimal_R2)) %>% 
  group_by(size) %>% summarise(mean_performance = mean(performance,na.rm = TRUE),
                               max_performance = max(performance,na.rm = TRUE),
                               min_performance =min(performance,na.rm = TRUE))
performance_top_matched$type = 'top'
performance_top_matched$color = '#684e94'
performance_top_matched$alpha = 0.8

mean(performance_top_matched$mean_performance)
sd(performance_top_matched$mean_performance)
mean(performance_top_matched$min_performance)
sd(performance_top_matched$min_performance)

mean(performance_random$mean_performance)
sd(performance_random$mean_performance)
mean(performance_random$min_performance)
sd(performance_random$min_performance)


sum(performance_top_matched$mean_performance)/sum(performance_random$mean_performance) - 1 
sum(performance_top_matched$min_performance)/sum(performance_random$min_performance) - 1


performance_df <- rbind(performance_random,performance_top_matched)
p3<-ggplot(performance_df)+
  geom_ribbon(aes(x=size,ymin=min_performance,ymax=max_performance,
                  alpha = alpha*0.4, fill = type))+
  geom_line(aes(x=size, y=min_performance, color = color)) +
  geom_line(aes(x=size, y=max_performance, color = color)) +
  
  geom_line(aes(x=size, y=mean_performance, color = color),size=0.6) +
  geom_point(aes(x=size, y=mean_performance, color = color,shape = type),size=1.5) +
  
  scale_color_identity()+#scale_fill_identity()+#scale_y_log10()+
  scale_x_continuous(breaks = seq(5, 27,2)) + # scale_y_log10()+
  scale_y_continuous(breaks = seq(0, 1,0.2)) +
  scale_alpha_identity()+
  scale_fill_manual(values  = c('#999999','#684e94'))+
  scale_shape_manual(values = c(15,16))+
  labs(# title = "Density",
    x = "Number of traits",
    y = "Network performance") + theme_classic() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

# axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()

P2 <-  ggarrange(p2,p3,ncol = 2,labels = c('b','c'))
P1P2 <- ggarrange(p1,P2,ncol = 1,labels = c('a','',''))
P1P2

ggsave('./_Results/Figures/Fig. 2 The dissimilarity between optimal reduced networks of different sizes to the full network, and the trait redundancy and performance of optimal reduced networks.jpg',
       width = 20 , height = 15, units = "cm",dpi = 600)
# ggsave('./_Results/Figures/Fig. 2 The dissimilarity between optimal reduced networks of different sizes to the full network, and the trait redundancy and performance of optimal reduced networks.pdf',
#        width = 20, height = 15, units = "cm",dpi = 600)

## Supplement: Trait redundancy bar plot ----
### random sampled
random_metrics <- read.csv('./_Results/Data/Random_sampled_network/metrics_of_random_1781855.csv')
random_metrics[c(7:33)] <- random_metrics[c(7:33)]/10000
best_matched_metrics <- read.csv('./_Results/Data/best_matching_networks/best_matched_network_metrics_0.05.csv')
best_matched_metrics[c(7:33)] <- best_matched_metrics[c(7:33)]/10000
random_metrics_summary <- random_metrics[c(3,7:33)] %>% 
  pivot_longer(cols = LA:SD,names_to = 'Trait',values_to = 'Value') %>% 
  mutate(Value =1- Value) %>%   mutate(Value =ifelse(Value<0,0,Value)) %>% 
  group_by(Trait) %>% summarise(value_mean = mean(Value,na.rm = TRUE),
                                value_sd = sd(Value,na.rm = TRUE))

random_metrics_summary$type = 'matched'
random_metrics_summary$color = '#684e94' #009E73 #CC79A7
colnames(random_metrics_summary)  <- c("Nodes","value_mean", "value_sd","type","color")
random_metrics_summary <- inner_join(random_metrics_summary,Vertices)
random_metrics_summary$Nodes <- factor(random_metrics_summary$Nodes,
                                       levels=random_metrics_summary$Nodes[order(random_metrics_summary$value_mean,decreasing = FALSE)])
random_metrics_summary$value_mean <- random_metrics_summary$value_mean -0.7

p1 <- ggplot(random_metrics_summary,
             aes(x=Nodes,y=value_mean,fill = Organ_color))+
  geom_bar(stat = 'identity',alpha = 0.7,color = '#000000',width = 0.6)+
  geom_errorbar(aes(ymin=value_mean-value_sd, ymax=value_mean+value_sd,color = Organ_color),
                width=.5)+
  scale_color_identity() +
  scale_fill_identity() +
  
  labs(# title = "Modularity",
    x = "Number of traits",
    y = "Trait redundancy") + theme_bw() + 
  scale_y_continuous(breaks = seq(0, 1,0.1),labels = seq(0.7,1.7,0.1)) +
  
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))


### best matched 

best_matched_metrics_summary <- best_matched_metrics[c(3,7:33)] %>% 
  pivot_longer(cols = LA:SD,names_to = 'Traits',values_to = 'Value') %>% 
  mutate(Value =1- Value) %>%   mutate(Value =ifelse(Value<0,0,Value)) %>% 
  group_by(Traits) %>% summarise(value_mean = mean(Value,na.rm = TRUE),
                                 value_sd = sd(Value,na.rm = TRUE))

best_matched_metrics_summary$type = 'matched'
best_matched_metrics_summary$color = '#684e94' #009E73 #CC79A7
colnames(best_matched_metrics_summary)  <- c("Nodes","value_mean", "value_sd","type","color")
best_matched_metrics_summary <- inner_join(best_matched_metrics_summary,Vertices)
best_matched_metrics_summary$Nodes <- factor(best_matched_metrics_summary$Nodes,
                                             levels=best_matched_metrics_summary$Nodes[order(best_matched_metrics_summary$value_mean,decreasing = FALSE)])
best_matched_metrics_summary$value_mean <- best_matched_metrics_summary$value_mean -0.7

p2 <- ggplot(best_matched_metrics_summary,
             aes(x=Nodes,y=value_mean,fill = Organ_color))+
  geom_bar(stat = 'identity',alpha = 0.7,color = '#000000',width = 0.6)+
  geom_errorbar(aes(ymin=value_mean-value_sd, ymax=value_mean+value_sd,color = Organ_color),
                width=.5)+
  scale_color_identity() +
  scale_fill_identity() +
  
  labs(# title = "Modularity",
    x = "Number of traits",
    y = "Trait redundancy") + theme_bw() + 
  scale_y_continuous(breaks = seq(0, 1,0.1),labels = seq(0.7,1.7,0.1)) +
  
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggarrange(p1,p2,ncol = 1,labels = c('a','b'))
ggsave('./_Results/Figures/Supplement Fig. 1 The averaged redundancy of traits in randomly selected networks and optimal reduced networks.jpg',
       width = 22 , height = 12, units = "cm",dpi = 600)

## Supplement: performance detail -----
performance_random <- read.csv('./_Results/Data/Random_sampled_network/Performance_of_random_1781855.csv') %>% 
  group_by(size) %>% summarise(mean_density_R2 = mean(density_R2,na.rm = TRUE),max_density_R2 = max(density_R2,na.rm = TRUE),min_density_R2 =min(density_R2,na.rm = TRUE),
                               mean_density_weighted_R2 = mean(density_weighted_R2,na.rm = TRUE),max_density_weighted_R2 = max(density_weighted_R2,na.rm = TRUE),min_density_weighted_R2 =min(density_weighted_R2,na.rm = TRUE),
                               mean_modularity_optimal_R2 = mean(modularity_optimal_R2,na.rm = TRUE),max_modularity_optimal_R2 = max(modularity_optimal_R2,na.rm = TRUE),min_modularity_optimal_R2 =min(modularity_optimal_R2,na.rm = TRUE))
performance_random$type = 'random'
performance_random$color = '#666666'
performance_random$alpha = 0.2


performance_top_matched <- read.csv('./_Results/Data/best_matching_networks/best_matched_network_performance_0.05.csv')
performance_top_matched <- performance_top_matched  %>% group_by(size) %>% 
  summarise(mean_density_R2 = mean(density_R2,na.rm = TRUE),max_density_R2 = max(density_R2,na.rm = TRUE),min_density_R2 =min(density_R2,na.rm = TRUE),
            mean_density_weighted_R2 = mean(density_weighted_R2,na.rm = TRUE),max_density_weighted_R2 = max(density_weighted_R2,na.rm = TRUE),min_density_weighted_R2 =min(density_weighted_R2,na.rm = TRUE),
            mean_modularity_optimal_R2 = mean(modularity_optimal_R2,na.rm = TRUE),max_modularity_optimal_R2 = max(modularity_optimal_R2,na.rm = TRUE),min_modularity_optimal_R2 =min(modularity_optimal_R2,na.rm = TRUE))

performance_top_matched$type = 'top'
performance_top_matched$color = '#684e94'
performance_top_matched$alpha = 0.6

summary(performance_random)

## density_R2 
mean(performance_random$mean_density_R2)
sd(performance_random$mean_density_R2)
mean(performance_random$min_density_R2)
sd(performance_random$min_density_R2)

mean(performance_top_matched$mean_density_R2)
sd(performance_top_matched$mean_density_R2)
mean(performance_top_matched$min_density_R2)
sd(performance_top_matched$min_density_R2)

# density_weighted_R2
mean(performance_random$mean_density_weighted_R2)
sd(performance_random$mean_density_weighted_R2)
mean(performance_random$min_density_weighted_R2)
sd(performance_random$min_density_weighted_R2)

mean(performance_top_matched$mean_density_weighted_R2)
sd(performance_top_matched$mean_density_weighted_R2)
mean(performance_top_matched$min_density_weighted_R2)
sd(performance_top_matched$min_density_weighted_R2)


# modularity_optimal_R2
mean(performance_random$mean_modularity_optimal_R2)
sd(performance_random$mean_modularity_optimal_R2)
mean(performance_random$min_modularity_optimal_R2)
sd(performance_random$min_modularity_optimal_R2)

mean(performance_top_matched$mean_modularity_optimal_R2)
sd(performance_top_matched$mean_modularity_optimal_R2)
mean(performance_top_matched$min_modularity_optimal_R2)
sd(performance_top_matched$min_modularity_optimal_R2)



performance_df <- rbind(performance_random,performance_top_matched)
p1<-ggplot(performance_df)+
  geom_line(aes(x=size, y=min_density_R2,color = color)) +
  geom_line(aes(x=size, y=max_density_R2, color = color)) +
  geom_line(aes(x=size, y=mean_density_R2, color = color),size=0.8) +
  geom_ribbon(aes(x=size,ymin=min_density_R2,ymax=max_density_R2,
                  alpha = alpha, fill = type))+
  scale_color_identity()+#scale_fill_identity()+#scale_y_log10()+
  scale_x_continuous(breaks = seq(5, 27,2)) + # scale_y_log10()+
  scale_y_continuous(breaks = seq(0, 1,0.1)) +
  scale_alpha_identity()+
  scale_fill_manual(values  = c('#666666','#684e94'))+
  labs(# title = "Density",
    x = "Number of traits",
    y = TeX("$R^2$ on Edge Density")) + theme_classic() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

# axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()

p2<-ggplot(performance_df)+
  geom_line(aes(x=size, y=min_density_weighted_R2,color = color)) +
  geom_line(aes(x=size, y=max_density_weighted_R2, color = color)) +
  geom_line(aes(x=size, y=mean_density_weighted_R2, color = color),size=0.8) +
  geom_ribbon(aes(x=size,ymin=min_density_weighted_R2,ymax=max_density_weighted_R2,
                  alpha = alpha, fill = type))+
  scale_color_identity()+#scale_fill_identity()+#scale_y_log10()+
  scale_x_continuous(breaks = seq(5, 27,2)) + # scale_y_log10()+
  scale_y_continuous(breaks = seq(0, 1,0.1)) +
  scale_alpha_identity()+
  scale_fill_manual(values  = c('#666666','#684e94'))+
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
  scale_color_identity()+#scale_fill_identity()+#scale_y_log10()+
  scale_x_continuous(breaks = seq(5, 27,2)) + # scale_y_log10()+
  scale_y_continuous(breaks = seq(0, 1,0.1)) +
  scale_alpha_identity()+
  scale_fill_manual(values  = c('#666666','#684e94'))+
  labs(# title = "Modularity",
    x = "Number of traits",
    y = TeX("$R^2$ on Modularity")) + theme_classic() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

ggarrange(p1,p2,p3,ncol = 3,labels = c('a','b','c'))

ggsave('./_Results/Figures/Supplement Fig. 2 Performance of network along size.jpg',
       width = 20 , height = 6, units = "cm",dpi = 600)
# ggsave('./_Results/Figures/Supplement Fig. 2 Performance of network along size.pdf',
#        width = 20 , height = 6, units = "cm",dpi = 600)


#*********************************************------------------
