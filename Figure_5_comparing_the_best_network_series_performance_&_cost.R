rm(list = ls())
library(latex2exp)
source('./0_Loading_variables.R')

# 4  marginal benefit of best performed network and cheapest network----
## Cost effective network series ------------------
best_performing_network_series_marginal_benefit <- read.csv('./_Results/Data/Network_series/best_performing_network_series_0.1_marginal_benefit_0.05.csv')
best_performing_start_number = best_performing_network_series_marginal_benefit$number[1]
best_performing_network_series_marginal_benefit$type = 'Best'

best_worthy_networks_series_marginal_benefit <- read.csv('./_Results/Data/Network_series/best_worthy_network_marginal_benefit_0.05.csv')
best_worthy_start_number = best_worthy_networks_series_marginal_benefit$number[1]
best_worthy_networks_series_marginal_benefit$type = 'Worty'



1-sum(best_worthy_networks_series_marginal_benefit$cost)/sum(best_performing_network_series_marginal_benefit$cost)
comparing_df <- rbind(best_worthy_networks_series_marginal_benefit,best_performing_network_series_marginal_benefit)


write.csv(comparing_df,'./_Results/Data/Network_series/comparing_df_cost_effective.csv',row.names = FALSE)
ggplot(comparing_df, aes(size, weight = - cost, fill = type)) +
  geom_hline(yintercept = seq(-1, 1, 0.2), color = '#cccccc',size = 0.1) +
  geom_bar(color = "#222222", width = .6, position = 'dodge') +
  
  geom_line(aes(x = size , y = performance, color = type))+
  geom_point(aes(x = size, y = performance, color = type),
             shape = 15)+
  
  # geom_line(aes(x = ifelse(type == 'Best',size-0.2,size+0.2) , y = performance, color = type))+
  # geom_point(aes(x = ifelse(type == 'Best',size-0.2,size+0.2), y = performance, color = type),
  #            shape = 15)+
  
  # geom_text(aes(x=ifelse(type == 'Best',size-0.25,size+0.25),
  #               y=ifelse(type == 'Best',performance+0.03,performance-0.03),label = round(performance,2),color = type), size = 2.5) +
  # geom_text(aes(x=size,y=0-0.02,label = round(cost,2),color = type), size = 2.5) +
  # geom_text(aes(label = newtrait,x=size,y=cost-0.02,color = type), size = 2.5) +
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = -1)+
  
  geom_text(aes(label = newtrait,x=size-0.3,hjust = 0,
                y=ifelse(type != 'Best',0.03,0.1),color = type),size = 2.5) +
  
  geom_text(aes(label = ifelse(type == 'Best',
                               paste0(paste0(PowerSetsBinary(best_performing_start_number),collapse = '_')),
                               paste0(paste0(PowerSetsBinary(best_worthy_start_number),collapse = '_'))),
                x=4.7,hjust = 0,
                y=ifelse(type != 'Best',0.85,0.95),color = type),size = 2.5) +
  geom_text(aes(x=size,y=-1.04,label = size),size = 2.5) +
  
  scale_fill_manual(values = c('#c9cba3','#990000'))+
  scale_color_manual(values = c('#c9cba3','#990000'))+
  
  scale_x_continuous(breaks = seq(5, 27,1)) + 
  scale_y_continuous(
    ~.*1,name = '      Network cost                              Network performance',
    breaks=seq(-1,1,0.2),labels = c(seq(1,0.2,-0.2),seq(0,1,0.2)))+ 
  labs(x = "Number of traits") + theme_classic() + 
  theme(panel.grid.minor = element_blank(),legend.position = "none",
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.line.x = element_blank())

ggsave('./_Results/Figures/Fig. 5 The cost and performance of the best performing network series and the most cost-efficient network series.jpg',
       width = 22,height =13, units = "cm",dpi = 600)
# ggsave('./_Results/Figures/Fig. 5 The cost and performance of the best performing network series and the most cost-efficient network series.pdf',
#        width = 22,height =13, units = "cm",dpi = 600)

## marginal benefit------------------
ggplot(comparing_df %>% filter(type != "Best"), aes(size, weight = - cost, fill = type)) +
  geom_line(aes(x = size, y = marginal_benefit, color = "#990000"))+
  geom_point(aes(x = size, y = marginal_benefit, color =  "#990000"),
             shape = 15)+
  scale_color_identity()+
  labs(x="Number of traits", y = "Marginal benefit")+
  scale_x_continuous(breaks = seq(5, 27,1)) + 
  theme_classic() + theme(legend.position = 'None')

ggsave('./_Results/Figures/Supplement Fig. 4 The marginal benefit of adding traits in the cost-efficient the most cost-efficient network.jpg',
       width = 22,height =13, units = "cm",dpi = 600)
# ggsave('./_Results/Figures/Supplement Fig. 4 The marginal benefit of adding traits in the cost-efficient the most cost-efficient network.pdf',
#        width = 22,height =13, units = "cm",dpi = 600)
