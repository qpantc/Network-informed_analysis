# Loading data or variables
# new_node_list <- c("Aarea","L15N", "LA","LCC.m","LCNR", "LD","LDM","LDMC","LL",
#                    "LNC.m","LPC.m","LSC", "LSLA","Lthick", "LWid","Max.H","Rarea",
#                    "RD","Rdia","RMI","RNC", "RSRL","RTD", "SD","SdM","SHV","SSD")
# new_node_list2 <- c("Aarea","L15N", "LA","LCC.m","LCNR", "LD","LDM","LDMC","LL",
#                    "LNC.m","LPC.m","LSC", "SLA","Lthick", "LWid","Height","Rarea",
#                    "RD","Rdia","RMI","RNC", "SRL","RTD", "SD","SdM","SHV","SSD")

source('./0_Functions_for_calculating_networks.R')
library(tseries)
n = 1000
x <- rnorm(n, mean=0, sd=1)
y <- rnorm(n, mean=0, sd=1)
z <- y * 1000

cor.test(x,y)
cor.test(x,z)

# Vertices #####
Vertices <- read.csv('./_data/Vertices_trait_info.csv')
Vertices_applicaiton <- subset(Vertices,Selection>1)

Vertices <- subset(Vertices,Selection<3)
Vertices$cost <- max(Vertices$TRY6_species )/Vertices$TRY6_species


## ------------------------------------------------------------------------------------------------------
## The following gap-filled plant trait data will be made available from Jens Kattge on request.
## ------------------------------------------------------------------------------------------------------


# Trait data normal examine and log transform ----------------------------------
# ecoregion_trait <- read.csv('../_data/Gap_filled_merged_by_ecoregion.csv')
# ecoregion_trait$GrowthForm[which(ecoregion_trait$GrowthForm=='herb/shrub')] = 'shrub'
# ecoregion_trait$GrowthForm[which(ecoregion_trait$GrowthForm=='herb/shrub/tree')] = 'shrub'
# ecoregion_trait$GrowthForm[which(ecoregion_trait$GrowthForm=='shrub/tree')] = 'tree'
# 
# 
# ## normal distribution test
# for (i in c(seq(2,26),32,33)) {
#   # print(colnames(ecoregion_trait)[i])
#   if (min(ecoregion_trait[,i],na.rm = T) <= 0) {
#     ecoregion_trait[,i] <- ecoregion_trait[,i] - min(ecoregion_trait[,i],na.rm = T)
#   }
#   transformed_x_nona <- ecoregion_trait[,i]
#   test<-jarque.bera.test(transformed_x_nona)
#   test_log<-jarque.bera.test(log_normalize(transformed_x_nona))
#   # print(paste(test$statistic,test_log$statistic))
#   difference<-(("teststatistik"=as.numeric(test$statistic)) - ("teststatistik_log"=as.numeric(test_log$statistic)))
#   if(difference > 0.0){	
#     ecoregion_trait[,i] <- log10(ecoregion_trait[,i]) # log transformation
#     ecoregion_trait[,i][is.infinite(ecoregion_trait[,i])] <- NA
#   } else {
#     print(colnames(ecoregion_trait)[i])
#     print(difference)
#   }
# 
#   # print(paste(mean(ecoregion_trait[,i]),sd(ecoregion_trait[,i])))
# }
# 
# ecoregion_trait <- ecoregion_trait %>% na.omit()
# 
# 
# mean(ecoregion_trait$SSD)
# # Trait data inner join with environmental data --------------------------------------
# ecoregion_env <- read.csv('../_data/environmental_factors_for_ecoregion.csv')
# colnames(ecoregion_env) <- c("ECO_ID",'Lon','Lat','AI','ET',
#                              'Bdod','Cec','Cfvo','Clay','Nitrogen','Ocd','pH','Sand','Silt','Soc',
#                              'MAT','Temp.warmest.quarter','Temp.coldest.quarter',
#                              'MAP','Precip.wettest.month','Precip.driest.month',
#                              'Precip.seasonality','Precip.wettest.Quarter','Precip.driest.Quarter',
#                              'Precip.warmest.Quarter','Precip.coldest.Quarter','Temp.range.monthly',
#                              'Isothermality','Temp.seasonality','MaxTemp.warmes.month',
#                              'MinT.coldest.month','Temp.range.annualy','Temp.wettest.Quarter',
#                              'Temp.driest.Quarter')
# 
# ecoregion_trait_with_env <- inner_join(ecoregion_trait,ecoregion_env)
# ecoregion_trait_with_env$gradient <- as.integer(as.factor(ecoregion_trait_with_env$ECO_ID))
# existed_ecoregion <- unique(ecoregion_trait_with_env[,c(30,31,37:71)])
# 
# colnames(ecoregion_trait_with_env)
# 
# # ecoregion_trait <- subset(ecoregion_trait,GrowthForm=='herb')
# # all_edges <- data_to_edges(ecoregion_trait_with_env[c(seq(2,26),32,33)])
# # write.csv(all_edges,'./tips/all_edges.csv',row.names = FALSE)
# 
# rm(list = c('ecoregion_trait','existed_ecoregion','ecoregion_env','i'))
