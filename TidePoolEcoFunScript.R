### New biogeochemistry analysis from Silbiger and Sorte 2018


### Load Libraries ####
library(tidyverse)
library(janitor)
library(ggfortify)

### Read in data from github repo ####

# Chemistry data
urlfile<-'https://raw.githubusercontent.com/njsilbiger/Biophysical_feedbacks_in_coastal_ecosystems/master/Data/CleanData/ChemDataByTimePoint.csv'
ChemData<-read.csv(urlfile)
#clean the headers
ChemData<-clean_names(ChemData)


# community composition data
#Mobile counts
mobileurl<-'https://raw.githubusercontent.com/njsilbiger/Biophysical_feedbacks_in_coastal_ecosystems/master/Data/RawData/MobileAll3.csv'
Mobile<-read.csv(mobileurl, stringsAsFactors=FALSE)
#benthic percent cover
sessileurl<-"https://raw.githubusercontent.com/njsilbiger/Biophysical_feedbacks_in_coastal_ecosystems/master/Data/RawData/SessileAll_edited.csv"
Sessile<-read.csv(sessileurl, stringsAsFactors=FALSE)

# Pull out the info for the subgroups
SessileGroups<- Sessile[1,] #pull out the first row
InvertGroups <- Mobile[1,]

#remove the first row and convert everything back to numeric
Sessile<-data.frame(Sessile[-1,])
Mobile<-data.frame(Mobile[-1,])

#convert everything back to numeric
for (i in 3: ncol(Mobile)){
  Mobile[,i]<-as.numeric(Mobile[,i])
}

for (i in 4: ncol(Sessile)){
  Sessile[, i]<-as.numeric(Sessile[,i])
}


#replace NA with 0
Sessile[is.na(Sessile)]<-0 
Mobile[is.na(Mobile)]<-0

# Make all the community data a relative percent
PercentSessile<-100*Sessile[4:ncol(Sessile)]/Sessile$Squares
#create % cover for the mobile inverts.
MobileCover<-Mobile

## Rescale the mobile inverts after measuring them to convert to percent cover
# these were calculated by outlining lots of individuals in ImageJ.  See supplemental methods from Silbiger and Sorte for details
MobileCover$Chlorostoma<-(Mobile$Chlorostoma*2.6)/100 # divide by 100 cm2 because that is how much one square is
MobileCover$Littorina <- (MobileCover$Littorina*0.28)/100 # Litorrina
MobileCover$Lottia.sp<- (MobileCover$Lottia.sp*1.84)/100 # limpets
MobileCover$Fissurella<- (MobileCover$Fissurella*1.84)/100 # limpets
MobileCover$Cyanoplax<- (MobileCover$Cyanoplax*1.19)/100 #chitons
MobileCover$Nuttalina<-(MobileCover$Nuttalina*1.19)/100 #chitons
MobileCover$Mopilia<-(MobileCover$Mopilia*1.19)/100 #chitons

MobileCover$Pagurus<- (MobileCover$Pagurus*2.6)/100 # most pagurus are in chlorostoma shells
# I did not measure Nucella, but they are slightly smaller than chlorostoma
MobileCover$Nucella<-(Mobile$Nucella*2)/100 # divide by 100 cm2 because that is how much one square is

# for the remaining groups
MobileCover[,which(InvertGroups=='Med')]<-(MobileCover[,which(InvertGroups=='Med')]*2)/100 # 2cm2
MobileCover[,which(InvertGroups=='Fish')]<-(MobileCover[,which(InvertGroups=='Fish')]*2)/100 # these aren't actually included in the analyses because there were so little
MobileCover[,which(InvertGroups=='Large')]<-(MobileCover[,which(InvertGroups=='Large')]*10)/100 # 10 cm2

#Add the Mobile cover to the sessile
PercentTotal<-cbind(PercentSessile,MobileCover[3:ncol(MobileCover)])
#normalize to the sum of the total cover (since it can be greater than 100%)
PercentTotal<- 100*PercentTotal/rowSums(PercentTotal)

## plot percent of pools by group for each site
##remove the site and square info from the groups and join them into one vector so its the same length as percent total
SessileGroups<-SessileGroups[-c(1:3)]
InvertGroups<- InvertGroups[-c(1:2)]
Groups<- cbind(SessileGroups, InvertGroups)

#divide by groups
AllAlgae<-rowSums(PercentTotal[,c(which(Groups=='Fleshy' | Groups=='Crust' | Groups=='Coralline'))]) #sum across the algae rows for total % cover
SurfGrass<- PercentTotal$Phyllospadix
FleshyAlgae<-rowSums(PercentTotal[,c(which(Groups=='Fleshy'))])
CorallineAlgae<-rowSums(PercentTotal[,c(which(Groups=='Coralline'))])
Crust<- PercentTotal$NonCorallinecCust
RockSand<- rowSums(PercentTotal[,c(which(Groups=='Rock'))])
Inverts<-rowSums(PercentTotal[,c(which(Groups=='Invert' | Groups=='Small' | Groups=='Med'| Groups=='Star'| Groups=='Large'))])#sum across the invert rows for total % cover
Fish<-rowSums(PercentTotal[,c(which(Groups=='Fish'))])
Mussel<-PercentTotal$Mytilus
OtherInverts<-Inverts-Mussel
# group them into one data frame
CoverbyGroups<-data.frame(Sessile$Site,RockSand,FleshyAlgae,CorallineAlgae,Crust,SurfGrass,Fish, Mussel, OtherInverts)
colnames(CoverbyGroups)[1]<-'site'

# add producer dominance where positive values are dominated by producers and negative by consumers
# the community composition metric
CoverbyGroups$ProducerDom<-CoverbyGroups$FleshyAlgae+CoverbyGroups$SurfGrass-(CoverbyGroups$Mussel +CoverbyGroups$OtherInverts)
# add back the poll ID
CoverbyGroups$pool_id<-Mobile$Pool
#change order of columns
CoverbyGroups<-CoverbyGroups[,c(1,11,2:10)]


### Calculate the difference beterrn pools and ocean for all timepoints of the chemistry data

# pull out all the ocean time points
Ocean<-ChemData[which(ChemData$pool_id=='Ocean'),]
# pull out the values that we are interested in and rename with _ocean
Ocean<-Ocean %>%
  select("site","time_point",'nn_umol_l', "po_umol_l", "nh4_umol_l", "p_h", "dic_umol_kg", "do_mg_l", "ta_u_e_kg")%>%
  set_names(
  names(.) %>%
  paste0("_ocean")) 
  colnames(Ocean)[1:2]<-c("site", "time_point")

# remove these from the pool dataset
PoolChem<-ChemData[-which(ChemData$pool_id=='Ocean'),]

# join the ocean data with the pool data to make it easier to subtract at time points
Chem_joined<-left_join(PoolChem, Ocean)

# subtract ocean data from pool data to see how different it is
Chem_joined$deltaNN<-with(Chem_joined,nn_umol_l-nn_umol_l_ocean)
Chem_joined$deltaNH4<-with(Chem_joined,nh4_umol_l-nh4_umol_l_ocean)
Chem_joined$deltaPO<-with(Chem_joined,po_umol_l-po_umol_l_ocean)
Chem_joined$deltapH<-with(Chem_joined,p_h-p_h_ocean)
#Chem_joined$deltaDIC<-with(Chem_joined,dic_umol_kg-dic_umol_kg_ocean)
#Chem_joined$deltaTA<-with(Chem_joined,ta_u_e_kg-ta_u_e_kg_ocean)
Chem_joined$deltaDO<-with(Chem_joined,do_mg_l-do_mg_l_ocean)
# make pool ID numeric
Chem_joined$pool_id<-as.numeric(Chem_joined$pool_id)

## make a summary dataframe of all the deltas
ChemSummary<-Chem_joined%>%
  select(starts_with("delta"), site, pool_id)%>%
  group_by(site, pool_id)%>%
  summarise_if(is.numeric, .funs = c('mean', 'var', 'max'))

# join the chem data with the percent cover data so that it is appropriately grouped by site and pool id
Chem_Comm<-left_join(ChemSummary, CoverbyGroups)

# last chem column
l<-ncol(Chem_Comm)-9
chem_pca<-prcomp(Chem_Comm[,c(3:l)],scale. = TRUE, center = TRUE)

# get variance explained by each PC
var_ex<-summary(chem_pca)$importance

ggplot(chem_pca$x)+
  geom_point(aes(x = PC1, y = PC2, col = Chem_Comm$ProducerDom))+
  scale_colour_gradient2(low = "red", mid = "grey",
                        high = "darkgreen", midpoint = 0, space = "Lab",
                        na.value = "grey50", guide = "colourbar", aesthetics = "colour")+
  geom_segment(data = 10*chem_pca$rotation, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'darkred')+
  geom_text(data = 10*chem_pca$rotation, aes(x=PC1, y=PC2, label=row.names(chem_pca$rotation)),color='darkred')+
  xlab(paste('PC1', 100*round(var_ex[2,1],2), "% explained"))+
  ylab(paste('PC2', 100*round(var_ex[2,2],2), "% explained"))+
  labs(color = "Producer Dominance")+
  theme_bw()



## Calculate change over time within each tide pool without accounting for the ocean
ChemData %>%
  group_by(site, pool_id, day_night) %>%
  arrange(day_night,time_point, site, pool_id) %>%
  mutate(diff = nn_umol_l - lag(nn_umol_l, default = first(nn_umol_l)))%>%
  View()
## CCA analysis with species and chemistry

var.cca<-cca(Chem_Comm[,24:31], Chem_Comm[,3:23])
mod<-lmer(chem_pca$x[,1]~Chem_Comm$ProducerDom+(1|Chem_Comm$site))
