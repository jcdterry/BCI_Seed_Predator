  ALLINPLOT<-ggtree(WellSampledTree, layout = 'circular', col='darkgrey') %<+% DF +
  scale_fill_manual(values = c( 'white' , 'Black' ),
                    labels= c('Not predated', 'Predated'),
                    name='')+
  theme(legend.position=c(0.1, 0.95), legend.background = element_blank())+
  geom_tippoint(aes(fill=Predated), shape=21)+
  geom_tiplab2(aes(image=WaspPic), geom='image', hjust = 0.1, offset=25, size=0.04 )+
  geom_tiplab2(aes(image=MothPic), geom='image', hjust = 0.1, offset=15, size=0.015)+
  geom_tiplab2(aes(image=BeetlePic), geom='image', hjust = 0.1, offset=5, size=0.015)+
  geom_tiplab2(size=2, hjust = -0.1, offset=30, fontface = "italic")

  
  ggsave(plot = ALLINPLOT, '../Figures/AllINPLOT.png', width = 10, height=10, units = 'in')
  

ggtree(WellSampledTree, layout = 'circular', col='darkgrey') %<+% DF +
  scale_fill_manual(values = c( 'white' , 'Black' ),
                    labels= c('Not predated', 'Predated'),
                    name='')+
  theme(legend.position=c(0.1, 0.9), legend.background = element_blank())+
  geom_tippoint(aes(fill=Predated), shape=21)+
  labs(caption='Phylogenetic distribution of observed predators amongst the "well-sampled"" species')