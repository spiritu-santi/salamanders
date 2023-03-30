taxa<-read.csv("level-4.csv")
taxa$Organism<-gsub("Cryptobranchus alleganiensis alleganiensis","Cryptobranchus alleganiensis",taxa$Organism)
taxa$Organism<-gsub("Cryptobranchus alleganiensis bishopi","Cryptobranchus alleganiensis",taxa$Organism)


Species<-unique(sort(taxa$Organism))


for(x in Species){
  tax<-taxa[taxa$Organism%in%x,][,2:232]
  bact<-as.data.frame(colMeans(tax))
  colnames(bact)<-x
  bact[,1]<-(bact[,1]/colSums(bact))*100
  if(exists("table_tax")){
    table_tax<-cbind(table_tax,bact)}else{
      table_tax<-bact
    }
}

include<-names(sort(rowSums(table_tax*100),decreasing = T)[1:40])
include<-include[-grep("[.]_",include)]
include<-include[-grep("unclassified",include)]
indx<-which(rownames(table_tax)%in%include)
table_tax<-table_tax[indx,]
Other<-100-colSums(table_tax[,-ncol(table_tax)])
table_tax<-rbind(table_tax,Other)
names<-row.names(table_tax)[-nrow(table_tax)]
names<-gsub("[.][.]",".",names)
names<-gsub("[.][.]",".",names)
names<-gsub("[.][.]",".",names)
names<-gsub("[.]$","",names)
names<-gsub(".*[.].*[.].*[.](.*)", "\\1", names)

row.names(table_tax)<-c(names,"Other")
table_tax$Taxa<-row.names(table_tax)


tablita_tax<-melt(table_tax)

colnames(tablita_tax)<-c("Taxa","Specie","Abundance")

colores<-c( "navy"          , "aliceblue",      "plum"    ,       "salmon"   ,      "darkgray"      , "tomato"     ,    "lawngreen" ,    
            "blue"          , "brown1"     ,    "seagreen2" ,     "lightsalmon3"  , "red"         ,   "firebrick"  ,    "darkmagenta" ,  
            "rosybrown1"    , "chocolate4"  ,   "burlywood4"    , "gray23"    ,     "darkred"      ,  "peru"      ,     "darkgoldenrod2",
            "maroon4"       , "dodgerblue"     ,"forestgreen"   , "lightyellow3" ,  "cyan"          , "darkgreen"    ,  "blueviolet"    ,
            "orange"        , "darkkhaki"  ,    "darkslategray" , "cornflowerblue", "salmon"   )


tablita_tax$Taxa<-with(tablita_tax,reorder(Taxa,Abundance))

tablita_tax$Specie<-gsub("Pseudoeurycea sp. n. Mozotal 1","Pseudoeurycea sp", tablita_tax$Specie)

tablita_tax$Specie_2<-factor(tablita_tax$Specie,levels = rev(c("Batrachuperus tibetanus",
                                                               "Andrias japonicus",
                                                               "Cryptobranchus alleganiensis",
                                                               "Salamandra salamandra",
                                                               "Taricha granulosa",
                                                               "Notophthalmus viridescens",
                                                               "Cynops pyrrhogaster",
                                                               "Lissotriton boscai",
                                                               "Lissotriton vulgaris",
                                                               "Lissotriton helveticus",
                                                               "Ichthyosaura alpestris",
                                                               "Triturus cristatus",
                                                               "Triturus marmoratus",
                                                               "Echinotriton andersoni",
                                                               "Ambystoma mexicanum",
                                                               "Ambystoma altamirani",
                                                               "Plethodon glutinosus",
                                                               "Plethodon cinereus",
                                                               "Pseudoeurycea sp",
                                                               "Pseudoeurycea rex",
                                                               "Pseudoeurycea lynchi",
                                                               "Aquiloeurycea cafetalera",
                                                               "Bolitoglossa lincolni",
                                                               "Bolitoglossa franklini",
                                                               "Parvimolge townsendi",
                                                               "Chiropterotriton nubilus",
                                                               "Batrachoseps attenuatus")))

abu<-ggplot() +
  geom_bar(data=tablita_tax, aes(x=Specie_2, y=Abundance, fill=Taxa),width = .7, position = "fill",stat = "identity")+
  labs(title="" ,x=NULL, y="Relative abundance (%)") + 
  theme_Publication()+
  theme(plot.title = element_text(size = rel(2)), axis.title.y = element_text(size = rel(1.5)), 
        axis.text = element_text(size=12)) +
  scale_fill_manual(values=colores) +
  theme(axis.text.x = element_text(angle = 90))+
  theme(legend.title = element_blank())+
  guides(fill = guide_legend(reverse=TRUE))+coord_flip()


print(abu)
