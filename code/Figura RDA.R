dferia <- readRDS("~/Desktop/RDA/sites.rds") %>% as_tibble()
Habcentroids1 <- readRDS("~/Desktop/RDA/centroids.rds") %>% as_tibble()
continuous_arrows <-  readRDS("~/Desktop/RDA/arrows.rds") %>% as_tibble()
dferia %>% rename("Class"=site) %>% mutate(Class="sample") -> dferia
Habcentroids1 %>% rename("Class"=Centroids) -> Habcentroids1
continuous_arrows %>% rename("Class"=class) -> continuous_arrows
mult = 5
continuous_arrows %>% filter(Class %in% all_of(c("sample","FamilyAmbystomatidae","FamilyCryptobranchidae","FamilyHynobiidae","FamilyPlethodontidae","FamilySalamandridae","bio1","bio2","bio16","bio18","bio19")))  -> continuous_arrows
pp <- ggplot(dferia,aes(x = CAP1, y = CAP2,fill=Habitat,color=Habitat,group=Habitat)) +
  geom_hline(yintercept=0,color="grey20") +
  geom_vline(xintercept=0,color="grey20") +
  geom_point(shape=19,stroke=0.3,size=2.5,alpha=0.95) +  
  stat_ellipse(size=1) +
  scale_color_manual(values=c(wesanderson::wes_palette("Darjeeling1")),name="Species") +
 #scale_fill_manual(values=c(wesanderson::wes_palette("Darjeeling1")),name="Family") +
  geom_segment(data = continuous_arrows,aes(x = 0, xend = mult * CAP1,y = 0, yend = mult * CAP2,group=Class,fill=NULL), arrow = arrow(length = unit(0.2, "cm")), colour = "black") + geom_text(data = continuous_arrows,
            aes(x= (mult + mult/5) * CAP1, y = (mult + mult/5) * CAP2, 
                label = Class,group=Class,fill=NULL), 
            size = 3,parse = TRUE,color="black") +
  xlim(-5, 5) + ylim(-5, 5) +
  labs(x="CAP1 (8.6%)",y="CAP2 (6.4%)") + theme(panel.background = element_blank(),panel.grid = element_blank(),axis.line = element_line()) + 
  NULL
pp

ggsave(filename = "MS/edit/RDA_2.pdf",plot = pp)







