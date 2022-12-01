library(tidyverse)
#subsampling vs known variant discovery
subs <- read.delim("~/Documents/MedGenetics/BioinfoAnalysis/SaturationAnalysis/subs.known.vars.txt", row.names = 1, header = T)
t(subs)

subs <- t(subs) %>% reshape2::melt()
subs$Var2 <- subs$Var2*100


ggplot(subs, aes(x=Var2, group=Var2, y=value ))+ geom_boxplot()+
  theme_bw()+
  theme(text = element_text(size=10),
        axis.text.x = element_text(size=15,angle=0, colour = "black"),axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=15, hjust=1, colour = "black"),axis.title.y = element_text(size=20))+
  ylab("Number of Known SNP discovery\n per sample")+xlab("% of Sequencing Reads")+
  scale_x_discrete(limits=c(20, 40, 60, 80, 100))






