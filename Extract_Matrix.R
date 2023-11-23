library(dplyr)
library(tidyr)
library(ggpubr)
#setwd("E:\\BreastCancer_Project")
#data <- read.csv("exp_mat.csv")
#sectar <- read.csv("Sec_Tar.csv")
#data <- data %>% separate(ID, c('Symbol', 'gene_id'))
#data2 <- data[data$Symbol %in% unique(sectar$Symbol),]
#row.names(data2) <- data2$Symbol
#data2 <- subset(data2, select = -Symbol )
#data2 <- subset(data2, select = -gene_id )
#write.csv(data.frame(ID=row.names(data2),data2), "Sectar_exp_mat.csv", quote = F, row.names = F)
data <- read.csv("Sectar_exp_mat.csv")
row.names(data) <- data$ID
data <- subset(data, select = -ID )
tdata <- as.data.frame(t(data))
pdf("corelation_plot.pdf",width=20,height=10,paper='special')
ggscatter(tdata, x = "ESR1", y = "CHRNA5", 
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "ESR1", ylab = "CHRNA5", color = "blue",
          fill = "gray", shape = 19,
          size = 2)
dev.off()
pdf("corelation_Mplots2.pdf",width=20,height=10,paper='special')
tdata %>% 
  gather(key = variable, value = values, ARHGAP11A:BMP5) %>% 
  ggplot(aes(ESR1, values)) + 
  geom_point() + 
  facet_grid(. ~ variable, scales = "free_x") + 
  geom_smooth(method = "lm", se = FALSE) 
dev.off()

tdata2 <- round(cor(tdata,use="pairwise.complete.obs"),2)
tdata3 <- as.data.frame(tdata2)
write.csv(data.frame(ID=row.names(tdata3),tdata3), "Corelation_mat.csv", quote = F, row.names = F)
library(ggcorrplot)
pdf("corelation_Matrix.pdf",width=20,height=10,paper='special')
ggcorrplot(tdata2,hc.order = TRUE,
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"))
dev.off()
