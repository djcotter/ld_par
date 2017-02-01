library(ggplot2)
library(ggthemes)
setwd("~/Projects/PAR_variation/")

subpop_data <- read.csv("data/data_combine/combined_subpops.txt", header = TRUE, sep = "\t")


subpop_data$Adj_position <- sapply(subpop_data$Position, function(x){(((x+(x+100000))/2)/100000)})

# col=sapply(subpop_data$Adj_position, function(x){if(x<=26.99){"red"}else if(x >= 881 & x <= 931){"blue"}else if(x>=1549){"red"}else{"black"}}), pch=20)

ALL_females_subset <- subset(subpop_data, subpop_data$Adj_position <= 100)

require('mgcv')
k_num = 18
p1 = ggplot(subpop_data, aes(x = Adj_position))
p1 = p1 + geom_smooth(aes(y = YRI, color = 'YRI'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = LWK, color = 'LWK'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p1 = p1 + geom_smooth(aes(y = GWD, color = 'GWD'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = MSL, color = 'MSL'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = ESN, color = 'ESN'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = ASW, color = 'ASW'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = ACB, color = 'ACB'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1  + coord_cartesian(xlim=c(0,150))
p1 = p1 + labs(list(x = 'Postion (kb)', y = expression(paste('Diversity (', pi, ')'))))
p1 = p1 + theme(axis.title=element_text(size=14), legend.title=element_text(face='bold'), legend.text=element_text(size=15))
p1

p2 = ggplot(subpop_data, aes(x = Adj_position))
p2 = p2 + geom_smooth(aes(y = CHB, color = 'CHB'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p2 = p2 + geom_smooth(aes(y = JPT, color = 'JPT'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p2 = p2 + geom_smooth(aes(y = CHS, color = 'CHS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p2 = p2 + geom_smooth(aes(y = CDX, color = 'CDX'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p2 = p2 + geom_smooth(aes(y = KHV, color = 'KHV'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p2 = p2  + coord_cartesian(xlim=c(0,150))
p2 = p2 + labs(list(x = 'Postion (kb)', y = expression(paste('Diversity (', pi, ')'))))
p2 = p2 + theme(axis.title=element_text(size=14), legend.title=element_text(face='bold'), legend.text=element_text(size=15))
p2

p3 = ggplot(subpop_data, aes(x = Adj_position))
p3 = p3 + geom_smooth(aes(y = CEU, color = 'CEU'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p3 = p3 + geom_smooth(aes(y = TSI, color = 'TSI'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p3 = p3 + geom_smooth(aes(y = FIN, color = 'FIN'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p3 = p3 + geom_smooth(aes(y = GBR, color = 'GBR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p3 = p3 + geom_smooth(aes(y = IBS, color = 'IBS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p3 = p3  + coord_cartesian(xlim=c(0,150))
p3 = p3 + labs(list(x = 'Postion (kb)', y = expression(paste('Diversity (', pi, ')'))))
p3 = p3 + theme(axis.title=element_text(size=14), legend.title=element_text(face='bold'), legend.text=element_text(size=15))
p3

p4 = ggplot(subpop_data, aes(x = Adj_position))
p4 = p4 + geom_smooth(aes(y = GIH, color = 'GIH'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p4 = p4 + geom_smooth(aes(y = PJL, color = 'PJL'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p4 = p4 + geom_smooth(aes(y = BEB, color = 'BEB'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p4 = p4 + geom_smooth(aes(y = STU, color = 'STU'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p4 = p4 + geom_smooth(aes(y = ITU, color = 'ITU'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p4 = p4  + coord_cartesian(xlim=c(0,150))
p4 = p4 + labs(list(x = 'Postion (kb)', y = expression(paste('Diversity (', pi, ')'))))
p4 = p4 + theme(axis.title=element_text(size=14), legend.title=element_text(face='bold'), legend.text=element_text(size=15))
p4

p5 = ggplot(subpop_data, aes(x = Adj_position))
p5 = p5 + geom_smooth(aes(y = MXL, color = 'MKL'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p5 = p5 + geom_smooth(aes(y = PUR, color = 'PUR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p5 = p5 + geom_smooth(aes(y = CLM, color = 'CLM'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p5 = p5 + geom_smooth(aes(y = PEL, color = 'PEL'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p5 = p5  + coord_cartesian(xlim=c(0,150))
p5 = p5 + labs(list(x = 'Postion (kb)', y = expression(paste('Diversity (', pi, ')'))))
p5 = p5 + theme(axis.title=element_text(size=14), legend.title=element_text(face='bold'), legend.text=element_text(size=15))
p5

colors = c("#ba6437",
            "#7f63b8",
            "#98a441",
            "#b94b75",
            "#50b47b")

p6 = ggplot(subpop_data, aes(x = Adj_position))
p6 = p6 + geom_smooth(aes(y = YRI, color = 'AFR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = LWK, color = 'AFR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p6 = p6 + geom_smooth(aes(y = GWD, color = 'AFR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = MSL, color = 'AFR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = ESN, color = 'AFR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = ASW, color = 'AFR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = ACB, color = 'AFR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = CHB, color = 'EAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = JPT, color = 'EAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p6 = p6 + geom_smooth(aes(y = CHS, color = 'EAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = CDX, color = 'EAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = KHV, color = 'EAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = CEU, color = 'EUR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = TSI, color = 'EUR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p6 = p6 + geom_smooth(aes(y = FIN, color = 'EUR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = GBR, color = 'EUR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = IBS, color = 'EUR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = GIH, color = 'SAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = PJL, color = 'SAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p6 = p6 + geom_smooth(aes(y = BEB, color = 'SAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = STU, color = 'SAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = ITU, color = 'SAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = MXL, color = 'AMR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = PUR, color = 'AMR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p6 = p6 + geom_smooth(aes(y = CLM, color = 'AMR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6 + geom_smooth(aes(y = PEL, color = 'AMR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p6 = p6  + coord_cartesian(xlim=c(0,150))

p6 = p6 + scale_color_manual(name = '', values = colors, breaks=c('AFR','SAS','AMR','EAS','EUR'))

p6 = p6 + labs(list(x = 'Postion (kb)', y = expression(paste('Diversity (', pi, ')'))))
p6 = p6 + theme(axis.title=element_text(size=14), legend.title=element_text(face='bold'), legend.text=element_text(size=15))
p6














# colors = c("#ba6437",
#            "#7f63b8",
#            "#98a441",
#            "#b94b75",
#            "#50b47b")
# 
# require('mgcv')
# k_num = 18
# p1 = ggplot(ALL_females, aes(x = position))
# #p1 = p1 + geom_smooth(aes(y = ALL_pi, color = 'ALL'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
# p1 = p1 + geom_smooth(aes(y = AFR_pi, color = 'AFR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
# p1 = p1 + geom_smooth(aes(y = SAS_pi, color = 'SAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
# p1 = p1 + geom_smooth(aes(y = AMR_pi, color = 'AMR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
# p1 = p1 + geom_smooth(aes(y = EUR_pi, color = 'EUR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
# p1 = p1 + geom_smooth(aes(y = EAS_pi, color = 'EAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
# 
# p1 = p1 + scale_color_manual(name = '', values = colors, breaks=c('AFR','SAS','AMR','EAS','EUR'))
# 
# p1 = p1  + coord_cartesian(xlim=c(0,150))
# 
# p1 = p1 + geom_point(aes(y = ALL_pi)) + labs(list(x = 'Postion (kb)', y = expression(paste('Diversity (', pi, ')'))))
# 
# p1 = p1 + theme(axis.title=element_text(size=14), legend.title=element_text(face='bold'), legend.text=element_text(size=15))
# 
# p1 + coord_cartesian(ylim = c(0,0.0025), xlim=c(0,150))
# 
# 
# ggsave(file="~/Research/thesis/tables_and_figures/defense_figures/superpop_diversity_along_boundary_all.png", width = 14, height = 5, units = 'in')
