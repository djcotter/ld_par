library(ggplot2)
library(ggthemes)

ALL_females <- read.delim("~/Research/thesis/superpop_100kb_window_data/ALL_females_100kb_no_diversity.txt", header=FALSE)
AFR_females <- read.delim("~/Research/thesis/superpop_100kb_window_data/AFR_females_100kb_no_diversity.txt", header=FALSE)
AMR_females <- read.delim("~/Research/thesis/superpop_100kb_window_data/AMR_females_100kb_no_diversity.txt", header=FALSE)
EUR_females <- read.delim("~/Research/thesis/superpop_100kb_window_data/EUR_females_100kb_no_diversity.txt", header=FALSE)
EAS_females <- read.delim("~/Research/thesis/superpop_100kb_window_data/EAS_females_100kb_no_diversity.txt", header=FALSE)
SAS_females <- read.delim("~/Research/thesis/superpop_100kb_window_data/SAS_females_100kb_no_diversity.txt", header=FALSE)

ALL_females <- ALL_females_100kb_no_diversity

ALL_females$position <- sapply(ALL_females$V2, function(x){(((x+(x+100000))/2)/100000)})

ALL_females$ALL_pi <- NA
ALL_females$ALL_pi <- ifelse(ALL_females$V5 == 'NA', NA, ALL_females$V4)

AFR_females$pi <- NA
AFR_females$pi <- ifelse(AFR_females$V5 == 'NA', NA, AFR_females$V4)

AMR_females$pi <- NA
AMR_females$pi <- ifelse(AMR_females$V5 == 'NA', NA, AMR_females$V4)

EUR_females$pi <- NA
EUR_females$pi <- ifelse(EUR_females$V5 == 'NA', NA, EUR_females$V4)

EAS_females$pi <- NA
EAS_females$pi <- ifelse(EAS_females$V5 == 'NA', NA, EAS_females$V4)

SAS_females$pi <- NA
SAS_females$pi <- ifelse(SAS_females$V5 == 'NA', NA, SAS_females$V4)

ALL_females$AFR_pi <- AFR_females$pi
ALL_females$AMR_pi <- AMR_females$pi
ALL_females$EUR_pi <- EUR_females$pi
ALL_females$EAS_pi <- EAS_females$pi
ALL_females$SAS_pi <- SAS_females$pi

plot(ALL_females$position, ALL_females$ALL_pi, col=sapply(ALL_females$position, function(x){if(x<=26.99){"red"}else if(x >= 881 & x <= 931){"blue"}else if(x>=1549){"red"}else{"black"}}), pch=20)

ALL_females_subset <- subset(ALL_females, ALL_females$position <= 100)

colors = c("#ba6437",
           "#7f63b8",
           "#98a441",
           "#b94b75",
           "#50b47b")

require('mgcv')
k_num = 18
p1 = ggplot(ALL_females, aes(x = position))
#p1 = p1 + geom_smooth(aes(y = ALL_pi, color = 'ALL'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = AFR_pi, color = 'AFR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = SAS_pi, color = 'SAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p1 = p1 + geom_smooth(aes(y = AMR_pi, color = 'AMR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = EUR_pi, color = 'EUR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = EAS_pi, color = 'EAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 

p1 = p1 + scale_color_manual(name = '', values = colors, breaks=c('AFR','SAS','AMR','EAS','EUR'))

p1 = p1  + coord_cartesian(xlim=c(0,150))

p1 = p1 + geom_point(aes(y = ALL_pi)) + labs(list(x = 'Postion (kb)', y = expression(paste('Diversity (', pi, ')'))))

p1 = p1 + theme(axis.title=element_text(size=14), legend.title=element_text(face='bold'), legend.text=element_text(size=15))

p1 + coord_cartesian(ylim = c(0,0.0025), xlim=c(0,150))

pdf(file = '~/Research/thesis/tables_and_figures/superpop_diversity_along_boundary.pdf', width = 12, height = 7)
p1
dev.off()

png(file = '~/Research/thesis/tables_and_figures/superpop_diversity_along_boundary.png', width = 12, height = 7, units = 'in', res = 120)
p1    
dev.off()










colors = c("#ba6437",
           "#7f63b8",
           "#98a441",
           "#b94b75",
           "#50b47b")

require('mgcv')
k_num = 18
p1 = ggplot(ALL_females, aes(x = position))
#p1 = p1 + geom_smooth(aes(y = ALL_pi, color = 'ALL'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = AFR_pi, color = 'AFR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = SAS_pi, color = 'SAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE)
p1 = p1 + geom_smooth(aes(y = AMR_pi, color = 'AMR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = EUR_pi, color = 'EUR'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 
p1 = p1 + geom_smooth(aes(y = EAS_pi, color = 'EAS'), method = 'gam', formula = y ~ s(x, k = k_num), size = 1, se = FALSE) 

p1 = p1 + scale_color_manual(name = '', values = colors, breaks=c('AFR','SAS','AMR','EAS','EUR'))

p1 = p1  + coord_cartesian(xlim=c(0,150))

p1 = p1 + geom_point(aes(y = ALL_pi)) + labs(list(x = 'Postion (kb)', y = expression(paste('Diversity (', pi, ')'))))

p1 = p1 + theme(axis.title=element_text(size=14), legend.title=element_text(face='bold'), legend.text=element_text(size=15))

p1 + coord_cartesian(ylim = c(0,0.0025), xlim=c(0,150))


ggsave(file="~/Research/thesis/tables_and_figures/defense_figures/superpop_diversity_along_boundary_all.png", width = 14, height = 5, units = 'in')
