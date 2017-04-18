ld_data <- YRI_females_avg_R2_100kb_windows_1Mb_ldBins
ld_data <- YRI_females_avg_R2_100kb_windows_200kb_ldBins

ld_data$position <- ((ld_data$X1 + ld_data$X2)/2)/1000000
ld_data_subset <- subset(ld_data, ld_data$position <= 30)
plot(ld_data$position, ld_data$X3, pch=20, ylim=c(0,1))
plot(ld_data_subset$position, ld_data_subset$X3, pch=20)
