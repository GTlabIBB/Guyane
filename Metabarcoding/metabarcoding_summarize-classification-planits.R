library(xlsx)

zotus <- read.table("classification_planits/zotus.sintax.txt", header=F, sep="\t", fill = TRUE)

zotustab <- read.table("output_pool_zotus/zotutab.txt", header=T, sep="\t", fill = TRUE, comment.char = "")
colnames(zotustab) <- c("otu", "x14M973", "x14M974", "x14M975")



zotus$V2[zotus$V2==""] <- "p:(),c:(),o:(),f:(),g:(),s:()"

id <- data.frame(do.call('rbind', strsplit(as.character(zotus$V1),';',fixed=TRUE)))[1]

classif_raw <- data.frame(do.call('rbind', strsplit(as.character(zotus$V2),',',fixed=TRUE)))


f <- data.frame(do.call('rbind', strsplit(gsub(")", "",classif_raw$X4),'(',fixed=TRUE)))
f$X1 <- gsub("f:", "", f$X1)
f$X2 <- gsub("f:", "", f$X2)
f$X2 <- as.numeric(f$X2)
g <- data.frame(do.call('rbind', strsplit(gsub(")", "",classif_raw$X5),'(',fixed=TRUE)))
g$X1 <- gsub("g:", "", g$X1)
g$X2 <- gsub("g:", "", g$X2)
g$X2 <- as.numeric(g$X2)
s <- data.frame(do.call('rbind', strsplit(gsub(")", "",classif_raw$X6),'(',fixed=TRUE)))
s$X1 <- gsub("s:", "", s$X1)
s$X2 <- gsub("s:", "", s$X2)
s$X2 <- as.numeric(s$X2)

classif <- cbind(id, f, g, s)
colnames(classif) <- c("otu", "family", "prob_family", "genus", "prob_genus", "species", "prob_species")


aggregated <- merge(classif, zotustab, by="otu")
aggregated <- as.data.frame(aggregated)

#filter out reads classified with probability <95%
aggregated_filtered <- subset(aggregated, (prob_species >= 0.95) & (species != "NA"))

aggregated_filtered_species <- aggregate(cbind(x14M973, x14M974, x14M975) ~ species, FUN = sum, aggregated_filtered)

#save as xlsx file
write.xlsx(aggregated_filtered_species, "classification_planits/zotus_aggregated_species95.xlsx", row.names=F)
write.xlsx(aggregated, "classification_planits/zotus_aggregated_probabilities.xlsx", row.names=F)
