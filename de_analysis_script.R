
### Andy's analysis script - summary from the R training course ###
### A clean script for future use###
### update 23/08/2023 ###



### ----- Load libraries ----- ###

library(ggplot2) # most popular plot library #

library(ggrepel) # tide the overlap text on plots #

library(reshape2) # transform data between wide and long formats # # use the melf function#


library(amap) # hierarchical clustering and principal component analysis #

library(clusterProfiler) # analyze and visualize functional profiles (GO and KEGG) of gene and gene clusters #

library(org.Mm.eg.db) 
# all of the known mouse annotations.
# The name means organism (org) mus musculus (Mm) database. 

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
# So, if you do this for human you will need to download and load a DIFFERENT database. For example, human is org.Hs.eg.db. 
# Here is the list of annotations by species: https://www.bioconductor.org/packages/release/data/annotation/

library(devEMF) # saving plots as EMF files #



### ----- Load data files ----- ###

## example: read files ##

# set working directory 
setwd("C:\\Users\\jg1u21\\OneDrive - University of Southampton\\Documents\\d10") 


# read table function 
EM = read.table("EM.csv", header = TRUE, row.names = 1, sep="\t")
DE = read.table("DE_Senes_vs_Prolif.csv", header=TRUE, row.names=1, sep= "\t")
ANNO = read.table("Human_Background_GRCh38.p13.csv", header=TRUE, row.names=1, sep= "\t")


# join all of the tables together to create a “Master” table 
master_temp = merge(EM, DE, by.x = 0, by.y = 0)
master = merge(master_temp, ANNO , by.x = 1, by.y = 0)


# change the names to make the row names become gene symbols 
row.names(master) = master[,"Row.names"]
row.names(master) = master[,"SYMBOL"]
names(master)[1] = "ENSEMBL"


# create EM SYMBOLS 
EM_SYMBOLS = master[,c(2:7)]


# remove rows from the master table that have NA.s 
master = na.omit(master)


# sort master table by p 
order(master[,"p.adj"], decreasing=FALSE)
sorted_order=order(master[,"p.adj"],decreasing=FALSE)
master = master[sorted_order,] # don't forget the comma 


# make a new column in your master table for mean expression 
master$mean = rowMeans(master[,2:7])


# add a column for -log10p to the master table 
master$mlog10p = -log10(master $p)


# add a column flagging significance to the master table 
master$sig = as.factor(master$p.adj<0.05&abs(master$log2fold)>1.0)


# create sig_genes 
master_sig = subset(master, sig==TRUE)
sig_genes = row.names(master_sig) # select significant genes names as a new list #


# make a scaled expression matrix 
EM_scaled=na.omit(data.frame(t(scale(t(EM_SYMBOLS)))))


# create EM sig table 
EM_SYMBOLS_sig = EM_SYMBOLS[sig_genes,]
EM_scaled_sig = EM_scaled[sig_genes,]


# save your symbols expression table to disk 
write.table(EM_SYMBOLS,file="EM_SYMBOLS.csv",sep="\t")
write.table(EM_scaled,file="EM_scaled.csv",sep="\t")
write.table(master,file="master.csv",sep="\t")
write.table(master_sig,file="master_sig.csv",sep="\t")
write.table(EM_SYMBOLS_sig,file="EM_SYMBOLS_sig.csv",sep="\t")
write.table(EM_scaled_sig,file="EM_scaled_sig.csv",sep="\t")



### ----- Make plots ----- ###

## --- VOLCANO PLOT --- ##

# parse
master_up = subset(master,p.adj < 0.05 & log2fold > 1)
master_down = subset(master,p.adj < 0.05 & log2fold < -1)
master_up_top5 = master_up[1:5,]
master_down_top5 = master_down[1:5,]

# make plot
ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  
  # adds the dots
  geom_point(colour = "black") +
  geom_point(data=master_up, colour = "red") +
  geom_point(data=master_down, colour = "blue") +
  
  # adds the text labels for top 5 genes
  geom_text_repel(data=master_up_top5, aes(label=SYMBOL), colour = "red") + # the SYMBOL is the same name inthe column #
  geom_text_repel(data=master_down_top5, aes(label=SYMBOL), colour = "blue") + 
  
  # adds the fancy lines
  geom_vline(xintercept=-1,linetype="dashed") +
  geom_vline(xintercept=1,linetype="dashed") +
  geom_hline(yintercept=-log10(0.05),linetype="dashed") +
  
  # adds the theme and axis titles
  theme_bw() +
  labs(title = "Volcano", x= "Log2 fold change", y= "-log10 p")


# save plot in png, pdf and emf format respectively
png("volcano.png", width=500,height=500)
print(ggp)
dev.off()

pdf("volcano.pdf", width=8,height=8)
print(ggp)
dev.off()

emf("volcano.emf", width=8,height=8)
print(ggp)
dev.off()



## --- MA PLOT --- ##

# make plot
ggp = ggplot(master, aes(x=log10(mean), y=log2fold)) + 
  
  # adds the dots
  geom_point(colour = "black") +
  geom_point(data=master_up, colour = "red") +
  geom_point(data=master_down, colour = "blue") +
  
  # adds the fancy lines
  geom_hline(yintercept=1,linetype="dashed") +
  geom_hline(yintercept=-1,linetype="dashed") +
  
  # adds the theme and axis titles
  theme_bw() +
  labs(title = "MA", x= "Log10 mean expression", y= "Log2 fold change")

# Note: I didn't add dot lables, becuase I don't personally like them in MA.

# save the plot
png("ma.png", width=500,height=500)
ggp
dev.off()



## --- Theme --- ##

my_theme = theme(
  panel.grid = element_blank(),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=1),
  plot.background = element_blank(), 
  legend.background = element_rect(fill="transparent", colour=NA),
  legend.key = element_rect(fill="transparent", colour=NA),
  plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial", face="bold"),
  title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, family="Arial", face="bold"),
  axis.text.y = element_text(size = 11, margin = margin(r = 5),hjust=1,vjust=0.5, family="Arial", face="bold",colour="black"),
  axis.text.x = element_text(size = 11, margin = margin(t = 5),hjust=0.5,vjust=1, family="Arial", face="bold",colour="black"), 
  axis.title.y = element_text(size = 12, margin = margin(r = 10),angle = 90,hjust=0.5,vjust=0.5, family="Arial", face="bold"),
  axis.title.x = element_text(size = 12, margin = margin(t = 10),hjust=0.5,vjust=1, family="Arial", face="bold"),
  legend.text=element_text(size=12, family="Arial", face="bold"),
  legend.title=element_blank(), 
  legend.key.size=unit(1,"line"),
  plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"),
  strip.text.x = element_text(size = 12, family="Arial", face="bold", vjust=1),
  panel.spacing = unit(1, "lines")
)



## --- Volcano With Legend --- ##

# make the plot
ggp = ggplot(master, aes(x=log2fold, y=mlog10p)) + 
  
  # add the dots
  geom_point(aes(colour = "a")) +
  geom_point(data=master_up, aes(colour = "b")) +
  geom_point(data=master_down, aes(colour = "c")) +
  
  # add the data labels
  geom_text_repel(data=master_up_top5, aes(label=SYMBOL, colour = "b"),  show.legend=FALSE) + 
  geom_text_repel(data=master_down_top5, aes(label=SYMBOL, colour = "c"), show.legend=FALSE) + 
  
  # choose colours and legend names
  scale_color_manual(values=c("black", "red", "blue"), labels = c("No change", "Up", "Down"), name = "") + 
  
  # add the dashed lines
  geom_vline(xintercept=-1,linetype="dashed") +
  geom_vline(xintercept=1,linetype="dashed") +
  geom_hline(yintercept=-log10(0.05),linetype="dashed") +
  
  # add the theme and axis names
  my_theme +
  labs(title = "Volcano Plot", x= "Log2 fold change", y= "-Log10 p")

# save
png("plots/volcano.png", width=500,height=500)
ggp
dev.off()


##-- MA With Legend --##

# make plot
ggp = ggplot(master, aes(x=log10(mean), y=log2fold)) + 
  
  # adds the dots
  geom_point(aes(colour = "a")) +
  geom_point(data=master_up, aes(colour = "b")) +
  geom_point(data=master_down, aes(colour = "c")) +
  
  # choose colours and legend names
  scale_color_manual(values=c("black", "red", "blue"), labels = c("No change", "Up", "Down"), name = "") + 
  
  # adds the fancy lines
  geom_hline(yintercept=1,linetype="dashed") +
  geom_hline(yintercept=-1,linetype="dashed") +
  
  # adds the theme and axis titles
  my_theme +
  labs(title = "MA", x= "Log10 mean expression", y= "Log2 fold change")

# Note: I didn't add dot lables, becuase I don't personally like them in MA.

# save
png("plots/ma.png", width=500,height=500)
ggp
dev.off()



## --- PCA --- ##

# read sample info
ss = read.table("sample_sheet.csv", header = TRUE, sep="\t")

# change levels if necessary
ss$SAMPLE_GROUP = factor(ss$SAMPLE_GROUP, levels = c("Prolif", "Senes"))

# run PCA
xx = prcomp(t(EM_scaled))
pca_coordinates = data.frame(xx$x)

# get % variation 
vars = apply(xx$x, 2, var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100
prop_y = round(vars["PC2"] / sum(vars),4) * 100

x_axis_label = paste("PC1 (" ,prop_x, "%)", sep="")
y_axis_label = paste("PC2 (" ,prop_y, "%)", sep="")

# plot
ggp = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = ss$SAMPLE_GROUP)) + 
  geom_point() + 
  geom_label_repel(aes(label=ss$SAMPLE), show.legend=FALSE) + 
  labs(title = "PCA", x= x_axis_label, y= y_axis_label) + 
  my_theme
ggp

# saves - note scales png size to % variance - cool.
png("PCA.png", width=600*(prop_x / prop_y),height=600)
ggp
dev.off()



## --- Expression Density --- ##


# melts EM SYMBOLS
EM_SYMBOLS.m = melt(EM_SYMBOLS)

# gets a groups column - so we can colour the plot
rows = nrow(EM_SYMBOLS)
groups = c(rep(rep("a", rows),3), rep(rep("b", rows),3)) # no need to change "3" 
EM_SYMBOLS.m$groups = groups

# plot
ggp = ggplot(EM_SYMBOLS.m, aes(x = log10(value+0.01), fill = groups)) + 
  geom_density(alpha = 0.75) + 
  facet_wrap(~variable, ncol=3) + 
  my_theme + 
  theme(strip.background = element_rect(fill="transparent", size=0), legend.position = "none") + 
  labs(x = "Expression (log10)", y = "Density")

# saves
png("Expression_Density.png", width=800,height=800)
ggp
dev.off()



## --- Single Gene Boxplots --- ##


# creates the new table
gene_data = EM_SYMBOLS["CCND2", ]
gene_data = data.frame(t(gene_data))
gene_data$sample_group = ss$SAMPLE_GROUP
names(gene_data) = c("CCND2","groups")

# makes the plot
ggp = ggplot(gene_data, aes(x=groups, y= CCND2)) + 
  geom_boxplot(aes(fill = groups))  + 
  my_theme

# saves
png("CCND2 boxplot.png", height = 400, width = 400)
print(ggp)
dev.off()



## --- Facet Boxplot --- ##


# get the gene list
candidate_genes = row.names(master[1:10,])
candidate_genes = c(candidate_genes)

# create the table
gene_data = EM_scaled[candidate_genes, ]
gene_data = data.frame(t(gene_data))
gene_data$groups = ss$SAMPLE_GROUP

# melt
gene_data.m = melt(gene_data, id.vars = "groups")

# makes the plot
ggp = ggplot(gene_data.m, aes(x=groups, y=value, fill = groups)) + 
  geom_boxplot()  + 
  my_theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~variable, ncol=5) + 
  labs(x = "", y = "Expression (z-score)")

# saves
png("facet_boxplot.png", height = 400, width = 400)
print(ggp)
dev.off()



## --- Multi Gene Boxplot --- ##


# get the gene list
candidate_genes = row.names(master[1:10,])
candidate_genes = c(candidate_genes)

# create the table
gene_data = EM_scaled[candidate_genes, ]
gene_data = data.frame(t(gene_data))
gene_data$groups = ss$SAMPLE_GROUP

# melt
gene_data.m = melt(gene_data, id.vars = "groups")

# makes the plot
ggp = ggplot(gene_data.m, aes(x=variable, y=value, fill = groups)) + 
  geom_boxplot()  + 
  my_theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Expression (z-score)")

png("multi_boxplot.png", height = 400, width = 400)
print(ggp)
dev.off()



## --- Heatmap --- ## 


# makes a matrix
hm.matrix = as.matrix(EM_scaled_sig)

# does the y clustering
y.dist = Dist(hm.matrix, method="spearman")
y.cluster = hclust(y.dist, method="average")
y.dd = as.dendrogram(y.cluster)
y.dd.reorder = reorder(y.dd,0,FUN="average")
y.order = order.dendrogram(y.dd.reorder)

hm.matrix_clustered = hm.matrix[y.order,]

# melt and plot
hm.matrix_clustered.m = melt(hm.matrix_clustered)

# colour palette
colours = c("blue","black","yellow")
palette = colorRampPalette(colours)(100)

# plot
ggp = ggplot(hm.matrix_clustered.m, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = palette) + 
  labs(x="", y="") + 
  theme(legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'), axis.text.y = element_blank(), axis.ticks=element_blank())

# saves
png("heatmap.png", height = 800, width = 400)
print(ggp)
dev.off()



## --- Heatmap Rug --- ## # Rarely used #

# Function that gets the default gg plot colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# uses the function to get the colours
colours = gg_color_hue(3)

# rug for sample groups
groups_data = as.matrix(as.numeric(ss$SAMPLE_GROUP))
groups_data = melt(groups_data)

# makes the plot
ggp = ggplot(groups_data, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile(linetype="blank") + 
  my_theme + 
  scale_fill_gradientn(colours = colours) + 
  labs(x = "Samples") +
  theme(plot.margin=unit(c(0,1,1,1), "cm"), axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())

ggp



## --- Pathway --- ##

# converts from ensembl Symbols to Entrez - needed for cluster profiler
sig_genes_entrez = 
  bitr(sig_genes, 
       fromType = "SYMBOL", 
       toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)

# gets the enrichment
pathway_results = enrichGO(gene = sig_genes_entrez$ENTREZID,OrgDb = org.Hs.eg.db, readable = T,ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

#Plot the enriched ontologies in various different ways:
ggp = barplot(pathway_results, showCategory=10)
png("pathway_barplot.png", height = 800, width = 800)
ggp
print(ggp)
dev.off()

ggp = dotplot(pathway_results, showCategory=10)
png("pathway_dotplot.png", height = 800, width = 800)
print(ggp)
dev.off()

ggp = goplot(pathway_results, showCategory = 10)
png("pathway_goplot.png", height = 800, width = 800)
print(ggp)
dev.off()

ggp = cnetplot(pathway_results, categorySize="pvalue")
png("pathway_cnet.png", height = 800, width = 800)
print(ggp)
dev.off()


# get the gene list for the top ontology

# extract the info from go enrich 
gene_sets = pathway_results$geneID
description = pathway_results$Description
p.adj = pathway_results$p.adjust

# combine into a data frame
ora_results = data.frame(cbind(gene_sets, description, p.adj))

enriched_gene_set = as.character(ora_results[1,1])
candidate_genes = unlist(strsplit(enriched_gene_set, "/")) 


# Pathway Boxplot

# create the table
gene_data = EM_scaled[candidate_genes, ]
gene_data = data.frame(t(gene_data))
gene_data$groups = ss$SAMPLE_GROUP

# melt
gene_data.m = melt(gene_data, id.vars = "groups")

# makes the plot
ggp = ggplot(gene_data.m, aes(x=variable, y=value, fill = groups)) + 
  geom_boxplot()  + 
  my_theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png("plots/pathway_multi_boxplot.png", height = 400, width = 2000)
print(ggp)
dev.off()



# Pathway Heatmap

# makes a matrix - selecting only the candidate genes
hm.matrix = as.matrix(EM_scaled[candidate_genes,])

# does the y clustering
y.dist = Dist(hm.matrix, method="spearman")
y.cluster = hclust(y.dist, method="average")
y.dd = as.dendrogram(y.cluster)
y.dd.reorder = reorder(y.dd,0,FUN="average")
y.order = order.dendrogram(y.dd.reorder)

# re-orders
hm.matrix_clustered = hm.matrix[y.order,]

# melt
hm.matrix_clustered = melt(hm.matrix_clustered)

# colour palette
colours = c("blue","black","yellow")
palette = colorRampPalette(colours)(100)

# plot
ggp = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = palette) + 
  labs(x="", y="") + 
  theme(legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'), axis.text.y = element_blank(), axis.ticks=element_blank())

png("plots/pathway_heatmap.png", height = 800, width = 400)
print(ggp)
dev.off()



### ----- Save plots ----- ###

png("plot_path.png", height=800, width=800)
print(ggp)
dev.off()


pdf("plot_path.pdf", width=8,height=8)
print(ggp)
dev.off()

emf("plot_path.emf", height=10, width=10)
print(ggp)
dev.off()

# Note: height / width should be ~ 80 times smaller for emf than png (cm not pixels).







