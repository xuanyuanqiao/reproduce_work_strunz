library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(cluster)
library(tidyr)

##load Dropseq data
data_dir<-‘c:/2019/GEO/wholelung’
list.files(data_dir)
data<-Read10X(data.dir=data_dir)
adata<-CreateSeuratObject(raw.data = data, project = 'WholeLung')

##load cell info and add them to adata
cellinfo<-read.csv("c:/2019/GEO/wholelung/cellinfo.csv", header = TRUE)

res.2<-cellinfo$res.2
names(res.2)<-rownames(adata@meta.data)
adata@meta.data$res.2<-res.2

identifier<-cellinfo$identifier
names(identifier)<-rownames(rm_cc@meta.data)
rm_cc@meta.data$identifier<-identifier

grouping<-cellinfo$grouping
names(grouping)<-rownames(adata@meta.data)
adata@meta.data$grouping<-grouping

cell.type<-cellinfo$cell.type
names(cell.type)<-rownames(adata@meta.data)
adata@meta.data$cell.type<-cell.type

metacelltype<-cellinfo$metacelltype
names(metacelltype)<-rownames(adata@meta.data)
adata@meta.data$metacelltype<-metacelltype

spline_cluster<-cellinfo$spline_cluster
names(spline_cluster)<-rownames(adata@meta.data)
adata@meta.data$spline_cluster<-spline_cluster

##Filtercells
adata<-FilterCells(adata,subset.names='nUMI',high.thresholds = 5000)
##Normalization 
adata<-NormalizeData(adata,normalization.method = 'LogNormalize',scale.factor = 10000)
##Scale data to regress out the effect of UMI
adata <- ScaleData(adata,vars.to.regress = 'nUMI')

##measure the effect of cell cycle genes
##load the cell cycle genes file
cc_file<-read.csv("c:/2019/data/Mus_musculus.csv")

##convert the Ensembl ID to gene names(generated in 'David' website)
names<-read.table("c:/2019/data/names.txt")
cc_file$geneID<-names
genes<-cc_file$geneID
s_genes<-cc_file %>% dplyr::filter(phase=="S") %>% pull("geneID")
s_genes<-s_genes[['V1']]
g2m_gene<-cc_file %>% dplyr::filter(phase=="G2/M") %>% pull("geneID")
g2m_gene<-g2m_gene[['V1']]

##check those genes are present in my dataset
which(x = unlist(x = genes) %in% rownames(x = adata@data))

##CellcycleScoring to score the cell cycle genes
adata<-CellCycleScoring(adata,g2m.genes = g2m_gene,s.genes = s_genes,set.ident = TRUE)

##regress out the effect of cell cycle genes(method one)
##control <- ScaleData(adata, vars.to.regress = c("S.Score", "G2M.Score"))
##control<-RunPCA(control,features = c(s.genes,g2m.genes))
##DimPlot(control)

##another way to remove cell cycle influence(method two)
adata@meta.data$cc.dif<-adata@meta.data$S.Score-adata@meta.data$G2M.Score
rm_cc<-ScaleData(adata,vars.to.regress = c('nUMI','cc.dif'))
rm_cc<-RunPCA(rm_cc,features = c(s.genes,g2m.genes))
DimPlot(rm_cc)

##find variable genes and select top 7000
rm_cc<-FindVariableGenes(rm_cc, top.genes = 7000)

##component analysis
rm_cc<-RunPCA(rm_cc,pc.genes = rm_cc@var.genes,pcs.compute = 50)
rm_cc<-FindClusters(rm_cc,dims.use = 1:50,resolution = 2)
DimPlot(rm_cc)

##umap
rm_cc<-RunUMAP(rm_cc,dims.use = 1:50,n_neighbors = 10)

#saveRDS(rm_cc,file='c:/2019/rm_cc.rds')

rm_cc<-readRDS('c:/2019/rm_cc.rds')

##visualization
DimPlot(rm_cc,reduction.use = "umap",group.by = "grouping",vector.friendly = T)
DimPlot(rm_cc,reduction.use = "umap",vector.friendly = T,do.label = TRUE)

##discovery of cell type identity marker genes
cluster0.markers<-FindMarkers(rm_cc,ident.1 = 0,min.pct = 0.25)
head(cluster0.markers,n=5)
##cluster0 is AT2 cells

##findallmarkers() to detect maker genes of all 40 clusters
rm_cc.markers <- FindAllMarkers(rm_cc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers<-rm_cc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
head(rm_cc.markers)
##save the marker.csv and manually assign the clusters to cell type and meta-cell type based on paper's supplementary file
write.csv(markers,"c:/2019/data/markers.csv")

#Draw heatmap to check the marker genes show the difference between the clusters(used in discussion part
top2 <- rm_cc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
p<-DoHeatmap(rm_cc, genes.use = top2$gene,slim.col.label=TRUE) 
p

##edit the umap plot
##draw by days (order the sequence of the legend)
order=c("PBS","d3","d7","d10","d14","d21","d28")
##do.return returns a ggplot subject which can used to modify the plot
umap_plot<-DimPlot(rm_cc,reduction.use = "umap",group.by = "grouping",do.return = TRUE,vector.friendly = T,plot.order = rev(order))
umap_plot+theme(legend.position = "right")

##draw by groups
umap_plot<-DimPlot(rm_cc,reduction.use = "umap",do.return = TRUE,do.label = TRUE)
umap_plot+theme(legend.position = "right")


##FeaturePlot(rm_cc,reduction.use = "umap", features.plot = c("Gzma","Lox","Cxcl15"))


## Generate Fig S1a (show the distribution of bleo treated and control group)
color <- c("blue","red")

##add the identity
rm_cc@meta.data$pbs <- rm_cc@meta.data$grouping == "PBS"
r <- sample(rownames(rm_cc@meta.data))
rm_cc@meta.data <- rm_cc@meta.data[r,]
rm_cc@dr$umap@cell.embeddings <- rm_cc@dr$umap@cell.embeddings[r,]

b<-DimPlot(rm_cc, reduction.use = "umap", group.by = "pbs", vector.friendly = T, cols.use = color, pt.size = 0.5,do.return = T)
##try to change the legend but failed...
b + scale_fill_discrete(breaks=c('FALSE','TRUE'),labels=c('Bleomycin','PBS'))


##generate Fig 1a(assign the celltype to the clusters)
##genes here are typed based on the genes used in the x-axis
genes <- c("Acoxl", "Fabp1", "Lgi3", "Slc34a2", "Gzma", "Arhgef38", "Aox3", "Lamp3", "Ppp1r14c", "Lyz1",
           "Chi3l1", "Itih4", "Cxcl15", "Acot1", "Enpep", "Sftpa1", "Atp6v1c2", "Slc4a5", "Ank3", "Tinag",
           "Pgam1", "Lox", "Ecm1", "Tpm2", "Cxcl14", "C1qc", "Thbs1", "C1qa", "C1qb", "Mmp19", "Mmp12",
         "Emilin2", "S100a4", "Serpina3n", "Gpnmb", "Fn1", "Msr1", "Tnc", "Spp1", "Arg1")


###celltypes are assigned to the cluster manually
##rename the ident of the cluster.
##because using a list of id is failed, so I rename the ident one by one....(later I found I could use lapply() function)
##"Unassigned" means that according to the marker.csv I generated, I couldn't define the identity of the cluster, 
## because either the DE genes didn't appear in the list the author provided, or the genes are expressed several cell types.
##names(id)<-levels(rm_cc)
rm_cc<-RenameIdent(rm_cc,old.ident.name = '5',new.ident.name  = "Ciliated cells")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '35',new.ident.name  = "Ciliated cell subet")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '32',new.ident.name  = "Goblet cells")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '15',new.ident.name  = "Club cells")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '30',new.ident.name  = "AT2 cells")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '0',new.ident.name  = "activated AT2")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '29',new.ident.name  = "AT2 cells")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '19',new.ident.name  = "Krt8+ ADI")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '27',new.ident.name  = "LECs")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '4',new.ident.name  = "VEC")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '34',new.ident.name  = "CECs")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '38',new.ident.name  = "Vcam1+ VECs subet")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '39',new.ident.name  = "Vcam1+ VECs")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '20',new.ident.name  = "Mesothelial cells")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '10',new.ident.name  = "Fibroblasts")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '25',new.ident.name  = "CD103+ DCs")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '26',new.ident.name  = "CD103+ DCs")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '14',new.ident.name  = "DCs")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '16',new.ident.name  = "CD103+ DCs subet")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '24',new.ident.name  = "Ccl17+ DCs")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '23',new.ident.name  = "Plasma cells")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '28',new.ident.name  = "NK cells")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '12',new.ident.name  = "B-lymphocytes")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '11',new.ident.name  = "B-lymphocytes")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '3',new.ident.name  = "T-lymphocytes")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '36',new.ident.name  = "T-lymphocytes")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '22',new.ident.name  = "macrophage")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '1',new.ident.name  = "M2 macrophage")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '8',new.ident.name  = "recruited monocytes")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '17',new.ident.name  = "recruited monocytes")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '6',new.ident.name  = "CD163+/CD11c- IMs")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '7',new.ident.name  = "AM (Bleo)")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '18',new.ident.name  = "AM (PBS)")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '2',new.ident.name  = 'AM (PBS)')
rm_cc<-RenameIdent(rm_cc,old.ident.name = '21',new.ident.name  = "Neutrophils")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '33',new.ident.name  = "Fn1+ macrophages")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '9',new.ident.name  = "Mki67+ Proliferating cells")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '37',new.ident.name  = "unassigned")
rm_cc<-RenameIdent(rm_cc,old.ident.name = '31',new.ident.name  = "unassigned")
##drat the umap plot
DimPlot(rm_cc, reduction.use = "umap", group.by = "ident", vector.friendly = T,do.label = TRUE)

##generate Fig S1g (Dot plot of the marker gene expression of each cluster)

celltypes <- rbind(c(5, "Ciliated cells"),
              c(35, "Ciliated cell subet"),
              c(32, "Goblet cells"),
              c(15, "Club cells"),
              c(30, "AT2 cells"),
              c(0, "activated AT2"),
              c(29, "AT2 cells"),
              c(19, "Krt8+ ADI"),
              c(27, "LECs"),
              c(4, "VEC"),
              c(34, "CECs"),
              c(38, "Vcam1+ VECs subet"),
              c(39, "Vcam1+ VECs"),
              c(20, "Mesothelial cells"),
              c(10, "Fibroblasts"),
              c(26, "CD103+ DCs"),
              c(25, "CD103+ DCs"),
              c(14, "DCs"),
              c(16, "CD103+ DCs subet"),
              c(24, "Ccl17+ DCs"),
              c(23, "Plasma cells"),
              c(28, "NK cells"),
              c(12, "B-lymphocytes"),
              c(11, "B-lymphocytes"),
              c(3, "T-lymphocytes"),
              c(36, "T-lymphocytes"),
              c(22, "macrophage"),
              c(1, "M2 macrophages"),
              c(13, "resolution macrophage"),
              c(8, "recruited monocytes"),
              c(17, "recruited monocytes"),
              c(6, "CD163+/CD11c- IMs"),
              c(7, "AM (Bleo)"),
              c(18, "AM (PBS)"),
              c(2, "AM (PBS)"),
              c(21, "Neutrophils"),
              c(33, "Fn1+ macrophages"),
              c(9, "Mki67+ Proliferating cells"))

celltypes <- data.frame(celltypes)
colnames(celltypes) <- c("res.2", "name")
celltypes$res.2 <- as.numeric(as.character(celltypes$res.2))

##create the cluster - celltype identity of rm_cc data
rm_cc@meta.data$tmp <- NA
lapply(1:nrow(celltypes), function(x) 
  rm_cc@meta.data$tmp[which(rm_cc@meta.data$res.2 == celltypes$res.2[x])] <<- as.character(paste(celltypes$res.2[x], celltypes$name[x], sep = " - ")))

##filter the cells that are assigned cell types
cells <- rownames(rm_cc@meta.data)[which(!is.na(rm_cc@meta.data$tmp))]
##only use cells that are assigned a celltype
control <- SubsetData(rm_cc, cells.use = cells)

##Defien a function：only genes expression level > threshold will be selected
PercentAbove <- function (x, threshold) length(x = x[x > threshold])/length(x = x)
##generate data table
data <- data.frame(FetchData(object = control, vars.all = genes))
data$cell <- rownames(x = data)
data$id <- factor(control@meta.data$tmp, levels = rev(paste(celltypes$res.2, celltypes$name, sep = " - ")))

data <- gather(data, key = genes.plot, value = expression, -c(cell, id))
data <- data %>% group_by(id, genes.plot) %>% summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression, threshold = 0))
data <- data %>% ungroup() %>% group_by(genes.plot) %>% mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale, max = 2.5, min = -2.5))
data$genes.plot <- factor(x = data$genes.plot, levels = rev(x = sub(pattern = "-", replacement = ".", x = genes)))
data$pct.exp[data$pct.exp < 0] <- NA
##plot the dotplot
p <- ggplot(data = data, mapping = aes(x = genes.plot, y = id)) +
  geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) + 
  scale_radius(range = c(0, 6)) + scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p

