# Load R libs ####
library(readxl)
library(igraph)

#Read receptor-ligand table ####
merged <- read_excel("Desktop/hw/1800016203-轩辕乔-大作业/reference_paper/41467_2020_17358_MOESM9_ESM.xlsx",
                     sheet = 2)
merged <- merged[-which(merged$cluster.lig %in% c("plasma_cells", "unassigned")),]
merged <- merged[-which(merged$cluster.rec %in% c("plasma_cells", "unassigned")),]

#Add metacelltype groups ####
assignments <- list(
  T_cells = c("t_cells"),
  B_cells = c("b_cells"),
  NK_cells = "nk_cells",
  Dendritic_cells = c("ccl17_dc", "cd103_dc", "dc"),
  Monocytes = c("nonclassical_mo", "classical_mo"),
  Granulocytes = "granulocytes",
  Macrophages = c("am", "fn1_mcpg", "mono_im", "im", "cd163_im", 14, 19, 29),
  Endothelial_cells = c("cec", "lec", "vec", "vcam1_vec"),
  Mesothelial_cells = c("mesothelial_cells"),
  Fibroblasts = c("fibroblasts"),
  Smooth_muscle_cells = "smooth_muscle_cells",
  Alveolar_epithelium = c("alv_epithelium"),
  Club_cells = "club_cells",
  Goblet_cells = "goblet_cells",
  Ciliated_cells = c("ciliated_cells")
)

lapply(names(assignments), function(x){
  print(x)
  tmp <- intersect(assignments[[x]], merged$cluster.lig)
  merged$metacelltype.lig[which(merged$cluster.lig %in% tmp)] <<- x
  tmp <- intersect(assignments[[x]], merged$cluster.rec)
  merged$metacelltype.rec[which(merged$cluster.rec %in% tmp)] <<- x
})

# Generate adjacency table of all pairs ####
tmp <- table(paste(merged$metacelltype.lig, merged$metacelltype.rec, sep = "|"))
tmp <- data.frame(do.call(rbind, 
                          lapply(names(tmp),function(x) strsplit(x, '|', fixed = T)[[1]])),
                  as.numeric(tmp))
colnames(tmp) <- c('receptor', 'ligand', 'freq')
tmp$receptor <- as.character(tmp$receptor)
tmp$ligand <- as.character(tmp$ligand)
celltypes <- setdiff(unique(c(tmp[,1], tmp[,2])), "NA")
matr <- matrix(0, length(celltypes), length(celltypes))
colnames(matr) <- rownames(matr) <- celltypes
lapply(1:nrow(tmp), function(x){
  try(matr[tmp[x,1], tmp[x,2]] <<- tmp[x, 3]  )
})
adj_all <- matr

# Generate plots ####
vertex_order <- c("Ciliated_cells", "Club_cells", "Goblet_cells", "Alveolar_epithelium",
                  "B_cells", "Dendritic_cells", "Granulocytes", "Macrophages", "Monocytes", 
                  "NK_cells", "T_cells", "Mesothelial_cells", "Smooth_muscle_cells",
                  "Fibroblasts", "Endothelial_cells")
library(igraph)
Graph <- as.directed(graph.adjacency(adj_all, weighted = TRUE, mode = "directed", diag = F))
#mode:directed or un directed(有指向或无指向)
layout_style<-layout_in_circle(Graph, order = vertex_order) #设置图的布局方式
V(Graph)$size <- degree(Graph) #节点大小与点中心度成正比，中心度即与该点相连的点的总数
V(Graph)$label.color <- 'black'#设置节点标记的颜色
E(Graph)$width <- E(Graph)$weight/20#根据频次列设置边宽度
E(Graph)$arrow.size=0.3 #设置箭头大小
#cell cell communication plot
plot(Graph, 
     layout=layout_style,
     vertex.color="#08d9d6",
     remove.loops=F, 
     edge.color="#ff2e63")

#Heatmap
library(RColorBrewer)
library(pheatmap)
pheatmap(adj_all, 
         main = "Cell Communication Heatmap",
         display_numbers = TRUE,number_color = "Black",
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = F,
         cutree_cols = 1, 
         show_rownames = T, 
         col=brewer.pal(n = 9, name = "Reds"))