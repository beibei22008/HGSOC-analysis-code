

allowWGCNAThreads()
pvalueFilter = 0.05  
qvalueFilter = 1  

expFile = "gene.txt"  
expro = read.table(expFile, sep = "\t", header = T, check.names = F)
dim(expro)  
expro = as.matrix(expro)      
rownames(expro) = expro[, 1]   
exp = expro[, 2:ncol(expro)]    
dimnames = list(rownames(exp), colnames(exp))   
expro = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames) 
expro = log2(expro+1)   
data = normalizeBetweenArrays(expro)    
data=expro[apply(data,1, sd)>0.8,]      
m.vars=apply(data,1,var)  
expro.upper=data[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.1))[4]),]  
dim(expro.upper)  
datExpr=as.data.frame(t(expro.upper))  
nGenes = ncol(datExpr)          
nSamples = nrow(datExpr)  
codata_all <- read.csv("codata_all.csv")
rownames(codata_all) <- codata_all[,1]
codata_all <- codata_all[,2:ncol(codata_all)]
nomo_all <- read.csv("nomo_all.csv")
rownames(nomo_all) <- nomo_all[,1]
nomo_all <- nomo_all[,2:ncol(nomo_all)]
rowname <- rownames(codata_all)
nomo_all <- nomo_all[rowname,]
nomo <- nomo_all[, 3:ncol(nomo_all)]
codata <- codata_all[, 3:ncol(codata_all)]
nomo_codata_all <- cbind(nomo, codata)
nomo_codata_all <- cbind(rownames(nomo_codata_all), nomo_codata_all)
colnames(nomo_codata_all)[1] <- "id"
write.csv(nomo_codata_all, file = "nomo_codata_all.csv", row.names = F)
trait_data <- read.csv("nomo_codata_all.csv", row.names = 1)
rowname <- intersect(rownames(trait_data), rownames(datExpr))
trait_data <- trait_data[rowname,]
datExpr <- datExpr[rowname,]
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average")
pdf("Sample_clustering.pdf")
plot(sampleTree, main = "Sample clustering"
     , sub="", xlab="",cex=0.8,font=2,lwd=1.5,cex.lab=1.5)
dev.off()
powers = c(seq(1,10,by = 1), seq(12, 20, by = 2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("Soft_Threshold(power).pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
net = blockwiseModules(datExpr, power = 5, maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AS-green-FPKM-TOM",
                       verbose = 3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
mergedColors = labels2colors(net$colors)
pdf("Cluster_Dendrogram.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
png("Cluster_Dendrogram.png")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
trait <- trait_data
modTraitCor = cor(MEsWW, trait, method  = "spearman")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
modules = MEsWW
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
png("Module_trait_ralationships_text.png", width = 1000, height = 618)
par(mar=c(10, 6, 3, 1))
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(trait), yLabels = names(MEsWW), cex.lab = 1,  yColorWidth=0.01, 
               xColorWidth = 0.03,textMatrix = textMatrix,
               ySymbols = colnames(modules), colorLabels = FALSE, colors = blueWhiteRed(50), 
               setStdMargins = FALSE, cex.text = 1, zlim = c(-.4,.4)
               , main = paste("Module-trait relationships"))
dev.off()
png("Module_trait_ralationships_without_text.png", width = 1000, height = 618)
par(mar=c(10, 6, 3, 1))
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(trait), yLabels = names(MEsWW), cex.lab = 1,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(modules), colorLabels = FALSE, colors = blueWhiteRed(50), 
               setStdMargins = FALSE, cex.text = 1, zlim = c(-.4,.4)
               , main = paste("Module-trait relationships"))
dev.off()

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
moduleTraitCor = cor(MEsWW, trait, use = "p")
moduleTraitCor = cor(MEsWW, trait, method = "spearman")
max(moduleTraitCor)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
modNames = substring(names(MEsWW), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(trait)
geneTraitSignificance = as.data.frame(cor(datExpr, trait, use = "p"))
max(geneModuleMembership)
max(geneTraitSignificance)
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")
setwd("D:\\newlassocox\\WGCNA\\WGCNAguohong\\scatter")
par(mai = c(1, 1, 1 ,1))
for(trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for(module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      outPdf=paste(trait, "_", module,".png",sep="")
      png(file=outPdf)
      #par(mfrow = c(1,1))
      par(mai = c(1, 1, 1 ,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste( module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         cex.main = 1.2, cex.lab = 2, cex.axis = 2, col = module, pch = 8,
                         font.lab = 2, font.axis = 2, lwd = 2
      )
      dev.off()
    }
  }
}

module_name <- modNames
module_name <- module_name[1:length(modNames)-1]
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)
png("Eigengene_adjacency_heatmap.png", width = 600, height = 370.8)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(4,4,4,2),yLabels = module_name,
                      xLabels = module_name,
                      plotDendrograms = F, xLabelsAngle = 60, setMargins = T)
dev.off()
for(mod in 1:nrow(table(moduleColors))){
  modules=names(table(moduleColors))[mod]
  probes=colnames(datExpr)
  inModule=(moduleColors==modules)
  moduleGenes=probes[inModule]
  write.table(moduleGenes,file = paste("WGCNA_",modules,".txt"),sep = "\t",row.names = F,col.names = F,quote = F)
}

module_names <- colnames(modules)
module_names <- module_names[1:length(module_names)-1] 
for(i in 1:length(module_names)){
  module_names[i] <- substr(module_names[i], 3,nchar(module_names[i]))
}
for(st in module_names){
  inputFile <- paste("WGCNA_ ",st," .txt", sep = "")
  rt = read.table(inputFile, sep = "\t", header = FALSE)
  genes = as.vector(rt[,1])
  entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
  entrezIDs <- as.character(entrezIDs)
  out = cbind(rt, entrezID = entrezIDs)
  colnames(out)[1] = "Gene"
  write.table(out, file = paste(st,"_id.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
}
for(st in module_name){
  input <- paste(st,"_id.txt", sep = "")
  rt = read.table(input, sep = "\t", header = TRUE, check.names = FALSE) 
  rt = rt[is.na(rt[,"entrezID"]) == FALSE,] 
  colnames(rt)[1] = "Gene"
  gene = rt$entrezID
  kk = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, qvalueCutoff = 1,
                ont = "all", readable = TRUE)
  GO = as.data.frame(kk)
  GO = GO[(GO$pvalue < pvalueFilter & GO$qvalue < qvalueFilter), ] 
  write.table(GO, file = paste(st,"_GO.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  showNum = 10
  colorSel = "qvalue"
  if(qvalueFilter > 0.05){
    colorSel = "pvalue"
  }
  bar <- barplot(kk, drop = TRUE, showCategory = showNum,cex.lab = 2,
                 split = "ONTOLOGY", color = colorSel, orient = "v",
  ) + facet_grid(ONTOLOGY~., scales = "free")
  png(paste("GO_", st, "barplot.png", sep = ""), width = 800, height = 480)
  plot(bar)
  dev.off()
  dot <- dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", 
                 split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scales = "free")
  png(paste("GO_", st, "bubplot.png", sep = ""), width = 800, height = 520)
  plot(dot)
  dev.off()
  if(st == "green"){
    rt = read.table("green_id.txt", sep = "\t", header = TRUE, check.names = FALSE) 
    rt = rt[is.na(rt[,"entrezID"]) == FALSE,]  
    colnames(rt)[1] = "Gene"
    gene = rt$entrezID
    kk = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, qvalueCutoff = 1,
                  ont = "all", readable = TRUE)
    GO = as.data.frame(kk)
    GO = GO[(GO$pvalue < pvalueFilter & GO$qvalue < qvalueFilter), ] 
    write.table(GO, file = paste(st,"_GO.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
    showNum = 10
    colorSel = "qvalue"
    if(qvalueFilter > 0.05){
      colorSel = "pvalue"
    }
    bar <- barplot(kk, drop = TRUE, showCategory = showNum,cex.lab = 2,
                   split = "ONTOLOGY", color = colorSel, orient = "v") + facet_grid(ONTOLOGY~., scales = "free")
    png(paste("GO_", st, "barplot.png", sep = ""), width = 1600, height = 480)
    plot(bar)
    dev.off()
    dot <- dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", 
                   split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scales = "free")
    png(paste("GO_", st, "bubplot.png", sep = ""), width = 1600, height = 520)
    plot(dot)
    dev.off()
  }
}

probes = colnames(datExpr)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder = order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)

for(st in module_name){
  input <- paste(st,"_id.txt", sep = "")
  rt = read.table(input, sep = "\t", header = TRUE, check.names = FALSE)
  rt = rt[is.na(rt[,"entrezID"]) == FALSE,] 
  colnames(rt)[1] = "Gene"
  gene = rt$entrezID
  kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
  KEGG = as.data.frame(kk)
  KEGG = KEGG[(KEGG$pvalue < pvalueFilter & KEGG$qvalue < qvalueFilter), ]
  write.table(KEGG, file = paste(st, "KEGG.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  colorSel = "qvalue"
  if(qvalueFilter > 0.5){
    colorSel = "pvalue"
  }
  showNum = 30
  if(nrow(KEGG) < 30){
    showNum = nrow(KEGG)
  }
  bar <- barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
  png(paste("KEGG_", st, "barplot.png", sep = ""), width = 700, height = 500)
  plot(bar)
  dev.off()
  bub <- dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", color = colorSel)
  png(paste("KEGG_", st, "bubplot.png", sep = ""), width = 700, height = 500)
  plot(bub)
  dev.off()
}

