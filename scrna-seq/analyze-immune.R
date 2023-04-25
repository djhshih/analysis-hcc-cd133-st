library(scuttle)
library(scran)
library(scater)
library(io)
library(ggplot2)

library(scDblFinder)
library(BiocParallel)

library(SingleR)
library(celldex)


mc.cores <- 8;
options(mc.cores = mc.cores);

source("R/common.R");

.like <- function(s) {
	grep(s, rownames(sce), value=TRUE, ignore.case=TRUE)
}

out.fn <- filename("hcc-cd133", path="immune", tag="immune");
pdf.fn <- insert(out.fn, ext="pdf");
rds.fn <- insert(out.fn, ext="rds");


sce <- qread("sce/immune.rds");
print(dim(sce))

rownames(sce) <- rowData(sce)$Symbol;

sce <- logNormCounts(sce);

sce <- scDblFinder(sce, samples="Sample", BPPARAM=MulticoreParam(mc.cores));
table(sce$Sample, sce$scDblFinder.class)

set.seed(1337);
cl <- quickCluster(sce);
table(cl)

ggplot(data.frame(colData(sce), cl), aes(x=cl, y=total)) +
	theme_classic() +
	geom_violin() +
	scale_y_log10()

ggplot(data.frame(colData(sce), cl), aes(x=cl, y=subsets_mito_percent)) +
	theme_classic() +
	geom_violin() +
	scale_y_log10()

colLabels(sce) <- cl;

dec <- modelGeneVar(sce);
with(dec, plot(mean, total, xlab="mean log-expr", ylab="variance"));
curve(metadata(dec)$trend(x), col="blue", add=TRUE);

sce <- runTSNE(sce);

cidx <- sample(1:ncol(sce));
plotTSNE(sce[, cidx], colour_by="label", text_by="label", point_alpha=0.3);

plotTSNE(sce[, cidx], colour_by="total", text_by="label", point_alpha=0.3);
plotTSNE(sce[, cidx], colour_by="sizeFactor", text_by="label");
plotTSNE(sce[, cidx], colour_by="subsets_mito_percent", text_by="label", point_alpha=0.3);
plotTSNE(sce[, cidx], colour_by="Sample", text_by="label", point_alpha=0.3);
plotTSNE(sce[, cidx], colour_by="scDblFinder.score", text_by="label", point_alpha=0.3);
plotTSNE(sce[, cidx], colour_by="scDblFinder.class", text_by="label", point_alpha=0.3);

sce <- runUMAP(sce);
plotUMAP(sce[, cidx], colour_by="label", text_by="label", point_alpha=0.3);
plotUMAP(sce[, cidx], colour_by="sizeFactor", text_by="label");
plotUMAP(sce[, cidx], colour_by="scDblFinder.class", text_by="label", point_alpha=0.3);

# identify cluster markers

markers <- scoreMarkers(sce);

grep("mt-", rownames(sce), value=TRUE)

markers.sel <- lapply(markers,
	function(m) {
		d <- m[order(m$mean.AUC, decreasing=TRUE)[1:50], c("self.average", "other.average", "mean.AUC")];
		#rownames(d) <- rowData(sce)$Symbol[match(rownames(d), rowData(sce)$ID)]
		d
	}
);

plotTSNE(sce, colour_by = "Ms4a4b")
plotTSNE(sce, colour_by = "Nkg7")

with(colData(sce), plot(total, scDblFinder.score))

plotUMAP(sce[, cidx], colour_by="scDblFinder.score", text_by="label", point_alpha=0.3);
plotUMAP(sce[, cidx], colour_by="sizeFactor", text_by="label", point_alpha=0.3);
plotUMAP(sce[, cidx], colour_by="subsets_mito_percent", text_by="label", point_alpha=0.3);

markers.sel[[4]]   # main B cell cluster
markers.sel[[8]]   # T cell cluster
markers.sel[[19]]  # T cell cluster
markers.sel[[11]]  # T cell cluster (activated?)
markers.sel[[2]]   # macrophages
markers.sel[[10]]  # macrophages
markers.sel[[21]]  # hepatocytes
markers.sel[[5]]   # macrophages?

markers.sel[[15]]  # dying cells?
markers.sel[[20]]  # dying cells?

markers.sel[[1]]   # T cells?
markers.sel[[7]]   # monocytes/macrophages?
markers.sel[[16]]  # B cells + monocytes?
markers.sel[[17]]  # T + B cells
markers.sel[[22]]  # B cells + monocytes?
markers.sel[[24]]  # heptocytes?

plotUMAP(sce, colour_by = "Cd79a", text_by = "label")   # B cells
plotUMAP(sce, colour_by = "Cd3g", text_by = "label")    # T / NK cells
plotUMAP(sce, colour_by = "Trac", text_by = "label")    # T cells
plotUMAP(sce, colour_by = "Cd4", text_by = "label")     # CD 4 T cells
plotUMAP(sce, colour_by = "Cd8a", text_by = "label")    # CD8 T cells
plotUMAP(sce, colour_by = "Cd14", text_by = "label")    # monocytes
plotUMAP(sce, colour_by = "Fcgr3", text_by = "label")   # monocytes / NK cells
plotUMAP(sce, colour_by = "Ncam1", text_by = "label")   # NK cells
plotUMAP(sce, colour_by = "Itgam", text_by = "label")   # monocytes
plotUMAP(sce, colour_by = "Lyz2", text_by = "label")    # monocytes/macrophages
plotUMAP(sce, colour_by = "Itgax", text_by = "label")   # dendritic cells

plotUMAP(sce, colour_by = "Mpeg1", text_by = "label")   # macrophages
plotUMAP(sce, colour_by = "Cybb", text_by = "label")    # macrophages

# B cell

plotTSNE(sce, colour_by = "Cd79a", text_by = "label")
plotTSNE(sce, colour_by = "Cd79b")
plotTSNE(sce, colour_by = "Ighd")
plotTSNE(sce, colour_by = "Ighm")

# T cell / NK cell

plotTSNE(sce, colour_by = "Cd3g")
plotTSNE(sce, colour_by = "Cd3d")
plotTSNE(sce, colour_by = "Cd3e")

# T cell

plotTSNE(sce, colour_by = "Trac")
plotTSNE(sce, colour_by = "Trbc1")
plotTSNE(sce, colour_by = "Cd4")
plotTSNE(sce, colour_by = "Cd8a")
plotTSNE(sce, colour_by = "Cd8b1")

# NK cell

plotTSNE(sce, colour_by = "Ncr1")

plotTSNE(sce, colour_by = "Klrc1")
plotTSNE(sce, colour_by = "Klrc2")

# Monocyte / NK cell

plotTSNE(sce, colour_by = "Fcgr3")

# Dendritic

plotTSNE(sce, colour_by = "Cd14")
plotTSNE(sce, colour_by = "Cd68")

# Dendritic / Granulocytes / Monocyte

plotTSNE(sce, colour_by = "Itgam")

plotTSNE(sce, colour_by = "Ceacam1")


plotHighestExprs(sce)

# ---

# remove dying clusters
cl.dying <- c("15", "20");
cidx.dying <- sce$label %in% cl.dying;

# remove likely doublets
cidx.doublet <- sce$scDblFinder.class == "doublet";

# remove rRNA contamination
# Gm42418 and AY036118 (overlaps with Rn45s)
# https://doi.org/10.1242/dev.183855
gidx.ribo <- rownames(sce) %in% c("Gm42418", "AY036118");

sce.f <- sce[!gidx.ribo, !(cidx.dying | cidx.doublet)];

dim(sce)
dim(sce.f)

dec.f <- modelGeneVar(sce.f);
with(dec.f, plot(mean, total, xlab="mean log-expr", ylab="variance"));
curve(metadata(dec.f)$trend(x), col="blue", add=TRUE);

# top highly variable genes
top.hvgs <- getTopHVGs(dec.f, prop=0.3);

# reduce data dimension using PCA
sce.f <- denoisePCA(sce.f, dec.f, subset.row=top.hvgs, min.rank=10);

set.seed(1337);
g <- buildSNNGraph(sce.f, use.dimred="PCA");
cl.f <- factor(igraph::cluster_walktrap(g)$membership)
table(cl.f)
colLabels(sce.f) <- cl.f;

ggplot(data.frame(colData(sce.f)), aes(x=label, y=total)) +
	theme_classic() +
	geom_violin() +
	scale_y_log10()

ggplot(data.frame(colData(sce.f)), aes(x=label, y=scDblFinder.score)) +
	theme_classic() +
	geom_violin() +
	scale_y_log10()

ggplot(data.frame(colData(sce.f)), aes(x=label, y=subsets_mito_percent)) +
	theme_classic() +
	geom_violin() +
	scale_y_log10()

plotTSNE(sce.f, colour_by="label", text_by="label", point_alpha=0.3);
plotTSNE(sce.f[, cidx.f], colour_by="Sample", text_by="label", point_alpha=0.3);
plotTSNE(sce.f, colour_by="scDblFinder.score", text_by="label", point_alpha=0.3);

# ---

# remove remaining doublet cluster
cidx.doublet2 <- cl.f %in% c("20", "26")
table(cidx.doublet2)

sce.f2 <- sce.f[, !cidx.doublet2];

# ---

dec.f2 <- modelGeneVar(sce.f2);
with(dec.f2, plot(mean, total, xlab="mean log-expr", ylab="variance"));
curve(metadata(dec.f2)$trend(x), col="blue", add=TRUE);

# top highly variable genes
top.hvgs <- getTopHVGs(dec.f2, prop=0.3);

# reduce data dimension using PCA
sce.f2 <- denoisePCA(sce.f2, dec.f2, subset.row=top.hvgs, min.rank=5);

set.seed(1337);
g <- buildSNNGraph(sce.f2, use.dimred="PCA");
cl.f <- factor(igraph::cluster_leiden(g, resolution_parameter=0.01)$membership);
table(cl.f)
colLabels(sce.f2) <- cl.f;

sce.f2 <- runTSNE(sce.f2);
#sce.f <- runTSNE(sce.f, dimred="PCA");

sce.f2 <- runUMAP(sce.f2);

cidx.f2 <- sample(1:ncol(sce.f2));
table(sce.f2$Sample)

plotTSNE(sce.f2, colour_by="label", text_by="label", point_alpha=0.3);
plotTSNE(sce.f2[, cidx.f2], colour_by="Sample", text_by="label", point_alpha=0.3);
plotTSNE(sce.f2, colour_by="scDblFinder.score", text_by="label", point_alpha=0.3);

plotUMAP(sce.f2, colour_by="label", text_by="label", point_alpha=0.3);
plotUMAP(sce.f2[, cidx.f], colour_by="Sample", text_by="label", point_alpha=0.3);
plotUMAP(sce.f2, colour_by="scDblFinder.score", text_by="label", point_alpha=0.3);

plotTSNE(sce.f2, colour_by = "Cd79a", text_by = "label")   # B cells
plotTSNE(sce.f2, colour_by = "Ighm", text_by = "label")    # B cells
plotTSNE(sce.f2, colour_by = "Igkc", text_by = "label")    # B cells
plotTSNE(sce.f2, colour_by = "Cd3g", text_by = "label")    # T / NK cells
plotTSNE(sce.f2, colour_by = "Trac", text_by = "label")    # T cells
plotTSNE(sce.f2, colour_by = "Cd14", text_by = "label")    # monocytes
plotTSNE(sce.f2, colour_by = "Fcgr3", text_by = "label")   # monocytes / NK cells
plotTSNE(sce.f2, colour_by = "Klrk1", text_by = "label")   # NK cells
plotTSNE(sce.f2, colour_by = "Ncam1", text_by = "label")   # NK cells
plotTSNE(sce.f2, colour_by = "Itgam", text_by = "label")   # monocytes / NK / neutrophils
plotTSNE(sce.f2, colour_by = "Lyz2", text_by = "label")    # monocytes/macrophages
plotTSNE(sce.f2, colour_by = "Itgax", text_by = "label")   # dendritic cells

plotTSNE(sce.f2, colour_by = "Aldob", text_by = "label")    # hepatocytes


markers.f2 <- scoreMarkers(sce.f2);

markers.f2.sel <- lapply(markers.f2,
	function(m) {
		m[
			order(m$mean.AUC, decreasing=TRUE)[1:50],
			c("self.average", "other.average", "mean.AUC")
		]
	}
);

markers.f2.sel.min <- lapply(markers.f2,
	function(m) {
		m[
			order(m$min.AUC, decreasing=TRUE)[1:50],
			c("self.average", "other.average", "min.AUC")
		]
	}
);


# Neat1+ Sorl1+ Lyst1+ activated neutrophils
markers.f2.sel[[1]]
markers.f2.sel.min[[1]]
plotTSNE(sce.f2, colour_by = "Malat1", text_by = "label")
plotTSNE(sce.f2, colour_by = "Neat1", text_by = "label")
plotTSNE(sce.f2, colour_by = "Sorl1", text_by = "label")
plotTSNE(sce.f2, colour_by = "Lyst", text_by = "label")
plotTSNE(sce.f2, colour_by = "Nlrp3", text_by = "label")  # inflammasome

markers.f2.sel[[6]]  # macrophages
plotTSNE(sce.f2, colour_by = "Lyz2", text_by = "label")

# a subcluster
plotTSNE(sce.f2, colour_by = "Ly86", text_by = "label")
plotTSNE(sce.f2, colour_by = "Cd74", text_by = "label")

markers.f2.sel[[5]]   # B cells

markers.f2.sel[[3]]   # T cells

markers.f2.sel[[11]]  # Mpeg1-high activated macrophages
plotTSNE(sce.f2, colour_by = "Mpeg1", text_by = "label")
plotTSNE(sce.f2, colour_by = "Psap", text_by = "label")
plotTSNE(sce.f2, colour_by = "Cybb", text_by = "label")

plotTSNE(sce.f2, colour_by = "Msr1", text_by = "label")
plotTSNE(sce.f2, colour_by = "Prf1", text_by = "label")

# markers.f2.sel[[12]]  # activated T cells
# markers.f2.sel.min[[12]]  # Lef1+ Tcf7+ Il2ra- Foxp3- T cells
# plotTSNE(sce.f2, colour_by = "Lef1", text_by = "label")
# plotTSNE(sce.f2, colour_by = "Tcf7", text_by = "label")
# plotTSNE(sce.f2, colour_by = "Il7r", text_by = "label")
# plotExpression(sce.f2, colour_by = "Lef1")
# plotTSNE(sce.f2, colour_by = "Foxp3", text_by = "label")
# plotTSNE(sce.f2, colour_by = "Cd4", text_by = "label")
# plotTSNE(sce.f2, colour_by = "Il2ra", text_by = "label")

markers.f2.sel[[2]]  # neutrophils
markers.f2.sel[[7]]  # neutrophils
markers.f2.sel[[8]]  # neutrophils
markers.f2.sel[[9]]  # neutrophils
plotTSNE(sce.f2, colour_by = "Srgn", text_by = "label")
plotTSNE(sce.f2, colour_by = "Il1b", text_by = "label")

markers.f2.sel[[10]]  # hepatocytes


markers.f2.sel[[4]] # B1 cells
plotTSNE(sce.f2, colour_by = "Mzb1", text_by = "label")

with(colData(sce.f2), prop.table(table(Sample, label), margin=1))


# ---

# reference-based cell type prediction

ref <- MouseRNAseqData();
pred <- SingleR(test = sce.f2, ref=ref, labels=ref$label.main);
table(pred$labels)
table(pred$labels, sce.f2$label);

sce.f2$label.ref.pred <- pred$labels;

# predicted labels largely agree with manual assignments
plotTSNE(sce.f2, colour_by = "label.ref.pred", text_by = "label")

# ---

sce.f2$label1 <- NA;
sce.f2$label1[sce.f2$label == "1"] <- "Neat1+ neutrophil";
sce.f2$label1[sce.f2$label == "2"] <- "neutrophil";
sce.f2$label1[sce.f2$label == "3"] <- "NK / T cell";
sce.f2$label1[sce.f2$label == "4"] <- "B1 cell";
sce.f2$label1[sce.f2$label == "5"] <- "B cell";
sce.f2$label1[sce.f2$label == "6"] <- "macrophage";
sce.f2$label1[sce.f2$label == "7"] <- "neutrophil";
sce.f2$label1[sce.f2$label == "8"] <- "neutrophil";
sce.f2$label1[sce.f2$label == "9"] <- "neutrophil";
sce.f2$label1[sce.f2$label == "10"] <- "hepatocyte";
sce.f2$label1[sce.f2$label == "11"] <- "Mpeg1-high macrophage";
sce.f2$label1 <- factor(sce.f2$label1);

table(pred$labels, sce.f2$label1);

qdraw(
	ggrastr::rasterize(
		plotTSNE(sce.f2, colour_by = "label1", text_by = "label1") +
		theme(legend.position="none") +
		coord_fixed()
	),
	width = 8, height = 8,
	file = insert(pdf.fn, c("tsne", "clusters"))
)

cidx.f2 <- sample(1:ncol(sce.f2));
qdraw(
	ggrastr::rasterize(
		plotTSNE(sce.f2[, cidx.f2], colour_by = "Sample", text_by = "label1") +
		coord_fixed()
	),
	width = 8, height = 8,
	file = insert(pdf.fn, c("tsne", "sample"))
)


.draw_gene <- function(gene) {
	qdraw(
		ggrastr::rasterize(
			plotTSNE(sce.f2, colour_by = gene) + coord_fixed()
		),
		width = 8, height = 8,
		file = insert(pdf.fn, c("gene", tolower(gene)))
	)
}

# Mpeg1 / Perforin-2
.draw_gene("Mpeg1")

# Msr1 / CD204
.draw_gene("Msr1")

# M2 macrophage
.draw_gene("Il10")

.draw_gene("Ighm")  # B cell
.draw_gene("Trac")  # T cell
.draw_gene("Klrk1") # NK cell
.draw_gene("Lyz2")  # macrophages
.draw_gene("Itgam") # macrophages / neutrophils
.draw_gene("Neat1")


# ---

with(colData(sce.f2), enrich_test(label1, Sample))


# depleted in Prom1-DTA samples:
# cluster 11 (Mpeg1-high macrophages)
# cluster 5 (B cells)
# cluster 3 (T cells)
# cluster 6 (macrophages)
enrich.d[enrich.d$adj.p < 0.05, ]


d.label1 <- proportions(sce.f2$Sample, sce.f2$label1);
d.label1$group <- sub("-IMM", "", d.label1$group);
d.label1$group <- sub("Prom1-", "", d.label1$group);
d.label1$value[d.label1$value == "Mpeg1-high macrophage"] <- "Mpeg1-high MP";

d.sub <- d.label1[!grepl("neutrophil", d.label1$value), ];

qdraw(
	ggplot(d.sub, aes(x=group, y=mean, ymin=lower, ymax=upper, fill=value)) +
		theme_classic() + theme(strip.background = element_blank()) +
		geom_col() + 
		geom_errorbar(width=0.3, show.legend=FALSE) +
		facet_wrap(~ value, nrow=1) +
		scale_y_continuous(n.breaks=3) +
		theme(legend.position="none") +
		scale_fill_nejm() +
		xlab("") + ylab("proportion")
	,
	width = 7, height = 3,
	file = insert(pdf.fn, c("cell-type-prop"))
);

d.neut <- d.label1[grepl("neutrophil", d.label1$value), ];

qdraw(
	ggplot(d.neut, aes(x=group, y=mean, ymin=lower, ymax=upper, fill=value)) +
		theme_classic() + theme(strip.background = element_blank()) +
		geom_col() + 
		geom_errorbar(width=0.3, show.legend=FALSE) +
		facet_wrap(~ value, nrow=1) +
		scale_y_continuous(n.breaks=3) +
		theme(legend.position="none") +
		scale_fill_nejm() +
		xlab("") + ylab("proportion")
	,
	width = 3, height = 3,
	file = insert(pdf.fn, c("cell-type-prop", "neutrophil"))
);

qwrite(sce.f2, "sce/immune_filtered.rds");

