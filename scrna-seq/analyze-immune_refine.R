library(scuttle)
library(scran)
library(scater)
library(io)
library(ggplot2)

library(scDblFinder)
library(BiocParallel)

library(SingleR)
library(celldex)

library(msigdbr)
library(dplyr)

fdr.cut <- 0.05;

mc.cores <- 8;
options(mc.cores = mc.cores);

source("R/common.R");

.like <- function(s) {
	grep(s, rownames(sce), value=TRUE, ignore.case=TRUE)
}

gsets.bp <- msigdbr(species, "C5", "BP");

.gene_set <- function(keyword) {
	species <- "Mus musculus";
	gsets.bp.b.cell <- filter(gsets.bp, grepl(keyword, gs_name));
	genes.sel <- intersect(gsets.bp.b.cell$gene_symbol, rownames(sce));
}

.markers <- function(sce, n=50) {
	if (is.null(n) || nrow(sce) < n) {
		n <- nrow(sce);
	}
	markers <- scoreMarkers(sce);
	markers.mean <- lapply(markers,
		function(m) {
			m[
				order(m$mean.AUC, decreasing=TRUE)[1:n],
				c("self.average", "other.average", "mean.AUC")
			]
		}
	);
	markers.min <- lapply(markers,
		function(m) {
			m[
				order(m$min.AUC, decreasing=TRUE)[1:n],
				c("self.average", "other.average", "min.AUC")
			]
		}
	);
	list(mean = markers.mean, min = markers.min)
}


out.fn <- filename("hcc-cd133", path="immune", tag="immune");
pdf.fn <- insert(out.fn, ext="pdf");
rds.fn <- insert(out.fn, ext="rds");

sce <- qread("sce/immune_filtered.rds");
print(dim(sce))

sce.neut <- sce[, grepl("neutrophil", sce$label1)];
table(sce.neut$label1)

sce.mp <- sce[, grepl("macrophage", sce$label1)];
table(sce.mp$label1)

sce.b <- sce[, grepl("B cell", sce$label1)];
table(sce.b$label1)

sce.t <- sce[, grepl("T cell", sce$label1)];
table(sce.t$label1)

plotHighestExprs(sce.b)

sce.b <- runTSNE(sce.b);
sce.neut <- runTSNE(sce.neut);
sce.mp <- runTSNE(sce.mp);
sce.t <- runTSNE(sce.t);

plotTSNE(sce.b, colour_by="Sample", point_alpha=0.3);
plotTSNE(sce.b, colour_by="total", point_alpha=0.3);
plotTSNE(sce.b, colour_by="sizeFactor", point_alpha=0.3);
plotTSNE(sce.b, colour_by="Cd74", point_alpha=0.3);

plotTSNE(sce.t, colour_by="Sample", point_alpha=0.3);

plotTSNE(sce.mp, colour_by="Sample", point_alpha=0.3);

# ---

dec.b <- modelGeneVar(sce.b);
with(dec.b, plot(mean, total, xlab="mean log-expr", ylab="variance"));
curve(metadata(dec.b)$trend(x), col="blue", add=TRUE);

genes.b <- .gene_set("B_CELL");

# reduce data dimension using PCA
sce.b <- denoisePCA(sce.b, dec.b, subset.row=genes.b);
dim(reducedDim(sce.b, "PCA"))

sce.b <- runUMAP(sce.b, dimred="PCA");

qdraw(
	ggrastr::rasterize(
		plotUMAP(sce.b, colour_by="Sample") +
			coord_fixed()
	),
	width = 5, height = 5
)

sce.b <- runTSNE(sce.b, dimred="PCA");
plotTSNE(sce.b, colour_by="Sample", point_alpha=0.3);

set.seed(1337);
g.b <- buildSNNGraph(sce.b, use.dimred="PCA");
cl.b <- factor(igraph::cluster_leiden(g.b,
	resolution_parameter=0.1)$membership);
table(cl.b)
colLabels(sce.b) <- cl.b;

qdraw(
	ggrastr::rasterize(
		plotUMAP(sce.b, colour_by="label") +
			coord_fixed()
	),
	width = 5, height = 5
)

with(colData(sce.b), prop.table(table(Sample, label), 1))

markers.b <- .markers(sce.b);
markers.b$mean
markers.b$min

plotUMAP(sce.b, colour_by="Ighm")
plotUMAP(sce.b, colour_by="Igkc")
plotUMAP(sce.b, colour_by="Ptpn6")

# no obvious clusters among B cells

# ---


dec.t <- modelGeneVar(sce.t);
with(dec.t, plot(mean, total, xlab="mean log-expr", ylab="variance"));
curve(metadata(dec.t)$trend(x), col="blue", add=TRUE);

genes.t <- .gene_set("T_CELL");

# reduce data dimension using PCA
sce.t <- denoisePCA(sce.t, dec.t, subset.row=genes.t);
dim(reducedDim(sce.b, "PCA"))

sce.t <- runUMAP(sce.t, dimred="PCA");

qdraw(
	ggrastr::rasterize(
		plotUMAP(sce.t[, sample(1:ncol(sce.t))], colour_by="Sample") +
			coord_fixed()
	),
	width = 5, height = 5
)

set.seed(1337);
g.t <- buildSNNGraph(sce.t, use.dimred="PCA");
cl.t <- factor(igraph::cluster_leiden(g.t,
	resolution_parameter=0.1)$membership);
table(cl.t)
colLabels(sce.t) <- cl.t;

qdraw(
	ggrastr::rasterize(
		plotUMAP(sce.t[, sample(1:ncol(sce.t))], colour_by="label") +
			coord_fixed()
	),
	width = 5, height = 5
)


plotUMAP(sce.t, colour_by="Cd4")
plotUMAP(sce.t, colour_by="Cd8a")
plotUMAP(sce.t, colour_by="Klrk1")


markers.t <- .markers(sce.t);
markers.t$mean
markers.t$min

# cluster 1: Cd4+ T cell (Lef1+, Tcf7+)
data.frame(markers.t$min[[1]])
plotUMAP(sce.t, colour_by="Lef1", text_by="label")
plotUMAP(sce.t, colour_by="Tcf7", text_by="label")

# cluster 2: Cd8+ T cell (Lef1+, Tcf7+)
data.frame(markers.t$min[[2]])
data.frame(markers.t$mean[[2]])
plotUMAP(sce.t, colour_by="Cd8b1", text_by="label")
plotUMAP(sce.t, colour_by="Ccr7", text_by="label")

# cluster 3: NK cell (Ncr1+, Prf1+)
data.frame(markers.t$min[[3]])
plotUMAP(sce.t, colour_by="Ncr1", text_by="label")
plotUMAP(sce.t, colour_by="Prf1", text_by="label")

# cluster 4: Cd8+ T cell, activated (Gzmk+, Ccl5+, Nkg7+)
data.frame(markers.t$min[[4]])
data.frame(markers.t$mean[[4]])
plotUMAP(sce.t, colour_by="Cd48", text_by="label")
plotUMAP(sce.t, colour_by="Ccl5", text_by="label")
plotUMAP(sce.t, colour_by="Nkg7", text_by="label")
plotUMAP(sce.t, colour_by="Gzmk", text_by="label")

# cluster 5: Cd4+ T cell
data.frame(markers.t$min[[5]])
plotUMAP(sce.t, colour_by="Capg", text_by="label")

# cluster 6: Cd8+ T cell, activated (Cd7+, Gzmb+, Klra5+, Gzmk+, Ccl5+, Nkg7+)
data.frame(markers.t$min[[6]])
plotUMAP(sce.t, colour_by="Cd7", text_by="label")
plotUMAP(sce.t, colour_by="Gzmb", text_by="label")
plotUMAP(sce.t, colour_by="Klra5", text_by="label")

# cluster 7: NK cells (Il4+, Socs2+ Gzmb+, Ccl5+)
data.frame(markers.t$min[[7]])
data.frame(markers.t$mean[[7]])
plotUMAP(sce.t, colour_by="Il4", text_by="label")
plotUMAP(sce.t, colour_by="Socs2", text_by="label")

# add other cells
cl.t.lab <- factor(cl.t,
	levels=1:7,
	labels=c(
		"Lef1+ Cd4+ T cell",
		"Lef1+ Cd8+ T cell",
		"Ncr1+ NK cell",
		"Gzmk+ Cd8+ T cell",
		"Capg+ Cd4+ T cell",
		"Klra5+ Gzmb+ Cd8+ T cell",
		"Il4+ Socs2+ NK cell"
	)
);

sce$label.t.cell <- "non T cell";
cidx <- match(colnames(sce.t), colnames(sce));
sce$label.t.cell[cidx] <- as.character(cl.t.lab);
table(sce$label.t.cell)

enrich.t <- with(colData(sce), enrich_test(label.t.cell, Sample));
enrich.t[enrich.t$q < fdr.cut, ]

# ---


