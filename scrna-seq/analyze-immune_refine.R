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

# ---

.like <- function(s) {
	grep(s, rownames(sce), value=TRUE, ignore.case=TRUE)
}

species <- "Mus musculus";
gsets.bp <- msigdbr(species, "C5", "BP");

.gene_set <- function(keyword) {
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

.draw_umap <- function(sce, width = 5, height = 5, file = NULL, ...) {
	cidx <- sample(1:ncol(sce));
	qdraw(
		ggrastr::rasterize(
			plotUMAP(sce[, cidx], ...) + coord_fixed()
		),		
		width = width, height = height, file = file
	);
}

# ---

in.fn <- "sce/immune_filtered.rds";

out.fn <- filename("hcc-cd133", path="immune", tag="immune");
pdf.fn <- insert(out.fn, ext="pdf");
rds.fn <- insert(out.fn, ext="rds");
csv.fn <- insert(out.fn, ext="csv");

sce <- qread(in.fn);
print(dim(sce))

sce.b <- sce[, grepl("B cell", sce$label1)];
table(sce.b$label1)

sce.t <- sce[, grepl("T cell", sce$label1)];
table(sce.t$label1)

sce.mp <- sce[, grepl("macrophage", sce$label1)];
table(sce.mp$label1)

sce.neut <- sce[, sce$label1 == "neutrophil"];
table(sce.neut$label1)

# ---

# sce.b <- runTSNE(sce.b);
# sce.t <- runTSNE(sce.t);
# sce.mp <- runTSNE(sce.mp);
# sce.neut <- runTSNE(sce.neut);

# ---

# Analyze B cells

plotHighestExprs(sce.b)

plotTSNE(sce.b, colour_by="Sample", point_alpha=0.3);
plotTSNE(sce.b, colour_by="total", point_alpha=0.3);
plotTSNE(sce.b, colour_by="sizeFactor", point_alpha=0.3);
plotTSNE(sce.b, colour_by="Cd74", point_alpha=0.3);


dec.b <- modelGeneVar(sce.b);
with(dec.b, plot(mean, total, xlab="mean log-expr", ylab="variance"));
curve(metadata(dec.b)$trend(x), col="blue", add=TRUE);

genes.b <- .gene_set("B_CELL");

# reduce data dimension using PCA
sce.b <- denoisePCA(sce.b, dec.b, subset.row=genes.b);
dim(reducedDim(sce.b, "PCA"))

sce.b <- runUMAP(sce.b, dimred="PCA");

.draw_umap(
	sce.b, colour_by="Sample",
	file = insert(pdf.fn, c("b-cell", "umap", "sample"))
)

set.seed(1337);
g.b <- buildSNNGraph(sce.b, use.dimred="PCA");
cl.b <- factor(igraph::cluster_leiden(g.b,
	resolution_parameter=0.1)$membership);
table(cl.b)
colLabels(sce.b) <- cl.b;

.draw_umap(
	sce.b, colour_by="label",
	file = insert(pdf.fn, c("b-cell", "umap", "clusters"))
)

with(colData(sce.b), prop.table(table(Sample, label), 1))

markers.b <- .markers(sce.b);
markers.b$mean
markers.b$min

plotUMAP(sce.b, colour_by="Ighm")
plotUMAP(sce.b, colour_by="Igkc")

.draw_umap(
	sce.b, colour_by="Ptpn6",
	file = insert(pdf.fn, c("b-cell", "umap", "ptpn6"))
)

# no obvious clusters among B cells

# ---

# Analyze T cells

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

.draw_umap(
	sce.t, colour_by="Sample",
	file = insert(pdf.fn, c("t-cell", "umap", "sample"))
);

set.seed(1337);
g.t <- buildSNNGraph(sce.t, use.dimred="PCA");
cl.t <- factor(igraph::cluster_leiden(g.t,
	resolution_parameter=0.1)$membership);
table(cl.t)
colLabels(sce.t) <- cl.t;

.draw_umap(
	sce.t, colour_by="label",
	file = insert(pdf.fn, c("t-cell", "umap", "clusters"))
);

.draw_umap(
	sce.t, colour_by="Cd4",
	file = insert(pdf.fn, c("t-cell", "umap", "cd4"))
);

.draw_umap(
	sce.t, colour_by="Cd8a",
	file = insert(pdf.fn, c("t-cell", "umap", "cd8a"))
);

.draw_umap(
	sce.t, colour_by="Cd8b1",
	file = insert(pdf.fn, c("t-cell", "umap", "cd8b1"))
);

.draw_umap(
	sce.t, colour_by="Klrk1",
	file = insert(pdf.fn, c("t-cell", "umap", "klrk1"))
);

markers.t <- .markers(sce.t);
markers.t$mean
markers.t$min

# cluster 1: Cd4+ T cell (Lef1+, Tcf7+)
data.frame(markers.t$min[[1]])
.draw_umap(
	sce.t, colour_by="Lef1",
	file = insert(pdf.fn, c("t-cell", "umap", "lef1"))
);
.draw_umap(
	sce.t, colour_by="Tcf7",
	file = insert(pdf.fn, c("t-cell", "umap", "tcf7"))
);

# cluster 2: Cd8+ T cell (Lef1+, Tcf7+)
data.frame(markers.t$min[[2]])
data.frame(markers.t$mean[[2]])
.draw_umap(
	sce.t, colour_by="Ccr7",
	file = insert(pdf.fn, c("t-cell", "umap", "ccr7"))
);

# cluster 3: NK cell (Ncr1+, Prf1+)
data.frame(markers.t$min[[3]])
.draw_umap(
	sce.t, colour_by="Ncr1",
	file = insert(pdf.fn, c("t-cell", "umap", "ncr1"))
);
.draw_umap(
	sce.t, colour_by="Prf1",
	file = insert(pdf.fn, c("t-cell", "umap", "prf1"))
);

# cluster 4: Cd8+ T cell, activated (Gzmk+, Ccl5+, Nkg7+)
data.frame(markers.t$min[[4]])
data.frame(markers.t$mean[[4]])
.draw_umap(
	sce.t, colour_by="Cd48",
	file = insert(pdf.fn, c("t-cell", "umap", "cd48"))
);
.draw_umap(
	sce.t, colour_by="Ccl5",
	file = insert(pdf.fn, c("t-cell", "umap", "ccl5"))
);
.draw_umap(
	sce.t, colour_by="Nkg7",
	file = insert(pdf.fn, c("t-cell", "umap", "nkg7"))
);
.draw_umap(
	sce.t, colour_by="Gzmk",
	file = insert(pdf.fn, c("t-cell", "umap", "gzmk"))
);

# cluster 5: Cd4+ T cell
data.frame(markers.t$min[[5]])
.draw_umap(
	sce.t, colour_by="Capg",
	file = insert(pdf.fn, c("t-cell", "umap", "capg"))
);

# cluster 6: Cd8+ T cell, activated (Cd7+, Gzmb+, Klra5+, Gzmk+, Ccl5+, Nkg7+)
data.frame(markers.t$min[[6]])
.draw_umap(
	sce.t, colour_by="Cd7",
	file = insert(pdf.fn, c("t-cell", "umap", "cd7"))
);
.draw_umap(
	sce.t, colour_by="Gzmb",
	file = insert(pdf.fn, c("t-cell", "umap", "gzmb"))
);
.draw_umap(
	sce.t, colour_by="Klra5",
	file = insert(pdf.fn, c("t-cell", "umap", "klra5"))
);

# cluster 7: NK cells (Il4+, Socs2+ Gzmb+, Ccl5+)
data.frame(markers.t$min[[7]])
data.frame(markers.t$mean[[7]])
.draw_umap(
	sce.t, colour_by="Il4",
	file = insert(pdf.fn, c("t-cell", "umap", "il4"))
);
.draw_umap(
	sce.t, colour_by="Socs2",
	file = insert(pdf.fn, c("t-cell", "umap", "socs2"))
);

# add other cells
cl.t.lab <- factor(cl.t,
	levels=1:7,
	labels=c(
		"Lef1+ Cd4+ T cell",
		"Lef1+ Cd8+ T cell",
		"Ncr1+ NK cell",
		"Gzmk+ Cd8+ T cell",
		"Capg+ Cd4+ T cell",
		"Klra5+ Cd8+ T cell",
		"Il4+ NK cell"
	)
);

sce.t$label <- cl.t.lab
.draw_umap(
	sce.t, colour_by="label",
	file = insert(pdf.fn, c("t-cell", "umap", "clusters"))
);

sce$label.t.cell <- "non-T cell";
cidx <- match(colnames(sce.t), colnames(sce));
sce$label.t.cell[cidx] <- as.character(cl.t.lab);
table(sce$label.t.cell)

enrich.t <- with(colData(sce), enrich_test(label.t.cell, Sample));
enrich.t[enrich.t$q < fdr.cut, ]

# ---

# Analyze macrophages

dec.mp <- modelGeneVar(sce.mp);
with(dec.mp, plot(mean, total, xlab="mean log-expr", ylab="variance"));
curve(metadata(dec.mp)$trend(x), col="blue", add=TRUE);

genes.mp <- .gene_set("MACROPHAGE");

# reduce data dimension using PCA
sce.mp <- denoisePCA(sce.mp, dec.mp, subset.row=genes.mp);
dim(reducedDim(sce.b, "PCA"))

sce.mp <- runUMAP(sce.mp, dimred="PCA");

.draw_umap(
	sce.mp, colour_by="Sample",
	file = insert(pdf.fn, c("mp", "umap", "sample"))
);

set.seed(1337);
g.mp <- buildSNNGraph(sce.mp, use.dimred="PCA");
cl.mp <- factor(igraph::cluster_leiden(g.mp,
 	resolution_parameter=0.03)$membership);
table(cl.mp)
colLabels(sce.mp) <- cl.mp;

.draw_umap(
	sce.mp, colour_by="label",
	file = insert(pdf.fn, c("mp", "umap", "clusters"))
);

.draw_umap(
	sce.mp, colour_by="Cd14",
	file = insert(pdf.fn, c("mp", "umap", "cd14"))
);

.draw_umap(
	sce.mp, colour_by="Fcgr3",
	file = insert(pdf.fn, c("mp", "umap", "fcgr3"))
);

markers.mp <- .markers(sce.mp);
markers.mp$mean
markers.mp$min

# cluster 1: S100a6-high macrophage
data.frame(markers.mp$min[[1]])
.draw_umap(
	sce.mp, colour_by="S100a6",
	file = insert(pdf.fn, c("mp", "umap", "s100a6"))
);
.draw_umap(
	sce.mp, colour_by="S100a4",
	file = insert(pdf.fn, c("mp", "umap", "s100a4"))
);
.draw_umap(
	sce.mp, colour_by="S100a11",
	file = insert(pdf.fn, c("mp", "umap", "s100a11"))
);
.draw_umap(
	sce.mp, colour_by="Msrb1",
	file = insert(pdf.fn, c("mp", "umap", "msrb1"))
);

# cluster 2: Cd74+ C1q+ macrophage
data.frame(markers.mp$min[[2]])
.draw_umap(
	sce.mp, colour_by="Cd74",
	file = insert(pdf.fn, c("mp", "umap", "cd74"))
);
.draw_umap(
	sce.mp, colour_by="C1qb",
	file = insert(pdf.fn, c("mp", "umap", "c1qb"))
);
.draw_umap(
	sce.mp, colour_by="C1qc",
	file = insert(pdf.fn, c("mp", "umap", "c1qc"))
);

# cluster 3: Mpeg1-high Msr1-high macrophage
data.frame(markers.mp$min[[3]])
.draw_umap(
	sce.mp, colour_by="Mpeg1",
	file = insert(pdf.fn, c("mp", "umap", "mpeg1"))
);
.draw_umap(
	sce.mp, colour_by="Msr1",
	file = insert(pdf.fn, c("mp", "umap", "msr1"))
);
.draw_umap(
	sce.mp, colour_by="Ptprc",
	file = insert(pdf.fn, c("mp", "umap", "ptprc"))
);
.draw_umap(
	sce.mp, colour_by="Itgam",
	file = insert(pdf.fn, c("mp", "umap", "itgam"))
);
.draw_umap(
	sce.mp, colour_by="Cybb",
	file = insert(pdf.fn, c("mp", "umap", "cybb"))
);

# add other cells
cl.mp.lab <- factor(cl.mp,
	levels=1:3,
	labels=c(
		"S100a6-high macrophage",
		"Cd74+ C1q+ macrophage",
		"Mpeg1-high macrophage"
	)
);

sce.mp$label <- cl.mp.lab;
.draw_umap(
	sce.mp, colour_by="label",
	file = insert(pdf.fn, c("mp", "umap", "clusters"))
);

sce$label.mp.cell <- "non-macrophage";
cidx <- match(colnames(sce.mp), colnames(sce));
sce$label.mp.cell[cidx] <- as.character(cl.mp.lab);
table(sce$label.mp.cell)

enrich.mp <- with(colData(sce), enrich_test(label.mp.cell, Sample));
enrich.mp[enrich.mp$q < fdr.cut, ]

# ---

# Analyze neutrophils

dec.neut <- modelGeneVar(sce.neut);
with(dec.neut, plot(mean, total, xlab="mean log-expr", ylab="variance"));
curve(metadata(dec.neut)$trend(x), col="blue", add=TRUE);

genes.neut <- .gene_set("NEUTROPHIL");

# reduce data dimension using PCA
sce.neut <- denoisePCA(sce.neut, dec.neut, subset.row=genes.neut);
dim(reducedDim(sce.b, "PCA"))

sce.neut <- runUMAP(sce.neut, dimred="PCA");

.draw_umap(
	sce.neut, colour_by="Sample",
	file = insert(pdf.fn, c("neut", "umap", "sample"))
);

set.seed(1337);
g.neut <- buildSNNGraph(sce.neut, use.dimred="PCA");
cl.neut <- factor(igraph::cluster_leiden(g.neut,
 	resolution_parameter=0.005)$membership);
table(cl.neut)
colLabels(sce.neut) <- cl.neut;

.draw_umap(
	sce.neut, colour_by="label",
	file = insert(pdf.fn, c("neut", "umap", "clusters"))
);

markers.neut <- .markers(sce.neut);
markers.neut$mean
markers.neut$min

# cluster 1: Ccl6+ Gngt2- Cxcl2+ neutrophil
data.frame(markers.neut$min[[1]])
.draw_umap(
	sce.neut, colour_by="S100a9",
	file = insert(pdf.fn, c("mp", "umap", "s100a9"))
);
.draw_umap(
	sce.neut, colour_by="S100a8",
	file = insert(pdf.fn, c("mp", "umap", "s100a8"))
);
.draw_umap(
	sce.neut, colour_by="S100a11",
	file = insert(pdf.fn, c("mp", "umap", "s100a11"))
);
.draw_umap(
	sce.neut, colour_by="Ccl6",
	file = insert(pdf.fn, c("mp", "umap", "ccl6"))
);

# cluster 2: Gngt2+ Cxcl2+ neutrophil
data.frame(markers.neut$min[[2]])
data.frame(markers.neut$mean[[2]])
.draw_umap(
	sce.neut, colour_by="Cxcl2",
	file = insert(pdf.fn, c("mp", "umap", "cxcl2"))
);
.draw_umap(
	sce.neut, colour_by="Il1b",
	file = insert(pdf.fn, c("mp", "umap", "il1b"))
);

# cluster 3: Gngt2+ Cxcl2- neutrophil
data.frame(markers.neut$min[[3]])
data.frame(markers.neut$mean[[3]])
.draw_umap(
	sce.neut, colour_by="Gngt2",
	file = insert(pdf.fn, c("mp", "umap", "gngt2"))
);
.draw_umap(
	sce.neut, colour_by="Ppia",
	file = insert(pdf.fn, c("mp", "umap", "ppia"))
);
.draw_umap(
	sce.neut, colour_by="Cst3",
	file = insert(pdf.fn, c("mp", "umap", "cst3"))
);

# add other cells
cl.neut.lab <- factor(cl.neut,
	levels=1:3,
	labels=c(
		"Gngt2- Cxcl2+ neutrophil",
		"Gngt2+ Cxcl2+ neutrophil",
		"Gngt2+ Cxcl2- neutrophil"
	)
);

sce.neut$label <- cl.neut.lab;
.draw_umap(
	sce.neut, colour_by="label",
	file = insert(pdf.fn, c("neut", "umap", "clusters"))
);

sce$label.neut.cell <- "non-neutrophil";
cidx <- match(colnames(sce.neut), colnames(sce));
sce$label.neut.cell[cidx] <- as.character(cl.neut.lab);
table(sce$label.neut.cell)

enrich.neut <- with(colData(sce), enrich_test(label.neut.cell, Sample));
enrich.neut

# ---

# Identify enriched or depleted cell subtypes in treatment vs. control

sce$label2 <- as.character(sce$label1);
cidx <- match(colnames(sce.neut), colnames(sce));
sce$label2[cidx] <- as.character(cl.neut.lab);
cidx <- match(colnames(sce.mp), colnames(sce));
sce$label2[cidx] <- as.character(cl.mp.lab);
cidx <- match(colnames(sce.t), colnames(sce));
sce$label2[cidx] <- as.character(cl.t.lab);

table(sce$label2)

enrich.all <- with(colData(sce), enrich_test(label2, Sample));
enrich.all <- dplyr::rename(enrich.all, group = value);

enrich.all[enrich.all$q < fdr.cut, ]

qwrite(enrich.all, insert(csv.fn, "enrich"));

# ---

# enrichment plots

d.label2 <- proportions(sce$Sample, sce$label2);
d.label2$group <- sub("-IMM", "", d.label2$group);
d.label2$group <- sub("Prom1-", "", d.label2$group);
d.label2$group <- relevel(factor(d.label2$group), "WT");

d.label2$cell_type <- gsub(
	"(.*)((NK cell)|(T cell)|(B cell)|(neutrophil)|(macrophage))", "\\2",
	d.label2$value
);

d.label2$cell_subtype <- trimws(gsub(
	"(.*)((NK cell)|(T cell)|(B cell)|(neutrophil)|(macrophage))", "\\1",
	d.label2$value
));

d.sub <- dplyr::filter(d.label2,
	! value %in% c("B1 cell", "hepatocyte")
);

d.subs <- split(d.sub, d.sub$cell_type);

cols.cell.type <- pal_nejm()(length(d.subs));

for (ctype in names(d.subs)) {
	k <- length(unique(d.subs[[ctype]]$cell_subtype));
	qdraw(
		ggplot(d.subs[[ctype]], aes(x=group, y=mean, ymin=lower, ymax=upper)) +
			labs(title = ctype) +
			theme_classic() + 
			geom_col(fill = cols.cell.type[k]) +
			geom_errorbar(width=0.3, show.legend=FALSE) +
			facet_grid(. ~ cell_subtype, scales="free") +
			scale_y_continuous(n.breaks=3) +
			scale_fill_nejm() +
			theme(
				legend.position = "none", 
				strip.background = element_blank(),
				strip.clip = "off"
			) +
			xlab("") + ylab("proportion")
		,
		width = 0.5 + k * 1, height = 3,
		file = insert(pdf.fn, c("cell-subtype-prop", sanitize_fn(ctype)))
	);
}

d.sub.sig <- left_join(d.sub, enrich.all, by=c("value"="group"));

qdraw(
	ggplot(d.sub.sig,
			aes(x = group, y = cell_subtype, size = mean, colour = cell_type, alpha = q)
		) +
		theme_classic() +
		theme(
			strip.background = element_blank(),
			strip.clip = "off",
		) +
		scale_size_continuous(name = "proportion", trans = "log10") +
		scale_colour_nejm(guide = FALSE) +
		scale_alpha_continuous(name = "FDR", trans = revlog_trans()) +
		facet_grid(cell_type ~ ., scales="free_y", space = "free", switch = "y") +
		geom_point() +
		xlab("treatment") + ylab("cell type")
	,
	width = 3, height = 5,
	file = insert(pdf.fn, c("cell-subtype-prop"))
);

# ---

# Generate overall TSNE and UMAP plots

genes.immune <- .gene_set("IMMUN");

genes <- unique(c(genes.immune, genes.b, genes.t, genes.mp, genes.neut));

dec <- modelGeneVar(sce);
with(dec, plot(mean, total, xlab="mean log-expr", ylab="variance"));
curve(metadata(dec)$trend(x), col="blue", add=TRUE);

# reduce data dimension using PCA
sce <- denoisePCA(sce, dec, subset.row=genes, name="PCAsub");
dim(reducedDim(sce, "PCAsub"))

sce <- runTSNE(sce, dimred="PCA", name="TSNEsub");
sce <- runUMAP(sce, dimred="PCA", name="UMAPsub");


qdraw(
	ggrastr::rasterize(
		plotReducedDim(sce, "TSNEsub", colour_by="label2", point_alpha = 0.3) +
			coord_fixed()
	),
	width = 8, height = 8,
	file = insert(pdf.fn, tag=c("tsne", "clusters2"))
)

qdraw(
	ggrastr::rasterize(
		plotReducedDim(sce, "UMAPsub", colour_by="label2") +
			coord_fixed()
	),
	width = 8, height = 8,
	file = insert(pdf.fn, tag=c("umap", "clusters2"))
)

sce.f <- sce[, ! sce$label2 %in% c("hepatocyte", "B1 cell")];

qdraw(
	ggrastr::rasterize(
		plotReducedDim(sce.f, "TSNEsub", colour_by="label2") +
			coord_fixed()
	),
	width = 8, height = 8,
	file = insert(pdf.fn, tag=c("tsne", "clusters2", "clean"))
)

qdraw(
	ggrastr::rasterize(
		plotReducedDim(sce.f, "UMAPsub", colour_by="label2") +
			coord_fixed()
	),
	width = 8, height = 8,
	file = insert(pdf.fn, tag=c("umap", "clusters2", "clean"))
)

qwrite(sce, in.fn);

