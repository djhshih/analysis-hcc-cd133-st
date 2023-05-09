# Plot a dot plot of expression level for genes in gene sets

library(io)
library(scran)
library(SingleCellExperiment)
library(reshape2)
library(dplyr)
library(ggplot2)

source("../R/common.R");
source("./plot_common.R");

sce <- qread("../sce/immune_filtered.rds");
markers <- qread("./mouse-signature_fixed.gmt"); # t-cell-signatures_azizi18.gmt

out.fn <- filename("immune", tag="mp");

# re-calculate the logcounts to overcome delayed array error
sce <- logNormCounts(sce);

sce.mp <- sce[, sce$label.mp.cell != "non-macrophage"];
dim(sce.mp) # [1] 32283  1455

# if include conditions
#sce.mp$label.mp.cell.condition <- paste0(sce.mp$label.mp.cell, sce.mp$Sample);
#sce.mp$label.mp.cell.condition <- sub("Prom1-DTA-IMM", 
#  "-DTA", sce.mp$label.mp.cell.condition);
#sce.mp$label.mp.cell.condition <- sub("Prom1-WT-IMM",
#  "-WT", sce.mp$label.mp.cell.condition);
#table(sce.mp$label.mp.cell.condition);

table(sce.mp$label.mp.cell);

gset.names <- c("M1 Macrophage Polarization", "M2 Macrophage Polarization");
gsets <- markers$data[gset.names];

clusters <- sce.mp$label.mp.cell;
#clusters <- sce.mp$label.mp.cell.condition; # if include condition

ds <- lapply(gsets,
	function(gset) summarize_gset(sce.mp, gset, clusters)
);


# combine the g
d <- Reduce(rbind, 
	mapply(
		function(d, gset) data.frame(d, gset=gset),
		ds, gset.names,
		SIMPLIFY = FALSE
	)
);


d$gset <- sub(" Macrophage Polarization", "", d$gset);
d$group <- sub(" macrophage", "", d$group);

output_plot <- plot_graph(d, title = "macrophage");

qdraw(output_plot, width = 4, height = 11.5,
  file = insert(out.fn, c("gene-expr", "m1-m2-gset-unfiltered"), ext="pdf"))


# filter low expression
exprs_cutoff <- 0.1;
d_filtered <- d[-which(d$p_expressed <= exprs_cutoff), ];
dim(d_filtered)

output_plot <- plot_graph(d_filtered, title = "macrophage-filtered");

qdraw(output_plot, width = 4, height = 8.5,
      file = insert(out.fn, c("gene-expr", "m1-m2-gset-filtered"), ext="pdf"))
