# Plot a dot plot of expression level for genes in gene sets

library(io)
library(scran)
library(SingleCellExperiment)
library(reshape2)
library(dplyr)
library(ggplot2)

source("../R/common.R");

sce <- qread("../sce/immune_filtered.rds");
markers <- qread("./t-cell-signatures_azizi18.gmt");

out.fn <- filename("immune", tag="mp");

# re-calculate the logcounts to overcome delayed array error
sce <- logNormCounts(sce);

sce.mp <- sce[, sce$label.mp.cell != "non-macrophage"];
dim(sce.mp)
table(sce.mp$label.mp.cell)

summarize_gset <- function(sce, gset, groups) {

	idx <- match(tolower(gset), tolower(rowData(sce)$Symbol));
	valid <- !is.na(idx);
	message("unmatched fraction: ", mean(valid))
	idx.valid <- idx[valid];

	expr <- logcounts(sce[idx.valid, ]);

	ms <- scoreMarkers(expr, groups);

	d <- do.call(rbind, mapply(
		function(m, k) {
			data.frame(
				group = k,
				gene = rownames(m),
				logfc = m[, "self.average"] - m[, "other.average"],
				p_expressed = m[, "self.detected"]
			)
		},
		ms, names(ms),
		SIMPLIFY = FALSE
	));
	rownames(d) <- NULL;

	d
}

gset.names <- c("M1 Macrophage Polarization", "M2 Macrophage Polarization");
gsets <- markers$data[gset.names];

clusters <- sce.mp$label.mp.cell;

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

qdraw(
	ggplot(d,
		aes(
			x = group, y = gene,
			colour = bound(logfc, c(-1.5, 1.5)),
			size = p_expressed * 100
		)
	) +
		theme_classic() +
		labs(title = "macrophage") +
		geom_point() +
		facet_grid(gset ~ ., scales="free_y", space="free", switch="y") +
		scale_colour_gradient2(breaks = (-1):1,
			low="royalblue4", 
			mid="grey90",
			high="orangered2"
		) +
		guides(
			colour = guide_colourbar("log FC"),
			size = guide_legend("% expressed")
		) +
		theme(
			strip.text.y.left = element_text(angle = 0, hjust=0),
			strip.background = element_blank(),
			axis.ticks.y = element_blank(),
			axis.ticks.x = element_blank(),
			axis.text.x = element_text(angle = 45, hjust=1)
		) +
		scale_size_continuous(range = c(-0.5, 3.5)) +
		scale_y_discrete(limits = rev) +
		xlab("") + ylab("")
	,
	width = 4, height = 10,
	file = insert(out.fn, c("gene-expr", "m1-m2-gset"), ext="pdf")
)

