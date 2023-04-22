library(DropletUtils)
library(scuttle)
library(scran)
library(filenamer)
library(io)
library(ggplot2)

min.umi.count <- 500;

minfo <- read10xMolInfo("./cellranger_all/Prom1-WT/count_Prom1-WT/outs/molecule_info.h5");

samples <- list.files("cellranger", full.names=FALSE);
names(samples) <- samples;

e.outs <- lapply(samples,
	function(s) {
		out.fn <- filename("empty-drops", path=c("dropletutils", s), ext="rds", date=NA);
		sce <- read10xCounts(file.path("cellranger", s, "raw_feature_bc_matrix.h5"));
		e.out <- emptyDrops(sce, lower=min.umi.count);
		qwrite(e.out, out.fn);
		e.out
	}
)

# ---

fdr.cut <- 0.01;

for (s in samples) {

	out.fn <- filename("empty-drops", path=c("dropletutils", s), date=NA);
	pdf.fn <- insert(out.fn, ext="pdf");

	e.out <- e.outs[[s]];
	is.cell <- e.out$FDR < fdr.cut & !is.na(e.out$FDR);
	n.cells <- sum(is.cell);
	print(n.cells)
	# parameter npts should be increased 
	# if there are any non-signicant & limited barcodes
	table(limited = e.out$Limited, significant = is.cell)

	# diagnosis plots

	qdraw(
		with(e.out,
			plot(Total, -LogProb, col=ifelse(is.cell, "red", "black"),
				xlab="total UMI count", ylab="- log prob",
				xlim = c(0, 20000), ylim=c(0, 20000)
			)
		),
		insert(pdf.fn, tag="logp-vs-umi")
	)

	qdraw(
		qplot(e.out$Total, geom="histogram", bins=80) +
			theme_classic() +
			scale_x_log10(n.breaks=10) +
			annotation_logticks() +
			xlab("total UMI count") + ylab("count")
		,
		width = 8,
		insert(pdf.fn, tag=c("hist", "umi-count"))
	)

	qdraw(
		qplot(e.out$Total, geom="histogram", bins=80) +
			theme_classic() +
			scale_x_log10(n.breaks=5, limits=c(100, 1e5)) +
			annotation_logticks() +
			xlab("total UMI count") + ylab("count")
		,
		width = 8,
		insert(pdf.fn, tag=c("hist", "umi-count", "xlim"))
	)

	graphics.off()

}

