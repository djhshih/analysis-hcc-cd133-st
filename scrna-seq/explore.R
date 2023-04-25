library(scuttle)
library(scran)
library(scater)
library(io)
library(ggplot2)

#sce <- qread("sce/Prom1-WT.rds");
sce <- qread("sce/Prom1-DTA.rds");
print(dim(sce))

rownames(sce) <- rowData(sce)$Symbol;

sce <- logNormCounts(sce);

cl <- quickCluster(sce);
cl <- factor(cl);
table(cl)

ggplot(data.frame(colData(sce), cl), aes(x=cl, y=prob_compromised)) +
	theme_classic() +
	geom_violin() +
	scale_y_log10()

colLabels(sce) <- cl;

sce <- runTSNE(sce);
plotTSNE(sce, colour_by="label", text_by="label");
plotTSNE(sce, colour_by="dropletutils_fdr", text_by="label");
plotTSNE(sce, colour_by="sizeFactor", text_by="label");
plotTSNE(sce, colour_by="subsets_mito_percent", text_by="label");

dec <- modelGeneVar(sce);
with(dec, plot(mean, total, xlab="mean log-expr", ylab="variance"));
curve(metadata(dec)$trend(x), col="blue", add=TRUE);

markers <- scoreMarkers(sce);

lapply(markers,
	function(m) {
		d <- m[order(m$mean.AUC, decreasing=TRUE)[1:50], c("self.average", "other.average", "mean.AUC")];
		#rownames(d) <- rowData(sce)$Symbol[match(rownames(d), rowData(sce)$ID)]
		d
	}
)

# Prom1-DTA

plotTSNE(sce, colour_by = "Cd55")


# Prom1-WT

plotTSNE(sce, colour_by = "Gm42418")
# Gm42418 gene overlaps rRNA element Rn45s and represent rRNA contamination
plotTSNE(sce, colour_by = "Alb")
plotTSNE(sce, colour_by = "Hp")
plotTSNE(sce, colour_by = "Apoa1")
plotTSNE(sce, colour_by = "Spp1")

plotTSNE(sce, colour_by = "Tm4sf4")
plotTSNE(sce, colour_by = "Vtn")

plotTSNE(sce, colour_by = "Clu")
plotTSNE(sce, colour_by = "Tesc")

plotTSNE(sce, colour_by = "Serpina1a")

plotTSNE(sce, colour_by = "Col1a2")
plotTSNE(sce, colour_by = "Col1a1")

