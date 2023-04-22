library(DropletUtils)
library(scuttle)
library(scran)
library(miQC)
library(io)
library(ggplot2)

samp <- "Prom1-WT"
#samp <- "Prom1-DTA"

#mito.posterior.cut <- 0.90;
#fdr.cut <- 0.01;
#tag <- NULL;

# mito.posterior.cut <- 1;
# mito.cut <- 5;
# fdr.cut <- 0.05;
# tag <- "relaxed";

# samp <- "Prom1-WT-IMM";
# samp <- "Prom1-DTA-IMM";

mito.posterior.cut <- 1;
mito.cut <- 5;
fdr.cut <- 0.01;
tag <- NULL;


minfo <- read10xMolInfo(sprintf("./cellranger_all/%s/count_%s/outs/molecule_info.h5", samp, samp));
sce <- read10xCounts(file.path("cellranger", samp, "raw_feature_bc_matrix.h5"));
cell.calls <- qread(file.path("./dropletutils/", samp, "empty-drops.rds"));

out.fn <- filename(samp, path="sce", date=NA);

print(dim(sce))

idx <- which(cell.calls$FDR < fdr.cut);
sce.f <- sce[, idx];
cell.calls.f <- cell.calls[idx, c("PValue", "Limited", "FDR")];
colnames(cell.calls.f) <- paste0("dropletutils_", tolower(colnames(cell.calls.f)));

print(dim(sce.f))

is.mito <- grep("mt-", rowData(sce)$Symbol);

qc.cells <- perCellQCMetrics(sce.f, subsets=list(mito=is.mito));
colData(sce.f) <- cbind(colData(sce.f), qc.cells, cell.calls.f);

with(colData(sce.f), plot(dropletutils_pvalue, detected, pch="."))
with(colData(sce.f), plot(dropletutils_pvalue, detected, xlim=c(1e-5, 0.05), log="x"))
#with(colData(sce.f), smoothScatter(dropletutils_pvalue, detected))
with(colData(sce.f), plot(dropletutils_fdr, detected, pch="."))
#with(colData(sce.f), smoothScatter(dropletutils_fdr, detected))

with(colData(sce.f), table(detected < 500))


summary(qc.cells$detected)
summary(qc.cells$subsets_mito_percent)

with(qc.cells, plot(detected, subsets_mito_percent))

sum.is.low <- isOutlier(qc.cells$sum, type="lower", log=TRUE);
table(sum.is.low)

model.miqc <- miQC::mixtureModel(sce.f, model_type = "linear");
miQC::plotModel(sce.f, model.miqc);
miQC::plotModel(sce.f, model.miqc) + ylim(0, 25);
#miQC::plotFiltering(sce.f, model.miqc);

hist(sce.f$subsets_mito_percent, breaks=1000);

with(colData(sce.f), table(subsets_mito_percent >= 10))
with(colData(sce.f), table(subsets_mito_percent >= 5))
with(colData(sce.f), prop.table(table(subsets_mito_percent >= 10)))
with(colData(sce.f), prop.table(table(subsets_mito_percent >= 5)))


if (mito.posterior.cut < 1) {
	sce.f2 <- miQC::filterCells(sce.f, model.miqc, posterior_cutoff = mito.posterior.cut,
		keep_all_below_boundary = TRUE,
		enforce_left_cutoff = TRUE,
	);
} else {
	sce.f2 <- sce.f[, colData(sce.f)$subsets_mito_percent < mito.cut]
}
print(1 - ncol(sce.f2) / ncol(sce.f))

with(colData(sce.f), smoothScatter(detected, subsets_mito_percent))

with(colData(sce.f2), smoothScatter(detected, subsets_mito_percent))

with(colData(sce.f2), hist(detected, breaks=1000))
abline(v=500, col="red")

with(colData(sce.f2), table(detected < 500))
with(colData(sce.f2), table(detected >= 500 & detected <= 3000))

summary(colData(sce.f)$subsets_mito_percent)
summary(colData(sce.f2)$subsets_mito_percent)

dim(sce.f2)

#sce.f3 <- quickPerCellQC(sce.f2, sub.fields="subsets_mito_percent");
#dim(sce.f3)

#qc.features <- perFeatureQCMetrics(sce.f2);
#summary(qc.features$mean)
#summary(qc.features$detected)

#counts.avg <- calculateAverage(sce.f2);
#summary(counts.avg)

qwrite(sce.f2, insert(out.fn, tag=tag, ext="rds"))

