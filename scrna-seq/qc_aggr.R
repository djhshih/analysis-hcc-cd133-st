library(scuttle)
library(scran)
library(io)
library(ggplot2)

# samples <- c("Prom1-DTA", "Prom1-WT");
# group <- "non-immune"

samples <- c("Prom1-DTA-IMM", "Prom1-WT-IMM");
group <- "immune"
tag <- NULL;

minfo <- read10xMolInfo(sprintf("./cellranger_all/%s/count_%s/outs/molecule_info.h5", samples[1], samples[1]));

# aggregated sce
sce <- read10xCounts(file.path("aggr", group, "outs", "count", "filtered_feature_bc_matrix.h5"));

# individually filtered sce
sce.fs <- lapply(samples,
	function(samp) {
		qread(file.path("./sce/", paste0(samp, ".rds")));
	}
);

out.fn <- filename(group, path="sce", date=NA);

print(dim(sce))

# cross-reference against filtered sce
# TODO
barcodes <- lapply(sce.fs, function(s) colData(s)$Barcode);
# barcodes in aggregated data have been modified
barcodes[[2]] <- sub("-1$", "-2", barcodes[[2]]);

# ensure that sample order is correct
idxs <- lapply(barcodes, function(b) match(b, colData(sce)$Barcode));
stopifnot(prop.table(table(!is.na(idxs[[1]])))[2] > 0.9)
stopifnot(prop.table(table(!is.na(idxs[[2]])))[2] > 0.9)

idxs <- lapply(barcodes, function(b) match(colData(sce)$Barcode, b));
table(!is.na(idxs[[1]]))
table(!is.na(idxs[[2]]))
valid <- !is.na(idxs[[1]]) | !is.na(idxs[[2]]);
table(valid)

sce.f <- sce[, valid];

dim(sce)
dim(sce.f)

# annotate sample information
sample.idx <- unlist(lapply(strsplit(sce.f$Barcode, "-", fixed=TRUE), function(ss) as.integer(ss[2])));
table(sample.idx)
colData(sce.f)$Sample <- samples[sample.idx];

is.mito <- grep("mt-", rowData(sce)$Symbol);

qc.cells <- perCellQCMetrics(sce.f, subsets=list(mito=is.mito));
colData(sce.f) <- cbind(colData(sce.f), qc.cells);

with(colData(sce.f), table(detected < 500))

summary(qc.cells$detected)
summary(qc.cells$subsets_mito_percent)

with(qc.cells, plot(detected, subsets_mito_percent))

sum.is.low <- isOutlier(qc.cells$sum, type="lower", log=TRUE);
table(sum.is.low)

summary(sce.f$subsets_mito_percent)
hist(sce.f$subsets_mito_percent, breaks=1000);

qwrite(sce.f, insert(out.fn, tag=tag, ext="rds"))

