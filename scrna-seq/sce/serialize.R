library(io)
library(scran)

# convert DelayedMatrix to dgCMatrix (compressed parse column matrix)
# so that we do not need to refer to a HD5F file for the data

#in.fn <- "immune.rds";
in.fn <- "non-immune.rds";

x <- qread(in.fn);

m <- counts(x);
counts(x) <- as(m, "dgCMatrix");

qwrite(x, in.fn);

