library(io)
library(scran)

#in.fn <- "immune.rds";
in.fn <- "non-immune.rds";

x <- qread(in.fn);

fp <- seed(counts(x))@filepath;

# assume that we are in a subdirectory of the project;
# therefore, we need to ascent to the parent directory
parent.new <- normalizePath(file.path(getwd(), ".."));

# match the deepest directory with the name of the project directory
# and replace the path
fp.new <- sub(sprintf(".*/%s", basename(parent.new)), parent.new, fp);

# fix the filepath 
seed(counts(x))@filepath <- fp.new;

counts(x)

qwrite(x, in.fn);

