#!/usr/bin/env Rscript

library(io)
library(scran)
library(argparser)
library(DelayedArray)

pr <- arg_parser("Fix path in DelayedMatrix");
pr <- add_argument(pr, "input", help="input RDS file");

argv <- parse_args(pr);

in.fn <- argv$input;

x <- qread(in.fn);

fix_path <- function(x) {

	fp <- seed(counts(x))@filepath;

	# assume that we are in a subdirectory of the project;
	# therefore, we need to ascent to the parent directory
	parent.new <- normalizePath(file.path(getwd(), ".."));

	# match the deepest directory with the name of the project directory
	# and replace the path
	fp.new <- sub(sprintf(".*/%s", basename(parent.new)), parent.new, fp);

	# fix the filepath 
	seed(counts(x))@filepath <- fp.new;

	x
}

x <- fix_path(x);

counts(x)

qwrite(x, in.fn, serialize=FALSE);

