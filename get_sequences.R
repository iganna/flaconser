library(optparse)
source('utils.R')


option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "table file with all cleaned sequences", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "table file with output palindromes", metavar = "character")
)

# Parsing command line options
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Assigning values to variables
file.in <- args$input
file.out <- args$output
min.len <- args$min_len

# Check if the necessary files are specified
if(is.null(file.in) ) {
  stop("The input file must be specified", call. = FALSE)
}


# ---- Main Analysis ----

# 
# file.in = '../tir/out.rds'
# file.out = '../tir/palindromes.rds'
# min.len = 10000

# ---- Read ----
x = readRDS(file.in)
# 
# for(irow in which(x$dir == '-')){
#   x$V9[irow] = nt2seq(revCompl(seq2nt(x$V9[irow])))
# }

seqs = x$V9
names(seqs) = x$name

writeFastaMy(seqs, file.out)




