library(optparse)
source('utils.R')


option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "table file with all cleaned sequences", metavar = "character")
)

# Parsing command line options
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Assigning values to variables
file.in <- args$input

# Check if the necessary files are specified
if(is.null(file.in) ) {
  stop("The input file must be specified", call. = FALSE)
}



# ---- Main Analysis ----
# ---- Read ----
x = readRDS(file.in)

min.len = 10000
idx = which((diff(as.numeric(as.factor(x$dir))) != 0) & (abs(diff(x$V4) < min.len)) & (diff(as.numeric(as.factor(x$V8))) == 0))
x$beg = 0
x$beg[idx] = idx
x$beg[idx+1] = idx

# ---- Find palindrome pairs ----









