library(optparse)
source('utils.R')


option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "table file with all cleaned sequences", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "table file with output palindromes", metavar = "character"),
  make_option(c("-r", "--residual"), type = "character", default = NULL,
              help = "table file with residual blast hits", metavar = "character"),
  make_option(c("-l", "--min_len"), type = "character", default = NULL,
              help = "table file with all cleaned sequences", metavar = "character")
)

# Parsing command line options
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Assigning values to variables
file.in <- args$input
file.out <- args$output
file.resid <- args$residual
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
x = x[,-9]
x = x[order(x$V4),]
x = x[order(x$V8),]

idx = which((diff(as.numeric(as.factor(x$dir))) != 0) & (abs(diff(x$V4) < min.len)) & (diff(as.numeric(as.factor(x$V8))) == 0))
idx = idx[x$dir[idx] == '+']

# Disentangle tabdem palindromes
while (T) {
  idx.problem = idx[which(diff(idx) == 1)]
  if(length(idx.problem) == 0) break
  print(idx.problem[1] + 1)
  idx = setdiff(idx, idx.problem[1] + 1)   # +1 is important!!!
}

x$pair = 0
x$pair[idx] = idx
x$pair[idx+1] = idx

# ---- Find palindrome pairs ----

y = x[x$pair != 0,]
ny = nrow(y)
idx.pal = seq(1, ny, 2)

z = data.frame(V1 = y$V8[idx.pal], beg =  y$V4[idx.pal], end = y$V5[idx.pal + 1], dir = y$dir[idx.pal], id = y$pair[idx.pal],
               stringsAsFactors = F)
z$len = z$end - z$beg + 1

if(sum(z$dir == '-') != 0) stop('Something is wrong with palindrome directions')


# ---- Save ----

saveRDS(z, file.out)
saveRDS(x[x$pair == 0,], file.resid)

# z[z$len > 1000,]$V1

# ---- Get sequences ----

