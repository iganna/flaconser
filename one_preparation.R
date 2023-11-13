library(optparse)
source('utils.R')


option_list <- list(
  make_option(c("-m", "--merged_file"), type = "character", default = NULL,
              help = "path to merged file", metavar = "character"),
  make_option(c("-t", "--target_file"), type = "character", default = NULL,
              help = "path to target file", metavar = "character"),
  make_option(c("-o", "--out_file"), type = "character", default = NULL,
              help = "output file with all long sequences", metavar = "character"),
  make_option(c("-c", "--seq_cover"), type = "numeric", default = 0.9,
              help = "sequence coverage", metavar = "numeric")
)

# Parsing command line options
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Assigning values to variables
file.merged <- args$merged_file
file.target <- args$target_file
file.target <- args$target_file
seq.cover <- args$seq_cover

# Check if the necessary files are specified
if(is.null(file.merged) || is.null(file.target)) {
  stop("Both merged and target files must be specified", call. = FALSE)
}

# Check ifcoverage are specified and valid
if(is.null(seq.cover) || seq.cover <= 0 || seq.cover > 1) {
  stop("Sequence coverage (seq.cover) must be specified and between 0 and 1", call. = FALSE)
}


# ---- Main Analysis ----

# file.merged = '../tir/merged.txt'
# file.target = '../candidates/tir.fasta'
# 
# seq.len = 140
# seq.cover = 0.9

seqs.target = readFastaMy(file.target)

x = read.table(file.merged, stringsAsFactors=F)
query.len = as.numeric(sapply(x$V1, function(s) strsplit(s, '\\|')[[1]][5]))
x = x[(x$V3 - x$V2 + 1) >= query.len * seq.cover,]

# Just take them all
idx =  x$V4 >  x$V5
tmp =  x$V4[idx]
x$V4[idx] =  x$V5[idx]
x$V5[idx] = tmp
x$dir = c('+', '-')[1 * idx + 1]
x$len = x$V5 - x$V4 + 1
x$V9 <- gsub("-", "", x$V9)
x$name = paste(x$V8, x$V4, x$V5, x$dir, x$len, sep = '|')

x = x[!duplicated(x$name),]

# ---- Merge with previous sequences ----



# ---- Analysis of seqeunces ----


# Check overlap


idx.unique = !duplicated(x$V9)
seqs = x$V9[idx.unique]
names(seqs) = x$name[idx.unique]

writeFastaMy(seqs, file.out)


# # For Fun
# 
# x$dir = (x$V4 > x$V5) * 1
# x = x[order(x$V4),]
# x = x[order(x$V8),]
# 
# min.len = 1000
# idx = which((diff(x$dir) != 0) & (abs(diff(x$V4) < min.len)) & (diff(as.numeric(as.factor(x$V8))) == 0))
# x$beg = 0
# x$beg[idx] = idx
# x$beg[idx+1] = idx

