library(optparse)
source('utils.R')

# Parsing command line options
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Assigning values to variables
file.merged <- args$merged_file
file.target <- args$target_file
seq.len <- args$seq_len
seq.cover <- args$seq_cover

# Check if the necessary files are specified
if(is.null(file.merged) || is.null(file.target)) {
  stop("Both merged and target files must be specified", call. = FALSE)
}

# Check if sequence length and coverage are specified and valid
if(is.null(seq.len) || seq.len <= 0) {
  stop("Sequence length (seq.len) must be specified and greater than 0", call. = FALSE)
}
if(is.null(seq.cover) || seq.cover <= 0 || seq.cover > 1) {
  stop("Sequence coverage (seq.cover) must be specified and between 0 and 1", call. = FALSE)
}


# ---- Main Analysis ----

file.merged = 'bl_tir_merged.txt'
file.target = ''

seq.len = 140
seq.cover = 0.9


x = read.table(file.merged, stringsAsFactors=F)
x = x[x$V7 >= seq.len * seq.cover,]

# Just take them all
pos.beg = x$V4
pos.end = x$V5
idx = pos.beg > pos.end
tmp = pos.beg[idx]
pos.beg[idx] = pos.end[idx]
pos.end[idx] = tmp
pos.dir = c('+', '-')[1 * idx + 1]

seq.names = paste(x$V8, pos.beg, pos.end, pos.dir, sep = '|')
idx.unique = !duplicated(seq.names)

seqs = x$V9[idx.unique]
names(seqs) = seq.names[idx.unique]

writeFastaMy(seqs, file.target)


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

