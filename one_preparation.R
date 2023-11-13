



library(optparse)
source('utils.R')


option_list <- list(
  make_option(c("-q", "--file_query"), type = "character", default = NULL,
              help = "path to target file for the BLAST", metavar = "character"),
  make_option(c("-m", "--file_merged"), type = "character", default = NULL,
              help = "path to merged file after the BLAST", metavar = "character"),
  make_option(c("-o", "--file_out"), type = "character", default = NULL,
              help = "table file with all sequences", metavar = "character"),
  make_option(c("-s", "--seq_cover"), type = "numeric", default = 0.9,
              help = "sequence coverage", metavar = "numeric")
)

# Parsing command line options
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Assigning values to variables
file.merged <- args$file_merged
file.query <- args$file_query
file.out <- args$file_out
seq.cover <- as.numeric(args$seq_cover)

# Check if the necessary files are specified
if(is.null(file.merged) || is.null(file.query) || is.null(file.out) ) {
  stop("Three files - merged, target and out - must be specified", call. = FALSE)
}

# Check if coverage are specified and valid
if(is.null(seq.cover) || seq.cover <= 0 || seq.cover > 1) {
  stop("Sequence coverage (seq.cover) must be specified and between 0 and 1", call. = FALSE)
}


# ---- Main Analysis ----

# file.merged = '../tir/m.txt'
# file.query = '../candidates/tir.fasta'
# file.out = '../tir/out.rds'
# seq.cover = 0.9

# seqs.target = readFastaMy(file.query)

x = read.table(file.merged, stringsAsFactors=F)
query.len = as.numeric(sapply(x$V1, function(s) strsplit(s, '\\|')[[1]][5]))
x = x[(x$V3 - x$V2 + 1) >= query.len * seq.cover,]
# Just take them all
idx =  x$V4 >  x$V5
if(sum(idx) > 0){
  tmp =  x$V4[idx]
  x$V4[idx] =  x$V5[idx]
  x$V5[idx] = tmp
}

x$dir = c('+', '-')[1 * idx + 1]
x$len = x$V5 - x$V4 + 1
x$V9 <- gsub("-", "", x$V9)
x$name = paste(x$V8, x$V4, x$V5, x$dir, x$len, sep = '|')

x = x[!duplicated(x$name),]

# ---- Merge with previous sequences ----

seq.names.prev = c('')
if(file.exists(file.out)){
  y = readRDS(file.out)
  seq.names.prev = y$names
  x = rbind(x, y[!(y$name %in% x$name),])
  rm(y)
}

# ---- Analysis of seqeunces ----

# Check overlaps

x = x[order(-x$V5),]
x = x[order(x$V4),]
x = x[order(x$V8),]
x = x[order(x$dir),]


while(T){
  nx = nrow(x)
  idx.cover = which((x$V4[-nx] <= x$V4[-1]) & (x$V5[-nx] >= x$V5[-1]) & 
                      (x$dir[-1] == x$dir[-nx]) & (x$V8[-1] == x$V8[-nx]))  
  if(length(idx.cover) != 0){
    x = x[-(idx.cover+1),]
  } else {
    break
  }
}


# partial coverage
idx.partial = which((x$V4[-nx] >= x$V4[-1]) & (x$V4[-nx] <= x$V5[-1]) & 
              (x$dir[-1] == x$dir[-nx]) & (x$V8[-1] == x$V8[-nx]))
if(length(idx.partial) > 0){
  stop('Partial overlap')
  idx.partial = rev(idx.partial)
  for(i in idx.partial){
    print(i) 
  }
}

# If no new sequences were found
if(length(setdiff(seq.names.prev, x$name)) == 0){
  
  
  stop('Final')
  
} else {
  
  # Save the table
  
  x = x[order(-x$V5),]
  x = x[order(x$V4),]
  x = x[order(x$V8),]
  saveRDS(x, file.out, compress = F)
  
  # Save something to search
  idx.unique = !duplicated(x$V9)
  seqs = x$V9[idx.unique]
  names(seqs) = x$name[idx.unique]
  
  writeFastaMy(seqs, file.query)
  
}

