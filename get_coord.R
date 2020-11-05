args <- commandArgs(trailingOnly = TRUE)
snp.info.file <- args[1]
out.file <- args[2]
begin.hg19 <- as.numeric(args[3])
end.hg19 <- as.numeric(args[4])

hg19pos <- read.delim(snp.info.file, stringsAsFactors = F, sep=" ", head=F)[,1]

hg19.begin.index <- which(hg19pos == begin.hg19)
if (length(hg19.begin.index) == 0) { #Find the position before the specified beginning that is closest 
  hg19pos.delta <- begin.hg19 - hg19pos
  hg19pos.delta[hg19pos.delta < 0 ] <- NA
  hg19.begin.index <- which.min(hg19pos.delta)
}
hg19.end.index <- which(hg19pos == end.hg19)
if (length(hg19.end.index) == 0) { #Find the position after the specified beginning that is closest 
  hg19pos.delta <- hg19pos - end.hg19
  hg19pos.delta[hg19pos.delta < 0 ] <- NA
  hg19.end.index <- which.min(hg19pos.delta)
}

cat(hg19pos[hg19.begin.index], hg19pos[hg19.end.index], sep="\t", file=out.file, append=F)
