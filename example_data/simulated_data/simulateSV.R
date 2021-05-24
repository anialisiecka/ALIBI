library(RSVSim)
library(BSgenome.Hsapiens.UCSC.hg38)

args <- commandArgs(trailingOnly = TRUE)

# number of SV of each type
n <- as.numeric(args[1])

# number of repeats
k <- as.numeric(args[2])

# genomic sequence name
gname="chr13_KI270842v1_alt"

# number of simulated sequences
m = 10


simSVn <- function(gname, m, n, k){
  genome = DNAStringSet(Hsapiens[[gname]])
  names(genome) = c("alt")
  outfile = paste(gname,"_",m,"_",n,"_",k,".fasta", sep="")
  for(i in 1:m){
    sim = simulateSV(output = NA, genome = genome, verbose = F,
                     bpSeqSize=6, percCopiedIns = 0.2, maxDups = 2,
                     sizeDels = 20, sizeIns = 20, sizeInvs = 200, sizeDups = 500,
                     dels = n,  ins = n, invs = n, dups = n
    )
    names(sim) = c(paste("alt",i, sep=""))
    writeXStringSet(sim, outfile, append = (i>1))
  }
}

simSVn(gname, m, n, k)
