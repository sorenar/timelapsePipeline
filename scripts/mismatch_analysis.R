#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

result_path= args[1]
exp_name= args[2]
chr_name= args[3]
chr_len= as.numeric(args[4])
#print(exp_name)
#print(paste0("/dfs3/samlab/sorenar/OsO-seq/20190905/filtering/",exp_name,"/sub/"))
setwd(paste0("/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/results/",result_path))

# setwd("/Users/sorenarahmanian/Documents/4SU_OsO/timelapse/HS/filtering/")
# 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsamtools")

library(Rsamtools)
library(gtools)
# library(seqinr)

#which <- IRangesList("chr8"=IRanges(1,145000000))
which <- GRanges(chr_name,IRanges(start=1,end=(chr_len-1)))

what <- c("qname","rname", "strand", "pos", "qwidth", "seq","flag","cigar","mrnm","qual")
param <- ScanBamParam(tag=c("NM","MD"), what=what,which=which)

bamfile <- system.file("extdata", paste0(exp_name,"_sorted.bam"), package="Rsamtools")
bam <- scanBam(paste0(exp_name,"_sorted.bam"), param=param)


# fastafile = system.file("genome.fa",package = "seqinR")
# fasta = read.fasta("genome.fa",as.string = TRUE,strip.desc = TRUE)

snp_data = read.delim("/dfs3/samlab/sorenar/OsO-seq/timelapsePipeline/ref/MEF_TCSNP.txt",header = FALSE,stringsAsFactors = FALSE)
snp_data$pos = paste0(snp_data[,1],"_",snp_data[,2])

.unlist <- function (x)
   {
     ## do.call(c, ...) coerces factor to integer, which is undesired
       x1 <- x[[1L]]
       if (is.factor(x1)) {
         structure(unlist(x), class = "factor", levels = levels(x1))
        } else {
           do.call(c, x)
               }
        }
bam <- unname(bam) # names not useful in unlisted result
elts <- setNames(bamWhat(param), bamWhat(param))
lst <- lapply(elts, function(elt) .unlist(lapply(bam, "[[", elt)))

bam_data = data.frame(do.call("DataFrame", lst))
bam_data$NM = as.character(bam[[1]][["tag"]][[1]])
bam_data$MD = as.character(bam[[1]][["tag"]][[2]])


bam_data$cat = 0
write.table("",paste0(exp_name,"_read_profile_",chr_name,".txt"),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")

for(i in 1:dim(bam_data)[1]){
        
        nm = as.numeric(bam_data$NM[i])
        if(nm > 0){
        md = bam_data$MD[i]
        cigar = bam_data$cigar[i]
        seq = bam_data$seq[i]
        qual = bam_data$qual[i]
        len = bam_data$qwidth[i]
        rname = bam_data$rname[i]
        pos = bam_data$pos[i]
        
        markers = gregexpr("[MSNID]",cigar)[[1]]
        if (length(markers) == 1){
                cigar_struct= data.frame(type = substr(cigar,markers[1],markers[1]),len = as.numeric(substr(cigar,1,markers[1]-1)))
        }else{
                cigar_struct= data.frame(type = apply(data.frame(markers),1,function(x) substr(cigar,x,x)),
                                         len = apply(cbind(c(0,markers[1:(length(markers)-1)])+1,markers-1),1,function(x) as.numeric(substr(cigar,x[1],x[2]))))
        }
        
        sum(cigar_struct$len[cigar_struct$type %in% c("M","I","S")])
        cigar_struct$seq_len=0
        cigar_struct$seq_len[cigar_struct$type %in% c("M","I","S")]= cigar_struct$len[cigar_struct$type %in% c("M","I","S")]
        cigar_struct$seq_pos = cumsum(cigar_struct$seq_len)
        cigar_struct$md_len = 0
        cigar_struct$md_len[cigar_struct$type %in% c("M")]= cigar_struct$len[cigar_struct$type %in% c("M")]
        cigar_struct$md_pos = cumsum(cigar_struct$md_len)
        cigar_struct$diff = cigar_struct$seq_pos - cigar_struct$md_pos
        cigar_struct$query_len=0
        cigar_struct$query_len[cigar_struct$type %in% c("D","N")]= cigar_struct$len[cigar_struct$type %in% c("D","N")]
        cigar_struct$qref_adj = cumsum(cigar_struct$query_len)        
        
        md_markers = gregexpr("[ATCG\\^]",md)[[1]]
        md_struct= data.frame(type = apply(data.frame(md_markers),1,function(x) substr(md,x,x)),
                              len = apply(cbind(c(0,md_markers[1:(length(md_markers)-1)])+1,md_markers-1),1,function(x) as.numeric(substr(md,x[1],x[2]))))
        
        md_struct$md_len = 0
        md_struct$md_len = md_struct$len + 1
        md_struct$md_len[md_struct$type == "^"] = md_struct$len[md_struct$type == "^"]
        md_struct = md_struct[!is.na(md_struct$len),]
        md_struct$md_pos = cumsum(md_struct$md_len)
        md_struct = md_struct[!md_struct$type == "^",]
        
        md_struct$miamatch = apply(md_struct,1,function(x) {
                ind = as.numeric(x[4])+cigar_struct$diff[cigar_struct$md_pos >= as.numeric(x[4])][1]
                substr(seq,ind,ind)
        })
        
        md_struct$qual = apply(md_struct,1,function(x) {
                ind = as.numeric(x[4])+cigar_struct$diff[cigar_struct$md_pos >= as.numeric(x[4])][1]
                if(ind > 2 && ind < (len-2)) {
                        asc(substr(qual,ind,ind))-33
                } else {
                        -1
                }
                
        })
        
        md_struct$pos = apply(md_struct,1,function(x) {
                ind = as.numeric(x[4])+cigar_struct$qref_adj[cigar_struct$md_pos >= as.numeric(x[4])][1]+pos
                paste0(rname,"_",ind)
        })
        
        md_struct = md_struct[md_struct$qual >= 30,]
        
        
        mm_TC = sum(md_struct$type == "T" & md_struct$miamatch == "C" & !(md_struct$pos %in% snp_data$pos))
        mm_tot = dim(md_struct)[1]
        # print(mm_TC)
        if(mm_TC > 2*mm_tot/3){
                bam_data$cat[i] = -1
        }else if(mm_TC == 0){
                bam_data$cat[i] = 0
        }else if(mm_TC >= 5){
                bam_data$cat[i] = -2
        }else{
                bam_data$cat[i] = mm_TC
        }
}
        write.table(bam_data[i,c("qname","cat")],paste0(exp_name,"_read_profile_",chr_name,".txt"),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t",append = TRUE)

}




