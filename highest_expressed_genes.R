### script to identify high expressed genes (features)

# load libs
source('functions.R')

# load tool locations
tool_path <- load_tool_location()

# get args
get_dir <- commandArgs(trailingOnly = TRUE)

# check input files
in_bam <- get_dir[1]
if(grep('.sort.bam$',in_bam) != 1){
  stop(
    paste0(in_bam, ' is not a BAM file!') 
  )
}
in_bed <- get_dir[2]
if(grep('.bed$',in_bed) != 1){
  stop(
    paste0(in_bed, ' is not a BED file!') 
  )
}

# check output dir
out_path <- parse_dir(get_dir[3])
if(!dir.exists(out_path)){
  dir.create(out_path)
}
if(!dir.exists(paste0(out_path,'highest_expressed_genes/'))){
  dir.create(paste0(out_path,'highest_expressed_genes/'))
}

# check optional threshold for "high" expression
t_high <- get_dir[4]
if(is.na(t_high) || is.null(t_high) || t_high == '' ){
  t_high <- 0.1
} 

### get highest expressed genes
# count
cat('Get highest expressed genes ...')
system(
  paste0(
    tool_path['bedtools'],' coverage',
    ' -s ', # strand specific
    ' -a ',in_bed,
    ' -b ',in_bam,
    ' > ',out_path,'highest_expressed_genes/counts_',gsub('.sort.bam','',basename(in_bam)),'.txt'
  )
)
cat('DONE \n')
# load
c_all <- read.table(
  paste0(out_path,'highest_expressed_genes/counts_',gsub('.sort.bam','',basename(in_bam)),'.txt'),
  header = F,
  stringsAsFactors = F,
  sep = "\t" 
)

### filter highest expressed 
# get unique genes
u_genes <- unique(c_all[,4])
# if only one exon/CDS per gene
if(length(u_genes)==nrow(c_all)){
  c_high <- c_all[order(c_all[,7]/c_all[,9],decreasing = T),]
  c_high <- c_high[c_high[,7]/c_high[,9] > t_high, 1:6] 
# else (eukaryotic workflow)
} else {
  c_sum <- aggregate(c_all[,c(7,9)], by = list(name=c_all[,4]), sum)
  c_sum <- c_sum[order(c_all[,2]/c_all[,3],decreasing = T),]
  c_sum <- c_sum[c_sum[,2]/c_sum[,3] > t_high, ]
  c_high <- c_all[is.element(c_all[,4],c_sum[,1]),1:6]
}
write.table(
  c_high, 
  paste0(out_path,'highest_expressed_genes/highest_expressed_genes.bed'),
  col.names = F,
  row.names = F,
  sep = "\t",
  quote = F
)


