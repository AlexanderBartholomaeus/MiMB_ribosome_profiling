### script count reads per base and generate coverage plot around start and stop

# load libs
library(foreach)
library(doParallel)
source('functions.R')

# load tool locations
tool_path <- load_tool_location()

# set max cores
max_cores <- get_max_cores() # from functions.R
# register max 8 cores because filtering process is I/O intense (mainly read and writing)
if(max_cores > 8) {
  max_cores <- 8
}
registerDoParallel(makeCluster(max_cores))

# get args
get_dir <- commandArgs(trailingOnly = TRUE)

# check input 
in_path <- parse_dir(get_dir[1])
if(!dir.exists(in_path)){
  stop(
    paste0(in_path, ' does not exist!') 
  )
}
# check if bam files
in_files <- list.files(
  path = in_path, 
  pattern = '.bam',
  ignore.case = TRUE
)
if(length(in_files) < 1){
  stop('No BAM files found!')
}

# check bed file
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
if(!dir.exists(paste0(out_path,'coverage_start_stop/'))){
  dir.create(paste0(out_path,'coverage_start_stop/'))
}
if(!dir.exists(paste0(out_path,'coverage_start_stop/counts/'))){
  dir.create(paste0(out_path,'coverage_start_stop/counts/'))
}

### count
cat('Count reads per base ... ')
for_res <- foreach(i = 1:length(in_files)) %dopar% {
  system(
    paste0(
      tool_path['bedtools'],' coverage',
      ' -d ', # count per base
      ' -s ', # strand specific
      ' -a ',in_bed,
      ' -b ',in_path,in_files[i],
      ' > ',out_path,'coverage_start_stop/counts/bedtools_',gsub('.bam','',basename(in_files[i])),'.txt'
    )
  )
}
cat('DONE \n')


### read counts and generate plot and file
cat('Normalize counts and create plot ... ')
for_res <- foreach(i = 1:length(in_files)) %dopar% {
  # read counts
  counts <- read.table(
    paste0(
      out_path,'coverage_start_stop/counts/bedtools_',gsub('.bam','',basename(in_files[i])),'.txt'
    ),
    header = F,
    stringsAsFactors = F,
    sep = "\t"
  )
  # get unique genes
  u_genes <- unique(counts[,4])
  # return normalized counts around start and stop
  c_list <- sapply(1:length(u_genes),function(j){
    sel <- counts[counts[,4]==u_genes[j],]
    sel <- sel[order(sel[,2]+sel[,7]),]
    cou <- sel[,8]
    dir <- sel[1,6]
    norm_cou <- norm_counts_start_stop(
      cou,
      dir,
      cds_length = 150)
  })
  # remove null
  u_genes2 <- u_genes[!unlist(lapply(c_list,function(j){is.null(j)}))]
  c_list <- c_list[!unlist(lapply(c_list,function(j){is.null(j)}))]
  # parse to matrix
  c_mat <- round(matrix(unlist(c_list), nrow = length(c_list), ncol = 400, byrow = T),digits=5)
  c_mat_mean <- round(apply(c_mat, 2, mean, na.rm=T),digits=5)
  colnames(c_mat) <- paste0('pos_',c(-50:-1,1:150,-150:-1,1:50))
  rownames(c_mat) <- u_genes2
  # write tables
  write.table(
    c_mat, 
    paste0(out_path,'coverage_start_stop/counts/all_genes_',gsub('.bam','',basename(in_files[i])),'.csv'),
    col.names = T,
    row.names = T,
    sep = ','
  )
  write.table(
    data.frame(position=c(-50:-1,1:150,-150:-1,1:50),counts=c_mat_mean),
    paste0(out_path,'coverage_start_stop/counts/plot_',gsub('.bam','',basename(in_files[i])),'.csv'),
    col.names = T,
    row.names = F,
    sep = ',',
    quote = T
  )
  # mean
  pdf(
    width = 7,
    height = 5,
    paste0(out_path,'coverage_start_stop/coverage_',gsub('.bam','',basename(in_files[i])),'.pdf')
  )
  plot(
    c(c_mat_mean[1:200],rep(NA,50),c_mat_mean[201:400]),
    type='l',
    xlab = 'Position [nt]',
    ylab = 'Counts [normalized]',
    las=1,
    xaxt = 'n',
    main = gsub('.bam','',basename(in_files[i]))
  )
  axis(
    1, 
    at = c(1,51,100,150,200,251,301,351,400,450), 
    label = c('-50','START','50','100','150','-150','-100','-50','STOP','+50'))
  dev.off()
}
cat('DONE \n')


