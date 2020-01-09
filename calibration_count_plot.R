### script to count reads and generate plot for calibration of the reads

# load libs
library(foreach)
library(doParallel)
source('functions.R')

# load parameters
conf <- load_mapping_config()

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

# check input dir
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
if(!dir.exists(paste0(out_path,'calibration/'))){
  dir.create(paste0(out_path,'calibration/'))
}
if(!dir.exists(paste0(out_path,'calibration/counts'))){
  dir.create(paste0(out_path,'calibration/counts'))
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
      ' > ',out_path,'calibration/counts/bedtools_',gsub('.bam','',basename(in_files[i])),'.txt'
    )
  )
}
cat('DONE \n')

### plot counts in frame
cat('Normalizing counts and generating plots ... ')
for_res <- foreach(i = 1:length(in_files)) %dopar% {
  # read counts
  counts <- read.table(
    paste0(
      out_path,'calibration/counts/bedtools_',gsub('.bam','',basename(in_files[i])),'.txt'
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
    # normalize counts
    norm_cou <- norm_counts_start_stop(
      cou,
      dir,
      cds_length = 150,
      normalize = F)
    # check which are not dividable by 3
    # if((length(cou)-100)%%3 != 0){
    #   #system(paste0('cat ',u_genes[j],'\n not dividable by 3.'))
    #   print(paste0(u_genes[j],' not dividable by 3.'))
    # }
    # frame counts
    frame_cou <- frame_counts(
      cou,
      dir
    )
    # return
    list(
      norm_cou,
      frame_cou
    )
  })
  # remove null
  c_list_norm <- c_list[seq(1,length(c_list),by=2)]
  u_genes2 <- u_genes[!unlist(lapply(c_list_norm,function(j){is.null(j)}))]
  c_list_norm <- c_list_norm[!unlist(lapply(c_list_norm,function(j){is.null(j)}))]
  c_list_frame <- c_list[seq(2,length(c_list),by=2)]
  
  # parse to matrix
  c_mat <- round(matrix(unlist(c_list_norm), nrow = length(c_list_norm), ncol = 400, byrow = T),digits=5)
  c_mat_mean <- round(apply(c_mat, 2, sum, na.rm=T),digits=5) # sum
  #c_mat_mean <- round(apply(c_mat, 2, mean, na.rm=T),digits=5) # mean (do normalization above)
  colnames(c_mat) <- paste0('pos_',c(-50:-1,1:150,-150:-1,1:50))
  rownames(c_mat) <- u_genes2
  c_frame_mat <- round(matrix(unlist(c_list_frame), nrow = length(c_list_frame), ncol = 3, byrow = T),digits=5)
  c_frame_sum <- round(apply(c_frame_mat, 2, sum, na.rm=T),digits=5)
  
  # write count tables
  write.table(
    c_mat, 
    paste0(out_path,'calibration/counts/all_genes_',gsub('.bam','',basename(in_files[i])),'.csv'),
    col.names = T,
    row.names = T,
    sep = ','
  )
  write.table(
    data.frame(position=c(-50:-1,1:150,-150:-1,1:50),counts=c_mat_mean),
    paste0(out_path,'calibration/counts/plot_',gsub('.bam','',basename(in_files[i])),'.csv'),
    col.names = T,
    row.names = F,
    sep = ',',
    quote = T
  )
  # plot for calibration 
  pdf(
    width = 7,
    height = 5,
    paste0(out_path,'calibration/calibration_',gsub('.bam','',basename(in_files[i])),'.pdf')
  )
  # plot count for frames
  par(mar=c(5.1,5.1,4.1,2.1))
  barplot(
    c_frame_sum,
    names.arg = c(0,1,2),
    col = c('darkgreen','orange','navyblue'),
    las = 1,
    xlab = 'Frame',
    cex.axis = 1.2,
    cex.names = 1.2,
    cex.lab = 1.2,
    main = gsub('.sort.bam','',basename(in_files[i])),
    ylab = '')
  title(ylab = 'Counts', line = 4, cex.lab = 1.2)
  
  # plot around the start and stop codon region
  if(length(grep('firstBase',in_files[i]))==1){
    colVec <- c(
      rep('black',35),
      rep(c('darkgreen','orange','navyblue'),55),
      rep('black',50),
      rep(c('darkgreen','orange','navyblue'),45),
      rep('black',65)
    )
  } else {
    colVec <- c(
      rep('black',65),
      rep(c('darkgreen','orange','navyblue'),45),
      rep('black',50),
      rep(c('darkgreen','orange','navyblue'),55),
      rep('black',35)
    )
  }
  par(mar=c(5.1,4.1,4.1,2.1))
  plot(
    c(c_mat_mean[1:200],rep(NA,50),c_mat_mean[201:400]),
    col = colVec,
    type='h',
    xlab = 'Position [nt]',
    ylab = 'Counts',
    las=1,
    xaxt = 'n',
    main = gsub('.sort.bam','',basename(in_files[i]))
  )
  axis(
    1, 
    at = c(1,51,100,150,200,251,301,351,400,450), 
    label = c('-50','START','50','100','150','-150','-100','-50','STOP','+50'))
  dev.off()
}
cat('DONE \n')
