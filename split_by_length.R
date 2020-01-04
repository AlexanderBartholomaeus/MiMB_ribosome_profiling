### script to split by read length (to be used for calibration)

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
  pattern = '.sort.bam',
  ignore.case = TRUE
)
if(length(in_files) < 1){
  stop('No BAM files found!')
}

# check output dir
out_path <- parse_dir(get_dir[2])
if(!dir.exists(out_path)){
  dir.create(out_path)
}
if(!dir.exists(paste0(out_path,'calibration/'))){
  dir.create(paste0(out_path,'calibration/'))
}
if(!dir.exists(paste0(out_path,'calibration/split_by_length'))){
  dir.create(paste0(out_path,'calibration/split_by_length'))
}

# check length to filter
len <- get_dir[3]
if(is.na(len) || is.null(len) || len == '' ){
  len <- c(24,30)
} else {
  if(grep('-',len) == 1){
    len <- unlist(strsplit(len, '-'))
    if(length(len) != 2){
      stop('Error while parsing read length. Please use " - " (minus) as delimiter! For example: 26-28.')
    }
    len <- as.numeric(len)
    len <- len[!is.na(len)]
    if(length(len) != 2){
      stop('Seems that you did not enter numeric values for the read length.')
    }
  } else {
    stop('Wrong read length delimiter. Please use " - " (minus) ! For example: 26-28.')
  }
}

# build vec for parallel exec
par_files <- rep(in_files, each = (len[2]-len[1]+1) )
par_read_length <- rep(len[1]:len[2], length(in_files))
# all length
#par_files <- rep(in_files)
#par_read_length <- rep('all', length(in_files))

### split by length and take first nucleotide only
cat('Splitting reads by length and modifying ... ')
for_res <- foreach(i = 1:length(par_files)) %dopar% {
  system(
    paste0(
      tool_path['samtools'],' view -h ',
      in_path,par_files[i],
      " | awk '$1 ~ /@/ || length($10) == ",par_read_length[i],"'",
      " | awk '{
if($1 ~ /@/)
  print $0;
else if($2 == 0) 
  print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t1M\t\" $7 \"\t\" $8 \"\t\" $9 \"\t\" substr($10,1,1) \"\t\" substr($11,1,1); 
else if($2 == 16) 
  print $1 \"\t\" $2 \"\t\" $3 \"\t\" ($4+length($10)-1) \"\t\" $5 \"\t1M\t\" $7 \"\t\" $8 \"\t\" $9 \"\t\" substr($10,length($10),1) \"\t\" substr($11,length($11),1);  
}'",
     ' | ',tool_path['samtools'],' view -S -b - > ',
     out_path,'calibration/split_by_length/',gsub(".sort.bam","",par_files[i]),'_',par_read_length[i],'_firstBase.sort.bam'
    )
  )
}
cat('DONE \n')

### create config file for later read calibration
cali_config <- data.frame(matrix(NA, nrow = (len[2]-len[1]+1), ncol = length(in_files)+1))
colnames(cali_config) <- c('length',gsub(".sort.bam","",in_files))
cali_config[,1] <- len[1]:len[2]
write.table(
  cali_config, 
  paste0(out_path,'calibration/calibration_config.csv'),
  col.names = T,
  row.names = F,
  sep = ','
)
