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

# check 5' or 3' calibration
firstLast <- get_dir[3]
if(firstLast != 'first' && firstLast != 'last' && firstLast != '5prime' && firstLast != '3prime'){
  # check if read_base is missing and length range was set
  if(grep("[0-9]{2}-[0-9]{2}",firstLast)==1){
    get_dir[4] <- get_dir[3]
  }
  # default switch to 5' calibration
  cat('No valid string to process first (5\') or last (3\') nucleotide detected. Please enter "first", "5prime", "last" or "3prime". Switching to first (5\') nucleotide processing.')
  firstLast <- 'first'
} else if(firstLast == '5prime'){
  firstLast <- 'first'
} else if(firstLast == '3prime'){
  firstLast <- 'last'
}

# check length to filter
len <- get_dir[4]
if(is.na(len) || is.null(len) || len == '' ){
  len <- c(24,30)
} else {
  if(length(grep('-',len)) == 1){
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
if(firstLast == 'first'){
  f_read <- "print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t1M\t\" $7 \"\t\" $8 \"\t\" $9 \"\t\" substr($10,1,1) \"\t\" substr($11,1,1); " 
  r_read <- "print $1 \"\t\" $2 \"\t\" $3 \"\t\" ($4+length($10)-1) \"\t\" $5 \"\t1M\t\" $7 \"\t\" $8 \"\t\" $9 \"\t\" substr($10,length($10),1) \"\t\" substr($11,length($11),1);  "
} else if(firstLast == 'last') {
  r_read <- "print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t1M\t\" $7 \"\t\" $8 \"\t\" $9 \"\t\" substr($10,1,1) \"\t\" substr($11,1,1); " 
  f_read <- "print $1 \"\t\" $2 \"\t\" $3 \"\t\" ($4+length($10)-1) \"\t\" $5 \"\t1M\t\" $7 \"\t\" $8 \"\t\" $9 \"\t\" substr($10,length($10),1) \"\t\" substr($11,length($11),1);  "
}
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
  ",f_read,"
else if($2 == 16) 
  ",r_read,"
}'",
     ' | ',tool_path['samtools'],' view -S -b - > ',
     out_path,'calibration/split_by_length/',gsub(".sort.bam","",par_files[i]),'_',par_read_length[i],'_',firstLast,'Base.sort.bam'
    )
  )
}
cat('DONE \n')

### create config file for later read calibration (5'prime)
cali_config <- data.frame(matrix(NA, nrow = (len[2]-len[1]+1), ncol = length(in_files)+1))
colnames(cali_config) <- c('length',gsub(".sort.bam","",in_files))
cali_config[,1] <- len[1]:len[2]
write.table(
  cali_config, 
  paste0(out_path,'calibration/calibration_5prime_config.csv'),
  col.names = T,
  row.names = F,
  sep = ','
)
write.table(
  cali_config, 
  paste0(out_path,'calibration/calibration_3prime_config.csv'),
  col.names = T,
  row.names = F,
  sep = ','
)
