### script to filter(=keep) reads of a certain length

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
# check if fasta files
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
if(!dir.exists(paste0(out_path,'filtered_length/'))){
  dir.create(paste0(out_path,'filtered_length/'))
}

# check length to filter
len <- get_dir[3]
if(is.na(len) || is.null(len) || len == '' ){
  stop('No read length given!')
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

### filter read length
for_res <- foreach(i = 1:length(in_files)) %dopar% {
  system(
    paste0(
      tool_path['samtools'],' view -h ',
      in_path,in_files[i],
      " | awk '$1 ~ /@/ || length($10) >= ",len[1],"'",
      " | awk '$1 ~ /@/ || length($10) <= ",len[2],"'",
      ' | ',tool_path['samtools'],' view -S -b - > ',
      out_path,'filtered_length/',gsub(".sort.bam","",in_files[i]),'_',len[1],'-',len[2],'.sort.bam'
    )
  )
}
