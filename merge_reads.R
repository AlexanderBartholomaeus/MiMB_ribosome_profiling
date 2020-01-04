### script to merge and sort calibrated reads

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

# check output dir
if(is.na(get_dir[2])){
  out_path <- in_path
} else {
  out_path <- parse_dir(get_dir[2])
}
if(!dir.exists(out_path)){
  dir.create(out_path)
}

# get base names
in_f <- unique(gsub("_[0-9]{2}_calibrated.bam","",in_files))

### merge
cat('Merging read files ... ')
#for_res <- foreach(i = 1:length(in_f)) %dopar% {
for(i in 1:length(in_f)) {
  # get files to merge
  f_merge <- paste0(in_path,grep(in_f[i], in_files, value = T))
  # merge
  # run
  system(
    paste0(
      tool_path['samtools'],' merge ',
      out_path,in_f[i],'_calibrated.bam ',
      paste(f_merge, collapse = ' ')
    )
  )
}  
cat('DONE \n')