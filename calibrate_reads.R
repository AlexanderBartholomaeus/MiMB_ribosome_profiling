### script to calibrate reads based on manual inspection

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

# check config for calibration
in_conf <- get_dir[2]
if(grep('.csv$',in_conf) != 1){
  stop(paste0(in_conf, ' is not a CSV file!'))
}
# read and check
cali_conf <- read.table(in_conf, header = T, stringsAsFactors = F, sep = ',')
if(nrow(cali_conf) < 1){
  stop(paste0(in_conf, ' has no rows! No read length  processed!'))
}
in_f <- gsub(".sort.bam","",in_files)
in_f2 <- in_f[is.element(in_f,colnames(cali_conf))]
cat(paste0(
  length(in_f), ' / ', length(in_f),' found in calibration config file! \n'
))

# check output dir
out_path <- parse_dir(get_dir[3])
if(!dir.exists(out_path)){
  dir.create(out_path)
}
if(!dir.exists(paste0(out_path,'calibration/'))){
  dir.create(paste0(out_path,'calibration/'))
}
if(!dir.exists(paste0(out_path,'calibration/calibrated/'))){
  dir.create(paste0(out_path,'calibration/calibrated/'))
}

# build vec for parallel exec
par_files <- rep(in_f2, each = nrow(cali_conf) )
par_read_length <- rep(cali_conf[,1], length(in_f2))

### split by length and take first nucleotide only
cat('Calibrating reads by length and determined nucleotide position ... ')
for_res <- foreach(i = 1:length(par_files)) %dopar% {
  # get nucleotide to take
  nt <- cali_conf[cali_conf[,1]==par_read_length[i],colnames(cali_conf)==par_files[i]]
  # check if NA or > 0
  if(!is.na(nt) && !is.null(nt) && nt > 0){
    # run
    system(
      paste0(
        tool_path['samtools'],' view -h ',
        in_path,par_files[i],'.sort.bam',
        " | awk '$1 ~ /@/ || length($10) == ",par_read_length[i],"'",
        " | awk '{
  if($1 ~ /@/)
    print $0;
  else if($2 == 0) 
    print $1 \"\t\" $2 \"\t\" $3 \"\t\" ($4-1+",nt,") \"\t\" $5 \"\t1M\t\" $7 \"\t\" $8 \"\t\" $9 \"\t\" substr($10,",nt,",1) \"\t\" substr($11,",nt,",1); 
  else if($2 == 16) 
    print $1 \"\t\" $2 \"\t\" $3 \"\t\" ($4+length($10)-",nt,") \"\t\" $5 \"\t1M\t\" $7 \"\t\" $8 \"\t\" $9 \"\t\" substr($10,(length($10)+1-",nt,"),1) \"\t\" substr($11,(length($11)+1-",nt,"),1);  
  }'",
        ' | ',tool_path['samtools'],' view -S -b - > ',
        out_path,'calibration/calibrated/',gsub(".sort.bam","",par_files[i]),'_',par_read_length[i],'_calibrated.bam'
      )
    )
  }
}
cat('DONE \n')
