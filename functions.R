### general functions

# function to summarize into reading frame 0,1,2
frame_counts = function(count, dir = '+', utr_length = 50){
  # check length
  if(length(count) > (2*utr_length)){
    # check if dividable by 3
    if( (length(count)-(2*utr_length))%%3 == 0){
      # if negative strand --> reverse
      if(dir == '-'){
        count <- rev(count)
      }
      # cut utr
      count <- count[(utr_length+1):(length(count)-utr_length)]
      # count frames
      frame <- c(
        sum(count[seq(1,length(count),by=3)]),
        sum(count[seq(2,length(count),by=3)]),
        sum(count[seq(3,length(count),by=3)])
      )
    } else {
      print('Not dividable by 3! No frame count performed!')
      NULL
    }
  }
}

# function for getting start and stop region normalized
norm_counts_start_stop = function(count, dir = '+', utr_length = 50, cds_length = 100, normalize = TRUE) {
  # check length
  if(length(count) > (cds_length + (2 * utr_length))){
    # if negative strand --> reverse
    if(dir == '-'){
      count <- rev(count)
    }
    if(normalize){
      # get mean 
      c_mean <- mean(count[(utr_length+1):(length(count)-utr_length)], na.rm = T)
      if(c_mean > 0){
        # normalize by mean
        count <- count / c_mean
        # return 
        c(count[1:(utr_length + cds_length)],count[(length(count) - cds_length - utr_length + 1):length(count)])
      } else {
        NULL
      }
    } else {
      # return 
      c(count[1:(utr_length + cds_length)],count[(length(count) - cds_length - utr_length + 1):length(count)])
    }
  }
}

# parse in dir
parse_dir = function(in_dir){
  # replace double //
  in_dir <- gsub("//","/",in_dir)
  # add / at end
  in_dir <- paste0(
   gsub("/$","",in_dir, perl = T),
   '/'
  )
  return(in_dir)
}

# get max cores to use from config file
get_max_cores = function(){
  # read
  max_cores <- read.table(
    'config/cores_max.csv', 
    header = F,
    stringsAsFactors = F
    )
  # return
  if(is.integer(max_cores[1,2])){
    return(max_cores[1,2])
  } else {
    cat(' Could not detect configured number of cores. Using 1 core only.\n')
    return(1)
  }
}

# load tool location
load_tool_location = function(){
  # load tool and location
  conf_t <- read.table('config/tools_location.csv', header = T, stringsAsFactors = F, sep = ",")
  # check if parameter_name 's are OK
  t_names <- c(
    'cutadapt',
    'bowtie',
    'samtools',
    'bedtools'
  )
  if(length(t_names) != sum(is.element(t_names,conf_t[,1])) ){
    stop('config/tools_location.csv not OK! tool name(s) are missing!')
  } else {
    # build list
    l_conf_t <- as.list(conf_t[,2])
    names(l_conf_t) <- conf_t[,1]
    return(l_conf_t)
  }
}

# check mapping_config.csv for validity
load_mapping_config = function(){
  # load conf
  conf_m <- read.table('config/mapping_config.csv', header = T, stringsAsFactors = F, sep = ",")
  # check if parameter_name 's are OK
  p_names <- c(
    'cutadapt_3prime_adapter',
    'cutadapt_5prime_adapter',
    'cutadapt_quality_trimming',
    'cutadapt_min_read_length',
    'bowtie_m',
    'bowtie_v'
  )
  if(length(p_names) != sum(is.element(p_names,conf_m[,1])) ){
    stop('config/mapping_config.csv not OK! parameter_name(s) are missing!')
  } else {
    # build list
    l_conf_m <- as.list(conf_m[,2])
    names(l_conf_m) <- conf_m[,1]
    return(l_conf_m)
  }
}