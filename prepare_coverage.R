### script to prepare annotation BED file to coverage start and stop region

# load libs
source('functions.R')

# load tool locations
tool_path <- load_tool_location()

# get args
get_dir <- commandArgs(trailingOnly = TRUE)

# check input file
in_bed <- get_dir[1]
if(grep('.bed$',in_bed) != 1){
  stop(
    paste0(in_bed, ' is not a BED file!') 
  )
}

# check output dir
# out_path <- parse_dir(get_dir[2])
# if(!dir.exists(out_path)){
#   dir.create(out_path)
# }
# if(!dir.exists(paste0(out_path,'start_stop_region/'))){
#   dir.create(paste0(out_path,'start_stop_region/'))
# }

# load
anno <- read.table(
  in_bed,
  stringsAsFactors = F,
  sep = "\t" 
)

# modify annotation
anno <- t(sapply(1:nrow(anno),function(j){
  # check for different exons
  anno_ <- anno[j,]
  idx <- which(anno[,4]==anno_[1,4])
  # if only one exon
  if(length(idx) == 1){
    anno_[1,2] <- anno_[1,2]-50
    anno_[1,3] <- anno_[1,3]+50
    return(anno_)
  # else check if first
  } else {
    if(anno_[1,2] == min(anno[idx,2])){
      anno_[1,2] <- anno_[1,2]-50
      return(anno_)
    } else {
      return(anno_)
    }
    if(anno_[1,2] == max(anno[idx,2])){
      anno_[1,3] <- anno_[1,3]+50
      return(anno_)
    } else {
      return(anno_)
    }
  }
}))
# write
write.table(
  anno, 
  paste0(
    gsub(".bed$","",in_bed),
    '_plus_50nt.bed'
  ),
  col.names = F,
  row.names = F,
  sep = "\t",
  quote = F
)



