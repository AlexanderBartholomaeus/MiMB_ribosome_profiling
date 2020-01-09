### script to calculate read length distribution on aligned reads (BAM files)

# load libs
library(foreach)
library(doParallel)
source('functions.R')

# load tool locations
tool_path <- load_tool_location()

# set max cores
max_cores <- get_max_cores() # from functions.R
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
  stop(
    paste0('No BAM files found!')
  )
}

# check output dir
out_path <- parse_dir(get_dir[2])
if(!dir.exists(out_path)){
  dir.create(out_path)
}
if(!dir.exists(paste0(out_path,'read_length_distribution/'))){
  dir.create(paste0(out_path,'read_length_distribution/'))
}

### count read length
cat('Count read length ... ')
for_res <- foreach(i = 1:length(in_files)) %dopar% {
  system(
    paste0(
      tool_path['samtools']," view -F 4 ",
      in_path,in_files[i],
      " | awk '{print length($10)}' | sort -n | uniq -c | awk '{sub(/^[ ]+/,\"\");print }' > ",
      out_path,"read_length_distribution/",gsub(".sort.bam","",in_files[i]),".txt"
    )
  )
}
cat('DONE \n')

### generate plot 
cat('Create plot and table of read length distribution ... \n')
pdf(
  paste0(out_path,"read_length_distribution/all_read_length.pdf"),
  width = 7,
  height = 5
)
length_dist <- data.frame('length' = c(20:50))
for(i in 1:length(in_files)){
  # read file
  l_dist <- read.table(
    paste0(out_path,"read_length_distribution/",gsub(".sort.bam","",in_files[i]),".txt"),
    header = F,
    sep = " ",
    stringsAsFactors = F
  )
  # remove whitespace
  length_dist <- merge(length_dist, l_dist, by.x = 1, by.y = 2, all.x = T)
  
  # file
  par(mar=c(5.1,6.1,4.1,2.1))
  b <- barplot(
    length_dist[,i+1],
    las=1, 
    main=in_files[i], 
    xlab = 'Read length [nt]',
    ylab = '',
    cex.names=1.2, 
    cex.lab=1.2, 
    cex.axis = 1.2) 
  title(ylab = 'Counts', line = 5, cex.lab = 1.2)
  axis(1, at = b[seq(1,31,by=5)], labels = c(seq(20,50,by=5)), cex.axis = 1.2)
}
dev.off()
colnames(length_dist)[2:ncol(length_dist)] <- in_files
write.table(
  length_dist,
  paste0(out_path,'read_length_distribution/all_read_length.csv'),
  sep = "\t",
  col.names = T,
  row.names = F
)
cat('DONE \n')

