### script to generate BED (6 column) from GTF/GFF

# load
source('functions.R')

# get args
get_dir <- commandArgs(trailingOnly = TRUE)

# check input file
in_gtf <- get_dir[1]
if(grep('.(gtf|gff[3]{0,1})$', in_gtf, ignore.case = T) != 1){
  stop(
    paste0(in_gtf, ' without gtf or gff file ending!') 
  )
}

# check output dir
out_path <- parse_dir(get_dir[2])
if(!dir.exists(out_path)){
  dir.create(out_path)
}

# read
tab <- read.table(in_gtf, sep='\t', stringsAsFactors = F, header=T)

# check columns
if(ncol(tab) < 9){
  stop('Less than 9 columns provided! Please provide correct GTF/GFF file.')
}

# create bed
bed <- cbind(tab[,c(1,4,5,9)],'.',tab[,7])

# correct to 0-based left position
bed[,2] <- bed[,2]-1 # correct for 0-base start in bed file

# check rows with NA or without name
no_name_idx <- unique(
  which(is.na(bed[,4]) == TRUE),
  which(bed[,4] == '')
)
if(length(no_name_idx) > 0){
  # set names
  bed[no_name_idx, 4] <- paste0(
    'featureName',
    1:length(no_name_idx)
  )
  # print message
  cat('Some rows do not contain names. Names set automatically to featureName1...n \n')
}

# write
write.table(
  bed,
  paste0(out_path,'/gtf_to_bed.bed'),
  row.names = F,
  col.names = F,
  sep="\t",
  quote=F
)
# print message
cat('BED written \n')
cat('DONE \n')

