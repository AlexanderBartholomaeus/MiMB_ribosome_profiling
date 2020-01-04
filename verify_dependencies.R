# read config file
conf <- read.table(
  file = 'config/tools_location.csv',
  header = T,
  stringsAsFactors = F,
  sep = ','
)

# test if tools are named correctly and appear
tools <- c(
  'cutadapt',
  'bowtie',
  'samtools'
  )
if(sum(is.element(
  tools,
  conf[,1])) != length(tools)) {
  stop(
    'The tools name in the first column of the tools_location.csv were renamed. 
    The pipeline needs the correct names to identify the tool. If you feel that
    you damaged the file, please check config/default/tools_location.csv for 
    another copy of this file.')
}

# test if tool can be found
for(i in 1:nrow(conf)){
  if(Sys.which(conf[i,2]) != '') {
    cat(
      paste0(
        conf[i,1]," ... OK \n"
      )
    )
  } else {
    cat(
      paste0(
        conf[i,1]," ... not found!!! Please check the location or install it! \n"
      )
    )
  }
}
