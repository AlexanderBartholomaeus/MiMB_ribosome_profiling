### script to pre-process and map/align sequencing reads

# load libs
source('functions.R')

# load parameters
conf <- load_mapping_config()

# load tool locations
tool_path <- load_tool_location()

# set max cores
max_cores <- get_max_cores() # from functions.R

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
  pattern = '.fastq',
  ignore.case = TRUE
)
if(length(in_files) < 1){
  stop(
    paste0('No FASTQ file(s) found!')
  )
}

# check output dir
out_path <- parse_dir(get_dir[2])
if(!dir.exists(out_path)){
  dir.create(out_path)
}
if(!dir.exists(paste0(out_path,'mapping/'))){
  dir.create(paste0(out_path,'mapping/'))
}
if(!dir.exists(paste0(out_path,'bowtie_index/'))){
  dir.create(paste0(out_path,'bowtie_index/'))
}
if(!dir.exists(paste0(out_path,'adapterRemoval/'))){
  dir.create(paste0(out_path,'adapterRemoval/'))
}

# check reference fasta input
fa_file <- get_dir[3]
if(!file.exists(fa_file)){
  stop('No FASTA file given!')
} else if(grep(".(fa|fasta|fna|mfa)$", fa_file, ignore.case = TRUE) != 1){
  stop('Given FASTA file NOT with fasta | fa | mfa | fna file end!')
}

### adapter removal
# activate 5' adapter removal if set
if(conf['cutadapt_5prime_adapter'] != ''){
  five_prime <- paste0(' -g ',conf['cutadapt_5prime_adapter'])
} else {
  five_prime <- ''
}
# remove adapters
cat('Adapter removal: \n')
for(i in 1:length(in_files)){
  # status print
  cat(paste0('  ',in_files[i]," ... \n"))
  # run
  system(
    paste0(
      tool_path['cutadapt'],
      ' --cores ',max_cores,
      ' -m ',conf['cutadapt_min_read_length'],
      ' -q ',conf['cutadapt_quality_trimming'],
      ' -a ',conf['cutadapt_3prime_adapter'],
      five_prime,
      ' -o ',out_path,'adapterRemoval/',in_files[i],
      ' ',in_path,in_files[i],
      ' > ',out_path,'adapterRemoval/',gsub(".fastq.*","",in_files[i]),'_info.txt'
    )
  )
  # status print
  cat('DONE \n')
}

### read mapping / alignment 
# create index if not present
cat('Creating index for mapping \n')
system(
  paste0(
    tool_path['bowtie'],'-build ',
    ' -q --threads ',max_cores,
    ' -f ',fa_file,
    ' ',out_path,'bowtie_index/index'
  )
)

# mapping
cat('Mapping sequencing reads: \n')
for(i in 1:length(in_files)){
  # status print
  cat(paste0('  ',in_files[i],' ... '))
  # run
  system(
    paste0(
      tool_path['bowtie'],
      ' -t',
      ' -p ',max_cores,
      ' -v ',conf['bowtie_v'],
      ' -m ',conf['bowtie_m'],
      ' -S ',out_path,'bowtie_index/index',
      ' ',out_path,'adapterRemoval/',in_files[i],
      ' 2> ',out_path,'mapping/',gsub(".fastq.*","",in_files[i]),'_info.txt',
      ' | samtools sort -O bam',
      ' -@ ',max_cores,
      ' -T ',out_path,'mapping/',gsub(".fastq.*","",in_files[i]),'.sort',
      ' -o ',out_path,'mapping/',gsub(".fastq.*","",in_files[i]),'.sort.bam'
    )
  )
  # status print
  cat('DONE \n')
}


