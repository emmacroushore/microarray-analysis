workdir=$(pwd)
cd $workdir
mkdir GSE164372  # see dataset https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164372

module load stack/2022.2-base_arch # See https://uiowa.atlassian.net/wiki/spaces/hpcdocs/pages/76513504/2022.2+stack for the list of modules loaded
module load sra-tools/3.0.3_gcc-9.5.0
cd $workdir/GSE164372

# select SRR files of interest
for (( i = 43; i <= 60; i++ ))
  do
  prefetch SRR133761$i
done

# convert SRR files to fastq.gz and download
for (( i = 43; i <= 60; i++ ))
  do
  parallel-fastq-dump --sra-id SRR133761$i --threads 8 --outdir out/ --split-files --gzip 
done
