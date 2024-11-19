#!/usr/bin/bash

#' this program is used to run salmon program 
#' 11/26/2020 Zhaozhi Li

#' @ usage: bash salmon.sh -f fastq_dir -i salmon_index
#' @ parameters: 
#' @ fastaq_dir: path to the folder that has all fastq files
#' @ salmon_index: path to the salmon index file 

#' @ salmon parameters
#' @ --validateMappings
#' Selective alignment
#' now the default mapping strategy (in version 1.0.0 forward)
#' it adopts a considerably more sensitive scheme that we have developed for finding the potential mapping 
#' loci of a read, and score potential mapping loci using the chaining algorithm introdcued in minimap2

#' @ -l
#' https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype
#' To allow Salmon to automatically infer the library type, simply provide -l A or --libType A to Salmon.

#' @ 11/26/2020
#' @ bash salmon.sh -f /space/mindds/1/ThurmanRNASeq -i /space/mindds/1/ThurmanRNASeq/gencode_salmon_index 
#' @ /space/mindds/1/ThurmanRNASeq
#bash salmon.sh -f /space/mindds/1/projects/ThurmanRNASeq/data/merged/20211129_Wheeler_Muscle_stem_cell_transcriptome -i /space/mindds/1/pipelines/rnaseq-pipeline/annotation/Mouse/GRCm39/salmonIndex_31nt &> /space/mindds/1/projects/ThurmanRNASeq/scripts/20211129_Wheeler_Muscle_stem_cell_transcriptome/logs/salmon_new.log


#' set up
#fastq_dir=''
#index_dir=''
print_usage() {
  printf "\nUsage: bash salmon.sh -f fastq_dir -i salmon_index
            args: 
                fastaq_dir: path to the folder that has all fastq files
                salmon_index: path to the salmon index file\n\n"
}

#' get options
#' ref: https://www.shellscript.sh/tips/getopts/
while getopts ':f:i:' flag; do
  case "${flag}" in
    f) fastq_dir=${OPTARG};;
    i) index_dir=${OPTARG};;
    *) print_usage
       exit 1 ;;
  esac
done 

echo ${fastq_dir}
echo ${index_dir}

fastqLabels=`find ${fastq_dir} -type f | grep 001.fastq.gz | awk '{gsub("_R[12]_001.fastq.gz", "", $0); print $0}' | sort | uniq` 
# fastqLabels=`ls $1 | grep .fastq.gz | awk '{gsub("_R[12].fastq.gz", "", $0); print $0}' | uniq` 
echo ${fastqLabels[@]}

i=1
for fastqL in $fastqLabels
do  
    # sample id
    fastqS=(${fastqL//\// }) 
    echo "processing sample $i : " ${fastqS[2]} "..."
    echo -e "cmd:  salmon quant -i ${index_dir}
                    -l A 
                    -1 ${fastqL}_R1_001.fastq.gz
                    -2 ${fastqL}_R2_001.fastq.gz
                    --validateMappings
                    --gcBias  
                    -o results/${fastqS[2]}"

/space/mindds/1/software/salmon-latest_linux_x86_64/bin/./salmon quant -i $index_dir -l A -1 $fastqL"_R1_001.fastq.gz" -2 $fastqL"_R2_001.fastq.gz" --validateMappings --gcBias -o transcripts_quants/${fastqS[2]}
    i=$((i + 1))
done

echo "done."

#salmon quant -i salmon_index -l A -1 reads1.fq -2 reads2.fq --validateMappings -o transcripts_quant
