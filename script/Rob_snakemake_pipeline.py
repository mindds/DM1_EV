import os
import os.path
from os import path
import sys
import getpass
import shutil
import time
import pandas as pd

## RNA-seq pipeline usage:
help = '''

RNA-seq pipeline

An example command to run the pipeline is:
    BASE=/space/mindds/1/piplines/rnaseq-pipeline/rnaseq-pipeline
    PATH_OUT=/tmp
    snakemake -s $BASE/snakemake_pipeline.py \
              --cores 2 \
              --directory=$PATH_OUT \
              --configfile="$BASE/TESTING/config_mac.yaml" 
              
Add the '-n -p' flags to do a dry run (and print the shell commands)

The two required input files are:
    - a yaml config file (passed to snakemake via the --configfile parameter)
    - a csv sample table (specified in the yaml config file via the "sample_table" : [path/tp/sampleTable.csv])

The output is written to a subdirectory of --directory=/path/to/output

An example of a valid config file is:
    "path_S3_in" : "s3://kitchen-mgh-data/Projects/TestData/RNAseq/Fastq"
    "path_S3_out" : "s3://kitchen-mgh-data/Projects/TestData/RNAseq/Processed"
    "sample_table" : "/Users/rk504/Pipelines/rnaseq-pipeline/TESTING/sampleTable.csv"
    "path_annotation" : "/Users/rk504/Dropbox_MGH/Annotations/Human"
    "path_tmp" : "/tmp/RNAseq"
    "rapid_run" : "True"
    "additional_args_salmon" : "--some other salmon args"
    "force" : "False"

And the sample table must be formatted as such:
    sampleID,read,fastqName
    test,R1,test_R1_001.fq.gz
    test,R2,test_R2_001.fq.gz

In the sample table it is imperative that:
    - the 'read' is specified exactly as either 'R1' or 'R2' or the pipeline will fail
    - the fastqName is present at the S3 URI specified in 'path_S3_in' in the config file
'''


## hard-coded parameters
path_bbmap = "/efs/bin/bbtools/bbmap"
exe_salmon = "/efs/bin/salmon/salmon-0.13.1_linux_x86_64/bin/salmon"
exe_samtools = "/efs/bin/samtools-1.9/samtools"
exe_bedtools = "/efs/bin/bedtools2/bin/bedtools"
path_pipelineScripts = "/efs/Pipelines/RNAseq_pipeline/scripts"
job_name_sep = "-"


## Read parameters from the config file
path_S3_in = config['path_S3_in']
path_S3_out = config['path_S3_out']
path_EFS_out = config['path_EFS_out']
path_ann = config["path_annotation"]

path_tmp = "/tmp"
if "path_tmp" in config:
    path_tmp = config["path_tmp"]

run_rapid = "False"
if "rapid_run" in config:
    run_rapid = config["rapid_run"]

additional_args_salmon = ""
if "additional_args_salmon" in config:
    additional_args_salmon = config["additional_args_salmon"]

forceRerun = "False"
if "force" in config:
    forceRerun = config["force"]


## Get a list of the samples and read numbers we need to process
samples = pd.read_csv(config["sample_table"], comment='#')
unique_sampleIDs_orig = samples['sampleID'].unique()
unique_readNumbers = samples['read'].unique()

if not (all(elem in unique_readNumbers for elem in ['R1','R2'])):
    sys.exit("FATAL ERROR: Read numbers must be specified as 'R1' and/or 'R2' - please address this in the sample_table")



## Get the system user
sys_username = getpass.getuser()


##
## Check for existing output
##
print("Processing samples:")
unique_sampleIDs = []
for s in unique_sampleIDs_orig:
    p0 = path_EFS_out + "/" + s + "/checkpoints/cleanup_and_sync.chk"
    p1 = path_EFS_out + "/" + s + "/salmon/quant.sf"
    p2 = path_EFS_out + "/" + s + "/salmon/quant.genes.sf"
    if path.exists(p0) and path.exists(p1) and path.exists(p2):
        if forceRerun == "True":
            unique_sampleIDs.append(s)
            print(s + ":\texists, queued (due to force=True)")
        else:
            print(s + ":\texists, skipping")
    else:
        unique_sampleIDs.append(s)
        print(s + ":\tqueued")



## Function to get the fastq filenames for all sample_R1/R2 combinations in the sample table
##  -> returns a list of --include= statements for downloading from the S3 bucket via 'aws s3 sync'
def get_fastq(wildcards):
    tmp = samples.loc[(samples['sampleID'] == wildcards.sampleID)  &  (samples['read'] == wildcards.read)]
    fileNames = tmp['fastqName'].values
    return fileNames


## Primary rule - determine which sub-rules below are executed
if not run_rapid == "True":
    rule all:
        input:
            expand("{sampleID}/{sampleID}.{read}_fastqc.html", sampleID=unique_sampleIDs, read=unique_readNumbers),
            expand("{sampleID}/trimmed/{sampleID}.{read}_fastqc.html", sampleID=unique_sampleIDs, read=unique_readNumbers),
            #expand("{sampleID}/salmon/{sampleID}.coverage.csv.gz", sampleID=unique_sampleIDs),
            expand("{sampleID}/checkpoints/cleanup_and_sync.chk", sampleID=unique_sampleIDs)
else:
    rule all:
        input:
            expand("{sampleID}/checkpoints/cleanup_and_sync.chk", sampleID=unique_sampleIDs)
            


rule cleanup_and_sync:
    input:
        "{sampleID}/checkpoints/sortBAM_pos.chk",
        #"{sampleID}/checkpoints/sortBAM_name.chk",
        "{sampleID}/checkpoints/calculate_stats.chk"
    output:
        "{sampleID}/checkpoints/cleanup_and_sync.chk"
    params:
        threads=1,
        usethreads=1,
        sampleID="{sampleID}",
        runtime="00:30:00",
        priority=10,
        name="{sampleID}"+job_name_sep+"clean_up"
    shell:
        '''
        rm -f {params.sampleID}/{params.sampleID}.R*.fastq.gz \
        && rm -f -r {params.sampleID}/rawFastq \
        && rm -f -r {params.sampleID}/trimmed \
        && rm -f {params.sampleID}/salmon/{params.sampleID}.bam \
        && aws s3 sync {params.sampleID} {path_S3_out}/{params.sampleID} --quiet \
        && rm -f {params.sampleID}/salmon/*.bam* \
        && mkdir -p {path_EFS_out} \
        && cp -r {params.sampleID} {path_EFS_out}/ \
        && touch {output} \
        && touch {path_EFS_out}/{output}
        '''

rule coverage:
    input:
        "{sampleID}/checkpoints/sortBAM_pos.chk"
    output:
        "{sampleID}/salmon/{sampleID}.coverage.csv.gz"
    params:
        threads=1,
        usethreads=1,
        sampleID="{sampleID}",
        runtime="02:00:00",
        priority=9,
        name="{sampleID}"+job_name_sep+"coverage"
    log:
        "{sampleID}/logs/coverage.log"
    shell:
        '''
        {exe_bedtools} coverage -sorted -d \
            -g {path_ann}/transcripts.bedtoolsGenome.txt \
            -a {path_ann}/transcripts.bed \
            -b {params.sampleID}/salmon/{params.sampleID}.sortedPos.bam \
        | {path_pipelineScripts}/coverage_slidingWindow.py - \
        | gzip -c \
        > {output}
        '''
## Todo: Add to coverage:
# cat sampleA/salmon/quant.sf | sed '1d' | awk -F '\t' '{if($5 >= 1){print $1}}' > tmp_tx.txt \

rule sortBAM_name:
    input:
        "{sampleID}/checkpoints/quantify_transcripts.chk"
    output:
        "{sampleID}/checkpoints/sortBAM_name.chk"
    params:
        threads=8,
        usethreads=8,
        sampleID="{sampleID}",
        runtime="01:00:00",
        priority=9,
        name="{sampleID}"+job_name_sep+"sortBAM_name"
    log:
        "{sampleID}/logs/sortBAM_name.log"
    shell:
        '''
        mkdir -p {path_tmp}/{params.sampleID} \
        && {exe_samtools} sort -n --threads {params.usethreads} -m 1G \
            -T {path_tmp}/{params.sampleID} \
            -o {params.sampleID}/salmon/{params.sampleID}.sortedName.bam \
            {params.sampleID}/salmon/{params.sampleID}.bam \
            >> {log} 2>&1 \
        && touch {output}
        '''

rule sortBAM_pos:
    input:
        "{sampleID}/checkpoints/quantify_transcripts.chk"
    output:
        "{sampleID}/checkpoints/sortBAM_pos.chk"
    params:
        threads=8,
        usethreads=8,
        sampleID="{sampleID}",
        runtime="01:00:00",
        priority=9,
        name="{sampleID}"+job_name_sep+"sortBAM_pos"
    log:
        "{sampleID}/logs/sortBAM_pos.log"
    shell:
        '''
        mkdir -p {path_tmp}/{params.sampleID} \
        && {exe_samtools} sort --threads {params.usethreads} -m 1G \
            -T {path_tmp}/{params.sampleID} \
            -o {params.sampleID}/salmon/{params.sampleID}.sortedPos.bam \
            {params.sampleID}/salmon/{params.sampleID}.bam \
            >> {log} 2>&1 \
        && {exe_samtools} index {params.sampleID}/salmon/{params.sampleID}.sortedPos.bam >> {log} 2>&1 \
        && touch {output}
        '''


rule calculate_stats:
    input:
        "{sampleID}/checkpoints/quantify_transcripts.chk"
    output:
        chk="{sampleID}/checkpoints/calculate_stats.chk",
        stats="{sampleID}/sampleStats.csv"
    params:
        threads=1,
        usethreads=1,
        sampleID="{sampleID}",
        runtime="00:10:00",
        priority=10,
        name="{sampleID}"+job_name_sep+"calculate_stats"
    log:
        "{sampleID}/logs/calculate_stats.log"
    shell:
        '''
        echo "What,ReadCount" > {output.stats} \
        && cat {params.sampleID}/logs/trim_adapters.log | grep -w "^Input:" \
            | awk -F " " '{{print "Input,"$2/2}}' >> {output.stats} \
        && cat {params.sampleID}/logs/trim_adapters.log | grep -w "^Result:" \
            | awk -F " " '{{print "TrimmedFiltered,"$2/2}}' >> {output.stats} \
        && cat {params.sampleID}/salmon/aux_info/meta_info.json | grep -w "num_mapped" \
            | sed 's/,//g' | awk -F " " '{{print "TranscriptomeMapped,"$2}}' >> {output.stats} \
        && touch {output.chk}
        '''


rule quantify_transcripts:
    input:
        R1="{sampleID}/checkpoints/trim_adapters_R1.chk",
        R2="{sampleID}/checkpoints/trim_adapters_R2.chk"
    output:
        "{sampleID}/checkpoints/quantify_transcripts.chk"
    params:
        threads=16,
        usethreads=16,
        sampleID="{sampleID}",
        runtime="02:00:00",
        priority=10,
        name="{sampleID}"+job_name_sep+"quantify_transcripts"
    log:
        "{sampleID}/logs/quantify_transcripts_salmon.log"
    run:
        shell('mkdir -p {path_tmp}/{params.sampleID}')
        shell('gunzip -c {path_ann}/annotation.gtf.gz > {path_tmp}/{params.sampleID}/annotation.gtf')
        shell('''
        {exe_salmon} quant \
            -i {path_ann}/salmonIndex_31nt -g {path_tmp}/{params.sampleID}/annotation.gtf \
            --threads={params.usethreads} --seqBias --gcBias -l A --numBootstraps=20 \
            {additional_args_salmon} \
            --validateMappings --allowDovetail --writeMappings={path_tmp}/{params.sampleID}/{params.sampleID}.sam \
            -o $PWD/{params.sampleID}/salmon \
            -1 $PWD/{params.sampleID}/trimmed/{params.sampleID}.R1.fastq.gz \
            -2 $PWD/{params.sampleID}/trimmed/{params.sampleID}.R2.fastq.gz \
            >>{log} 2>>{log} 
        ''')
        shell('touch {params.sampleID}/checkpoints/quantify_transcripts_salmon.chk')
        shell('''
        {exe_samtools} view -bS --threads {params.usethreads} \
            -o {params.sampleID}/salmon/{params.sampleID}.bam \
            {path_tmp}/{params.sampleID}/{params.sampleID}.sam 
        ''')
        shell('rm {path_tmp}/{params.sampleID}/{params.sampleID}.sam')
        shell('touch {params.sampleID}/checkpoints/quantify_transcripts.chk')


rule fastQC_trimmed:
    input:
        "{sampleID}/checkpoints/trim_adapters_{read}.chk"
    output:
        "{sampleID}/trimmed/{sampleID}.{read}_fastqc.html"
    params:
        threads=4,
        usethreads=4,
        sampleID="{sampleID}",
        read="{read}",
        runtime="01:00:00",
        priority=1,
        name="{sampleID}_{read}"+job_name_sep+"fastQC_trimmed"
    log:
        "{sampleID}/logs/fastQC_trimmed_{read}.log"
    shell:
        '''
        docker run -u `id -u {sys_username}` -v $PWD/{params.sampleID}:/base rkitchen/fastqc \
            -o /base/trimmed/ -t {params.usethreads} /base/trimmed/{params.sampleID}.{params.read}.fastq.gz \
            >> {log} 2>&1
        '''

rule trim_adapters:
    input:
        R1="{sampleID}/checkpoints/combine_fastqs_R1.chk",
        R2="{sampleID}/checkpoints/combine_fastqs_R2.chk"
    output:
        R1="{sampleID}/checkpoints/trim_adapters_R1.chk",
        R2="{sampleID}/checkpoints/trim_adapters_R2.chk"
    params:
        threads=16,
        usethreads=16,
        memory="28g",
        sampleID="{sampleID}",
        runtime="02:00:00",
        priority=10,
        name="{sampleID}"+job_name_sep+"trim_adapters"
    log:
        "{sampleID}/logs/trim_adapters.log"
    run:
        shell('''
        {path_bbmap}/bbduk.sh \
            -Xmx{params.memory} \
            ref={path_bbmap}/resources/adapters.fa \
            in=$PWD/{params.sampleID}/{params.sampleID}.R1.fastq.gz \
            in2=$PWD/{params.sampleID}/{params.sampleID}.R2.fastq.gz \
            out=$PWD/{params.sampleID}/trimmed/{params.sampleID}.R1.fastq.gz \
            out2=$PWD/{params.sampleID}/trimmed/{params.sampleID}.R2.fastq.gz \
            ktrim=r k=21 mink=11 tbo tpe hdist=2 minlen=31 threads={params.usethreads} \
            qtrim=r trimq=10 maq=10 \
            entropy=0.3 entropywindow=50 entropyk=5 \
            >>{log} 2>>{log}
        ''')
        shell("touch {output.R1}")
        shell("touch {output.R2}")


rule fastQC_raw:
    input:
        "{sampleID}/checkpoints/combine_fastqs_{read}.chk"
    output:
        "{sampleID}/{sampleID}.{read}_fastqc.html"
    params:
        threads=4,
        usethreads=4,
        sampleID="{sampleID}",
        read="{read}",
        runtime="01:00:00",
        priority=1,
        name="{sampleID}_{read}"+job_name_sep+"fastQC_raw"
    log:
        "{sampleID}/logs/fastQC_untrimmed_{read}.log"
    shell:
        '''
        docker run -u `id -u {sys_username}` -v $PWD/{params.sampleID}:/base rkitchen/fastqc \
            -o /base/ -t {params.usethreads} /base/{params.sampleID}.{params.read}.fastq.gz \
            >> {log} 2>&1
        '''


rule combine_fastqs:
    params:
        sampleID="{sampleID}",
        read="{read}",
        files=get_fastq,
        threads=2,
        usethreads=1,
        runtime="01:00:00",
        priority=10,
        name="{sampleID}_{read}"+job_name_sep+"combine_fastqs"
    output:
        "{sampleID}/checkpoints/combine_fastqs_{read}.chk"
    #log:
    #    "{sampleID}/logs/combine_fastqs_{read}.log"
    run:
        shell("mkdir -p {params.sampleID}/logs")
        shell("mkdir -p {params.sampleID}/checkpoints")
        shell("rm -f -r {params.sampleID}/{params.sampleID}.{params.read}.fastq.gz")
        for f in params.files:
            shell('aws s3 cp "{path_S3_in}/' + f + '" - >> {params.sampleID}/{params.sampleID}.{params.read}.fastq.gz')
        shell("touch {output}")


