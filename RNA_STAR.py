
"""RNA-Seq Pipeline for Humans and Mouse"""

import sys
import os
import subprocess
import click


TEMP = '/ifs/work/leukgen/home/yellapav/star_align/tmp'


def human_b37_e75():
    global GENOME,GTF,GDIR,SALMON_INDEX,STAR_FUSION_LIB,REFFLAT_FILE,RIBOSOME_LIST
    GENOME = '/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/homo_sapiens.fasta'
    GENOME = '/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta'
    GTF = '/home/yellapav/local/resources/Homo_sapiens.GRCh37.75.gtf'
    GDIR = '/home/yellapav/local/resources/star_h37d5_ens75'
    SALMON_INDEX = '/ifs/res/leukgen/home/yellapav/resources/salmon_index/Homo_sapiens.GRCh37.75.cdna.all'
    STAR_FUSION_LIB='/ifs/e63data/armstronglab/noushin/local/Hg19_CTAT_resource_lib'
    REFFLAT_FILE="/ifs/res/leukgen/home/yellapav/resources/misc/ensembl75_b37_refFlat.txt"
    RIBOSOME_LIST="/ifs/res/leukgen/home/yellapav/resources/misc/ensembl75_b37_rrna.intervals"
    print("Using human build 37 and ensembl 75 annatonations") 
    return


def mouse_b38_e75():
    global GENOME,GTF,GDIR,SALMON_INDEX,STAR_FUSION_LIB,REFFLAT_FILE,RIBOSOME_LIST
    GENOME = '/ifs/res/leukgen/home/yellapav/resources/GRCm38_e75/Mus_musculus.GRCm38.75.dna.primary_assembly.fa'
    GTF = '/ifs/res/leukgen/home/yellapav/resources/GRCm38_e75/Mus_musculus.GRCm38.75.gtf'
    GDIR = '/ifs/res/leukgen/home/yellapav/resources/GRCm38_e75/star_GRCm38_e75'
    SALMON_INDEX = '/ifs/res/leukgen/home/yellapav/resources/GRCm38_e75/Mus_musculus_GRCm38_75_salmon_index'
    STAR_FUSION_LIB='/ifs/res/leukgen/home/yellapav/resources/star_fusion/mm38_ens75_take2'
    REFFLAT_FILE="/ifs/res/leukgen/home/yellapav/resources/misc/mm_ensembl75_b37_refFlat.txt"
    RIBOSOME_LIST="/ifs/res/leukgen/home/yellapav/resources/misc/mm_ensembl75_b37_rrna.intervals"
    print("Using mouse build 38 and ensembl 75 annatonations")
    return





### RG Tags if we choose to have picard mark Duplicates or GATK run.
### Currently these are dummy values 
RPL="Illumina"
RLB="K1234"
ID="Kdgfgfhdd"



### Add RG Tags
@click.command()
@click.option('--bam', required=True, help='Bam File')
@click.option('--sample_name', required=True, help='Sample Name')
def add_rgid(bam, sample_name):
    """Add read-groups to bam."""
    cmd0 = [
        "java", " -Xmx4g", " -jar", " /opt/common/CentOS_6-dev/picard/v1.140/picard.jar",
        " AddOrReplaceReadGroups", " I=", bam ," O=", sample_name ,".rg.bam",
        " SO=coordinate", " RGID=", ID, " RGLB=", RLB, " RGPL=", RPL, " RGPU=machine",
        " RGSM=", sample_name,
        " VALIDATION_STRINGENCY=SILENT" 
        ]

    print("".join(cmd0))
    cmd0_run=("".join(cmd0))
    subprocess.run(cmd0_run,check=True,shell=True)



### STAR Alignment
@click.command()
@click.option('--fastq1', required=True, help='Fastq Read 1')
@click.option('--fastq2', required=True, help='Fastq Read 2')
@click.option('--species', required=True, help='Type human or mouse')
def star_align(fastq1, fastq2, species):
    """Align using STAR."""

    if species=="human":
        human_b37_e75()
    elif species=="mouse":
        mouse_b38_e75()
    else:
        print("Species should be either human or mouse")

    cmd1 = [
        "STAR", "--genomeDir", GDIR,
        "--sjdbGTFfile", GTF,
        "--readFilesIn", fastq1, fastq2,
        "--outFilterMismatchNmax 10",
        "--seedSearchStartLmax 30",
        "--chimSegmentMin 15",
        "--chimJunctionOverhangMin 15",
        "--outFilterType BySJout",
        "--runThreadN 8",
        "--outReadsUnmapped Fastx",
        "--genomeLoad NoSharedMemory",
        "--outSAMstrandField intronMotif",
        "--outSAMtype BAM SortedByCoordinate",
        "--outSAMmode Full"
        ]

    cmd1_run=(" ".join(cmd1))
    subprocess.run(cmd1,check=True)



### Picard Alignment Summary
@click.command()
@click.option('--bam', required=True, help='Bam File')
@click.option('--sample_name', required=True, help='Sample Name')
@click.option('--species', required=True, help='Type human or mouse')
def alignment_metrics(bam, sample_name, species):
    """Get alingment metrics."""

    if species=="human":
        human_b37_e75()
    elif species=="mouse":
        mouse_b38_e75()
    else:
        print("Species should be either human or mouse")

    cmd2 = [
        "java", " -jar", " /opt/common/CentOS_6-dev/picard/v1.140/picard.jar",
        " CollectAlignmentSummaryMetrics INPUT=",bam," OUTPUT=",sample_name,".picard.alignment_summary_metrics REFERENCE_SEQUENCE=",GENOME," TMP_DIR=",TEMP,
        " VALIDATION_STRINGENCY=LENIENT"
        ]

    print("".join(cmd2))
    cmd2_run=("".join(cmd2))
    subprocess.run(cmd2_run,check=True,shell=True)



### Samtools Flagstat
@click.command()
@click.option('--bam', required=True, help='Bam File')
@click.option('--sample_name', required=True, help='Sample Name')
def samtools_flagstat(bam, sample_name):
    """Samtools flagstat metrics."""
    cmd3 = [
        "samtools","flagstat", bam," > ", sample_name,".samStats"
        ]

    print("".join(cmd3))
    subprocess.run(cmd3,check=True)



### HTSEq which could be skipped since these counts can be recapitulted from STAR
@click.command()
@click.option('--bam', required=True, help='Bam File')
@click.option('--sample_name', required=True, help='Sample Name')
@click.option('--species', required=True, help='Type human or mouse')
def counts_htseq(bam, sample_name, species):
    """Generate gene level counts using htseq."""

    if species=="human":
        human_b37_e75()
    elif species=="mouse":
        mouse_b38_e75()
    else:
        print("Species should be either human or mouse")


    cmd4 = [
        "samtools view ", bam ," > ", sample_name,".sam;",
        "htseq-count -m intersection-strict -s no ", sample_name,".sam ", GTF," > ", sample_name,".counts; ",
        "rm ", sample_name,".sam"
        ]
    cmd4_run=("".join(cmd4))
    print(cmd4_run)
    subprocess.run(cmd4_run,check=True,shell=True)



### Salmon
@click.command()
@click.option('--fastq1', required=True, help='Fastq Read 1')
@click.option('--fastq2', required=True, help='Fastq Read 2')
@click.option('--output', required=True, help='Output file name')
@click.option('--threads', default="1", help='Number of threads to use')
@click.option('--species', required=True, help='Type human or mouse')
def salmon(fastq1, fastq2, threads, output, species):
    """Generate transcript level expression."""
    new_fastq1=fastq1.replace(".fastq", ".fastq.gz")
    new_fastq1=new_fastq1.replace(",", " ")
    new_fastq1=new_fastq1.rstrip()
    new_fastq2=fastq2.replace(".fastq", ".fastq.gz")
    new_fastq2=new_fastq2.replace(",", " ")
    new_fastq2=new_fastq2.rstrip()

    if species=="human":
        human_b37_e75()
    elif species=="mouse":
        mouse_b38_e75()
    else:
        print("Species should be either human or mouse")



    cmd5 = [
        "gzip *.fastq;",
        "salmon quant -i ", SALMON_INDEX," -l A ",
        "-1 ", new_fastq1," ",
        "-2 ", new_fastq2," ",
        "-p ", threads, " -o quants/", output
        ]

    cmd5_run=("".join(cmd5))
    print(cmd5_run)
    subprocess.run(cmd5_run, check=True, shell=True)




### Picard Collect RNA Metrics 
@click.command()
@click.option('--bam', required=True, help='Bam File')
@click.option('--sample_name', required=True, help='Sample Name')
@click.option('--species', required=True, help='Type human or mouse')
def collect_rna_metrics(bam, sample_name, species):
    """Get RNA metrics."""

    if species=="human":
        human_b37_e75()
    elif species=="mouse":
        mouse_b38_e75()
    else:
        print("Species should be either human or mouse")


    cmd6 = [
        "java -jar /opt/common/CentOS_6-dev/picard/v1.140/picard.jar CollectRnaSeqMetrics REF_FLAT=", REFFLAT_FILE, 
        " RIBOSOMAL_INTERVALS=", RIBOSOME_LIST," STRAND_SPECIFICITY=NONE REFERENCE_SEQUENCE=", GENOME,
        " INPUT=", bam, " OUTPUT=",sample_name, ".picRNAMetrics ",
        " CHART_OUTPUT=", sample_name,".picRNAMetrics.pdf VALIDATION_STRINGENCY=SILENT TMP_DIR=", TEMP
        ]
    cmd6_run=("".join(cmd6))
    print(cmd6_run)
    subprocess.run(cmd6_run,check=True,shell=True)




### STAR Fusion
@click.command()
@click.option('--fastq1', required=True, help='Fastq Read 1')
@click.option('--fastq2', required=True, help='Fastq Read 2')
@click.option('--sample_name', required=True, help='Sample Name')
@click.option('--species', required=True, help='Type human or mouse')
def star_fusion(fastq1, fastq2, sample_name, species):
    """Run STAR-Fusion."""

    OUTDIR="_STAR_FUSION"
    OUTDIR= sample_name + OUTDIR

    if species=="human":
        human_b37_e75()
    elif species=="mouse":
        mouse_b38_e75()
    else:
        print("Species should be either human or mouse")


    cmd7 = [
        "source /home/yellapav/.bashrc_fusion; ",
        "/ifs/e63data/armstronglab/noushin/local/STAR-Fusion_v0.5.4/STAR-Fusion ",
        "--genome_lib_dir ", STAR_FUSION_LIB, 
        " --right_fq ", fastq1, " --left_fq ", fastq2, " --output_dir ",OUTDIR 
        ]

    cmd7_run=("".join(cmd7))
    print(cmd7_run)
    subprocess.run(cmd7_run,check=True, shell=True)




@click.group()
def commands():
    """RNA-Seq pipeline has the following options."""
    pass


commands.add_command(add_rgid)
commands.add_command(star_align)
commands.add_command(alignment_metrics)
commands.add_command(samtools_flagstat)
commands.add_command(counts_htseq)
commands.add_command(salmon)
commands.add_command(collect_rna_metrics)
commands.add_command(star_fusion)

if __name__ == '__main__':
    commands()
