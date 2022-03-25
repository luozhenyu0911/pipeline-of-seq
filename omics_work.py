#coding=utf-8
from __future__ import print_function 
from __future__ import division
import sys
import subprocess
import argparse
example_text = '''example:
    python this.py -w “hisat2, bowtie2 or bwa” -i genome_index -f1 fastq1 -f2 fastq2 -g gff -t thread -p prefix
'''
parser = argparse.ArgumentParser(description="the script for omics work, Note: pair data",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--whichOfomics','-w',type=str,help="Which omics process(Choose one of the three: “hisat2, bowtie2, bwa”)",required= True,metavar='')
parser.add_argument('--genome_index','-i',type=str,help="thread",required= True,metavar='')
parser.add_argument('--fastq1', '-f1',type= str,help="fastq1",required= True,metavar='')
parser.add_argument('--fastq2', '-f2',type= str,help="fastq2",required= True,metavar='')
parser.add_argument('--gff', '-g',type= str,help="gff",required= True,metavar='')
parser.add_argument('--thread', '-t',type= str,help="thread",required= True,metavar='')
parser.add_argument('--prefix', '-p',type= str,help="the prefix of output file",required= True,metavar='')
args = parser.parse_args()

#记得习惯性空一格
whichOfomics = args.whichOfomics
fq1 = args.fastq1
fq2 = args.fastq2
genome_index = args.genome_index
gff = args.gff
thread = args.thread
prefix = args.prefix

if str(whichOfomics) == "hisat2":
    subprocess.call(["hisat2","-p",thread,"--dta","-x",genome_index,"-1",fq1,"-2",fq2,"-S",prefix+".sam"])
    subprocess.call(["samtools","sort","-@",thread,"-O","BAM","-o",prefix+".raw.bam",prefix+".sam"])
    subprocess.call(["rm",prefix+".sam"])
    subprocess.call(["samtools","index","-@",thread,prefix+".raw.bam"])
    subprocess.call(["samtools","view","-@",thread,"-q","20","-O","BAM","-o",prefix+".uniq.bam",prefix+".raw.bam"])
    subprocess.call(["samtools","index","-@",thread,prefix+".uniq.bam"])

    data_info=open('data_info.txt',"a")
    data_info.write("RNA-seq")
    data_info.write(prefix)
    subprocess.call(["samtools","view","-@",thread,"-c",prefix+".raw.bam"],stdout=data_info)
    subprocess.call(["samtools","view","-@",thread,"-F","4","-c",prefix+".raw.bam"],stdout=data_info)
    subprocess.call(["samtools","view","-@",thread,"-c",prefix+".uniq.bam"],stdout=data_info)
    data_info.close()
    subprocess.call(["stringtie","-e","-p",thread,"-G",gff,"-o",prefix+".gtf","-l",prefix,"-A",prefix+"_genefpkm.txt",prefix+".uniq.bam"])
    subprocess.call(["bamCoverage","-p",thread,"-b",prefix+".uniq.bam","-o",prefix+".uniq.bw","--centerReads","-bs","20","--smoothLength","60","--normalizeUsing","RPKM"])

elif str(whichOfomics) == "bowtie2":
    subprocess.call(["bowtie2","-p",thread,"-N","1","-x",genome_index,"-1",fq1,"-2",fq2,"-S",prefix+"raw.sam"])
    subprocess.call(["samtools","sort","-@",thread,"-O","BAM","-o",prefix+"raw.bam",prefix+"raw.sam"])
    subprocess.call(["samtools","view","-@",thread,"-q","20","-O","BAM","-o",prefix+"tmp.bam",prefix+"raw.sam"])
    subprocess.call(["samtools","sort","-O","BAM","-@",thread,"-n","-o",prefix+"sort_tmp.bam",prefix+"tmp.bam"])
    # 降重的不同方法
    # subprocess.call(["java","-jar","/gss1/biosoft/picard.jar","MarkDuplicates","I="+sys.argv[3]+"sort_tmp.bam","O="+sys.argv[3]+"rmdup.bam","M=marked_dup_metrics.txt","REMOVE_DUPLICATES=true"])

    subprocess.call(["samtools","fixmate","-@",thread,"-m",prefix+"sort_tmp.bam",prefix+"fixmate.bam"])
    subprocess.call(["samtools","sort","-@",thread,"-O","BAM","-o",prefix+"fixmatesort.bam",prefix+"fixmate.bam"])
    subprocess.call(["samtools","markdup","-@",thread,"-r","-O","BAM",prefix+"fixmatesort.bam",prefix+"rmdup.bam"])
    subprocess.call(["samtools","index","-@",thread,prefix+"rmdup.bam"])
    
    data_info=open('data_info.txt',"a")
    data_info.write("bowtie2 for seq")
    data_info.write(prefix)
    subprocess.call(["samtools","view","-@",thread,"-c",prefix+"raw.sam"])
    subprocess.call(["samtools","view","-@",thread,"-F","4","-c",prefix+"raw.sam"])
    subprocess.call(["samtools","view","-@",thread,"-c",prefix+"rmdup.bam"])
    data_info.close()
    subprocess.call(["bamCoverage","-p",thread,"-b",prefix+"rmdup.bam","-o",prefix+"rmdup.bw","--centerReads","--binSize","20","--smoothLength","60","--normalizeUsing","RPKM"])
    subprocess.call(["rm",prefix+"raw.sam"])
    subprocess.call(["rm",prefix+"tmp.bam"])
    subprocess.call(["rm",prefix+"sort_tmp.bam"])
    subprocess.call(["rm",prefix+"fixmate.bam"])
    subprocess.call(["rm",prefix+"fixmatesort.bam"])

elif str(whichOfomics) == "bwa":
    subprocess.call(["bwa","mem","-t","20","-M",genome_index,fq1,fq2,"-o",prefix+"raw.sam"])
    subprocess.call(["samtools","sort","-@",thread,"-O","BAM","-o",prefix+"raw.bam",prefix+"raw.sam"])
    subprocess.call(["samtools","view","-@",thread,"-q","20","-O","BAM","-o",prefix+"tmp.bam",prefix+"raw.sam"])
    subprocess.call(["samtools","sort","-O","BAM","-@",thread,"-n","-o",prefix+"sort_tmp.bam",prefix+"tmp.bam"])
    # 降重的不同方法
    # subprocess.call(["java","-jar","/gss1/biosoft/picard.jar","MarkDuplicates","I="+sys.argv[3]+"sort_tmp.bam","O="+sys.argv[3]+"rmdup.bam","M=marked_dup_metrics.txt","REMOVE_DUPLICATES=true"])

    subprocess.call(["samtools","fixmate","-@",thread,"-m",prefix+"sort_tmp.bam",prefix+"fixmate.bam"])
    subprocess.call(["samtools","sort","-@",thread,"-O","BAM","-o",prefix+"fixmatesort.bam",prefix+"fixmate.bam"])
    subprocess.call(["samtools","markdup","-@",thread,"-r","-O","BAM",prefix+"fixmatesort.bam",prefix+"rmdup.bam"])
    subprocess.call(["samtools","index","-@",thread,prefix+"rmdup.bam"])
    
    data_info=open('data_info.txt',"a")
    data_info.write("bowtie2 for seq")
    data_info.write(prefix)
    subprocess.call(["samtools","view","-@",thread,"-c",prefix+"raw.sam"])
    subprocess.call(["samtools","view","-@",thread,"-F","4","-c",prefix+"raw.sam"])
    subprocess.call(["samtools","view","-@",thread,"-c",prefix+"rmdup.bam"])
    data_info.close()
    subprocess.call(["bamCoverage","-p",thread,"-b",prefix+"rmdup.bam","-o",prefix+"rmdup.bw","--centerReads","--binSize","20","--smoothLength","60","--normalizeUsing","RPKM"])
    subprocess.call(["rm",prefix+"raw.sam"])
    subprocess.call(["rm",prefix+"tmp.bam"])
    subprocess.call(["rm",prefix+"sort_tmp.bam"])
    subprocess.call(["rm",prefix+"fixmate.bam"])
    subprocess.call(["rm",prefix+"fixmatesort.bam"])
