(bioinf) franza@ip-172-31-6-96:/$ cd /home/franza
#Per rimuovere una directory vuota:
(bioinf) franza@ip-172-31-6-96:~/my_project$ rmdir scripts (es)
# cancellare l'intera directory allign con tutto ciò che contiene:
 (bioinf) franza@ip-172-31-6-96:~/my_project$ rm -r nome directory
 
PS C:\Users\valen> cd C:\Users\valen\Desktop\tec_seq
PS C:\Users\valen\Desktop\tec_seq> ssh -i studentkey.pem franza@15.161.204.72
franza@ip-172-31-6-96:~$ source .bashrc
franza@ip-172-31-6-96:~$ conda
franza@ip-172-31-6-96:~$ conda activate bioinf
# creazione directory del progetto
(bioinf) franza@ip-172-31-6-96:~$ mkdir my_pro
# creazione directory per le sequenze scaricate 
(bioinf) franza@ip-172-31-6-96:~$ mkdir data
# scaricato le sequenze da SRA-Explorer 
(bioinf) franza@ip-172-31-6-96:~/data$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR150/082/SRR15093882/SRR15093882_1.fastq.gz 
                                            ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR150/082/SRR15093882/SRR15093882_2.fastq.gz
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ ls
SRR15093882_1.fastq.gz  SRR15093882_2.fastq.gz

(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ fastqc -t 2 -o qc *.fastq.gz
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ multiqc . -n SRR15093882_qc_report.html
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ ls

(bioinf) franza@ip-172-31-6-96:~/my_pro/data/qc$ ls
SRR15093882_1.fastq.gz  SRR15093882_2.fastq.gz  SRR15093882_qc_report.html  SRR15093882_qc_report_data  qc
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ pwd
/home/franza/my_pro/data


#shell pc per scaricare il file nella cartella tec_seq per visualizzare multi QC
PS C:\Users\valen\Desktop\tec_seq>  scp -i studentkey.pem franza@15.161.204.72:/home/franza/my_pro/data/SRR15093882_qc_report.html .
SRR15093882_qc_report.html                                                            100% 1111KB   2.8MB/s   00:00
PS C:\Users\valen\Desktop\tec_seq>
#analisi multi QC : in sezione Adapter Content ho visualizzato gli adattori per efettuare il Trimming 
#illumina small rna 3' adapter : TGGAATTCTCGGGTGCCAAGG (read1)
#illumina universal adapter : AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC (read2)

# Mapping - Allineamento whole genome

(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ mkdir trimmed
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ fastp --in1 SRR15093882_1.fastq.gz \
--in2  SRR15093882_2.fastq.gz  \
--out1 trimmed/SRR15093882_1.trimmed.fq.gz \
--out2 trimmed/SRR15093882_2.trimmed.fq.gz \
--adapter_sequence TGGAATTCTCGGGTGCCAAGG \
--adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

output: Read1 before filtering:
total reads: 3454249
total bases: 521591599
Q20 bases: 502846868(96.4062%)
Q30 bases: 469674098(90.0463%)

Read2 before filtering:
total reads: 3454249
total bases: 521591599
Q20 bases: 495312553(94.9618%)
Q30 bases: 465469497(89.2402%)

Read1 after filtering:
total reads: 3308572
total bases: 495865865
Q20 bases: 483022283(97.4099%)
Q30 bases: 452049625(91.1637%)

Read2 after filtering:
total reads: 3308572
total bases: 495960950
Q20 bases: 479572448(96.6956%)
Q30 bases: 452080853(91.1525%)

Filtering result:
reads passed filter: 6617144
reads failed due to low quality: 260112
reads failed due to too many N: 5852
reads failed due to too short: 25390
reads with adapter trimmed: 142735
bases trimmed due to adapters: 9522195

Duplication rate: 73.1825%

Insert size peak (evaluated by paired-end reads): 237

JSON report: fastp.json
HTML report: fastp.html

fastp --in1 SRR15093882_1.fastq.gz --in2 SRR15093882_2.fastq.gz --out1 trimmed/SRR15093882_1.trimmed.fq.gz --out2 trimmed/SRR15093882_2.trimmed.fq.gz --adapter_sequence TGGAATTCTCGGGTGCCAAGG --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
fastp v0.23.2, time used: 21 seconds

(bioinf) franza@ip-172-31-6-96:~/my_pro/data/trimmed$ ls
SRR15093882_1.trimmed.fq.gz  SRR15093882_2.trimmed.fq.gz
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/trimmed$ mkdir qc
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/trimmed$ fastqc -t 2 -o qc *.fq.gz
application/gzip
application/gzip
Started analysis of SRR15093882_1.trimmed.fq.gz
Started analysis of SRR15093882_2.trimmed.fq.gz
Approx 5% complete for SRR15093882_1.trimmed.fq.gz
Approx 5% complete for SRR15093882_2.trimmed.fq.gz
Approx 10% complete for SRR15093882_1.trimmed.fq.gz
Approx 10% complete for SRR15093882_2.trimmed.fq.gz
Approx 15% complete for SRR15093882_1.trimmed.fq.gz
Approx 15% complete for SRR15093882_2.trimmed.fq.gz
Approx 20% complete for SRR15093882_1.trimmed.fq.gz
Approx 20% complete for SRR15093882_2.trimmed.fq.gz
Approx 25% complete for SRR15093882_1.trimmed.fq.gz
Approx 25% complete for SRR15093882_2.trimmed.fq.gz
Approx 30% complete for SRR15093882_1.trimmed.fq.gz
Approx 30% complete for SRR15093882_2.trimmed.fq.gz
Approx 35% complete for SRR15093882_1.trimmed.fq.gz
Approx 35% complete for SRR15093882_2.trimmed.fq.gz
Approx 40% complete for SRR15093882_1.trimmed.fq.gz
Approx 40% complete for SRR15093882_2.trimmed.fq.gz
Approx 45% complete for SRR15093882_1.trimmed.fq.gz
..................................................
Approx 75% complete for SRR15093882_1.trimmed.fq.gz
Approx 75% complete for SRR15093882_2.trimmed.fq.gz
Approx 80% complete for SRR15093882_1.trimmed.fq.gz
Approx 80% complete for SRR15093882_2.trimmed.fq.gz
Approx 85% complete for SRR15093882_1.trimmed.fq.gz
Approx 85% complete for SRR15093882_2.trimmed.fq.gz
Approx 90% complete for SRR15093882_1.trimmed.fq.gz
Approx 90% complete for SRR15093882_2.trimmed.fq.gz
Approx 95% complete for SRR15093882_1.trimmed.fq.gz
Approx 95% complete for SRR15093882_2.trimmed.fq.gz
Analysis complete for SRR15093882_1.trimmed.fq.gz
Analysis complete for SRR15093882_2.trimmed.fq.gz
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/trimmed$ cd qc/
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/trimmed/qc$ ls
SRR15093882_1.trimmed_fastqc.html  SRR15093882_2.trimmed_fastqc.html
SRR15093882_1.trimmed_fastqc.zip   SRR15093882_2.trimmed_fastqc.zip
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/trimmed/qc$ multiqc . -n SRR15093882_trimmed_qc_report.html
/home/franza/miniconda3/envs/bioinf/lib/python2.7/site-packages/multiqc-1.0.dev0-py2.7.egg/multiqc/utils/config.py:44: YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
  configs = yaml.load(f)
/home/franza/miniconda3/envs/bioinf/lib/python2.7/site-packages/multiqc-1.0.dev0-py2.7.egg/multiqc/utils/config.py:50: YAMLLoadWarning: calling yaml.load() without Loader=... is deprecated, as the default Loader is unsafe. Please read https://msg.pyyaml.org/load for full details.
  sp = yaml.load(f)
[INFO   ]         multiqc : This is MultiQC v1.0.dev0
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching '.'
[INFO   ]          fastqc : Found 2 reports
[INFO   ]         multiqc : Report      : SRR15093882_trimmed_qc_report.html
[INFO   ]         multiqc : Data        : SRR15093882_trimmed_qc_report_data
[INFO   ]         multiqc : MultiQC complete
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/trimmed/qc$ ls
SRR15093882_1.trimmed_fastqc.html  SRR15093882_2.trimmed_fastqc.html  SRR15093882_trimmed_qc_report.html
SRR15093882_1.trimmed_fastqc.zip   SRR15093882_2.trimmed_fastqc.zip   SRR15093882_trimmed_qc_report_data

#cercare percorso
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/trimmed/qc$ pwd
/home/franza/my_pro/data/trimmed/qc

#andare in shell per scaricare il file 

PS C:\Users\valen\Desktop\tec_seq>  scp -i studentkey.pem franza@15.161.204.72:/home/franza/my_pro/data/trimmed/qc/SRR15093882_trimmed_qc_report.html .
SRR15093882_trimmed_qc_report.html                                                    100% 1109KB   2.9MB/s   00:00
PS C:\Users\valen\Desktop\tec_seq>

# Allineamento sequenze
(bioinf) franza@ip-172-31-6-96:~/my_pro/data mkdir map_reads
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ bwa mem \
  -t 4 \
  -M \
  -o map_reads/SRR15093882.sam \
  /resources/GRCh38.primary_assembly.genome.fa \
  trimmed/SRR15093882_1.trimmed.fq.gz \
  trimmed/SRR15093882_2.trimmed.fq.gz
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 266890 sequences (40000006 bp)...
[M::process] read 266936 sequences (40000260 bp)...
...................................................
[M::mem_process_seqs] Processed 212374 reads in 35.275 CPU sec, 8.796 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 4 -M -o map_reads/SRR15093882.sam /resources/GRCh38.primary_assembly.genome.fa trimmed/SRR15093882_1.trimmed.fq.gz trimmed/SRR15093882_2.trimmed.fq.gz
[main] Real time: 251.235 sec; CPU: 1000.243 sec
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ cd map_reads
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads$ ls
SRR15093882.sam
# Voglio contare quante letture hanno allineamenti secondari --> out: 12.509 letture che hanno almeno un allineamento secondario
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads$ grep "SA" SRR15093882.sam | cut -f1 | sort | uniq | wc -l
12509

#convertito file, ordinato e indicizzato
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads$ samtools sort -T SRR15093882.temp -o SRR15093882.sorted.bam SRR15093882.sam && samtools index SRR15093882.sorted.bam
[bam_sort_core] merging from 3 files and 1 in-memory blocks...
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads$ ls
SRR15093882.sam  SRR15093882.sorted.bam  SRR15093882.sorted.bam.bai

#Estrazione statistiche di allineamento
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads$ samtools stats SRR15093882.sorted.bam > qc/SRR15093882.samtools.stats
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads$ samtools flagstat SRR15093882.sorted.bam > qc/SRR15093882.samtools.flagstat
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads$ samtools depth SRR15093882.sorted.bam > qc/SRR15093882.samtools.depth
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads$ cd qc/
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads/qc$ ls
SRR15093882.samtools.depth  SRR15093882.samtools.flagstat  SRR15093882.samtools.stats

#visualizzare i file
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads/qc$ less SRR15093882.samtools.flagstat
6631762 + 0 in total (QC-passed reads + QC-failed reads)
6617144 + 0 primary
14618 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
6581418 + 0 mapped (99.24% : N/A)
6566800 + 0 primary mapped (99.24% : N/A)
6617144 + 0 paired in sequencing
3308572 + 0 read1
3308572 + 0 read2
6410972 + 0 properly paired (96.88% : N/A)
6566088 + 0 with itself and mate mapped
712 + 0 singletons (0.01% : N/A)
38240 + 0 with mate mapped to a different chr
34337 + 0 with mate mapped to a different chr (mapQ>=5)
SRR15093882.samtools.flagstat (END)

(bioinf) franza@ip-172-31-6-96:~/my_pro/data/map_reads/qc$ less SRR15093882.samtools.stats
# This file was produced by samtools stats (1.13+htslib-1.13) and can be plotted using plot-bamstats
# This file contains statistics for all reads.
# The command line was:  stats SRR15093882.sorted.bam
# CHK, Checksum [2]Read Names   [3]Sequences    [4]Qualities
# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)
CHK     f8ad76d3        589a2546        3d12a0e4
# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN      raw total sequences:    6617144 # excluding supplementary and secondary reads
SN      filtered sequences:     0
SN      sequences:      6617144
SN      is sorted:      1
SN      1st fragments:  3308572
SN      last fragments: 3308572
SN      reads mapped:   6566800
SN      reads mapped and paired:        6566088 # paired-end technology bit set + both mates mapped
SN      reads unmapped: 50344
SN      reads properly paired:  6410972 # proper-pair bit set
SN      reads paired:   6617144 # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      142626  # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 14618
SN      supplementary alignments:       0
SN      total length:   991826815       # ignores clipping
SN      total first fragment length:    495865865       # ignores clipping
SN      total last fragment length:     495960950       # ignores clipping
SN      bases mapped:   988722624       # ignores clipping
SN      bases mapped (cigar):   947801143       # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     2546165 # from NM fields
SN      error rate:     2.686392e-03    # mismatches / bases mapped (cigar)
SN      average length: 149
SN      average first fragment length:  150
SN      average last fragment length:   150
SN      maximum length: 151
SRR15093882.samtools.stats
.......................

#Allineamento alternativo: singolo cromosoma (chr 5)
# Trovare link chr5 --> https://genome.ucsc.edu/cgi-bin/hgGateway--> Downloads --> Dec. 2013 (GRCh38/hg38)
--> Sequence data by chromosome --> chr5.fa.gz link (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr5.fa.gz)

(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr5.fa.gz
--2025-07-03 09:02:58--  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr5.fa.gz
Resolving hgdownload.soe.ucsc.edu (hgdownload.soe.ucsc.edu)... 128.114.119.163
Connecting to hgdownload.soe.ucsc.edu (hgdownload.soe.ucsc.edu)|128.114.119.163|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 58835176 (56M) [application/x-gzip]
Saving to: ‘chr5.fa.gz’

chr5.fa.gz                    100%[=================================================>]  56.11M  16.9MB/s    in 3.6s

2025-07-03 09:03:02 (15.4 MB/s) - ‘chr5.fa.gz’ saved [58835176/58835176]
# Decomprimo e indicizzo 
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ gunzip chr5.fa.gz
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ bwa index chr5.fa
[bwa_index] Pack FASTA... 1.36 sec
[bwa_index] Construct BWT for the packed sequence...
.................................................
[bwt_gen] Finished constructing BWT in 97 iterations.
[bwa_index] 94.09 seconds elapse.
[bwa_index] Update BWT... 0.97 sec
[bwa_index] Pack forward-only FASTA... 0.82 sec
[bwa_index] Construct SA from BWT and Occ... 49.91 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index chr5.fa
[main] Real time: 147.769 sec; CPU: 147.154 sec
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ ls
SRR15093882_1.fastq.gz      SRR15093882_qc_report_data  chr5.fa      chr5.fa.bwt  fastp.html  qc
SRR15093882_2.fastq.gz      all_chr5                    chr5.fa.amb  chr5.fa.pac  fastp.json  trimmed
SRR15093882_qc_report.html  chr1                        chr5.fa.ann  chr5.fa.sa   map_reads
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ bwa mem \
  -t 4 \
  -M \
  -o all_chr5/SRR15093882_chr5.sam \
  chr5.fa \
  trimmed/SRR15093882_1.trimmed.fq.gz \
  trimmed/SRR15093882_2.trimmed.fq.gz
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 266890 sequences (40000006 bp)...
[M::process] read 266936 sequences (40000260 bp)...
................................................................................
[M::mem_process_seqs] Processed 212374 reads in 29.883 CPU sec, 7.440 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 4 -M -o all_chr5/SRR15093882_chr5.sam chr5.fa trimmed/SRR15093882_1.trimmed.fq.gz trimmed/SRR15093882_2.trimmed.fq.gz
[main] Real time: 222.926 sec; CPU: 900.891 sec

#Convertire il file .sam in .bam, ordinarlo e indicizzarlo 
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ ls
SRR15093882_chr5.sam  qc_chr5
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools view -b -o SRR15093882_chr5.bam SRR15093882_chr5.sam
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools sort -T SRR15093882_chr5.temp -o SRR15093882_chr5.sorted.bam SRR15093882_chr5.bam && samtools index SRR15093882_chr5.sorted.bam
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
#Elaborazione statistica 
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ ls
SRR15093882_chr5.bam  SRR15093882_chr5.sam  SRR15093882_chr5.sorted.bam  SRR15093882_chr5.sorted.bam.bai  qc_chr5
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools stats SRR15093882_chr5.sorted.bam > qc_chr5/SRR15093882.samtools.stats
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools flagstat SRR15093882_chr5.sorted.bam > qc_chr5/SRR15093882.samtools.flagstat
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools depth SRR15093882_chr5.sorted.bam > qc_chr5/SRR15093882.samtools.depth
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ cd qc_chr5/
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5/qc_chr5$ ls
SRR15093882.samtools.depth  SRR15093882.samtools.flagstat  SRR15093882.samtools.stats

 #Visualizzazione analisi 
 (bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5/qc_chr5$ less SRR15093882.samtools.flagstat
6619814 + 0 in total (QC-passed reads + QC-failed reads)
6617144 + 0 primary
2670 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
297044 + 0 mapped (4.49% : N/A)
294374 + 0 primary mapped (4.45% : N/A)
6617144 + 0 paired in sequencing
3308572 + 0 read1
3308572 + 0 read2
275790 + 0 properly paired (4.17% : N/A)
285552 + 0 with itself and mate mapped
8822 + 0 singletons (0.13% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
SRR15093882.samtools.flagstat (END)
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5/qc_chr5$ less SRR15093882.samtools.stats
# This file was produced by samtools stats (1.13+htslib-1.13) and can be plotted using plot-bamstats
# This file contains statistics for all reads.
# The command line was:  stats SRR15093882_chr5.sorted.bam
# CHK, Checksum [2]Read Names   [3]Sequences    [4]Qualities
# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)
CHK     fe4150a9        7bc7492a        87462d3a
# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN      raw total sequences:    6617144 # excluding supplementary and secondary reads
SN      filtered sequences:     0
SN      sequences:      6617144
SN      is sorted:      1
SN      1st fragments:  3308572
SN      last fragments: 3308572
SN      reads mapped:   294374
SN      reads mapped and paired:        285552  # paired-end technology bit set + both mates mapped
SN      reads unmapped: 6322770
SN      reads properly paired:  275790  # proper-pair bit set
SN      reads paired:   6617144 # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      8261    # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 2670
SN      supplementary alignments:       0
SN      total length:   991826815       # ignores clipping
SN      total first fragment length:    495865865       # ignores clipping
SN      total last fragment length:     495960950       # ignores clipping
SN      bases mapped:   44108077        # ignores clipping
SN      bases mapped (cigar):   29498993        # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     2462507 # from NM fields
SN      error rate:     8.347767e-02    # mismatches / bases mapped (cigar)
SN      average length: 149
SN      average first fragment length:  150
SN      average last fragment length:   150
SN      maximum length: 151
SRR15093882.samtools.stats
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5/qc_chr5$ head SRR15093882.samtools.depth
chr5    190388  2
chr5    190389  2
chr5    190390  2
chr5    190391  2
chr5    190392  2
chr5    190393  2
chr5    190394  2
chr5    190395  2
chr5    190396  2
chr5    190397  2

#Stimare per ciascuno degli autosomi la percentuale di basi coperte con almeno 20 reads
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5/qc_chr5$ awk --help
Usage: awk [POSIX or GNU style options] -f progfile [--] file ...
Usage: awk [POSIX or GNU style options] [--] 'program' file ...
POSIX options:          GNU long options: (standard)
        -f progfile             --file=progfile
        -F fs                   --field-separator=fs
        -v var=val              --assign=var=val
Short options:          GNU long options: (extensions)
        -b                      --characters-as-bytes
        -c                      --traditional
        -C                      --copyright
        -d[file]                --dump-variables[=file]
        -D[file]                --debug[=file]
        -e 'program-text'       --source='program-text'
        -E file                 --exec=file
        -g                      --gen-pot
        -h                      --help
        -i includefile          --include=includefile
        -I                      --trace
        -l library              --load=library
        -L[fatal|invalid|no-ext]        --lint[=fatal|invalid|no-ext]
        -M                      --bignum
        -N                      --use-lc-numeric
        -n                      --non-decimal-data
        -o[file]                --pretty-print[=file]
        -O                      --optimize
        -p[file]                --profile[=file]
        -P                      --posix
        -r                      --re-interval
        -s                      --no-optimize
        -S                      --sandbox
        -t                      --lint-old
        -V                      --version

To report bugs, use the `gawkbug' program.
For full instructions, see the node `Bugs' in `gawk.info'
which is section `Reporting Problems and Bugs' in the
printed version.  This same information may be found at
https://www.gnu.org/software/gawk/manual/html_node/Bugs.html.
PLEASE do NOT try to report bugs by posting in comp.lang.awk,
or by using a web forum such as Stack Overflow.

gawk is a pattern scanning and processing language.
By default it reads standard input and writes standard output.


#Stimare per ciascuno degli autosomi la percentuale di basi coperte con almeno 20 reads
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/ mkdir coverage_analysis
(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ for chr in {1..22}; do
    samtools depth -r chr${chr} map_reads/SRR15093882.sorted.bam > coverage_analysis/chr${chr}_coverage.txt
    TOTAL_BASES=$(wc -l < coverage_analysis/chr${chr}_coverage.txt)
    COVERED_20X=$(awk '$3>=20' coverage_analysis/chr${chr}_coverage.txt | wc -l)
    if [ $TOTAL_BASES -gt 0 ]; then
        PERCENTAGE=$(echo "scale=2; $COVERED_20X * 100 / $TOTAL_BASES" | bc -l)
        echo "chr${chr}: ${PERCENTAGE}% coverage ≥20x"
        fi
done
chr1: 69.21% coverage ≥20x
chr2: 68.04% coverage ≥20x
chr3: 76.18% coverage ≥20x
chr4: 75.92% coverage ≥20x
chr5: 56.03% coverage ≥20x
chr6: 58.80% coverage ≥20x
chr7: 70.23% coverage ≥20x
chr8: 68.77% coverage ≥20x
chr9: 68.91% coverage ≥20x
chr10: 50.97% coverage ≥20x
chr11: 68.37% coverage ≥20x
chr12: 85.06% coverage ≥20x
chr13: 60.44% coverage ≥20x
chr14: 68.71% coverage ≥20x
chr15: 61.71% coverage ≥20x
chr16: 51.08% coverage ≥20x
chr17: 77.35% coverage ≥20x
chr18: 74.90% coverage ≥20x
chr19: 75.78% coverage ≥20x
chr20: 90.87% coverage ≥20x
chr21: 57.69% coverage ≥20x
chr22: 57.25% coverage ≥20x

(bioinf) franza@ip-172-31-6-96:~/my_pro/data$ for chr in {1..22}; do
    samtools depth -a -r chr${chr} map_reads/SRR15093882.sorted.bam > coverage_analysis/chr${chr}_coverage.txt
    TOTAL_BASES=$(wc -l < coverage_analysis/chr${chr}_coverage.txt)
    COVERED_20X=$(awk '$3>=20' coverage_analysis/chr${chr}_coverage.txt | wc -l)
    if [ $TOTAL_BASES -gt 0 ]; then
        PERCENTAGE=$(echo "scale=2; $COVERED_20X * 100 / $TOTAL_BASES" | bc -l)
        echo "chr${chr}: ${PERCENTAGE}% coverage ≥20x"
                fi
done
chr1: 0% coverage ≥20x
chr2: 0% coverage ≥20x
chr3: 0% coverage ≥20x
chr4: 0% coverage ≥20x
chr5: 0% coverage ≥20x
chr6: 0% coverage ≥20x
chr7: 0% coverage ≥20x
chr8: 0% coverage ≥20x
chr9: 0% coverage ≥20x
chr10: 0% coverage ≥20x
chr11: .01% coverage ≥20x
chr12: 0% coverage ≥20x
chr13: 0% coverage ≥20x
chr14: 0% coverage ≥20x
chr15: .01% coverage ≥20x
chr16: 0% coverage ≥20x
chr17: .03% coverage ≥20x
chr18: 0% coverage ≥20x
chr19: .01% coverage ≥20x
chr20: .01% coverage ≥20x
chr21: .01% coverage ≥20x
chr22: 0% coverage ≥20x


#Selezionare 4 geni appartenenti al cromosoma scelto, tramite la lista in https://www.oncokb.org/cancer-genes, ed estrarre le rispettive coperture. 
#https://www.oncokb.org/cancer-genes --> search gene --> summary --> Ensembl Gene	ENSG..(GRCh37/GRCh38) --> GRCh38 --> Location (copiare e incolla la posizione)

(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ mkdir genes
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools depth \
  -r chr5:35852695-35879603  \
  SRR15093882_chr5.sorted.bam \
> genes/SRR15093882.samtools.depth.IL7R
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools coverage -r chr5:35852695-35879603 SRR15093882_chr5.sorted.bam
#rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq
chr5    35852695        35879603        0       0       0       0       0       0
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools depth \
  -r chr5:171387116-171411810  \
  SRR15093882_chr5.sorted.bam \
> genes/SRR15093882.samtools.depth.NPM1
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools coverage -r chr5:171387116-171411810 SRR15093882_chr5.sorted.bam
#rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq
chr5    171387116       171411810       42      259     1.0488  0.250091        31.4    60
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools depth \
  -r chr5:150053291-150113372  \
  SRR15093882_chr5.sorted.bam \
> genes/SRR15093882.samtools.depth.CSF1R
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools coverage -r chr5:150053291-150113372 SRR15093882_chr5.sorted.bam
#rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq
chr5    150053291       150113372       3       169     0.281282        0.00334543      35.3    20
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools depth \
  -r chr5:138465479-138469303  \
  SRR15093882_chr5.sorted.bam \
> genes/SRR15093882.samtools.depth.EGR1
(bioinf) franza@ip-172-31-6-96:~/my_pro/data/all_chr5$ samtools coverage -r chr5:138465479-138469303 SRR15093882_chr5.sorted.bam
#rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq
chr5    138465479       138469303       19805   265     6.9281  629.644 33.8    59.9


#Caricare i file di intersse nel pc 
PS C:\Users\valen\Desktop\tec_seq> scp -i studentkey.pem franza@15.161.204.72:/home/franza/my_pro/data/map_reads/qc/SRR15093882.samtools.flagstat .
SRR15093882.samtools.flagstat                                                                                             100%  522     6.1KB/s   00:00
PS C:\Users\valen\Desktop\tec_seq> scp -i studentkey.pem franza@15.161.204.72:/home/franza/my_pro/data/all_chr5/qc_chr5/SRR15093882.samtools.flagstat flagstat_chr5.txt
SRR15093882.samtools.flagstat                                                                                             100%  507     5.7KB/s   00:00
PS C:\Users\valen\Desktop\tec_s