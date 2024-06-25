##Genome indexing:
samtools faidx merge.fa
cut -f1-2 merge.fa.fai >merge.fa.size

##generate the fragments.csv
cooler digest -o DpnII_digetstion_fragment.bed merge.fa.size merge.fa DpnII
echo "chrom,start,end,fragment_length,fragment_id" > DpnII.vd.fragments.csv  #fragments.csv with columns [chrom, start, end, fragment_length, fragment_id]
awk -v OFS="," '{pirnt $1,$2,$3,$3-$2,NR}' DpnII_digetstion_fragment.bed >> DpnII.vd.fragments.csv

for i in 20231109-NPL2300560-P6-PAO58602-sup.pass10.fastq.gz 20231109-NPL2300560-P6-PAO58602-sup.pass12.fastq.gz 20231109-NPL2300560-P6-PAO58602-sup.pass3.fastq.gz 20231109-NPL2300560-P6-PAO58602-sup.pass5.fastq.gz 20231109-NPL2300560-P6-PAO58602-sup.pass6.fastq.gz 20231109-NPL2300560-P6-PAO58602-sup.pass7.fastq.gz 20231109-NPL2300560-P6-PAO58602-sup.pass8.fastq.gz 20231109-NPL2300560-P6-PAO58602-sup.pass9.fastq.gz 20231115-NPL2300560-P4-PAQ52955-sup.pass11.fastq.gz 20231115-NPL2300560-P4-PAQ52955-sup.pass13.fastq.gz 20231115-NPL2300560-P4-PAQ52955-sup.pass1.fastq.gz 20231115-NPL2300560-P4-PAQ52955-sup.pass2.fastq.gz 20231115-NPL2300560-P4-PAQ52955-sup.pass4.fastq.gz 20231115-NPL2300560-P4-PAQ52955-sup.pass.fastq.gz
do
set -evx
##mapping
minimap2 -cx map-ont -t 20 merge.fa $i > $i.paf
##Filtering
awk -v OFS="\t" '($12>0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$23}' $i.paf > $i.paf.filter
##To run the pipeline type (https://github.com/zhengdafangyuan/HiPore-C)
time sh ./porec_pipeline/Fragment_Annotation_linux.sh ./porec_pipeline/Read_Fragment_Annotation.py 02_paf2bam merge.fa $i $i.paf.filter DpnII.vd.fragments.csv ./porec_pipeline/minimap2_subreads_remapping.sh ./porec_pipeline/multichrom_check.sh
##Generate pairwise contact matrix (https://github.com/zhengdafangyuan/HiPore-C)
python ./porec_pipeline/Generate_Pairwise_contact_juicematrix.py -p ./02_paf2bam/tmpdir/$i/Read_Align_Fragment_RvdF.csv -o ./02_paf2bam/tmpdir/$i -s 0 -t 20 -c 1000000

##Converts the contact matrix file to BAM file format, as LACHESIS requires input files in BAM format.

##split bam and generate bin.bam
python BamModifier.py -i ./02_paf2bam/tmpdir/$i/contact_matrix.txt.bam -o ./02_paf2bam/tmpdir/$i/contact_matrix.txt.bam.bin.bam -b merge.fa.fai_bin.bed

done

Lachesis SMT.ini

##SMT.Bins.fasta (FASTA file split into bins of 200 kb)
##############
Sequence ID Example:
>h1tg000001l_1
>h1tg000001l_2
>h1tg000001l_3
>h1tg000001l_4
>h1tg000001l_5
>h1tg000001l_6
>h1tg000001l_7
###############

cat SMT.ini 
SPECIES = SMT
OUTPUT_DIR = result
DRAFT_ASSEMBLY_FASTA = ./SMT.Bins.fasta
SAM_DIR = ./bin_bam
SAM_FILES = sup.pass0.bin.bam   sup.pass12.bin.bam  sup.pass2.bin.bam  sup.pass5.bin.bam  sup.pass8.bin.bam sup.pass10.bin.bam  sup.pass13.bin.bam  sup.pass3.bin.bam  sup.pass6.bin.bam  sup.pass9.bin.bam sup.pass11.bin.bam  sup.pass1.bin.bam   sup.pass4.bin.bam  sup.pass7.bin.bam
RE_SITE_SEQ = GATC
USE_REFERENCE = 0
SIM_BIN_SIZE = 0
REF_ASSEMBLY_FASTA = 0
BLAST_FILE_HEAD = /xxx/xxx/xxx
DO_CLUSTERING = 1
DO_ORDERING = 1
DO_REPORTING = 1
OVERWRITE_GLM = 0
OVERWRITE_CLMS = 0
CLUSTER_N = 7
CLUSTER_CONTIGS_WITH_CENS = -1
CLUSTER_MIN_RE_SITES = 100
CLUSTER_MAX_LINK_DENSITY = 2.5
CLUSTER_NONINFORMATIVE_RATIO = 1.4
CLUSTER_DRAW_HEATMAP = 0 
CLUSTER_DRAW_DOTPLOT = 0 
ORDER_MIN_N_RES_IN_TRUNK = 60 
ORDER_MIN_N_RES_IN_SHREDS = 60 
ORDER_DRAW_DOTPLOTS = 0 
REPORT_EXCLUDED_GROUPS = 1 
REPORT_QUALITY_FILTER = 1 
REPORT_DRAW_HEATMAP = 0
