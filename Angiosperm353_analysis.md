# Angiosperm353_Paper

### This is the code that I used to analyze my Angiosperm-353 data in _Tricalysia_.


1. Trimmomatic - trim raw reads before running HybPiper

```
#loop for running Trimmomatic on all fastq files - be in directory with fastq files

for f1 in *R1_001.fastq.gz
do
f2=${f1%%R1_001.fastq.gz}"R2_001.fastq.gz"
java -jar /opt/modules/biology/trimmomatic/0.39/bin/trimmomatic-0.39.jar PE -phred33 $f1 $f2 $f1"_output"_R1_P.fq.gz $f1"_output"_R1_S.fq.gz $f1"_output"_R2_P.fq.gz $f1"_output"_R2_S.fq.gz ILLUMINACLIP:/opt/modules/biology/trimmomatic/0.39/bin/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36;
done

# moved cleaned reads to new folder - tricalysia 
# moved raw reads into rawreads folder - 20201223_miseqrun

mv *fastq.gz  ~/scratch/hybresults2/20201223_miseqrun

# uncompress zipped cleaned reads (make sure its just the cleaned reads in dir)
gunzip *.gz 

# put single reads into one file
while read name;
do cat "$name"_R1_S.fq "$name"_R2_S.fq > "$name"_12_S.fq;
done < namelist_tri2.txt;

```

2. Create Target file for HybPiper: I used NewTargets (https://github.com/chrisjackson-pellicle/NewTargets) and manually added _Tricalysia_ representatives for each gene from previous HybPiper run (with regular target file).

```
#NewTargets
#activate NewTarget environment
conda activate NewTargets

#change dir to NewTargets
#make select_file

#Rubiaceae
python filter_mega353.py mega353.fasta select_file_rubiaceae.txt -filtered_target_file rubiaceae_targetfile.fasta -report_filename rubiaceae_report.csv

conda deactivate
```

3. Run HybPiper with modified Target file

```
module load hybpiper
module load bwa

#reads_first.py
while read name;
do /opt/modules/biology/hybpiper/1.3.1/bin/reads_first.py -b rubiaceae_tricalysia_targetfile.fasta -r $name*_R*_P.fq --prefix $name --bwa --unpaired $name*_12_S.fq
done < namelist_tri2.txt

#get_seq_lengths.py
python /opt/modules/biology/hybpiper/1.3.1/bin/get_seq_lengths.py rubiaceae_tricalysia_targetfile.fasta namelist_tri2.txt dna > Seq_Length_Hybpiper_tri_mega353plus.txt

# hybpiper_stats.py
python /opt/modules/biology/hybpiper/1.3.1/bin/hybpiper_stats.py Seq_Length_Hybpiper_tri_mega353plus.txt namelist_tri2.txt > test_stats_Hybpiper_tri_mega353plus.txt

#intronerate.py
while read name;
do python /opt/modules/biology/hybpiper/1.3.1/bin/intronerate.py --prefix $name
done < namelist_tri2.txt

#paralog_investigator.py
while read i
do
echo $i
python /opt/modules/biology/hybpiper/1.3.1/bin/paralog_investigator.py $i
done < namelist_tri2.txt

#retrieve - SUPERCONTIGS
python /opt/modules/biology/hybpiper/1.3.1/bin/retrieve_sequences.py rubiaceae_tricalysia_targetfile.fasta . supercontig > supercontigs_tri_megaplus.log

# find empty files and delete
find . -size 0 -delete

```
4. Putative Paralog Detection (PPD) Analysis
```
# PPD Analysis

module load picard
module load gatk/4.1.7.0 
module load mafft
module load java/1.9.0 
module load bwa
module load samtools
module load python/2.7.10
```

PART 1: generates contigs with degenerate sequence does where heterozygotes are detected

Step1: Concatenate all the supercontigs into a single file. All supercontigs are stored in two new files named supercontig and exon. This script is modified from Mossmatters github "Alleles from HybSeq Data".
```
# Please type in the following arguments after Step1.sh in order: 
# * The directory of HybPiper.   ./
# * full name of namelist.  /mnt/ceph/masantos/scratch/tricalysia/namelist_tri5.txt

bash /mnt/ceph/masantos/putative_paralog/Step1.sh ./ /mnt/ceph/masantos/scratch/tricalysia/namelist_tri5.txt

bash /mnt/ceph/masantos/putative_paralog/Step1.sh ./ /mnt/ceph/masantos/scratch/tricalysia/namelist_tri5_2.txt
```

Step2: Generate the degenerated seuqences using IPUAC codes. This script is modified from Mossmatters github "Alleles from HybSeq Data".

 Please type in the following arguments after Step2.sh in order: 
1. BWA -k matching base length (if raw reads are 150 bp, we recommand use 100 bp here)
2. ploidy number of a species (2, 3, 4 ...)
3. a direcotory containing consensus sequences from Step1
4. a directory containing all raw reads file after QC control (keep raw reads in fastq.gz format)
5. output directory (to keep the degerated sequences file, which will be used for Step3)
6. full name of namelist
 
 ```
 bash /mnt/ceph/masantos/putative_paralog/Step2.sh 100 2 /mnt/ceph/masantos/scratch/tricalysia/supercontig /mnt/ceph/masantos/rawreads PP_Step2_output /mnt/ceph/masantos/scratch/tricalysia/namelist_tri.txt

 ```
PART 2: trims alignment and detectes the putative paralogs
This script only works when you use python/2.7.10
 
 ```
python /mnt/ceph/masantos/putative_paralog/PPD.py -ifa /mnt/ceph/masantos/scratch/tricalysia/PPD/PP_Step2_output -ina /mnt/ceph/masantos/scratch/tricalysia/namelist_tri6.txt -iref /mnt/ceph/masantos/scratch/tricalysia/Seq_Length_Hybpiper_tri_mega353plus.txt -io outgroups2.txt -o ./ -t supercontig -he 0.05 -gt 0.51 -hs 0.5 -nh 1 -w 10 -mi 10 -mo 8
 
 ```
 5. Remove samples that have long branches that are attributed with low data

```
# removing long branches  = "bad samples"

python3 /Users/maliasantos/Desktop/AMAS-master/amas/AMAS.py remove -x 2020_366 2020_156 2020_151 2019_961 2020_365 2019_912 2019_957 2020_147 2020_367 2020_364 2020_166 2019_956 2020_158 2019_921 2019_836 -d dna -f fasta -i *fasta -g removed_long_branches

```
6. After PPD analysis, we still had some hypervariable regions in our alignments. We manually edited and removed these sections by eye in AliView. 

7. IQ-tree: Maximum-Liklihood Analysis

```

# CONCAT species tree : All genes
iqtree2 -p MASedit -m TESTMERGE -T auto --prefix PPD_tree_part_7Feb22 -B 1000


# GENE TREES - on manually edited genes from PPD analysis (R3)
iqtree2 -S MASedit -m TESTMERGE -T auto --prefix PPD_R3_tree_genetrees_MASedit_7Feb22
```
8. Running ASTRAL to account for ILS

```
# RUNNING ASTRAL
java -jar astral.5.7.8.jar -i /mnt/ceph/masantos/scratch/tricalysia/result/supercontig/s8_rm_paralogs/PPD_R3_tree_genetrees_MASedit_7Feb22 -o PPD_R3_MASedit_ASTRALtree_7Feb22.tre
```