##################################################################01_READ_PROCESSING##################################################################
#####FILTER RAW READS (FASTP)
fastp --in1 MnAnt02_R1_001.fastq.gz --in2 MnAnt02_R2_001.fastq.gz --out1 MnAnt02_R1.trimmed.fastq.gz --out2 MnAnt02_R2.trimmed.fastq.gz -l 75 -q 15 -h MnAnt02.html &> MnAnt02.log
fastp --in1 MnAnt02_R1_001.fastq.gz --in2 MnAnt02_R2_001.fastq.gz --out1 MnAnt02_R1.trimmed.fastq.gz --out2 MnAnt02_R2.trimmed.fastq.gz -l 100 -q 30 --adapter_sequence GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2 GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --trim_poly_g -h MnAnt02.html &> MnAnt02.log

#####HIGH QUALITY READS MAP TO REFERENCE GENOME (BWA)
bwa mem -t 44 /home/celeminamaro/Genomes/Humpback_Whale/GCA_041834305.1_ASM4183430v1_genomic.fna MnAnt02_R1.trimmed.fastq.gz MnAnt02_R2.trimmed.fastq.gz > MnAnt02.sam

#####CONVERT SAM TO BAM (SAMTOOLS)
samtools view -@ 4 -b MnAnt02.sam > MnAnt02.bam

#####SORT BAM (SAMTOOLS)
samtools sort -@ 4 -o MnAnt02_sorted.bam MnAnt02.bam

#####ADD READ GROUPS (PICARD)
java -jar /home/celeminamaro/picard.jar AddOrReplaceReadGroups I=MnAnt02_sorted.bam O=MnAnt02_sorted_RG.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=MnAnt02
  
#####REMOVE PCR AND SEQUENCING DUPLICATES (PICARD)
java -jar /home/celeminamaro/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=MnAnt02_sorted_RG.bam O=MnAnt02_sorted_RG_MD.bam M=MnAnt02_marked_dup_metrics.txt

####REALIGNMENT AROUND INDELS (GATK3)
#1st:Create a target list of intervals to be realigned
java -jar /home/celeminamaro/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /home/celeminamaro/Genomes/Humpback_Whale/GCA_041834305.1_ASM4183430v1_genomic.fna -I MnAnt02_sorted_RG_MD.bam -o MnAnt02_target_intervals.list

#2nd:Realign
java -jar /home/celeminamaro/GenomeAnalysisTK.jar -T IndelRealigner -R /home/celeminamaro/Genomes/Humpback_Whale/GCA_041834305.1_ASM4183430v1_genomic.fna -I MnAnt02_sorted_RG_MD.bam -targetIntervals MnAnt02_target_intervals.list -o MnAnt02_sorted_RG_MD_realigned.bam

####SEX CHROMOSOMES FILTERING (SAMTOOLS)
samtools view -@ 4 -b -h -L mapped_scaffolds.txt MnAnt02_sorted_RG_MD_realigned.bam -U MnAnt02_sorted_RG_MD_realigned_filtered.bam

####REPETITIVE REGIONS FILTERING (SAMTOOLS)
samtools view -@ 4 -b -h -L /home/celeminamaro/Genomes/Humpback_Whale/GCA_041834305.1_ASM4183430v1_genomic_repetitive_regions.bed MnAnt02_sorted_RG_MD_realigned_filtered.bam -U MnAnt02_sorted_RG_MD_realigned_filtered2.bam

#####COVERAGE FILTERING 
#1st:Get coverage statistics (ANGSD)
angsd -i MnAnt02_sorted_RG_MD_realigned_filtered2.bam -doDepth 1 -out MnAnt02_pre_filtering -doCounts 1 -minMapQ 30 -minQ 20 -maxDepth 1000

#2nd:Calculate intermediate coverage (Rstudio)
depth<- read.table(file.choose(), header=F, sep = "\t", dec = ",")
total_number<-sum(depth$V2)
Percentage <- depth$V2 * 100 /total_number
weighted.mean(depth$V1, Percentage)

#3rd:Identify regions of too low/large covergae (GATK3)
java -jar /home/celeminamaro/GenomeAnalysisTK.jar -T CallableLoci -R /home/celeminamaro/Genomes/Humpback_Whale/GCA_041834305.1_ASM4183430v1_genomic.fna -I MnAnt02_sorted_RG_MD_realigned_filtered2.bam -maxDepth 12 --minBaseQuality 20 --minMappingQuality 30 --minDepth 2 --summary MnAnt02.summary -o MnAnt02.Callable.bed

#4th:Edit callable file
sed '/CALLABLE/!d' MnAnt02.Callable.bed > MnAnt02_Callable_edited.bed

#5th:Keep only callable regions of the genome (SAMTOOLS)
samtools view -@ 4 -b -h -L MnAnt02_Callable_edited.bed  MnAnt02_sorted_RG_MD_realigned_filtered2.bam > MnAnt02_sorted_RG_MD_realigned_filtered3.bam

#5th:Calculate final coverage (ANGSD)
angsd -i MnAnt02_sorted_RG_MD_realigned_filtered3.bam -doDepth 1 -doCounts 1 -minMapQ 30 -minQ 20 -maxDepth 1000 -out MnAnt02_post_coverage_filtering

#####DOWNSAMPLE TO 5X (PICARD)
#1st:Downsample
java -jar /home/celeminamaro/picard.jar DownsampleSam I=/work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/02_Final_nuclear_bam_files/MnAnt02_sorted_RG_MD_realigned_filtered3.bam O=MnAnt02_sorted_RG_MD_realigned_filtered3_downsampled.bam P=0.85

#2nd:Corroborate downsampling
angsd -i MnAnt02_sorted_RG_MD_realigned_filtered3_downsampled.bam -doDepth 1 -doCounts 1 -minMapQ 30 -minQ 20 -maxDepth 1000 -out MnAnt02_post_downsampling_coverage
######################################################################################################################################################

#######################################################################02_SATC########################################################################
#1st:Get idxstats from the _sorted_RG_MD.bam files (SAMTOOLS)
samtools idxstats MnAnt02_sorted_RG_MD.bam > MnAnt02_sorted_RG_MD.idxstats

#2nd:Identify sex by looking at coverage (SATC)
Rscript --vanilla /home/celeminamaro/SATC/satc.R -i sample_list.txt -o HW_downsampled_median --useMedian TRUE
######################################################################################################################################################

######################################################################03_NGSRELATE###################################################################
#1st:Estimate genotype likelihoods (ANGSD)
angsd -GL 1 -out Humpback_whale_relate -nThreads 10 -doGlf 3 -doMajorMinor 1 -minMaf 0.05 -minInd 15 -doMaf 3 -SNP_pval 1e-6 -bam /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/Final_nuclear_bam_files/HW_26_bam_list.txt

#2nd:Get allele frequencies
zcat Humpback_whale_relate.mafs.gz | cut -f5 |sed 1d > freq

#3rd:Estimate relatedness statistics (ANGSD)
./ngsRelate -g Humpback_whale_relate.glf.gz -n 26 -f freq -O NGSRelate_26_humpbacks

#4th:Plot (RSTUDIO)
######################################################################################################################################################

#################################################################05_GENOTYPE_LIKELIHOODS##############################################################
#1st:Estimate genotype likelihoods (ANGSD)
angsd -GL 1 -out 26_Humpback_whale -nThreads 10 -doGeno 4 -doBcf 1 -doPost 1 -doGlf 2 -doMajorMinor 1 -minMaf 0.05 -minInd 19 -doMaf 3 -SNP_pval 1e-6 -bam /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/Final_nuclear_bam_files/HW_26_bam_list.txt

#2nd:Estimate LD (NGSLD)
zcat 26_Humpback_whale.beagle.gz | awk '{print $1, $2}' > 26_Humpback_whale_SNPs_pos.txt
sed 's//\t/g' 26_Humpback_whale_SNPs_pos.txt > newfile.txt
awk '{print $1, $2}' newfile.txt > test.txt
sed 's//\t/g' test.txt > good.txt
mv good.txt 26_Humpback_whale_SNPs_pos.txt
cd /work/celeminamaro/celeminamaro/celeminamaro/bam_files/Genotype_likelihoods/ngsTools/ngsLD/
./ngsLD --geno /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/05_Genotype_likelihoods/26_Humpback_whales/26_Humpback_whale.beagle.gz --probs --n_ind 26 --n_sites 5695001 --n_threads 4 --max_kb_dist 20 --min_maf 0.05 --seed 1 --posH /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/05_Genotype_likelihoods/26_Humpback_whales/26_Humpback_whale_SNPs_pos.txt --out /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/Genotype_likelihoods/26_Humpback_whales/26_Humpback_whale_LD_20Kb

#3rd:Identify unlinked SNPs (NGSLD)
cd /work/celeminamaro/celeminamaro/celeminamaro/bam_files/Genotype_likelihoods/ngsTools/ngsLD/scripts/prune_graph/target/release
./prune_graph --in /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/05_Genotype_likelihoods/26_Humpback_whales/26_Humpback_whale_LD_20Kb --weight-field column_4 --weight-filter "column_4 >= 0.5" --out /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/Genotype_likelihoods/26_Humpback_whales/26_Humpback_whale_LD_unlinked.pos

#4th:Remove linked SNPs from the beagle file
awk '{gsub(":", "_");print}' 26_Humpback_whale_LD_unlinked.pos > 26_Humpback_whale_LD_unlinked_edited.pos
cat <(zcat 26_Humpback_whale.beagle.gz | head -n1) <(zcat 26_Humpback_whale.beagle.gz | awk 'NR==FNR{a[$1]=$1OFS$2;next}{$1=a[$1];print}' OFS='\t' 26_Humpback_whale_LD_unlinked_edited.pos - | grep -Pv "^\t") | gzip > 26_Humpback_whale_unlinked.beagle.gz

#5th:Filter related individuals from .beagle file
#24 Humpback whales
zcat /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/05_Genotype_likelihoods/Downsampled/26_Humpback_whales/26_Humpback_whale_downsampled_unlinked.beagle.gz | awk -F'\t' '{for (i=1; i<=81; i++) if (i != 10 && i != 11 && i != 12 && i != 76 && i != 77 && i != 78) printf "%s%s", $i, (i==81?"\n":"\t")}' > output.beagle
zcat /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/05_Genotype_likelihoods/Downsampled/26_Humpback_whales/26_Humpback_whale_downsampled_unlinked.beagle.gz | awk '{print $81}' > last_column.txt
paste output.beagle /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/05_Genotype_likelihoods/Downsampled/22_Humpback_whales/last_column.txt > 24_Humpback_whale_downsampled_unlinked.beagle
gzip 24_Humpback_whale_downsampled_unlinked.beagle

#22 Humpback whales
zcat /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/05_Genotype_likelihoods/Downsampled/26_Humpback_whales/26_Humpback_whale_downsampled_unlinked.beagle.gz | awk -F'\t' '{for (i=1; i<=81; i++) if (i != 10 && i != 11 && i != 12 && i != 25 && i != 26 && i != 27 && i != 55 && i != 56 && i != 57 && i != 76 && i != 77 && i != 78) printf "%s%s", $i, (i==81?"\n":"\t")}' > output.beagle
zcat /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/05_Genotype_likelihoods/Downsampled/26_Humpback_whales/26_Humpback_whale_downsampled_unlinked.beagle.gz | awk '{print $81}' > last_column.txt
paste output.beagle last_column.txt > 22_Humpback_whale_downsampled_unlinked.beagle
gzip 22_Humpback_whale_downsampled_unlinked.beagle
######################################################################################################################################################

#######################################################################05_PCANGSD#####################################################################
#1st:Estimate PCA for 26 HWs (PCANGSD)
pcangsd -b 66_Humpback_whale_downsampled_unlinked.beagle.gz -t 64 -o 66_Humpback_whales

#2nd:Estimate PCA for 24 HWs (PCANGSD)
pcangsd -b 24_Humpback_whale_unlinked.beagle.gz -t 64 -o 24_Humpback_whales

#3rd:Estimate PCA for 22 HWs (PCANGSD)
pcangsd -b 22_Humpback_whale_unlinked.beagle.gz -t 64 -o 22_Humpback_whales

#4th:Plot (RSTUDIO)
######################################################################################################################################################

#######################################################################06_NGSADMIX####################################################################
#1st: Estimate Admixture for 26 HWs (NGADMIX)
./create_NGS_admix_jobs_26_HW.sh

#2nd: Estimate Admixture for 22 HWs (NGADMIX)
./create_NGS_admix_jobs_22_HW.sh

#3rd:Plot (RSTUDIO)
######################################################################################################################################################

####################################################################07_HETEROZYGOSITY#################################################################
#1st:Estimate per individual saf file (ANGSD)
angsd -GL 1 -ref /home/celeminamaro/Genomes/Humpback_Whale/GCA_041834305.1_ASM4183430v1_genomic.fna -anc /home/celeminamaro/Genomes/Humpback_Whale/GCA_041834305.1_ASM4183430v1_genomic.fna -nThreads $SLURM_CPUS_PER_TASK -dosaf 1 -i /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/Final_nuclear_bam_files/MnAnt02_sorted_RG_MD_realigned_filtered3.bam -out DoSaf_MnAnt02

#2nd: Estimate per individual folded SFS (ANGSD)
realSFS DoSaf_MnAnt02.saf.idx -fold 1 -P $SLURM_CPUS_PER_TASK > MnAnt02_folded.sfs

#3rd: Calculate Heterozygosity (RSTUDIO)
sample_MnAnt02<-scan("MnAnt02_folded.sfs")
sample_MnAnt02_heterozygosity<-data.frame(sample_MnAnt02[2]/sum(sample_MnAnt02))

#4th:Plot (RSTUDIO)
######################################################################################################################################################

######################################################################08_MITOGENOME###################################################################
#1st:Map high quality reads to mitogenome (BWA)
bwa mem -t 44 /home/celeminamaro/Genomes/Humpback_Whale/Megaptera_novaengliae_mitogenome.fasta /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/Read_processing/MnAnt02_R1.trimmed.fastq.gz /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/Read_processing/MnAnt02_R2.trimmed.fastq.gz > MnAnt02.sam

#2nd:Convert sam to bam (SAMTOOLS)
samtools view -@ 4 -b MnAnt02.sam > MnAnt02.bam

#3rd:Sort bam (SAMTOOLS)
samtools sort -@ 4 -o MnAnt02_sorted.bam MnAnt02.bam

#4th:Add read groups (PICARD)
java -jar /home/celeminamaro/picard.jar AddOrReplaceReadGroups I=MnAnt02_sorted.bam O=MnAnt02_sorted_RG.bam RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=MnAnt02

#5th:Remove PCR and sequencing duplicates (PICARD)
java -jar /home/celeminamaro/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=MnAnt02_sorted_RG.bam O=MnAnt02_sorted_RG_MD.bam M=MnAnt02_marked_dup_metrics.txt

#6th:Get index file (SAMTOOLS)
samtools index MnAnt02_sorted_RG_MD.bam

#7th:Convert bam to fasta (ANGSD)
angsd -i MnAnt02_sorted_RG_MD.bam -ref /home/celeminamaro/Genomes/Humpback_Whale/Megaptera_novaengliae_mitogenome.fasta -remove_bads 1 -uniqueOnly 1 -baq 1 -C 50 -minMapQ 15 -only_proper_pairs 0 -minQ 15 -doCounts 1 -setMinDepth 6 -doDepth 1 -doFasta 2 -out MnAnt02

#8th:Haplotype network (POPART)

#9th:Phylogenetic tree with fin whale outgroup (W-IQTREE)
######################################################################################################################################################

######################################################################08_GENOTYPES####################################################################
#1st:Edit 26_Humpback_whale.geno file to put / between alleles
awk '{ for (i=1; i<=NF; i++) $i = substr($i,1,1) "/" substr($i,2,1); print }' 26_Humpback_whale_downsampled.geno > intermediate.geno

#2nd:Remove two first columns
awk '{for (i=3; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' intermediate.geno > intermediate2.geno

#3rd:Copy the firs two columns of 26_Humpback_whale.geno file to intermediate2.geno file
awk '{print $1, $2}' 26_Humpback_whale_downsampled.geno > two_first_columns.geno
paste two_first_columns.geno intermediate2.geno > intermediate3.geno

#4th:Copy the header row to intermediate3.geno file
cat header_geno_file.txt intermediate3.geno > intermediate4.geno

#5th:Convert white space in tab-separated
awk '{ gsub(/ +/, "\t"); print }' intermediate4.geno > intermediate5.geno
mv intermediate5.geno 26_Humpback_whale_downsampled_edited.geno

#6th:Convert geno in vcf (genoToVCF)
python genoToVCF.py -g /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/08_Genotypes/Downsampled/26_Humpback_whale_downsampled_edited.geno -f phased -r /home/celeminamaro/Genomes/Humpback_Whale/GCA_041834305.1_ASM4183430v1_genomic_CAPITAL_LETTERS.fna -o /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/08_Genotypes/Downsampled/26_Humpback_whale_downsampled.vcf

#7th: Annotate vcf (bcftools)
bcftools +fill-tags 26_Humpback_whale_downsampled.vcf > 26_Humpback_whale_downsampled_annotated.vcf

#8th: Filter vcf (bcftools)
vcftools --vcf 26_Humpback_whale_downsampled_annotated.vcf --min-alleles 2 --mac 2 --max-missing 1.0 --hwe 0.001 --recode --out 26_Humpback_whale_downsampled_annotated_filtered.vcf
gzip 26_Humpback_whale_downsampled_annotated_filtered.vcf
#####################################################################09_INBREEDING####################################################################
#1st:Convert vcf to plink (PLINK)
plink --allow-extra-chr --vcf /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/08_Genotypes/Downsampled/26_Humpback_whale_downsampled_annotated_filtered.vcf.gz --make-bed --out 26_Humpback_whale_downsampled_annotated_filtered

#2nd:Estimate ROH (PLINK)
plink --bfile 26_Humpback_whale_downsampled_annotated_filtered --allow-extra-chr --homozyg-snp 50 --homozyg-kb 100 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 3 --homozyg-window-missing 10 --homozyg-window-threshold 0.05 --out 26_Humpback_whale_annotated_filtered

#3rd:Estimate F (VCFTOOLS)
vcftools --gzvcf /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/08_Genotypes/Downsampled/26_Humpback_whale_downsampled_annotated_filtered.vcf.gz --het --out 26_Humpback_whale_downsample_het

#3rd:Plot ROH and F (RSTUDIO)
######################################################################################################################################################

##########################################################################10_SMC++####################################################################
#1st:Convert vcf into smc files (SMC++)
./ALL_cl_smc++_vcf2smc.sh
#2nd:Estimate population size (SMC++)
./ALL_cl_smc++_estimate.sh
#3rd:Get csv file (SMC++)
./ALL_cl_smc++_plot.sh
#4th:Plot Ne trough time (RSTUDIO)
######################################################################################################################################################

###########################################################################11_GONE####################################################################
#1st:Keep only 22 samples (VCFTOOLS)
vcftools --gzvcf 26_Humpback_whale_downsampled_annotated_filtered.vcf.gz --keep /work/celeminamaro/celeminamaro/celeminamaro/Humpback_raw_data/11_GONE/Downsampled/ALL/ALL_samples_to_keep.txt --recode --out output
#2nd:Convert vcf to .ped and .map (PLINK)
plink --allow-extra-chr --vcf --vcf Humpback_whale_downsampled_annotated_filtered_ALL.vcf --recode --out Humpback_whale_downsampled_annotated_filtered_ALL
#3rd:Edit .map file to put correct chromosome numbers
sed -E 's/\bCM085695\.1\b/1/g; s/\bCM085696\.1\b/2/g; s/\bCM085697\.1\b/3/g; s/\bCM085698\.1\b/4/g; s/\bCM085699\.1\b/5/g; s/\bCM085700\.1\b/6/g; s/\bCM085701\.1\b/7/g; s/\bCM085702\.1\b/8/g; s/\bCM085703\.1\b/9/g; s/\bCM085704\.1\b/10/g; s/\bCM085705\.1\b/11/g; s/\bCM085706\.1\b/12/g; s/\bCM085707\.1\b/13/g; s/\bCM085708\.1\b/14/g; s/\bCM085709\.1\b/15/g; s/\bCM085710\.1\b/16/g; s/\bCM085711\.1\b/17/g; s/\bCM085712\.1\b/18/g; s/\bCM085713\.1\b/19/g; s/\bCM085714\.1\b/20/g; s/\bCM085715\.1\b/21/g' Humpback_whale_downsampled_annotated_filtered_ALL.map > Humpback_whale_downsampled_annotated_filtered_ALL_edited.map
#4th:Run 100 replicates (GONE)
./GONE_ALL.sh
#5th:Plot Ne trough time (RSTUDIO)
######################################################################################################################################################
