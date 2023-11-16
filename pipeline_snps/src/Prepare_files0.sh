#!/bin/bash

# split alignments; remove MT (NC_040875.1) and Z chromosome (NC_031799.1) form mail file and make an extra alignment for Z chr
# grep -vP '^#' alignments/GCF_001522545.3_Parus_major1.1_assembly_report.txt | sed 's/Linkage Group/Linkage_Group/g' | sed 's/Primary Assembly/Primary_Assembly/g' | awk -v OFS='\t' '{print $7,"0",$9}' > alignments/chr_length
cat ../WGBS_Snakemake_Bismark/genome/chr_length | grep -v "NC_031799.1" | grep -v "NC_040875.1" | awk -v OFS='\t' '{print $1,"0",$2}' > genome/chr_keep.bed

# merge chromosomes per sample
picard -Xmx2g MergeVcfs \
          I=snps/W_ML_Pool4.Auto_SNPs_DPfilter_GT_Name.vcf \
          I=snps/W_ML_Pool4.Z_SNPs_DPfilter_GT_Name.vcf \
          O=test.vcf


# Format/Annotate vcf for plink
bcftools norm -O u -m -any snps/SNPs_merged.vcf.gz | bcftools norm -Ou -f ../wgbs_snakemake_reseq/genome/reference.fa | bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' -o snps/SNPs_merged.bcf


### PLINK
conda activate vcftools

## 1 INDIVIDUALS
# get vcf file into PLINK
plink --bcf snps/SNPs_Individuals.bcf --keep-allele-order --make-bed --const-fid --allow-extra-chr --chr-set 32 --out plink/individuals_raw/ERC_wgbs

# filter snps in GENEABLE
plink --bfile plink/individuals_raw/ERC_wgbs --allow-extra-chr --chr-set 32 --recode transpose --out plink/individuals_raw/ERC_wgbs_trans

# only keep good snps
plink --bfile plink/individuals_raw/ERC_wgbs --extract plink/individuals_.9_QC_ok/QC_good_snps -allow-extra-chr --chr-set 32 --make-bed -out plink/individuals_.9_QC_ok/ERC_wgbs


### Fst
## Fst plink
plink --bfile plink/individuals_.9_QC_ok/ERC_wgbs --within plink/individuals_QC_ok/groups -allow-extra-chr --chr-set 32 --fst --out plink/individuals_.9_QC_ok/Fst_ERC_wgbs

plink --bfile plink/individuals_.9_QC_ok/ERC_wgbs --within plink/individuals_QC_ok/groups -allow-extra-chr --chr-set 32 --hardy --out plink/individuals_.9_QC_ok/Hardy_ERC_wgbs

# fst vcftool (slifing window)
# plink --bfile plink/individuals_QC_ok/ERC_wgbs --allow-extra-chr --chr-set 32 --recode vcf-iid --out plink/individuals_QC_ok/vcf_ERC_wgbs
# cat plink/individuals_QC_ok/groups | awk '{if($3=="E") print $2}' > plink/individuals_QC_ok/early
# cat plink/individuals_QC_ok/groups | awk '{if($3=="L") print $2}' > plink/individuals_QC_ok/late
# vcftools --vcf plink/individuals_QC_ok/vcf_ERC_wgbs.vcf --weir-fst-pop plink/individuals_QC_ok/early --weir-fst-pop plink/individuals_QC_ok/late --fst-window-size 200000 --fst-window-step 50000 --out plink/individuals_QC_ok/Fst_window_ERC_wgbs


## Fst Arlequin
plink --bfile plink/individuals_.9_QC_ok/ERC_wgbs --allow-extra-chr --chr-set 32 --pheno plink/individuals_.9_QC_ok/groups_num --pheno-merge --update-sex plink/individuals_.9_QC_ok/sex --recode --out plink/individuals_.9_QC_ok_Arlequin/ERC_wgbs

PGDSpider2-cli  -Xmx5g -Xms512M -inputfile plink/individuals_.9_QC_ok_Arlequin/ERC_wgbs.ped -inputformat PED -outputfile plink/individuals_.9_QC_ok_Arlequin/ERC_wgbs.arp -outputformat ARLEQUIN -spid plink/individuals_.9_QC_ok_Arlequin/ERC_wgbs_arlequin.spid

# adding structure to the end of PGD conversion file (got it from Arlequin GUI-version), PGD doesn't add structure. Line is the group.
cat plink/individuals_QC_ok_Arlequin/struc >> plink/individuals_.9_QC_ok_Arlequin/ERC_wgbs.arp
./LaunchArlecore.sh


## Fst bayeScan
plink --bfile plink/individuals_.9_QC_ok/ERC_wgbs --allow-extra-chr --chr-set 32 --pheno plink/individuals_.9_QC_ok/groups_num --pheno-merge --update-sex plink/individuals_.9_QC_ok/sex --recode A --out plink/individuals_.9_QC_ok_Bayescan/ERC_wgbs

# use R packages adegenet and dartR to get snp data into BayeScan format
# R script Fst.R
# activate conda enviroment bayescan

# run BayeScan
/home/nioo/melaniel/Software/BayeScan2.1/binaries/BayeScan2.1_linux64bits plink/individuals_.9_QC_ok_Bayescan/ERC_wgbs.bayescan.txt -od plink/individuals_.9_QC_ok_Bayescan -threads 60 -pr_odds 1000

##  prior odds for neutral model; 1 would mean that for every locus, one assumes a priori that the model inlcuding selection is as likely as teh neutral model. This can lead to false positives when testing a large number of markers. Higgher prior odds will tend to eleminate false positives, but at the cost of reducing the power to detect any marker under selection. 10 seems reasonable for the identification of condidate loci within a few hunderds of markers, whereas values up to 10 000 are generally used in the context of genome wide association studies with millions of SNPs (idetification of top candidates)

# try lower prior oddds:
/home/nioo/melaniel/Software/BayeScan2.1/binaries/BayeScan2.1_linux64bits plink/individuals_.9_QC_ok_Bayescan/ERC_wgbs.bayescan.txt -od plink/individuals_.9_QC_ok_Bayescan_pr_odds100 -threads 60 -pr_odds 100

# try lower prior oddds + thin 50 and n 25000: later; server too busy
/home/nioo/melaniel/Software/BayeScan2.1/binaries/BayeScan2.1_linux64bits plink/individuals_.9_QC_ok_Bayescan/ERC_wgbs.bayescan.txt -od plink/individuals_.9_QC_ok_Bayescan_pr_odds100_more_samples -threads 60 -pr_odds 100 -n 25000 -thin 50


### RAiSD

plink --bfile plink/individuals_.9_QC_ok/ERC_wgbs --allow-extra-chr --chr-set 32 --pheno plink/individuals_.9_QC_ok/groups_num --pheno-merge --update-sex plink/individuals_.9_QC_ok/sex --recode vcf-iid --out plink/individuals_.9_QC_ok_RAiSD/ERC_wgbs

# run separately for early and late females
# early
cat plink/individuals_.9_QC_ok/groups_num | awk '{if($3=="1") print $2}' > plink/individuals_.9_QC_ok/early_females

for i in `cat plink/individuals_.9_QC_ok/early_females`; do echo "awk '{if(\$2==\"$i\"){print \$1,\$2}}' plink/individuals_.9_QC_ok/ERC_wgbs.fam";done | bash > plink/individuals_.9_QC_ok/early_females.out

# get snp data for early females
plink --bfile plink/individuals_.9_QC_ok/ERC_wgbs --keep plink/individuals_.9_QC_ok/early_females.out --allow-extra-chr --chr-set 32 --pheno plink/individuals_.9_QC_ok/groups_num --pheno-merge --update-sex plink/individuals_.9_QC_ok/sex --recode vcf-iid -out plink/individuals_.9_QC_ok_RAiSD_Early/ERC_wgbs

# run RAiSD
cd plink/individuals_.9_QC_ok_RAiSD_Early/
RAiSD -n ERC_wgbs_raisd -I ERC_wgbs.vcf -s -t -R

# format summary output
grep 'Set N' plink/individuals_.9_QC_ok_RAiSD_Early/RAiSD_Info.ERC_wgbs_raisd > plink/individuals_.9_QC_ok_RAiSD_Early/RAiSD_Info_formatted.ERC_wgbs_raisd


# late
cat plink/individuals_.9_QC_ok/groups_num | awk '{if($3=="2") print $2}' > plink/individuals_.9_QC_ok/late_females

for i in `cat plink/individuals_.9_QC_ok/early_females`; do echo "awk '{if(\$2==\"$i\"){print \$1,\$2}}' plink/individuals_.9_QC_ok/ERC_wgbs.fam";done | bash > plink/individuals_.9_QC_ok/early_females.out

# get snp data for late females
plink --bfile plink/individuals_.9_QC_ok/ERC_wgbs --keep plink/individuals_.9_QC_ok/early_females.out --allow-extra-chr --chr-set 32 --pheno plink/individuals_.9_QC_ok/groups_num --pheno-merge --update-sex plink/individuals_.9_QC_ok/sex --recode vcf-iid -out plink/individuals_.9_QC_ok_RAiSD_Early/ERC_wgbs

# run RAiSD
cd plink/individuals_.9_QC_ok_RAiSD_Early/
RAiSD -n ERC_wgbs_raisd -I ERC_wgbs.vcf -s -t -R

# format summary output
grep 'Set N' plink/individuals_.9_QC_ok_RAiSD_Late/RAiSD_Info.ERC_wgbs_raisd > plink/individuals_.9_QC_ok_RAiSD_Late/RAiSD_Info_formatted.ERC_wgbs_raisd


### RDA
plink --bfile plink/individuals_.9_QC_ok/ERC_wgbs -allow-extra-chr --chr-set 32 --recodeA --out plink/individuals_.9_QC_ok_RDA/ERC_wgbs


### get gene annotations for SNPs

# WGBS SNPs
cat temp/snp_locations | sort -k1,1 -k2,2n > temp/snp_locations.bed
intersectBed -a temp/snp_locations.bed -b temp/annotations.bed -wo > temp/annotated_Snps

# SNPs from Verhagen, I, Gienapp, P, Laine, VN, et al. Genetic and phenotypic responses to genomic selection for timing of breeding in a wild songbird. Funct Ecol. 2019; 33: 1708–1721. https://doi.org/10.1111/1365-2435.13360
cat temp/SNP.chip_SNP.positions.txt | sort -k1,1 -k2,2n > temp/SNP.chip_SNP.positions.bed
intersectBed -a temp/SNP.chip_SNP.positions.bed -b temp/annotations.bed -wo > temp/SNPchip.annotated_Snps
&& $4 > 16892000 && $4 < 16893000

# for RAIDS regions
cat temp/RAiDS_regions.txt | sort -k1,1 -k2,2n > temp/RAiDS_regions.bed
intersectBed -a temp/RAiDS_regions.bed -b temp/annotations.bed -wo > temp/annotated_RAiDS_regions
