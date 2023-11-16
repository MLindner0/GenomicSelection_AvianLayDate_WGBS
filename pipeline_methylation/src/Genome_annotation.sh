#!/bin/bash

### genome annotation

# get chr length file:
sed 's/Linkage Group/Linkage_Group/g' genome/GCF_001522545.3_Parus_major1.1_assembly_report.txt | sed 's/Primary Assembly/Primary_Assembly/g' | grep -vP '^#' | awk -v OFS='\t' '{print $7,$9}' > genome/chr_length
sed 's/Linkage Group/Linkage_Group/g' genome/GCF_001522545.3_Parus_major1.1_assembly_report.txt | sed 's/Primary Assembly/Primary_Assembly/g' | grep -vP '^#' | awk -v OFS='\t' '{print $1,$3,$7,$9}' | sed 's/na/Sc/g' > genome/chr_length_names

# get annotations:
grep -vP '^#' annotation/genes.gff3 | awk -v OFS='\t' '{ split($9,a,";"); split(a[1],b,"="); split(a[2],c,"="); print $1,$4,$5,$7,b[2],c[2],"gene"}' > annotation/genes.out
grep -vP '^#' annotation/downstream.laine.t.gff3 | awk -v OFS='\t' '{ split($9,a,";"); split(a[1],b,"="); split(a[2],c,"="); print $1,$4,$5,$7,b[2],c[2],"downstream10k"}' > annotation/downstream10k.out
grep -vP '^#' annotation/upstream.laine.t.gff3 | awk -v OFS='\t' '{ split($9,a,";"); split(a[1],b,"="); split(a[2],c,"="); print $1,$4,$5,$7,b[2],c[2],"upstream10k"}' > annotation/upstream10k.out

grep -vP '^#' annotation/promoters_2k.t.gff3 | awk -v OFS='\t' '{ split($9,a,";"); split(a[2],b,"="); split(a[3],c,"="); print $1,$4,$5,$7,b[2],c[2],"promoter"}' > annotation/promoter.out
grep -vP '^#' annotation/TSS.laine.t.gff3 | awk -v OFS='\t' '{ split($9,a,";"); split(a[2],b,"="); split(a[3],c,"="); print $1,$4,$5,$7,b[2],c[2],"tss"}' > annotation/tss.out

# remove wrong annotations
cat annotation/downstream10k.out | awk -v OFS='\t' '{print $0, $3-$2}' | awk -v OFS='\t' '{ if(($8 <= 0)) { print } }' > annotation/remove
cat annotation/downstream10k.out | awk -v OFS='\t' '{print $0, $3-$2}' | awk -v OFS='\t' '{ if(($8 > 0)) { print } }' | awk -v OFS='\t' '{$8=""; print $0}' > annotation/downstream10k_F.out
cat annotation/upstream10k.out | awk -v OFS='\t' '{print $0, $3-$2}' | awk -v OFS='\t' '{ if(($8 <= 0)) { print } }' >> annotation/remove
cat annotation/upstream10k.out | awk -v OFS='\t' '{print $0, $3-$2}' | awk -v OFS='\t' '{ if(($8 > 0)) { print } }' | awk -v OFS='\t' '{$8=""; print $0}' > annotation/upstream10k_F.out
cat annotation/promoter.out | awk -v OFS='\t' '{print $0, $3-$2}' | awk -v OFS='\t' '{ if(($8 <= 0)) { print } }' >> annotation/remove
cat annotation/promoter.out | awk -v OFS='\t' '{print $0, $3-$2}' | awk -v OFS='\t' '{ if(($8 > 0)) { print } }' | awk -v OFS='\t' '{$8=""; print $0}' > annotation/promoter_F.out
cat annotation/tss.out | awk -v OFS='\t' '{print $0, $3-$2}' | awk -v OFS='\t' '{ if(($8 <= 0)) { print } }' >> annotation/remove
cat annotation/tss.out | awk -v OFS='\t' '{print $0, $3-$2}' | awk -v OFS='\t' '{ if(($8 > 0)) { print } }' | awk -v OFS='\t' '{$8=""; print $0}' > annotation/tss_F.out

# combine annotations:
cat annotation/genes.out annotation/downstream10k_F.out annotation/upstream10k_F.out annotation/promoter_F.out annotation/tss_F.out | sort -k1,1 -k2,2n > annotation/annotations.bed

# get and sort all CpG sites in wgbs data
cat temp/CpGs_all | sort -k1,1 -k2,2n > temp/CpGs.bed
# cat temp/CpGs_all_pools  | sort -k1,1 -k2,2n > temp/CpGs_pools.bed

# get overlap between CpGs and annotations
intersectBed -a temp/CpGs.bed -b annotation/annotations.bed -wo > temp/annotated_CpGs
# intersectBed -a temp/CpGs_pools.bed -b annotation/annotations.bed -wo > temp/annotated_CpGs_pools

cat temp/annotated_CpGs | awk -v OFS='\t' '{print $1"_"$2,$1,$2,$10,$8,$9,$7}' > temp/annotated_CpGs.out
# cat temp/annotated_CpGs_pools | awk -v OFS='\t' '{print $1"_"$2,$1,$2,$10,$8,$9,$7}' > temp/annotated_CpGs_pools.out

# get annotations for CGIs
cat temp/CGIs | grep -vP '^#' | awk -v OFS='\t' '{print $1,$4,$5}' | sort -k1,1 -k2,2n > temp/CGIs.bed
intersectBed -a temp/CGIs.bed -b annotation/annotations.bed -wo > temp/annotated_CGIs
cat temp/annotated_CGIs | awk -v OFS='\t' '{print $1"_"$2"_"$3,$1,$2,$3,$10,$8,$9,$7}' > temp/annotated_CGIs.out

# get overlap between CGIs and CpGs to count number of CpGs per CGI
cat temp/CGIs | grep -vP '^#' | awk -v OFS='\t' '{split($9,a,"="); print $1,$4,$5,$7,a[2]}' | sort -k1,1 -k2,2n > temp/CGIs.0.bed
intersectBed -a temp/CGIs.0.bed -b temp/CpGs.bed -wo > temp/CpGs_in_CGI
