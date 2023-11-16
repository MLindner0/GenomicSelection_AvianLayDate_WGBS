#!/bin/bash


### --------------------------------------------------------------------------------
###
### Prepare files & Notes
### Project: WGBS Snp/SV calling
###

### prepare files for snakemake pipeline and notes for making Snakefile


## run 2020
# creat symlink for files (in directory of raw data)
for f in *.gz; do ln -s /mnt/nfs/bioinfdata/primary/AnE/Parus_major/raw_data/whole_genome_sequencing/bisulfite_sequencing/novoseq/selection_lines2018_redo18062020/data/$f /home/NIOO.INT/melaniel/projects/WGBS_Snakemake_Bismark/raw_data/reseq2020/$f; done
# remae files
conda activate rename
rename 's/_........_L00[0-9]//' *.gz #remove adapter sequences and lane
rename 's/_001//' *.gz # remove _001
conda deactivate

## run 2018
# creat symlink for files (in directory of raw data)
for f in *fastq; do ln -s /mnt/nfs/bioinfdata/primary/AnE/Parus_major/raw_data/whole_genome_sequencing/bisulfite_sequencing/novoseq/selection_lines2018/$f /home/NIOO.INT/melaniel/projects/WGBS_Snakemake_Bismark/raw_data/seq2018/$f; done
# remae files
conda activate rename
rename 's/_........_L00[0-9]//' *.fastq #remove adapter sequences and lane
rename 's/_001//' *.fastq # remove _001
conda deactivate

## run 2018
# creat symlink for files (in directory of raw data)
for f in *fastq.gz; do ln -s /mnt/nfs/bioinfdata/primary/AnE/Parus_major/raw_data/whole_genome_sequencing/bisulfite_sequencing/novoseq/selection_lines2018_redo10112020/data/$f /home/NIOO.INT/melaniel/projects/WGBS_Snakemake_Bismark/raw_data/reseq2020_SP/$f; done
# remae files
conda activate rename
rename 's/_......_L00M//' *.fastq.gz #remove adapter sequences and lane
rename 's/_001//' *.fastq.gz # remove _001
# fix library names
rename 's/M1/F3_E_BD_27211/' *.fastq.gz
rename 's/M2/F3_E_BD_27230/' *.fastq.gz
rename 's/M3/F3_E_BD_27369/' *.fastq.gz
rename 's/M4/F3_E_BD_27415/' *.fastq.gz
rename 's/M5/F3_L_BD_27286/' *.fastq.gz
rename 's/M6/F3_L_BD_27303/' *.fastq.gz
rename 's/M7/F3_L_BD_27344/' *.fastq.gz
rename 's/M8/F3_L_BD_27348/' *.fastq.gz
rename 's/P5/F1_E_Pool5/' *.fastq.gz
rename 's/P6/F1_L_Pool6/' *.fastq.gz
rename 's/P7/F3_E_Pool7/' *.fastq.gz
rename 's/P8/F3_L_Pool8/' *.fastq.gz
conda deactivate


# make env file .yml
conda activate BS_pipeline
conda env export > src/env_BS_pipeline.yml
conda deactivate

conda activate QC-NGS
conda env export > src/env_QC-NGS.yml
conda deactivate

# run snakemake
snakemake -n --use-conda
snakemake -j 32 --use-conda
snakemake -j 8 --use-conda

snakemake -n --use-conda out/merged/alignment/F3_E_BD_27272.deduplicated_merged.bam
snakemake -j 2 --use-conda out/merged/alignment/F3_E_BD_27272.deduplicated_merged.bam
snakemake -n --use-conda out/merged/alignment/F3_E_BD_27438.deduplicated_merged.bam
snakemake -j 2 --use-conda out/merged/alignment/F3_E_BD_27438.deduplicated_merged.bam

# get dag
snakemake --dag | sed -e '1,4d' | dot -Tsvg > dag.svg

## run fastqc and multi-qc
conda activate BS_pipeline
fastqc -o trimmed_data/QC/reseq2020_SP/reports/trimmed/ -t 24 trimmed_data/reseq2020_SP/*.fq.gz
fastqc -o trimmed_data/QC/reseq2020/reports/trimmed/ -t 10 trimmed_data/reseq2020/*.fq.gz
fastqc -o trimmed_data/QC/seq2018/reports/trimmed/ -t 10 trimmed_data/seq2018/*.fq
cd
multiqc -n report_all -o trimmed_data/QC/reseq2020_SP trimmed_data/QC/reseq2020_SP/reports/trimmed/
multiqc -n trimmed_data/QC/reseq2020/report_all -o trimmed_data/QC/reseq2020/reports/
multiqc -n trimmed_data/QC/seq2018/report_all -o trimmed_data/QC/seq2018/reports/
conda deactivate


# gz trimmed files (seq2018)
for f in *.fq; do pigz -p20 $f; done


### make dictionary for read groups
# creat symlink for files per lane (in directory of raw data)
for f in *R1_001.fastq.gz; do ln -s /mnt/nfs/bioinfdata/primary/AnE/Parus_major/raw_data/whole_genome_sequencing/bisulfite_sequencing/novoseq/selection_lines2018_redo18062020/data/$f /home/NIOO.INT/melaniel/projects/WGBS_Snakemake_Bismark/raw_data/reseq2020_full/$f; done
for f in *R1_001.fastq; do ln -s /mnt/nfs/bioinfdata/primary/AnE/Parus_major/raw_data/whole_genome_sequencing/bisulfite_sequencing/novoseq/selection_lines2018/$f /home/NIOO.INT/melaniel/projects/WGBS_Snakemake_Bismark/raw_data/seq2018_full/$f; done
for f in *fastq.gz; do ln -s /mnt/nfs/bioinfdata/primary/AnE/Parus_major/raw_data/whole_genome_sequencing/bisulfite_sequencing/novoseq/selection_lines2018_redo10112020/data/$f /home/NIOO.INT/melaniel/projects/WGBS_Snakemake_Bismark/raw_data/reseq2020_SP_full/$f; done
conda activate rename
# fix library names for third run
rename 's/M1/F3_E_BD_27211/' *.fastq.gz
rename 's/M2/F3_E_BD_27230/' *.fastq.gz
rename 's/M3/F3_E_BD_27369/' *.fastq.gz
rename 's/M4/F3_E_BD_27415/' *.fastq.gz
rename 's/M5/F3_L_BD_27286/' *.fastq.gz
rename 's/M6/F3_L_BD_27303/' *.fastq.gz
rename 's/M7/F3_L_BD_27344/' *.fastq.gz
rename 's/M8/F3_L_BD_27348/' *.fastq.gz
rename 's/P5/F1_E_Pool5/' *.fastq.gz
rename 's/P6/F1_L_Pool6/' *.fastq.gz
rename 's/P7/F3_E_Pool7/' *.fastq.gz
rename 's/P8/F3_L_Pool8/' *.fastq.gz
conda deactivate
