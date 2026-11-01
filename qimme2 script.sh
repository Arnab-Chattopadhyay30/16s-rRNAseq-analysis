mkdir emp-paired-end-sequences

# Import data

qiime tools import \
  --type EMPPairedEndSequences \
  --input-path emp-paired-end-sequences \
  --output-path emp-paired-end-sequences.qza

unzip -q casava-18-paired-end-demultiplexed.zip

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path casava-18-paired-end-demultiplexed \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

unzip -q pe-64.zip


#!/bin/bash

set -euo pipefail

# ==========================================
# GUT MICROBIOME AUTOMATED ANALYSIS PIPELINE
# ==========================================

# Select analysis:
# amplicon = QIIME 2 + DADA2 + SILVA
# shotgun  = MEGAHIT + SmORFinder

MODE="shotgun"

# Directories
RAW="/mnt/d/gut_microbiome"
CLEAN="$RAW/clean_reads"
RESULTS="$RAW/results"
ASSEMBLY="$RESULTS/assembly"
SMORF="$RESULTS/smorf"

# Threads
THREADS=4

# QIIME classifier
CLASSIFIER="$RAW/SILVA_144_SSURef_NR99_uniform_classifier_full-length.qza"

mkdir -p "$CLEAN" "$RESULTS" "$ASSEMBLY" "$SMORF"

echo "=================================="
echo "GUT MICROBIOME PIPELINE STARTED"
echo "Analysis mode: $MODE"
echo "=================================="

# ==========================================
# STEP 1: QUALITY CONTROL USING FASTP
# ==========================================

shopt -s nullglob

FILES=("$RAW"/*_1.fastq "$RAW"/*_1.fastq.gz)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: No forward FASTQ files found!"
    exit 1
fi

for read1 in "${FILES[@]}"
do

    filename=$(basename "$read1")
    sample="${filename%_1.fastq.gz}"
    sample="${sample%_1.fastq}"

    read2="${read1/_1.fastq/_2.fastq}"

    if [[ ! -f "$read2" ]]; then
        echo "Missing reverse read: $sample"
        exit 1
    fi

    echo "FASTP: Processing $sample"

    fastp \
        -i "$read1" \
        -I "$read2" \
        -o "$CLEAN/${sample}_1.fastq.gz" \
        -O "$CLEAN/${sample}_2.fastq.gz" \
        --detect_adapter_for_pe \
        -q 20 \
        -l 50 \
        -w "$THREADS" \
        -h "$RESULTS/${sample}_fastp.html" \
        -j "$RESULTS/${sample}_fastp.json"

done

echo "FASTP QUALITY CONTROL COMPLETED"

# ==========================================
# STEP 2A: AMPLICON ANALYSIS USING QIIME 2
# ==========================================

if [[ "$MODE" == "amplicon" ]]; then

    cd "$RESULTS"

    echo "Creating QIIME 2 manifest..."

    printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > manifest.tsv

    for read1 in "$CLEAN"/*_1.fastq.gz
    do
        filename=$(basename "$read1")
        sample="${filename%_1.fastq.gz}"
        read2="$CLEAN/${sample}_2.fastq.gz"

        printf "%s\t%s\t%s\n" \
            "$sample" "$read1" "$read2" >> manifest.tsv
    done

    echo "Importing reads into QIIME 2..."

    conda run -n qiime2-amplicon-2026.10 qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path manifest.tsv \
        --input-format PairedEndFastqManifestPhred33V2 \
        --output-path demux-paired.qza

    echo "Summarizing read quality..."

    conda run -n qiime2-amplicon-2026.10 qiime demux summarize \
        --i-data demux-paired.qza \
        --o-visualization demux-summary.qzv

    echo "Running DADA2..."

    conda run -n qiime2-amplicon-2026.10 qiime dada2 denoise-paired \
        --i-demultiplexed-seqs demux-paired.qza \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-n-threads "$THREADS" \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza

    echo "Summarizing feature table..."

    conda run -n qiime2-amplicon-2026.10 qiime feature-table summarize \
        --i-table table.qza \
        --o-visualization visual-summary.qzv

    echo "Classifying taxonomy..."

    conda run -n qiime2-amplicon-2026.10 qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER" \
        --i-reads rep-seqs.qza \
        --o-classification taxonomy.qza

    echo "Exporting representative sequences..."

    conda run -n qiime2-amplicon-2026.10 qiime tools export \
        --input-path rep-seqs.qza \
        --output-path exported-rep-seqs

    echo "AMPLICON ANALYSIS COMPLETED"

fi

# ==========================================
# STEP 2B: SHOTGUN METAGENOMIC ANALYSIS
# ==========================================

if [[ "$MODE" == "shotgun" ]]; then

    for read1 in "$CLEAN"/*_1.fastq.gz
    do

        filename=$(basename "$read1")
        sample="${filename%_1.fastq.gz}"

        read2="$CLEAN/${sample}_2.fastq.gz"

        echo "=================================="
        echo "MEGAHIT ASSEMBLY: $sample"
        echo "=================================="

        megahit \
            -1 "$read1" \
            -2 "$read2" \
            -o "$ASSEMBLY/${sample}_output" \
            -t "$THREADS"

        echo "SMORF PREDICTION: $sample"

        conda run -n smorfinder smorf meta \
            "$ASSEMBLY/${sample}_output/final.contigs.fa" \
            -o "$SMORF/${sample}_smorf_output"

    done

    echo "SHOTGUN METAGENOMIC ANALYSIS COMPLETED"

fi

echo "=================================="
echo "ALL PIPELINE STEPS COMPLETED"
echo "=================================="
