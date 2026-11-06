qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path demux-paired.qza

#Inspect read quality

qiime demux summarize \
  --i-data demux-paired.qza \
  --o-visualization demux-summary.qzv

#Denoise the paired-end reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza# Export your representative sequences as FASTA


 qiime feature-table summarize \
  --i-table table.qza \
  --o-feature-frequencies feature-frequencies.qza \
  --o-sample-frequencies sample-frequencies.qza \
  --o-summary visual-summary.qzv

qiime feature-classifier classify-sklearn \
  --i-classifier SILVA_144_SSURef_NR99_uniform_classifier_full-length.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza



qiime tools export \
  --input-path rep-seqs.qza \
  --output-path picrust2_input/sequences



biom convert \
  -i picrust2_input/table/feature-table.biom \
  -o picrust2_input/feature-table.tsv \
  --to-tsv

