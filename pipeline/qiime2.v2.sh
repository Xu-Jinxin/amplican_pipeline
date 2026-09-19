#!/bin/bash
export PYTHONWARNINGS="ignore::UserWarning"

conda activate qiime2-2026.1

# 去除引物
forward=GTGYCAGCMGCCGCGGTAA
reverse=GGACTACNVGGGTWTCTAAT
database="gtdb"  # 可以选择 "silva" 或 "greengenes"

while getopts "f:r:d:" opt; do
  case $opt in
    f) forward="$OPTARG" ;;
    r) reverse="$OPTARG" ;;
    d) database="$OPTARG" ;;
  esac
done

# 导入原始数据
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path paired-end-demux.qza

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f ${forward} \
  --p-front-r ${reverse} \
  --o-trimmed-sequences demux-noprimers.qza \
  --p-cores 6 \
  --quiet

# 去噪与特征序列提取
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-noprimers.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-trunc-q 20 \
  --p-n-threads 6 \
  --o-table table.qza \
  --o-representative-sequences representative-sequences.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza

# 物种鉴定
if [ $database == "silva" ]; then
    echo "Using SILVA database for classification."
    database_path=/home/mahoon/database/qiime2/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza
elif [ $database == "greengenes" ]; then
    echo "Using Greengenes database for classification."
    database_path=/home/mahoon/database/qiime2/2024.09.backbone.full-length.nb.qza
elif [  $database == "gtdb" ]; then
    echo "Using GTDB database for classification."
    database_path=/share/data01/project/xujinxin/database/qiime2/gtdb_classifier_r232/ssu_v4_515F_806R_classifier.qza
fi

# 物种分类
qiime feature-classifier classify-sklearn \
  --i-classifier $database_path \
  --p-n-jobs 6 \
  --i-reads representative-sequences.qza \
  --o-classification taxonomy.qza

# 过滤掉叶绿体和线粒体
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude "Chloroplast,Mitochondria" \
  --o-filtered-table table-filtered.qza

# 转置特征表; 合并 ASVs 为物种水平；导出 ASVs table 表
qiime feature-table transpose \
  --i-table table-filtered.qza \
  --o-transposed-feature-table transposed-table-filtered.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --m-input-file transposed-table-filtered.qza \
  --o-visualization taxo-table.qzv

qiime tools export \
  --input-path taxo-table.qzv \
  --output-path taxo-table-exported

cut -f 1-3 taxo-table-exported/metadata.tsv | sed "2d" | \
  awk 'BEGIN{FS=OFS="\t"}NR==1{print "ASVID",$0}NR>1{print "ASV" NR-1, $0}'  > seq2asvid.tsv
cut -f 4- taxo-table-exported/metadata.tsv | sed "2d" | \
  awk 'BEGIN{FS=OFS="\t"}NR==1{print "Taxon",$0}NR>1{print "ASV" NR-1, $0}'  > asv_table.tsv
rm -rf filtered-exported taxo-table-exported taxo-table.qzv

# 生成界、门、纲、目、科、属、种特征表
declare -A levels=(
  [2]="phylum"
  [3]="class"
  [4]="order"
  [5]="family"
  [6]="genus"
  [7]="species"
)

for level in "${!levels[@]}"; do
  name="${levels[$level]}"
  # 按照不同分类水平（门、纲、目、科、属、种）合并特征表
  qiime taxa collapse \
    --i-table table-filtered.qza \
    --i-taxonomy taxonomy.qza \
    --p-level ${level} \
    --o-collapsed-table collapsed-${name}.qza

  qiime tools export \
    --input-path collapsed-${name}.qza \
    --output-path table-filtered-exported

  biom convert --to-tsv \
    -i table-filtered-exported/feature-table.biom \
    -o ${name}_table.tsv

  sed -i "1d" ${name}_table.tsv
  sed -i "1s/#OTU ID/Taxon/" ${name}_table.tsv

  rm -rf table-filtered-exported
done

# 构建系统发育树
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences representative-sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# 导出系统发育树
qiime tools export \
  --input-path rooted-tree.qza \
  --output-path ./

awk -F'\t' 'NR>1{print "s/\\<" $2 "\\>/" $1 "/g"}' seq2asvid.tsv | \
  sed -f - tree.nwk > tree.nwk.temp; mv tree.nwk.temp tree.nwk

mkdir -p 02_qiime2_artifacts; mv *.qza manifest.tsv 02_qiime2_artifacts
mkdir -p 03_exported_tables; mv *.tsv tree.nwk 03_exported_tables
