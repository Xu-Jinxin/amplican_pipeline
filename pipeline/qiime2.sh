#!/bin/bash

conda activate qiime2-2024.5
export PYTHONWARNINGS="ignore::UserWarning"

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

# 去除引物
forward=ACTCCTACGGGAGGCAGCA
reverse=GGACTACHVGGGTWTCTAAT

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f ${forward} \
  --p-front-r ${reverse} \
  --o-trimmed-sequences demux-noprimers.qza \
  --p-cores 4 \
  --quiet

# 合并双端序列
qiime vsearch merge-pairs \
  --i-demultiplexed-seqs demux-noprimers.qza \
  --o-merged-sequences demux-merged.qza \
  --o-unmerged-sequences demux-unmerged.qza \
  --p-threads 4 \
  --quiet

# 过滤低质量序列
qiime quality-filter q-score \
  --i-demux demux-merged.qza \
  --o-filtered-sequences demux-merged-filtered.qza \
  --o-filter-stats demux-merged-filter-stats.qza

# 去噪与特征序列提取
qiime dada2 denoise-paired \
 --i-demultiplexed-seqs demux-noprimers.qza \
 --p-trunc-len-f 0 \
 --p-trunc-len-r 0 \
 --p-trunc-q 20 \
 --p-n-threads 12 \
 --o-table table.qza \
 --o-representative-sequences representative-sequences.qza \
 --o-denoising-stats denoising-stats.qza 

# 物种鉴定
database="greengenes"  # 可以选择 "silva" 或 "greengenes"

if [ $database == "silva" ]; then
    echo "Using SILVA database for classification."
    database_path=/home/mahoon/database/qiime2/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza
elif [ $database == "greengenes" ]; then
    echo "Using Greengenes database for classification."
    database_path=/home/mahoon/database/qiime2/2024.09.backbone.full-length.nb.qza
fi

# 物种分类
qiime feature-classifier classify-sklearn \
  --i-classifier $database_path \
  --p-n-jobs 12 \
  --i-reads representative-sequences.qza \
  --o-classification taxonomy.qza

# 过滤掉叶绿体和线粒体
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude "Chloroplast,Mitochondria" \
  --o-filtered-table table-filtered.qza

# 查找特征表中每个样本的总序列数
qiime tools export \
  --input-path table-filtered.qza \
  --output-path filtered-exported

# 最小样本量设置为每个样本序列数的95%，以保留大多数数据同时去除极端值
min_depth=$(biom summarize-table -i filtered-exported/feature-table.biom | \
  grep "Counts/sample detail" -A 1 | tail -1 | sed "s/,//" | awk '{print int($2*0.95)}')

# 稀释特征表以标准化样本间的测序深度
qiime feature-table rarefy \
  --i-table table-filtered.qza \
  --p-sampling-depth ${min_depth} \
  --o-rarefied-table table-filtered-rarefied.qza

# 转置特征表; 合并 ASVs 为物种水平；导出 ASVs table 表
qiime feature-table transpose \
  --i-table table-filtered-rarefied.qza \
  --o-transposed-feature-table transposed-table-filtered-rarefied.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --m-input-file transposed-table-filtered-rarefied.qza \
  --o-visualization taxo-table-rarefied.qzv

qiime tools export \
  --input-path taxo-table-rarefied.qzv \
  --output-path taxo-table-rarefied-exported

cut -f 1-3 taxo-table-rarefied-exported/metadata.tsv | sed "2d" | \
 awk 'BEGIN{FS=OFS="\t"}NR==1{print "ASVID",$0}NR>1{print "ASV" NR-1, $0}'  > seq2asvid.tsv

cut -f 4- taxo-table-rarefied-exported/metadata.tsv | sed "2d" | \
 awk 'BEGIN{FS=OFS="\t"}NR==1{print "Taxon",$0}NR>1{print "ASV" NR-1, $0}'  > asv_rarefied_table.tsv

rm -rf filtered-exported taxo-table-rarefied-exported taxo-table-rarefied.qzv

# 按照不同分类水平（门、纲、目、科、属、种）合并特征表
declare -A levels=(
  [2]="phylum"
  [3]="class"
  [4]="order"
  [5]="family"
  [6]="genus"
  [7]="species"
)

# 循环处理每个分类水平
for level in "${!levels[@]}"; do
  name="${levels[$level]}"
  
  qiime taxa collapse \
    --i-table table-filtered-rarefied.qza \
    --i-taxonomy taxonomy.qza \
    --p-level ${level} \
    --o-collapsed-table collapsed-${name}.qza

  qiime tools export \
  --input-path collapsed-${name}.qza \
  --output-path table-filtered-rarefied-exported

  biom convert --to-tsv \
    -i table-filtered-rarefied-exported/feature-table.biom \
    -o ${name}_rarefied_table.tsv

  sed -i "1d" ${name}_rarefied_table.tsv
  sed -i "1s/#OTU ID/Taxon/" ${name}_rarefied_table.tsv

  rm -rf table-filtered-rarefied-exported
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
