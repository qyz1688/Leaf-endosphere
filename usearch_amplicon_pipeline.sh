#!/usr/bin/env bash
set -euo pipefail

# USEARCH-based amplicon pipeline for bacterial 16S and fungal ITS data.
# Usage:
#   bash usearch_amplicon_pipeline.sh bacteria
#   bash usearch_amplicon_pipeline.sh fungi

MARKER=${1:-bacteria}   # bacteria or fungi
USEARCH=usearch
THREADS=30

SEQ_DIR=seq
META=metadata.txt

# Modify these database paths before running.
RDP_DB=/path/to/rdp_16s_v18.fa
UNITE_DB=/path/to/unite.fa

if [[ "${MARKER}" == "bacteria" ]]; then
    TAX_DB=${RDP_DB}
    STRIP_LEFT=25
    STRIP_RIGHT=24
    RARE_DEPTH=27000
elif [[ "${MARKER}" == "fungi" ]]; then
    TAX_DB=${UNITE_DB}
    STRIP_LEFT=22
    STRIP_RIGHT=20
    RARE_DEPTH=4000
else
    echo "ERROR: MARKER must be bacteria or fungi"
    exit 1
fi

MAXEE_RATE=0.01
MIN_UNIQUE_SIZE=8
MIN_ASV_SIZE=10
MAP_ID=0.97
SINTAX_CUTOFF=0.6

mkdir -p temp result/raw result/alpha result/beta result/log

echo "Step 1. Merge paired-end reads"
while read -r SAMPLE REST; do
    [[ "${SAMPLE}" == "SampleID" ]] && continue
    ${USEARCH} -fastq_mergepairs ${SEQ_DIR}/${SAMPLE}_R1.raw.fastq \
        -reverse ${SEQ_DIR}/${SAMPLE}_R2.raw.fastq \
        -fastqout temp/${SAMPLE}.merged.fq \
        -relabel ${SAMPLE}. \
        > result/log/${SAMPLE}.merge.log 2>&1
done < ${META}

cat temp/*.merged.fq > temp/all.merged.fq

echo "Step 2. Primer trimming and quality filtering"
${USEARCH} -fastq_filter temp/all.merged.fq \
    -fastq_stripleft ${STRIP_LEFT} \
    -fastq_stripright ${STRIP_RIGHT} \
    -fastq_maxee_rate ${MAXEE_RATE} \
    -fastaout temp/filtered.fa \
    > result/log/quality_filter.log 2>&1

echo "Step 3. Dereplication"
${USEARCH} -fastx_uniques temp/filtered.fa \
    -fastaout temp/uniques.fa \
    -sizeout \
    -relabel Uni_ \
    -minuniquesize ${MIN_UNIQUE_SIZE} \
    > result/log/dereplication.log 2>&1

echo "Step 4. ASV inference with UNOISE3"
${USEARCH} -unoise3 temp/uniques.fa \
    -zotus temp/zotus.fa \
    -minsize ${MIN_ASV_SIZE} \
    > result/log/unoise3.log 2>&1

sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/asvs.fa

echo "Step 5. Reference-based chimera checking"
${USEARCH} -uchime_ref temp/asvs.fa \
    -db ${TAX_DB} \
    -nonchimeras result/raw/asvs.fa \
    > result/log/chimera_ref.log 2>&1

echo "Step 6. Construct ASV table"
${USEARCH} -otutab temp/filtered.fa \
    -otus result/raw/asvs.fa \
    -otutabout result/raw/asv_table_raw.txt \
    -id ${MAP_ID} \
    -threads ${THREADS} \
    > result/log/asv_table.log 2>&1

echo "Step 7. Taxonomic classification with SINTAX"
${USEARCH} -sintax result/raw/asvs.fa \
    -db ${TAX_DB} \
    -sintax_cutoff ${SINTAX_CUTOFF} \
    -tabbedout result/raw/asvs.sintax \
    > result/log/sintax.log 2>&1

echo "Step 8. Filter non-target sequences"
python3 filter_non_target_asvs.py ${MARKER} \
    result/raw/asv_table_raw.txt \
    result/raw/asvs.sintax \
    result/asv_table.txt \
    result/asv_keep_ids.txt \
    result/raw/asv_discarded.txt

${USEARCH} -fastx_getseqs result/raw/asvs.fa \
    -labels result/asv_keep_ids.txt \
    -fastaout result/asvs.fa \
    > result/log/get_filtered_asvs.log 2>&1

echo "Step 9. Rarefy ASV table"
${USEARCH} -otutab_subsample result/asv_table.txt \
    -sample_size ${RARE_DEPTH} \
    -output result/asv_table_rare.txt \
    > result/log/rarefaction.log 2>&1

echo "Step 10. Alpha- and beta-diversity"
${USEARCH} -alpha_div result/asv_table_rare.txt \
    -output result/alpha/alpha_diversity.txt \
    > result/log/alpha_diversity.log 2>&1

${USEARCH} -alpha_div_rare result/asv_table_rare.txt \
    -output result/alpha/alpha_rarefaction.txt \
    -method without_replacement \
    > result/log/alpha_rarefaction.log 2>&1

${USEARCH} -cluster_agg result/asvs.fa \
    -treeout result/asvs.tree \
    > result/log/tree.log 2>&1 || true

${USEARCH} -beta_div result/asv_table_rare.txt \
    -filename_prefix result/beta/ \
    > result/log/beta_diversity.log 2>&1

echo "Step 11. Format taxonomy table"
cut -f 1,4 result/raw/asvs.sintax \
    | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
    > result/taxonomy2.txt || true

awk 'BEGIN{OFS=FS="\t"}{
    delete a;
    a["k"]="Unassigned"; a["p"]="Unassigned"; a["c"]="Unassigned";
    a["o"]="Unassigned"; a["f"]="Unassigned"; a["g"]="Unassigned"; a["s"]="Unassigned";
    split($2,x,";");
    for(i in x){split(x[i],b,"__"); if(length(b[1])>0 && length(b[2])>0) a[b[1]]=b[2];}
    print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];
}' result/taxonomy2.txt \
    | sed '1 s/^/ASVID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
    > result/taxonomy.txt || true

echo "Pipeline finished."
