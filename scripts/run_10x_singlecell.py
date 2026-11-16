#!/usr/bin/env python3
"""Run source-derived calling and filtering stages for one sample."""
PLATFORM='10x'
LEVEL='single-cell'
FILTERS={'samtools': {'10x_pooled': 'FMT/GT!="0/0" && QUAL>20 && ( (INFO/DP4[2]+INFO/DP4[3]) / (INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3]) ) > 0.05 && ( (TYPE="snp" && INFO/DP>10) || (TYPE!="snp" && INFO/DP>15) )', '10x_celltype': 'FMT/GT!="0/0" && QUAL>20 && ( (INFO/DP4[2]+INFO/DP4[3]) / (INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3]) ) > 0.05 && ((TYPE="snp" && INFO/DP>10) || (TYPE!="snp" && INFO/DP>15))', 'smartseq': 'FMT/GT!="0/0" && QUAL>20 && (INFO/DP4[2]+INFO/DP4[3])/(INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3])>0.05 && ( (TYPE="snp" && INFO/DP>10) || (TYPE!="snp" && INFO/DP>15) )'}, 'freebayes': {'10x_pooled': 'FMT/GT!="0/0" && QUAL>20 && INFO/AF[0]>=0.05 && ( (TYPE="snp" && INFO/DP[0]>=10) || (TYPE!="snp" && INFO/DP[0]>=15) )', '10x_celltype': 'FMT/GT!="0/0" && QUAL>20 && INFO/AF[0]>=0.05 && ( (TYPE="snp" && INFO/DP[0]>10) || (TYPE!="snp" && INFO/DP[0]>15) )', 'smartseq': 'FMT/GT!="0/0" && QUAL>20 && INFO/AF[0]>0.05 && ( (TYPE="snp" && INFO/DP[0]>10) || (TYPE!="snp" && INFO/DP[0]>15) )'}, 'varscan2': {'10x_pooled': 'FMT/GT!="0/0" && FMT/DP>10 && (FMT/AD / (FMT/RD + FMT/AD)) > 0.05 && ( (TYPE="snp" && FMT/DP>10) || (TYPE!="snp" && FMT/DP>15) )', '10x_celltype': 'FMT/GT!="0/0" && (FMT/AD / (FMT/RD + FMT/AD)) > 0.05 && ( (TYPE="snp" && FMT/DP>10) || (TYPE!="snp" && FMT/DP>15) )', 'smartseq': 'FMT/GT!="0/0" && FMT/DP>10 && (FMT/AD / (FMT/RD + FMT/AD)) > 0.05'}, 'strelka2': {'10x_pooled': 'ALT!="." && FMT/GT!="0/0" && FMT/GT!="0|0" && QUAL>20 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP > 10) || (TYPE!="snp" && FMT/DPI > 15) )', '10x_celltype': 'ALT!="." && FMT/GT!="0/0" && FMT/GT!="0|0" && QUAL>20 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP > 10) || (TYPE!="snp" && FMT/DPI > 15) )', 'smartseq': 'ALT!="." && FMT/GT!="0/0" && FMT/GT!="0|0" && QUAL>20 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP > 10) || (TYPE!="snp" && FMT/DPI > 15) )'}, 'monovar': {'10x_pooled': 'FILTER="PASS" && FMT/GT!="0/0" && FMT/GT!="./." && (FMT/AD[0:0] + FMT/AD[0:1]) > 0 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP>10) || (TYPE!="snp" && FMT/DP>15) )', '10x_celltype': 'FILTER="PASS" && FMT/GT!="0/0" && FMT/GT!="./." && (FMT/AD[0:0] + FMT/AD[0:1]) > 0 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP>10) || (TYPE!="snp" && FMT/DP>15) )', 'smartseq': 'FMT/GT!="0/0" && FMT/GT!="./." && (FMT/AD[0:0] + FMT/AD[0:1]) > 0 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP>10) || (TYPE!="snp" && FMT/DP>15) )'}, 'cellsnp': {'smartseq': 'FMT/GT!="0/0" && FMT/GT!="0|0" && FMT/GT!="." && FMT/GT!="./." && FMT/GT!=".|." && FMT/DP>10 && FMT/AD/FMT/DP>0.05'}}
SOURCES={}
SOURCES['scsnv_call']=r'''
###===== .sh file =====##
#!/bin/bash
export PATH=/home/tantingting/scmut/tools:$PATH

samples=( "P2_cSCC" "P4_cSCC" "P5_cSCC" "P6_cSCC" "P7_cSCC" "P9_cSCC" "P10_cSCC")

ref_dir="/home/tantingting/scmut/ref"
bwa_ref="${ref_dir}/bwa_ref"
anno_dir="${ref_dir}/scsnv_anno"

scsnv_bin="/home/tantingting/scmut/tools/scsnv/build/src/scsnv"
fastq_root="/mnt/storage/tantingting/scmut/fastp/1_cSCC_scRNA"
out_root="/mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/scsnv"

barcode_file="/home/tantingting/scmut/tools/scsnv/data/737K-august-2016.txt"

for sample in "${samples[@]}"; do
    echo ">>> Processing $sample"

    fastq_dir="${fastq_root}/${sample}"
    out_dir="${out_root}/${sample}"
    #mkdir -p "$out_dir"

    ## 1. count
    memusg -t "$scsnv_bin" count -l V2 -o "$out_dir/barcode" -k "$barcode_file" "$fastq_dir"
    echo "$sample: 1.count done"

    ## 2. map
    memusg -t "$scsnv_bin" map -l V2 -i "$bwa_ref/scsnv_ref" -g "$bwa_ref/genome.fa" -b "$out_dir/barcode" \
        -t 16 --bam-write 4 -q 4 -o "$out_dir/" "$fastq_dir"
    echo "$sample: 2.map done"

    ## 3. collapse
    memusg -t "$scsnv_bin" collapse -l V2 -r "$bwa_ref/genome.fa" -i "$bwa_ref/scsnv_ref" \
        -b "$out_dir/barcode_counts.txt.gz" --threads 16 --bam-write 8 -o "$out_dir/" "$out_dir/merged.bam"
    echo "$sample: 3.collapse done"

    ## 4. cells
    memusg -t scsnvmisc cells --skip-mt -o "$out_dir" "$out_dir/summary.h5"
    echo "$sample: 4.cells done"

    ## 5. pileup
    memusg -t "$scsnv_bin" pileup -l V2 -t 8 -x 8 -i "$bwa_ref/scsnv_ref" -r "$bwa_ref/genome.fa" \
         -p "$out_dir/passed_barcodes.txt.gz" -o "$out_dir/pileup" "$out_dir/collapsed.bam"
    echo "$sample: 5.pileup done"

    ## 6. annotate
    memusg -t scsnvmisc annotate -d "$anno_dir/1000Genome_change.txt" -e "$anno_dir/TABLE1_hg38_v2.txt" \
        "$out_dir/pileup"
    echo "$sample: 6.annotate done"

    memusg -t scsnv2mtx "$out_dir/pileup_annotated.h5" "$out_dir"
    echo "$sample: 7.scsnv2mtx done"

    echo ">>> $sample finished."
    echo "----------------------------------"
done

echo "All samples done."
'''
SOURCES['scsnv_filter']=r'''
samples=("P4_cSCC" "P5_cSCC" "P6_cSCC" "P7_cSCC" "P9_cSCC" "P10_cSCC")
root="/mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/scsnv"

for s in "${samples[@]}"; do
    in_csv="${root}/${s}/snvs.csv"
    out_prefix="${root}/${s}/snv.dp10_af005"

    echo ">>> Processing $s"

    awk -F',' 'BEGIN{
        OFS=","
    }
    NR==1{
        print "chrom","pos","ref","alt","ref_count","alt_count","strand","strand_AF","total_barcodes"
        next
    }
    {
        chrom=$1; pos=$2; ref=$3; alt=$4; refc=$5+0; altc=$8+0; strand=$15; af=$20+0; barcodes=$22+0; dp=refc+altc
        if (chrom !~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/) next
        if (dp <= 10) next
        if (af <= 0.05) next
        print chrom,pos,ref,alt,refc,altc,strand,af,barcodes
    }' "$in_csv" > "${out_prefix}.filtered.csv"

    awk -F',' 'BEGIN{
        OFS="\t"
        print "##fileformat=VCFv4.2"
        print "##source=scsnv_snv.csv_filtered_by_awk"
        print "##contig=<ID=chr1>"
        print "##contig=<ID=chr2>"
        print "##contig=<ID=chr3>"
        print "##contig=<ID=chr4>"
        print "##contig=<ID=chr5>"
        print "##contig=<ID=chr6>"
        print "##contig=<ID=chr7>"
        print "##contig=<ID=chr8>"
        print "##contig=<ID=chr9>"
        print "##contig=<ID=chr10>"
        print "##contig=<ID=chr11>"
        print "##contig=<ID=chr12>"
        print "##contig=<ID=chr13>"
        print "##contig=<ID=chr14>"
        print "##contig=<ID=chr15>"
        print "##contig=<ID=chr16>"
        print "##contig=<ID=chr17>"
        print "##contig=<ID=chr18>"
        print "##contig=<ID=chr19>"
        print "##contig=<ID=chr20>"
        print "##contig=<ID=chr21>"
        print "##contig=<ID=chr22>"
        print "##contig=<ID=chrX>"
        print "##contig=<ID=chrY>"
        print "##contig=<ID=chrM>"
        print "##INFO=<ID=REF_COUNT,Number=1,Type=Integer,Description=\"Reference allele count from scSNV CSV\">"
        print "##INFO=<ID=ALT_COUNT,Number=1,Type=Integer,Description=\"Alternative allele count from scSNV CSV\">"
        print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth calculated as ref_count + alt_count\">"
        print "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Strand field from scSNV CSV\">"
        print "##INFO=<ID=STRAND_AF,Number=1,Type=Float,Description=\"strand_AF from scSNV CSV\">"
        print "##INFO=<ID=TOTAL_BARCODES,Number=1,Type=Integer,Description=\"Total barcodes supporting site\">"
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", "SAMPLE"
    }
    NR==1{next}
    {
        chrom=$1; pos=$2; ref=$3; alt=$4; refc=$5+0; altc=$6+0; strand=$7; af=$8+0; barcodes=$9+0; dp=refc+altc
        info="REF_COUNT=" refc ";ALT_COUNT=" altc ";DP=" dp ";STRAND=" strand ";STRAND_AF=" af ";TOTAL_BARCODES=" barcodes
        gt="0/1"
        print chrom,pos,".",ref,alt,".","PASS",info,"GT",gt
    }' "${out_prefix}.filtered.csv" > "${out_prefix}.vcf"

    bcftools sort "${out_prefix}.vcf" -o "${out_prefix}.sorted.vcf"
    mv "${out_prefix}.sorted.vcf" "${out_prefix}.vcf"
    bgzip -f -c "${out_prefix}.vcf" > "${out_prefix}.vcf.gz"
    tabix -f -p vcf "${out_prefix}.vcf.gz"

    echo "$s done"
done

# Convert zero-based positions to one-based coordinates using pos + 1.
BASE="/mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/scsnv"
samples=(P4 P5 P6 P7 P9 P10)

for s in "${samples[@]}"; do
    sample="${s}_cSCC"
    in_vcf="${BASE}/${sample}/snv.dp10_af005.vcf"
    out_vcf="${BASE}/${sample}/snv.dp10_af005.pos1.vcf"

    echo ">>> Processing ${sample}" 

    awk 'BEGIN{OFS="\t"} 
         /^#/ {print; next} 
         {$2 = $2 + 1; print}' "$in_vcf" > "$out_vcf"

        # 2. Compress the VCF with bgzip.
    bgzip -f "$out_vcf"

    # 3. Create a tabix index.
    tabix -p vcf  "${out_vcf}.gz"

    echo "[OK] Output: $out_vcf" 

done

echo "=== Done ===" 

# Calculate the updated intersection.

'''
SOURCES['scsnv_pairs']=r'''
###===== python =====##
#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
from scipy.io import mmread

def main():
    parser = argparse.ArgumentParser(description="Extract filtered scSNV cell-mutation pairs and per-cell mutation counts.")
    parser.add_argument("--sample", required=True, help="Sample name, e.g. P5_cSCC")
    parser.add_argument("--indir", required=True, help="scSNV output directory")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--filtered_csv", default="snv.dp10_af005.filtered.csv", help="Filtered SNV csv file name")
    args = parser.parse_args()

    sample = args.sample
    indir = args.indir
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    barcode_file = os.path.join(indir, "barcodes.txt")
    snv_file = os.path.join(indir, "snvs.csv")
    alt_file = os.path.join(indir, "alts.mtx")
    filtered_file = os.path.join(indir, args.filtered_csv)
    print(f"[INFO] sample: {sample}")
    print(f"[INFO] input dir: {indir}")

    # 1. Read barcodes and remove optional -1 suffixes.
    barcodes = pd.read_csv(barcode_file, header=None, names=["cell"])
    barcodes["cell"] = barcodes["cell"].astype(str).str.replace("-1$", "", regex=True)
    barcodes["cell_index0"] = np.arange(barcodes.shape[0])
    print(f"[INFO] n_cells: {barcodes.shape[0]}")

    # 2. Read all SNV records; their row order matches alts.mtx.
    snvs = pd.read_csv(snv_file)
    snvs["snv_index0"] = np.arange(snvs.shape[0])
    snvs["mut"] = snvs["chrom"].astype(str) + ":" + snvs["pos"].astype(str) + ":" + snvs["ref"].astype(str) + ":" + snvs["alt"].astype(str)
    print(f"[INFO] n_all_snvs: {snvs.shape[0]}")

    # 3. Read the retained SNV list for subsequent variant-key filtering.
    filtered = pd.read_csv(filtered_file)
    filtered["mut"] = filtered["chrom"].astype(str) + ":" + filtered["pos"].astype(str) + ":" + filtered["ref"].astype(str) + ":" + filtered["alt"].astype(str)
    keep_mut = set(filtered["mut"]) # Convert retained variant keys to a set.
    print(f"[INFO] n_filtered_snvs: {len(keep_mut)}")

    # 4. Read alts.mtx and retain positive alternate-read cell-SNV records.
    print("[INFO] reading alts.mtx ...")
    alt_mtx = mmread(alt_file).tocoo()
    print(f"[INFO] alt matrix shape: {alt_mtx.shape}") # Check that the matrix and annotation dimensions agree.
    if alt_mtx.shape[0] != snvs.shape[0]:
        raise ValueError("Number of SNVs in alts.mtx does not match snvs.csv.")
    if alt_mtx.shape[1] != barcodes.shape[0]:
        raise ValueError("Number of cells in alts.mtx does not match barcodes.txt.")
    df = pd.DataFrame({"snv_index0": alt_mtx.row.astype(np.int64), "cell_index0": alt_mtx.col.astype(np.int64), "alt_count": alt_mtx.data.astype(np.int32)})# Store the SNV row index, cell column index and alternate-read count.
    df = df[df["alt_count"] > 0].copy()
    print(f"[INFO] cell-SNV records with alt_count > 0 before filtering: {df.shape[0]}")

    # 5. Join variant and cell identifiers using matrix indices.
    df = df.merge(barcodes[["cell_index0", "cell"]], on="cell_index0", how="left")
    df = df.merge(snvs[["snv_index0", "mut"]], on="snv_index0", how="left")

    # 6. Restrict records to variant keys retained in the filtered CSV.
    df = df[df["mut"].isin(keep_mut)].copy()
    print(f"[INFO] cell-SNV records after filtered SNV list: {df.shape[0]}")
    # 7. Write the two-column cell-mutation table.
    cell_mut = df[["cell", "mut"]].drop_duplicates()
    out_cell_mut = os.path.join(outdir, f"{sample}_scSNV_filtered_cell_mut.tsv.gz")
    cell_mut.to_csv(out_cell_mut, sep="\t", index=False, compression="gzip")
    print(f"[INFO] saved cell-mut table: {out_cell_mut}")
    # 8. Count retained variants for every cell, including zero-variant cells.
    cell_count = cell_mut.groupby("cell", as_index=False).agg(n_variants=("mut", "nunique"))
    all_cells = barcodes[["cell"]].copy()
    cell_count = all_cells.merge(cell_count, on="cell", how="left")
    cell_count["n_variants"] = cell_count["n_variants"].fillna(0).astype(int)
    cell_count.insert(0, "sample", sample)
    cell_count.insert(1, "tool", "scSNV")
    out_count = os.path.join(outdir, f"{sample}_scSNV_filtered_cell_variant_count.tsv")
    cell_count.to_csv(out_count, sep="\t", index=False)
    print(f"[INFO] saved cell count table: {out_count}")
    print("[DONE]")
if __name__ == "__main__":
    main()
'''
SOURCES['scsnv_valid']=r'''
import os
import pandas as pd

BASE_DIR = "/mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/scsnv"
BAM_DIR = "/mnt/storage/tantingting/scmut/bam/1_cSCC_scRNA"
samples = ["P4_cSCC", "P5_cSCC", "P6_cSCC", "P7_cSCC", "P9_cSCC", "P10_cSCC"]

all_count = []

for sample in samples:
    print(f"\n===== {sample} =====")
    cell_mut_file = f"{BASE_DIR}/{sample}/cell_snv/{sample}_scSNV_filtered_cell_mut.tsv.gz"
    reads_file = f"{BAM_DIR}/{sample}/outs/{sample}_barcode_reads_map.csv"
    cell_mut = pd.read_csv(cell_mut_file, sep="\t")
    reads_df = pd.read_csv(reads_file).rename(columns={"id": "reads"})
    # Normalize barcode suffixes before matching cell identifiers.
    cell_mut["cell"] = cell_mut["cell"].astype(str).str.replace("-1$", "", regex=True)
    reads_df["cell"] = reads_df["cell"].astype(str).str.replace("-1$", "", regex=True)
    # Cells in the read-count table define the valid Cell Ranger cell set.
    all_cells = reads_df[["cell", "reads"]].drop_duplicates()
    valid_cells = set(all_cells["cell"])
    print(f"[1] scSNV raw cell-mut pairs: {len(cell_mut)}")
    print(f"[1] scSNV raw SNVs: {cell_mut['mut'].nunique()}")
    print(f"[1] scSNV raw cells: {cell_mut['cell'].nunique()}")
    print(f"[2] CellRanger valid cells: {len(valid_cells)}")
    # Remove scSNV cells absent from the Cell Ranger read-count table before recurrence filtering.
    cell_mut_valid = cell_mut[cell_mut["cell"].isin(valid_cells)].drop_duplicates()
    removed_cells = set(cell_mut["cell"]) - valid_cells
    print(f"[3] removed scSNV-only cells: {len(removed_cells)}")
    print(f"[3] after valid-cell filter pairs: {len(cell_mut_valid)}")
    print(f"[3] after valid-cell filter SNVs: {cell_mut_valid['mut'].nunique()}")
    print(f"[3] after valid-cell filter cells: {cell_mut_valid['cell'].nunique()}")
    # Within valid cells, retain SNVs supported by at least two cells.
    snv_ncells = ( cell_mut_valid
        .groupby("mut")["cell"]
        .nunique()
        .reset_index(name="n_cells") )
    keep = set(snv_ncells.loc[snv_ncells["n_cells"] >= 2, "mut"])
    cell_mut_min2 = cell_mut_valid[cell_mut_valid["mut"].isin(keep)].drop_duplicates()
    print(f"[4] SNVs before min2cells: {cell_mut_valid['mut'].nunique()}")
    print(f"[4] SNVs after min2cells: {cell_mut_min2['mut'].nunique()}")
    print(f"[4] cell-mut pairs after min2cells: {len(cell_mut_min2)}")
    print(f"[4] cells with >=1 retained SNV: {cell_mut_min2['cell'].nunique()}")
    # Write the filtered cell-mutation records.
    out_long = f"{BASE_DIR}/{sample}/cell_snv/{sample}_scSNV_cell_mut_validcell_min2c.tsv.gz"
    cell_mut_min2.to_csv(out_long, sep="\t", index=False, compression="gzip")
    # Recalculate the retained variant count for every cell.
    cell_count = ( cell_mut_min2
        .groupby("cell", as_index=False)
        .agg(n_variants=("mut", "nunique")) )
    # Retain cells shared by scSNV and Cell Ranger.
    count_df = cell_count.merge(all_cells, on="cell", how="inner")
    count_df["n_variants"] = count_df["n_variants"].astype(int)
    count_df.insert(0, "sample", sample)
    count_df.insert(1, "tool", "scSNV")
    count_df["cell_id"] = count_df["sample"] + "_" + count_df["cell"].astype(str)
    print(f"[5] final output cells: {count_df['cell'].nunique()}")
    print(f"[5] final zero-variant cells: {(count_df['n_variants'] == 0).sum()}")
    print(f"[5] final cells with variants: {(count_df['n_variants'] > 0).sum()}")
    print(f"[5] reads NA count: {count_df['reads'].isna().sum()}")
    out_count = f"{BASE_DIR}/{sample}/cell_snv/{sample}_scSNV_cell_variant_count_validcell_min2c.tsv"
    count_df.to_csv(out_count, sep="\t", index=False)
    all_count.append(count_df)
    print(f"[OUT] long:  {out_long}")
    print(f"[OUT] count: {out_count}")

df_all = pd.concat(all_count, ignore_index=True)

out_all = f"{BASE_DIR}/scSNV_variant_count_vs_reads_input_validcell_min2c.tsv"
df_all.to_csv(out_all, sep="\t", index=False)

print("\n===== ALL DONE =====")
print(df_all.head())
print("total rows:", len(df_all))
print("unique cell_id:", df_all["cell_id"].nunique())
print("duplicated cell_id:", df_all["cell_id"].duplicated().sum())
print("reads NA:", df_all["reads"].isna().sum())
print("zero-variant cells:", (df_all["n_variants"] == 0).sum())
print("output:", out_all)

#===== P4_cSCC =====
#[1] scSNV raw cell-mut pairs: 3808239
#[1] scSNV raw SNVs: 131814
#[1] scSNV raw cells: 9095
# [2] CellRanger valid cells: 12795
# [3] removed scSNV-only cells: 38
# [3] after valid-cell filter pairs: 3805591
# [3] after valid-cell filter SNVs: 131814
# [3] after valid-cell filter cells: 9057
# [4] SNVs before min2cells: 131814
# [4] SNVs after min2cells: 131786
# [4] cell-mut pairs after min2cells: 3805563
# [4] cells with >=1 retained SNV: 9057
# [5] final output cells: 9057
# [5] final zero-variant cells: 0
# [5] final cells with variants: 9057
# [5] reads NA count: 0
# [OUT] long:  /mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/scsnv/P4_cSCC/cell_snv/P4_cSCC_scSNV_cell_mut_validcell_min2c.tsv.gz
# [OUT] count: /mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/scsnv/P4_cSCC/cell_snv/P4_cSCC_scSNV_cell_variant_count_validcell_min2c.tsv

# >>> print(df_all.head())
# sample   tool              cell  n_variants  reads                   cell_id
# 0  P4_cSCC  scSNV  AAACCTGCAGTTCATG         565  28952  P4_cSCC_AAACCTGCAGTTCATG
# 1  P4_cSCC  scSNV  AAACCTGCATCGGTTA         158   9990  P4_cSCC_AAACCTGCATCGGTTA
# 2  P4_cSCC  scSNV  AAAGATGAGGTCATCT         491  31052  P4_cSCC_AAAGATGAGGTCATCT
# 3  P4_cSCC  scSNV  AAAGCAAAGACCTTTG         291   6846  P4_cSCC_AAAGCAAAGACCTTTG
# 4  P4_cSCC  scSNV  AAATGCCCACGCATCG         453  16700  P4_cSCC_AAATGCCCACGCATCG
# >>> print("total rows:", len(df_all))
# total rows: 23372
# >>> print("unique cell_id:", df_all["cell_id"].nunique())
# unique cell_id: 23372
# >>> print("duplicated cell_id:", df_all["cell_id"].duplicated().sum())
# duplicated cell_id: 0
# >>> print("reads NA:", df_all["reads"].isna().sum())
# reads NA: 0
# >>> print("zero-variant cells:", (df_all["n_variants"] == 0).sum())
# zero-variant cells: 0
'''
SOURCES['monopogen_filter']=r'''
# Filter and process the variant calls.
#!/bin/bash
REF=/home/tantingting/scmut/ref/bwa_ref/genome.fa
MONO_BASE=/mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/monopo
samples=(P4_cSCC P5_cSCC P6_cSCC P7_cSCC P9_cSCC P10_cSCC)

for s in "${samples[@]}"; do
    echo "===== Processing $s ====="

    DIR=${MONO_BASE}/${s}/germline

    GL_ALL=${DIR}/${s}.gl.all.vcf.gz
    PHASED_ALL=${DIR}/${s}.phased.all.vcf.gz
    KEEP_KEYS=${DIR}/${s}.keep.keys.txt
    FINAL=${DIR}/${s}.phased.depth_filtered.vcf.gz

    # 1. Combine germline and phased calls across chromosomes.
    bcftools concat -a -Oz -o "$GL_ALL" $(find "$DIR" -maxdepth 1 -name "chr*.gl.vcf.gz" | sort -V)
    bcftools index -t "$GL_ALL"

    bcftools concat -a -Oz -o "$PHASED_ALL" $(find "$DIR" -maxdepth 1 -name "chr*.phased.vcf.gz" | sort -V)
    bcftools index -t "$PHASED_ALL"

    echo "[INFO] concat done for $s"

    # 2. Apply depth filtering to germline calls first.
    # SNP INFO/DP>10；INDEL INFO/DP>15
    GL_PASS=${DIR}/${s}.gl.depth_pass.vcf.gz
    bcftools view -i '
      (TYPE="snp"  && INFO/DP>10) ||
      (TYPE!="snp" && INFO/DP>15)
    ' "$GL_ALL" -Oz -o "$GL_PASS"
    bcftools index -t "$GL_PASS"

    # 3. Extract keys from filtered germline calls.
    bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "$GL_PASS" | sort -u > "$KEEP_KEYS"

    echo "[INFO] GL depth-pass sites: $(wc -l < "$KEEP_KEYS")"

    # 4. Filter phased calls using retained keys.
    bcftools view -T <(awk -F: 'BEGIN{OFS="\t"} {print $1,$2,$2,$3,$4}' "$KEEP_KEYS") "$PHASED_ALL" | \
    bcftools view -i '
      GT!="0/0" &&
      GT!="0|0" &&
      GT!="./." &&
      GT!=".|."
    ' -Oz -o "$FINAL"

    bcftools index -t "$FINAL"

    echo "[INFO] Final filtered phased VCF: $FINAL"
    echo "[INFO] Final kept sites: $(bcftools view -H "$FINAL" | wc -l)"

done

'''
SOURCES['monopogen_single_filter']=r'''
library(dplyr)
library(data.table)
library(stringr)

# =========================
# 0. Configure input paths and thresholds.
# =========================
MONO_BASE <- "/mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/monopo"
samples <- c("P4_cSCC", "P5_cSCC", "P6_cSCC", "P7_cSCC", "P9_cSCC", "P10_cSCC")
chroms <- paste0("chr", 1:22)

min_cells_with_any_reads <- 5
min_cells_with_alt_reads <- 2

# =========================
# 1. Parse reference and alternate read counts.
# =========================
parse_ref_alt_vec <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "."] <- "0/0"
  ref <- as.numeric(sub("/.*", "", x))
  alt <- as.numeric(sub(".*/", "", x))
  data.frame(ref = ref, alt = alt)
}

# =========================
# 2. Process one sample and one chromosome.
# =========================
process_one_chr <- function(sample, chr) {
  message("Processing ", sample, " ", chr, " ...")
  BASE_DIR <- file.path(MONO_BASE, sample, "somatic")
  OUT_DIR  <- file.path(MONO_BASE, sample, "somatic_combine")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  meta_file <- file.path(BASE_DIR, paste0(chr, ".putativeSNVs.csv"))
  mat_file  <- file.path(BASE_DIR, paste0(chr, ".SNV_mat.RDS"))

  if (!file.exists(meta_file)) {
    warning("Missing meta file: ", meta_file)
    return(NULL)
  }
  if (!file.exists(mat_file)) {
    warning("Missing mat file: ", mat_file)
    return(NULL)
  }

  meta <- fread(meta_file, data.table = FALSE)
  meta_filter <- meta %>%
    filter( Depth_total > 10, BAF_alt > 0.05, SVM_pos_score > 0.1 )

  message(chr, " meta total: ", nrow(meta), " | after meta filter: ", nrow(meta_filter))

  snv_pass <- paste0( meta_filter$chr, ":", meta_filter$pos, ":", meta_filter$Ref_allele, ":",
    meta_filter$Alt_allele )

  mat <- readRDS(mat_file)
  cell_cols <- colnames(mat)[19:ncol(mat)]
  common_snv <- intersect(snv_pass, rownames(mat))
  message(chr, " overlap SNVs: ", length(common_snv))
  mat_sub <- mat[common_snv, cell_cols, drop = FALSE]

  # Summarize cell support for each variant.
  variant_stat <- lapply(seq_len(nrow(mat_sub)), function(i) {
    x <- unlist(mat_sub[i, ], use.names = FALSE)
    pa <- parse_ref_alt_vec(x)

    data.frame( sample = sample, chr = chr, snv = rownames(mat_sub)[i],
      n_cells_with_any_reads = sum((pa$ref + pa$alt) > 0), n_cells_with_alt_reads = sum(pa$alt > 0),
      total_ref_reads = sum(pa$ref), total_alt_reads = sum(pa$alt), stringsAsFactors = FALSE )
  }) %>%
    bind_rows()

  variant_pass <- variant_stat %>%
    filter(
      n_cells_with_any_reads >= min_cells_with_any_reads, # Apply the configured minimum cell-support threshold.
      n_cells_with_alt_reads >= min_cells_with_alt_reads )

  message(chr, " after cell-level filter: ", nrow(variant_pass))
  message(chr, " cells: ", length(cell_cols))
  mat_pass <- mat_sub[variant_pass$snv, , drop = FALSE]

  # Count retained variants per cell.
  cell_count <- sapply(cell_cols, function(cc) {
    x <- as.character(mat_pass[, cc])
    alt <- as.numeric(sub(".*/", "", x))
    alt[is.na(alt)] <- 0
    sum(alt > 0) # Count sites with at least one alternate read in the cell.
  })

  cell_count_df <- data.frame( sample = sample, chr = chr, cell = names(cell_count),
    n_variants = as.integer(cell_count), stringsAsFactors = FALSE )

  # For cell-SNV records, retain alternate read counts > 0.
  long_df <- lapply(seq_len(nrow(mat_pass)), function(i) {
    x <- as.character(mat_pass[i, ])
    pa <- parse_ref_alt_vec(x)

    data.frame( sample = sample, chr = chr, snv = rownames(mat_pass)[i], cell = colnames(mat_pass),
      ref_reads = pa$ref, alt_reads = pa$alt, stringsAsFactors = FALSE
    ) %>%
    filter(alt_reads > 0) %>%
    mutate( cell_id = paste0(sample, "_", cell), total_reads = ref_reads + alt_reads,
      baf_cell = alt_reads / total_reads )
  }) %>%
    bind_rows()

  # Release large intermediate objects.
  rm(mat, mat_sub, mat_pass)
  gc()

  list( variant_stat = variant_stat, variant_pass = variant_pass, cell_count_df = cell_count_df,
    long_df = long_df )
}

# =========================
# 3. Process all chromosomes for one sample.
# =========================
process_one_sample <- function(sample) {
  message("==============================")
  message("Start sample: ", sample)
  message("==============================")

  OUT_DIR <- file.path(MONO_BASE, sample, "somatic_combine")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  res_list <- lapply(chroms, function(chr) {
    process_one_chr(sample, chr)
  })

  res_list <- res_list[!sapply(res_list, is.null)]

  if (length(res_list) == 0) {
    warning(sample, ": no valid chromosome results.")
    return(NULL)
  }

  variant_stat_all <- bind_rows(lapply(res_list, `[[`, "variant_stat"))
  variant_pass_all <- bind_rows(lapply(res_list, `[[`, "variant_pass"))
  cell_count_by_chr <- bind_rows(lapply(res_list, `[[`, "cell_count_df"))
  cell_snv_long_all <- bind_rows(lapply(res_list, `[[`, "long_df"))

  cell_variant_count <- cell_count_by_chr %>%
    group_by(sample, cell) %>%
    summarise( n_variants = sum(n_variants), .groups = "drop"
    ) %>%
    mutate( cell_id = paste0(sample, "_", cell) )

  # Write outputs.
  write.table( variant_stat_all,
    file = file.path(OUT_DIR, paste0(sample, "_Monopogen_variant_cell_support_stat.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE )

  write.table( variant_pass_all,
    file = file.path(OUT_DIR, paste0(sample, "_Monopogen_filtered_variant_pass.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE )

  write.table( cell_count_by_chr,
    file = file.path(OUT_DIR, paste0(sample, "_Monopogen_filtered_cell_variant_count_by_chr.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE )

  write.table( cell_variant_count,
    file = file.path(OUT_DIR, paste0(sample, "_Monopogen_filtered_cell_variant_count.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE )

  write.table( cell_snv_long_all,
    file = file.path(OUT_DIR, paste0(sample, "_Monopogen_filtered_cell_snv_long.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE )

  message("Done sample: ", sample)
  message("Filtered variants: ", nrow(variant_pass_all))
  message("Cell-SNV long pairs: ", nrow(cell_snv_long_all))
  message("Cells: ", nrow(cell_variant_count))
  message("Output: ", OUT_DIR)

  list( variant_stat_all = variant_stat_all, variant_pass_all = variant_pass_all,
    cell_count_by_chr = cell_count_by_chr, cell_variant_count = cell_variant_count,
    cell_snv_long_all = cell_snv_long_all )
}

# =========================
# 4. Iterate over the configured P4-P10 samples.
# =========================
all_sample_results <- lapply(samples, process_one_sample)
names(all_sample_results) <- samples

all_sample_results <- all_sample_results[!sapply(all_sample_results, is.null)]

# =========================
# 5. Combine the per-sample tables.
# =========================
ALL_OUT_DIR <- file.path(MONO_BASE, "all_samples_somatic_combine")
dir.create(ALL_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

all_cell_variant_count <- bind_rows(lapply(all_sample_results, `[[`, "cell_variant_count"))
all_cell_snv_long <- bind_rows(lapply(all_sample_results, `[[`, "cell_snv_long_all"))
all_variant_pass <- bind_rows(lapply(all_sample_results, `[[`, "variant_pass_all"))

write.table( all_cell_variant_count,
  file = file.path(ALL_OUT_DIR, "All_samples_Monopogen_filtered_cell_variant_count.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE )

write.table( all_cell_snv_long,
  file = file.path(ALL_OUT_DIR, "All_samples_Monopogen_filtered_cell_snv_long.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE )

write.table( all_variant_pass,
  file = file.path(ALL_OUT_DIR, "All_samples_Monopogen_filtered_variant_pass.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE )

message("All samples done!")
message("All cells: ", nrow(all_cell_variant_count))
message("All cell-SNV long pairs: ", nrow(all_cell_snv_long))
message("All filtered variants: ", nrow(all_variant_pass))
message("All output: ", ALL_OUT_DIR)
'''
SOURCES['cellsnp_single_filter']=r'''
import gzip
import os
import csv
from collections import defaultdict

samples = ["P4_cSCC", "P5_cSCC", "P6_cSCC", "P7_cSCC", "P9_cSCC", "P10_cSCC"]
CELL_BASE = "/mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/cellsnp"
OUT_BASE = "/mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/cellsnp/cell_variant_count"
os.makedirs(OUT_BASE, exist_ok=True)

def parse_info(info):
    d = {}
    for x in info.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
    return d

def get_pass_sites_from_base(base_vcf):
    """
    First-stage filtering criteria.
    chr1-22
    INFO/DP > 10
    INFO/AD / INFO/DP >= 0.05
    """
    pass_sites = set()
    with gzip.open(base_vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            arr = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt, qual, filt, info = arr[:8]
            if chrom not in set(map(str, range(1, 23))):
                continue
            info_d = parse_info(info)
            try:
                dp = int(info_d.get("DP", 0))
                ad = int(info_d.get("AD", 0))
            except ValueError:
                continue
            if dp <= 10:
                continue
            af = ad / dp if dp > 0 else 0
            if af <= 0.05:
                continue
            site_id = f"{chrom}:{pos}:{ref}:{alt}"
            pass_sites.add(site_id)
    return pass_sites

# Exclude reference genotypes even if alternate read counts are positive.
def count_cell_variants(sample):
    sample_dir = os.path.join(CELL_BASE, f"{sample}_1a")
    base_vcf = os.path.join(sample_dir, "cellSNP.base.vcf.gz")
    cell_vcf = os.path.join(sample_dir, "cellSNP.cells.vcf.gz")
    if not os.path.exists(base_vcf):
        print(f"[WARN] missing base VCF: {base_vcf}")
        return
    if not os.path.exists(cell_vcf):
        print(f"[WARN] missing cell VCF: {cell_vcf}")
        return
    print(f"\n===== {sample} =====")
    pass_sites = get_pass_sites_from_base(base_vcf)
    print(f"{sample}: first-round retained sites = {len(pass_sites)}")
    cell_barcodes = []
    cell_variant_count = defaultdict(int)
    cell_variant_ids = defaultdict(list)
    n_total_sites = 0
    n_round1_sites = 0
    n_round2_sites = 0
    with gzip.open(cell_vcf, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.split("\t")
                cell_barcodes = header[9:] # Store the cell identifiers.
                continue
            arr = line.split("\t")
            n_total_sites += 1
            chrom, pos, _id, ref, alt = arr[0], arr[1], arr[2], arr[3], arr[4]
            if chrom not in set(map(str, range(1, 23))):
                continue
            site_id = f"{chrom}:{pos}:{ref}:{alt}"
            if site_id not in pass_sites: # Count only sites retained by the first filtering stage.
                continue
            n_round1_sites += 1

            fmt = arr[8].split(":") # Locate the GT, AD and DP fields.
            try:
                gt_idx = fmt.index("GT")
                ad_idx = fmt.index("AD")
                dp_idx = fmt.index("DP")
            except ValueError:
                continue
            sample_fields = arr[9:]
            cells_with_reads = 0
            cells_with_alt = 0
            alt_cells = []
            for i, g in enumerate(sample_fields):
                if g == "." or g.startswith(".:"):
                    continue
                parts = g.split(":")
                if len(parts) <= max(gt_idx, ad_idx, dp_idx):
                    continue # Skip records with incomplete GT, AD or DP fields.
                try:
                    gt = parts[gt_idx]
                    ad = int(parts[ad_idx])
                    dp = int(parts[dp_idx])
                except ValueError:
                    continue
                if dp > 0:
                    cells_with_reads += 1
                # Count a cell-site mutation only when GT is non-reference and AD > 0.
                if ad > 0 and gt not in {"0/0", "0|0", "./.", "."}:
                    cells_with_alt += 1
                    alt_cells.append(i)
            if cells_with_reads >= 5 and cells_with_alt >= 2:
                n_round2_sites += 1
                for i in alt_cells:
                    cell = cell_barcodes[i]
                    cell_variant_count[cell] += 1
                    cell_variant_ids[cell].append(site_id)
    out_file = os.path.join( OUT_BASE, f"{sample}_cellSNP_filtered_cell_variant_count.tsv" )
    with open(out_file, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sample", "cell", "n_variants"])
        for cell in cell_barcodes:
            writer.writerow([ sample, cell, cell_variant_count.get(cell, 0)
            ])
    detail_file = os.path.join(OUT_BASE,f"{sample}_cellSNP_filtered_cell_variant_detail.tsv")
    with open(detail_file, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["sample", "cell", "variant_id"])
        for cell in cell_barcodes:
            for variant_id in cell_variant_ids.get(cell, []):
                writer.writerow([sample,cell,variant_id])
    print(f"{sample}: total cell VCF sites scanned = {n_total_sites}")
    print(f"{sample}: matched first-round sites = {n_round1_sites}")
    print(f"{sample}: final second-round sites = {n_round2_sites}")
    print(f"{sample}: output = {out_file}")

for sample in samples:
    count_cell_variants(sample)
'''
SOURCES['scallele_single_filter']=r'''
DIR="/mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/scallele"
REF="/home/tantingting/scmut/ref/bwa_ref/genome.fa"

samples=("P4_cSCC" "P5_cSCC" "P6_cSCC" "P7_cSCC" "P9_cSCC" "P10_cSCC")

CHR_LIST="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

ALL_LONG="${DIR}/scAllele_all_samples_cell_snv_long.tsv"
ALL_COUNT="${DIR}/scAllele_all_samples_cell_variant_count.tsv"

echo -e "tool\tsample\tcell\tsnv" > "$ALL_LONG"
echo -e "sample\tcell\tn_variants\tcell_id" > "$ALL_COUNT"

for s in "${samples[@]}"; do
  echo "========== Processing $s =========="

  INDIR="${DIR}/${s}"
  OUTDIR="${DIR}/${s}/PASS_norm_chr1_22"
  mkdir -p "$OUTDIR"

  SAMPLE_LONG="${OUTDIR}/${s}_cell_snv_long.tsv"
  SAMPLE_COUNT="${OUTDIR}/${s}_cell_variant_count.tsv"

  echo -e "tool\tsample\tcell\tsnv" > "$SAMPLE_LONG"
  echo -e "sample\tcell\tn_variants\tcell_id" > "$SAMPLE_COUNT"

  for RAW in "${INDIR}"/*.vcf; do
    [[ -e "$RAW" ]] || continue

    base=$(basename "$RAW" .vcf)
    cell=${base#${s}_map_}

    echo "--- Processing $cell ---"

    TMP_BODY="${OUTDIR}/${base}.body.vcf"
    TMP_FULL="${OUTDIR}/${base}.full.vcf"
    TMP_GZ="${OUTDIR}/${base}.full.vcf.gz"
    FINAL="${OUTDIR}/${base}_PASS_norm_chr1_22.vcf.gz"

    # 1. Retain the ten standard VCF columns and enforce tab separation.
    # Skip malformed records that do not have ten VCF columns.
    awk 'BEGIN{OFS="\t"} !/^#/ && NF==10 {
      print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10
    }' "$RAW" > "$TMP_BODY"

    # 2. Build a standard VCF header with a tab-separated final header line.
    {
      echo "##fileformat=VCFv4.2"
      echo "##source=scAllele"
      echo '##FILTER=<ID=PASS,Description="All filters passed">'
      echo '##FILTER=<ID=LowQual,Description="Low quality variant">'
      echo '##INFO=<ID=varEnd,Number=1,Type=Integer,Description="Variant end position">'
      echo '##INFO=<ID=varGroup,Number=1,Type=String,Description="Variant group">'
      echo '##INFO=<ID=varType,Number=1,Type=String,Description="Variant type">'
      echo '##INFO=<ID=varLength,Number=1,Type=Integer,Description="Variant length">'
      echo '##INFO=<ID=TandemRep,Number=1,Type=Integer,Description="Tandem repeat count">'
      echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
      echo '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">'
      echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">'
      echo '##FORMAT=<ID=AC,Number=1,Type=Integer,Description="Alt count">'
      echo '##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Ref count">'
      echo '##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allele balance">'
      echo '##FORMAT=<ID=RCL,Number=1,Type=Integer,Description="Read clipping length">'
      printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    } > "$TMP_FULL"

    cat "$TMP_BODY" >> "$TMP_FULL"

    # Create a header-only VCF when there are no variant records.
    bgzip -f -c "$TMP_FULL" > "$TMP_GZ"
    tabix -f -p vcf "$TMP_GZ"

    # 3. PASS + chr1-22 + norm
    bcftools view -r "$CHR_LIST" -f PASS "$TMP_GZ" -Ou | bcftools norm -f "$REF" -m -both -Oz -o "$FINAL"

    bcftools index -t "$FINAL"

    n=$(bcftools view -H "$FINAL" | wc -l)

    echo -e "${s}\t${cell}\t${n}\t${s}_${cell}" >> "$SAMPLE_COUNT"
    echo -e "${s}\t${cell}\t${n}\t${s}_${cell}" >> "$ALL_COUNT"

    if [[ "$n" -gt 0 ]]; then
      bcftools query -f "scAllele\t${s}\t${cell}\t%CHROM:%POS:%REF:%ALT\n" "$FINAL" >> "$SAMPLE_LONG"

      bcftools query -f "scAllele\t${s}\t${cell}\t%CHROM:%POS:%REF:%ALT\n" "$FINAL" >> "$ALL_LONG"
    fi

    rm -f "$TMP_BODY" "$TMP_FULL" "$TMP_GZ" "$TMP_GZ.tbi"

  done

  echo "[DONE] $s"
done

echo "========== ALL DONE =========="
echo "[INFO] all count: $ALL_COUNT"
echo "[INFO] all long:  $ALL_LONG"
'''
SOURCES['scallele_recurrence']=r'''
# Apply the additional recurrence filter.
library(data.table)
library(dplyr)

DIR <- "/mnt/storage/tantingting/scmut/vcf/1_cSCC_scRNA/scallele"

in_long <- file.path(DIR, "scAllele_all_samples_cell_snv_long.tsv")
in_count <- file.path(DIR, "scAllele_all_samples_cell_variant_count.tsv")

out_long <- file.path(DIR, "scAllele_all_samples_cell_snv_long_min2cells.tsv")
out_count <- file.path(DIR, "scAllele_all_samples_cell_variant_count_min2cells.tsv")

long <- fread(in_long) %>%
  distinct(tool, sample, cell, snv)

cell_all <- fread(in_count) %>%
  select(sample, cell, cell_id) %>%
  distinct(sample, cell, .keep_all = TRUE)

snv_keep <- long %>%
  group_by(sample, snv) %>%
  summarise(n_cells = n_distinct(cell), .groups = "drop") %>%
  filter(n_cells >= 2)

long_filt <- long %>%
  inner_join(snv_keep[, c("sample", "snv")], by = c("sample", "snv"))

count_filt_tmp <- long_filt %>%
  group_by(sample, cell) %>%
  summarise(n_variants = n_distinct(snv), .groups = "drop")

count_filt <- cell_all %>%
  left_join(count_filt_tmp, by = c("sample", "cell")) %>%
  mutate(n_variants = ifelse(is.na(n_variants), 0L, as.integer(n_variants)), cell_id = paste0(sample, "_", cell)) %>%
  select(sample, cell, n_variants, cell_id)

fwrite(long_filt, out_long, sep = "\t")
fwrite(count_filt, out_count, sep = "\t")

print(long %>%
  group_by(sample) %>%
  summarise(raw_snv = n_distinct(snv), raw_pairs = n(), raw_cells = n_distinct(cell)))
print(long_filt %>%
  group_by(sample) %>%
  summarise(filt_snv = n_distinct(snv), filt_pairs = n(), filt_cells = n_distinct(cell)))
print(count_filt %>%
  group_by(sample) %>%
  summarise(count_cells = n()))
'''
SOURCES['scallele_smart_filter']=r'''
#!/usr/bin/env bash

set -uo pipefail

REF="/home/tantingting/scmut/ref/bwa_ref/genome.fa"
MANIFEST="/mnt/storage/tantingting/scmut_smart/bam/star_qc_unique_ge10k_scrna_sample_wes.tsv"

TNBC_DIR="/mnt/storage/tantingting/scmut_smart/vcf/1_TNBC_scRNA/scallele"
BRCA_DIR="/mnt/storage/tantingting/scmut_smart/vcf/2_BRCA_scRNA/scallele"

CHR_LIST="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrM,chrX,chrY"

# Combine outputs from both datasets.
ALL_LONG="/mnt/storage/tantingting/scmut_smart/vcf/filter_stat/scAllele_Smartseq2_all_cell_snv_long.tsv"
ALL_COUNT="/mnt/storage/tantingting/scmut_smart/vcf/filter_stat/scAllele_Smartseq2_all_cell_variant_count.tsv"
MISSING="/mnt/storage/tantingting/scmut_smart/vcf/scAllele_Smartseq2_missing_vcf.tsv"

echo -e "tool\tcohort\tsample\tcell\tsnv" > "$ALL_LONG"
echo -e "cohort\tsample\tcell\tn_variants\tcell_id\tWES_SRR" > "$ALL_COUNT"
echo -e "scRNA_SRR\tsample\tWES_SRR\treason" > "$MISSING"

# Initialize the combined output files for each dataset.
for DIR in "$TNBC_DIR" "$BRCA_DIR"; do
  echo -e "tool\tcohort\tsample\tcell\tsnv" > "${DIR}/scAllele_all_cells_snv_long.tsv"
  echo -e "cohort\tsample\tcell\tn_variants\tcell_id\tWES_SRR" > "${DIR}/scAllele_all_cells_variant_count.tsv"
done

# Read the scRNA_SRR, Sample and WES_SRR manifest columns.
awk 'NR>1 && NF>=3 {print $1"\t"$2"\t"$3}' "$MANIFEST" |
while IFS=$'\t' read -r SCRNA_SRR SAMPLE WES_SRR; do

  # Identify TNBC or BRCA from the VCF directory.
  if [[ -f "${TNBC_DIR}/${SCRNA_SRR}.vcf" ]]; then
    DIR="$TNBC_DIR"
    COHORT="TNBC"
  elif [[ -f "${BRCA_DIR}/${SCRNA_SRR}.vcf" ]]; then
    DIR="$BRCA_DIR"
    COHORT="BRCA"
  else
    echo "[WARNING] VCF not found: ${SCRNA_SRR}"
    echo -e "${SCRNA_SRR}\t${SAMPLE}\t${WES_SRR}\tVCF_not_found" >> "$MISSING"
    continue
  fi

  echo "========== ${COHORT} | ${SAMPLE} | ${SCRNA_SRR} =========="

  RAW="${DIR}/${SCRNA_SRR}.vcf"
  OUTDIR="${DIR}/PASS_norm_chr1_22"
  mkdir -p "$OUTDIR"

  TMP_BODY="${OUTDIR}/${SCRNA_SRR}.body.vcf"
  TMP_FULL="${OUTDIR}/${SCRNA_SRR}.full.vcf"
  TMP_GZ="${OUTDIR}/${SCRNA_SRR}.full.vcf.gz"
  FINAL="${OUTDIR}/${SCRNA_SRR}_PASS.vcf.gz"

  DATASET_LONG="${DIR}/scAllele_all_cells_snv_long.tsv"
  DATASET_COUNT="${DIR}/scAllele_all_cells_variant_count.tsv"

  # Keep the ten standard VCF columns and enforce tab separation.
  awk 'BEGIN{OFS="\t"} !/^#/ && NF==10 {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10
  }' "$RAW" > "$TMP_BODY"

  # Reconstruct a standard VCF header.
  {
    echo "##fileformat=VCFv4.2"
    echo "##source=scAllele"
    echo '##FILTER=<ID=PASS,Description="All filters passed">'
    echo '##FILTER=<ID=LowQual,Description="Low quality variant">'
    echo '##INFO=<ID=varEnd,Number=1,Type=Integer,Description="Variant end position">'
    echo '##INFO=<ID=varGroup,Number=1,Type=String,Description="Variant group">'
    echo '##INFO=<ID=varType,Number=1,Type=String,Description="Variant type">'
    echo '##INFO=<ID=varLength,Number=1,Type=Integer,Description="Variant length">'
    echo '##INFO=<ID=TandemRep,Number=1,Type=Integer,Description="Tandem repeat count">'
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    echo '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">'
    echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">'
    echo '##FORMAT=<ID=AC,Number=1,Type=Integer,Description="Alt count">'
    echo '##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Ref count">'
    echo '##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allele balance">'
    echo '##FORMAT=<ID=RCL,Number=1,Type=Integer,Description="Read clipping length">'
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" "$SCRNA_SRR"
  } > "$TMP_FULL"

  cat "$TMP_BODY" >> "$TMP_FULL"
  bgzip -f -c "$TMP_FULL" > "$TMP_GZ"
  tabix -f -p vcf "$TMP_GZ"

  # Retain PASS calls on chromosomes 1-22, X, Y and M, then normalize.
     bcftools view -r "$CHR_LIST" -f PASS "$TMP_GZ" -Ou |
     bcftools norm -f "$REF" -m -both -Oz -o "$FINAL"
    bcftools index -f -t "$FINAL"
    N=$(bcftools view -H "$FINAL" | wc -l)

  CELL_ID="${SAMPLE}_${SCRNA_SRR}"

  # Write dataset-specific and combined summary tables.
  COUNT_LINE="${COHORT}\t${SAMPLE}\t${SCRNA_SRR}\t${N}\t${CELL_ID}\t${WES_SRR}"
  echo -e "$COUNT_LINE" >> "$DATASET_COUNT"
  echo -e "$COUNT_LINE" >> "$ALL_COUNT"

  if [[ "$N" -gt 0 ]]; then
    bcftools query -f "scAllele\t${COHORT}\t${SAMPLE}\t${SCRNA_SRR}\t%CHROM:%POS:%REF:%ALT\n" \
      "$FINAL" >> "$DATASET_LONG"

    bcftools query -f "scAllele\t${COHORT}\t${SAMPLE}\t${SCRNA_SRR}\t%CHROM:%POS:%REF:%ALT\n" \
      "$FINAL" >> "$ALL_LONG"
  fi

  rm -f "$TMP_BODY" "$TMP_FULL" "$TMP_GZ" "$TMP_GZ.tbi"

  echo "[DONE] ${SCRNA_SRR}: ${N} variants"
done

echo "========== ALL DONE =========="
echo "[INFO] TNBC count: ${TNBC_DIR}/scAllele_all_cells_variant_count.tsv"
echo "[INFO] TNBC long:  ${TNBC_DIR}/scAllele_all_cells_snv_long.tsv"
echo "[INFO] BRCA count: ${BRCA_DIR}/scAllele_all_cells_variant_count.tsv"
echo "[INFO] BRCA long:  ${BRCA_DIR}/scAllele_all_cells_snv_long.tsv"
echo "[INFO] ALL count:  ${ALL_COUNT}"
echo "[INFO] ALL long:   ${ALL_LONG}"
echo "[INFO] Missing:    ${MISSING}"
'''
SOURCES['scomatic_filter']=r'''
# Split SComatic calls into cell-type-specific tables.
import sys

infile = sys.argv[1]

out1 = open("Myeloid.tsv", "w")
out2 = open("T.tsv", "w")
out3 = open("Epithelial.tsv", "w")

for line in open(infile):
    line = line.rstrip()

    # Retain metadata comment lines in all three output files.
    if line.startswith("##"):
        out1.write(line + "\n")
        out2.write(line + "\n")
        out3.write(line + "\n")
        continue

    # Write the header.
    if line.startswith("#CHROM"):
        header = line.split()
        out1.write("\t".join(header[:25] + ["Myeloid"]) + "\n")
        out2.write("\t".join(header[:25] + ["T"]) + "\n")
        out3.write("\t".join(header[:25] + ["Epithelial"]) + "\n")
        continue

    x = line.split()

    chrom = x[0]
    start = x[1]
    end = x[2]
    ref = x[3]
    alt = x[4]
    filt = x[5]
    celltypes = x[6]
    ct_filters = x[23] # Retain this field without splitting its commas.

    # Remove records without alternate alleles or assigned cell types.
    if alt == "." or celltypes == ".":
        continue

    filt_low = filt.lower()

    # Remove noise and low-complexity upstream/downstream records.
    if ("noisy_site" in filt_low) or ("lc_upstream" in filt_low) or ("noise" in filt_low) or ("lc_downstream" in filt_low):
      continue

    chrom2 = chrom.replace("chr", "")
    if chrom2 not in [str(i) for i in range(1, 23)] + ["X", "Y", "M"]:
        continue

    # These fields correspond to individual cell types.
    alts = x[4].split(",")
    dps = x[10].split(",") # Depth is stored in the eleventh column.
    ncs = x[11].split(",")
    bcs = x[12].split(",")
    ccs = x[13].split(",")
    vafs = x[14].split(",")
    ccfs = x[15].split(",")
    bcps = x[16].split(",")
    ccps = x[17].split(",")
    cts = x[6].split(",")

    # Replicate singleton fields so their lengths match the cell-type list.
    if len(alts) == 1 and len(cts) > 1:
        alts = alts * len(cts)

    if len(dps) == 1 and len(cts) > 1:
        dps = dps * len(cts)

    if len(ncs) == 1 and len(cts) > 1:
        ncs = ncs * len(cts)

    if len(bcs) == 1 and len(cts) > 1:
        bcs = bcs * len(cts)

    if len(ccs) == 1 and len(cts) > 1:
        ccs = ccs * len(cts)

    if len(vafs) == 1 and len(cts) > 1:
        vafs = vafs * len(cts)

    if len(ccfs) == 1 and len(cts) > 1:
        ccfs = ccfs * len(cts)

    if len(bcps) == 1 and len(cts) > 1:
        bcps = bcps * len(cts)

    if len(ccps) == 1 and len(cts) > 1:
        ccps = ccps * len(cts)

    for i in range(len(cts)):
        ct = cts[i]

        # Use the depth for the current cell type.
        dp = dps[i]
        if dp == "." or dp == "NA":
            continue
        if float(dp) <= 10: # Exclude records with DP <= 10.
            continue

        this_alt = alts[i]
        this_nc = ncs[i]
        this_bc = bcs[i]
        this_cc = ccs[i]
        this_vaf = vafs[i]
        this_ccf = ccfs[i]
        this_bcp = bcps[i]
        this_ccp = ccps[i]

        # Split multiple alternate alleles within a cell type, such as G|T.
        alt_list = this_alt.split("|")
        bc_list = this_bc.split("|")
        cc_list = this_cc.split("|")
        vaf_list = this_vaf.split("|")
        ccf_list = this_ccf.split("|")
        bcp_list = this_bcp.split("|")
        ccp_list = this_ccp.split("|")

        for j in range(len(alt_list)):
            new = x[:]   # Copy the original row.

            # Update fields for the current cell type and allele.
            new[4] = alt_list[j]                 # ALT
            new[6] = ct                          # Cell_types
            new[9] = "1"
            new[10] = dp                         # Dp
            new[11] = this_nc                    # Nc
            new[12] = bc_list[j] if j < len(bc_list) else this_bc
            new[13] = cc_list[j] if j < len(cc_list) else this_cc
            new[14] = vaf_list[j] if j < len(vaf_list) else this_vaf
            new[15] = ccf_list[j] if j < len(ccf_list) else this_ccf
            new[16] = bcp_list[j] if j < len(bcp_list) else this_bcp
            new[17] = ccp_list[j] if j < len(ccp_list) else this_ccp

            # Retain the first 25 columns and the current cell-type column.
            if ct == "Myeloid":
                keep = new[:25] + [new[25]]
                out1.write("\t".join(keep) + "\n")
            elif ct == "T":
                keep = new[:25] + [new[26]]
                out2.write("\t".join(keep) + "\n")
            elif ct == "Epithelial":
                keep = new[:25] + [new[27]]
                out3.write("\t".join(keep) + "\n")

out1.close()
out2.close()
out3.close()
'''
import argparse
import contextlib
import csv
import datetime
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

TOOLS = {
    ('10x', 'pooled'): ['gatk','freebayes','samtools','varscan2','strelka2','monopogen','scsnv','cellsnp'],
    ('10x', 'cell-type'): ['scomatic','cellsnp','samtools','freebayes','varscan2','strelka2','monovar'],
    ('10x', 'single-cell'): ['monopogen','scsnv','scallele','cellsnp'],
    ('smartseq', 'single-cell'): ['gatk','samtools','freebayes','varscan2','strelka2','cellsnp','scallele','monovar']
}
CHROMS = ','.join('chr'+str(i) for i in range(1,23))+',chrX,chrY,chrM'


def quote(value):
    return shlex.quote(str(value))


def required(record, key):
    value = record.get(key)
    unset=isinstance(value,str) and 'REPLACE_' in value
    if isinstance(value,list):unset=any(isinstance(item,str) and 'REPLACE_' in item for item in value)
    if value is None or value == '' or unset:
        raise ValueError('Configure '+key)
    return value


def validate_id(value):
    if not re.fullmatch(r'[A-Za-z0-9_.-]+',value) or value in ('.','..'):
        raise ValueError('Invalid sample or unit identifier: '+value)
    return value


def source_text(name):
    return SOURCES[name]


def job(name, environment, command, inputs=(), outputs=(), values=None):
    return dict(name=name,environment=environment,command=command,
                inputs=list(map(str,inputs)),outputs=list(map(str,outputs)),values=values or {})


def python_job(name, environment, script, arguments, inputs=(), outputs=(), values=None):
    command="python - "+' '.join(quote(a) for a in arguments)+" <<'SCMUT_PY'\n"+script+"\nSCMUT_PY"
    return job(name,environment,command,inputs,outputs,values)


def run_jobs(jobs, cfg, root, dry_run):
    if dry_run:
        for step in jobs:
            print('\n['+step['name']+'] '+cfg.get('environments',{}).get(step['environment'],step['environment']))
            if step['values']:
                print(' '.join(k+'='+quote(v) for k,v in step['values'].items()))
            print(step['command'])
        return
    conda = required(cfg,'conda')
    envs = required(cfg,'environments')
    produced = {p for step in jobs for p in step['outputs']}
    missing = sorted({p for step in jobs for p in step['inputs'] if p not in produced and not Path(p).exists()})
    if missing:
        raise ValueError('Missing input files:\n'+'\n'.join(missing))
    for step in jobs:
        required(envs,step['environment'])
    if not Path(conda).is_file():
        raise ValueError('Conda executable not found: '+conda)
    root.mkdir(parents=True,exist_ok=True)
    with (root/'.workflow.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        stamp=datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f')
        logs=root/'logs'/stamp
        logs.mkdir(parents=True)
        records=[]
        for index,step in enumerate(jobs,1):
            log=logs/(f'{index:03d}-'+step['name']+'.log')
            environment=os.environ.copy()
            environment.update({k:str(v) for k,v in step['values'].items()})
            environment['SCMUT_OUTPUT_ROOT']=str(root)
            env_name=envs[step['environment']]
            selector='--prefix' if '/' in env_name else '--name'
            args=[conda,'run','--no-capture-output',selector,env_name,'bash','-e','-o','pipefail','-c',step['command']]
            print('Running '+step['name']+' in '+env_name,flush=True)
            with log.open('w') as stream:
                subprocess.run(args,env=environment,stdout=stream,stderr=subprocess.STDOUT,check=True)
            absent=[p for p in step['outputs'] if not Path(p).is_file()]
            if absent:
                raise RuntimeError('Expected outputs missing after '+step['name']+': '+', '.join(absent))
            records.append(dict(step=step['name'],environment=env_name,log=str(log),inputs=step['inputs'],outputs=step['outputs'],command=step['command'],values=step['values']))
            temporary=logs/'completed.json.tmp'
            temporary.write_text(json.dumps(records,indent=2)+'\n')
            temporary.replace(logs/'completed.json')


def fixed_filter(tool, platform, level):
    profile='smartseq' if platform=='smartseq' else '10x_'+('pooled' if level=='pooled' else 'celltype')
    return FILTERS[tool][profile]


def filtered_job(tool, ctx, cfg, platform, level, joint_sample=None):
    expression=fixed_filter(tool,platform,level)
    select=''
    if joint_sample:
        select=' -s '+quote(joint_sample)
    prepare=''
    if not ctx['RAW'].endswith('.gz'):
        prepare='bcftools view "$RAW" -Oz -o "$OUT/raw.indexed.vcf.gz"\nbcftools index -f -t "$OUT/raw.indexed.vcf.gz"\nRAW="$OUT/raw.indexed.vcf.gz"\n'
    view='bcftools view'+select+' -r "$CHROMS" "$RAW" -Ou'
    if tool=='strelka2':
        view=view.replace(' -r "$CHROMS"',' -i \'ALT!="." && FILTER="PASS"\' -r "$CHROMS"',1)
    norm='bcftools norm -f "$REF" -m -both'
    filt='bcftools view -i '+quote(expression)
    stages=[view,filt+' -Ou',norm+' -Oz -o "$FINAL"'] if tool=='monovar' else [view,norm+' -Ou',filt+' -Oz -o "$FINAL"']
    command=prepare+' | '.join(stages)+'\nbcftools index -f -t "$FINAL"'

    return job(tool+'-filter','scmut',command,[ctx['RAW']],[ctx['FINAL'],ctx['FINAL']+'.tbi'],ctx)


def gatk_jobs(ctx, platform):
    guarded=platform=='smartseq'
    snp=[('QD < 2.0','QD2'),('QUAL < 30.0','QUAL30'),('SOR > 3.0','SOR3'),('FS > 60.0','FS60'),('MQ < 40.0','MQ40')]
    snp += [("(vc.hasAttribute('MQRankSum') && MQRankSum < -12.5)" if guarded else 'MQRankSum < -12.5','MQRS12.5'),
            ("(vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)" if guarded else 'ReadPosRankSum < -8.0','RPRS8')]
    indel=[('QD < 2.0','QD2'),('QUAL < 30.0','QUAL30'),('FS > 200.0','FS200'),
           ("(vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)" if guarded else 'ReadPosRankSum < -20.0','RPRS20')]
    def flags(items):
        return ' '.join('-filter '+quote(expr)+' --filter-name '+quote(label) for expr,label in items)
    call='''mkdir -p "$OUT"
gatk SplitNCigarReads -R "$REF" -I "$BAM" -O "$OUT/split.bam" --create-output-bam-index true
gatk HaplotypeCaller -R "$REF" -I "$OUT/split.bam" -O "$RAW" \\
  --native-pair-hmm-threads 16 --standard-min-confidence-threshold-for-calling 30.0
bcftools index -f -t "$RAW"'''
    filtering='''bcftools norm -f "$REF" -m -both "$RAW" -Oz -o "$OUT/norm.vcf.gz"
bcftools index -f -t "$OUT/norm.vcf.gz"
gatk SelectVariants -R "$REF" -V "$OUT/norm.vcf.gz" -select-type SNP -O "$OUT/snp.vcf.gz"
gatk SelectVariants -R "$REF" -V "$OUT/norm.vcf.gz" -select-type INDEL -O "$OUT/indel.vcf.gz"
gatk VariantFiltration -R "$REF" -V "$OUT/snp.vcf.gz" '''+flags(snp)+''' -O "$OUT/snp.filtered.vcf.gz"
gatk VariantFiltration -R "$REF" -V "$OUT/indel.vcf.gz" '''+flags(indel)+''' -O "$OUT/indel.filtered.vcf.gz"
'''
    if guarded:
        filtering+='''bcftools view -f PASS -i 'FMT/DP[0]>10' "$OUT/snp.filtered.vcf.gz" -Oz -o "$OUT/snp.pass.vcf.gz"
bcftools view -f PASS -i 'FMT/DP[0]>15' "$OUT/indel.filtered.vcf.gz" -Oz -o "$OUT/indel.pass.vcf.gz"
'''
    else:
        filtering+='''bcftools view -i 'FMT/DP[0]>10' "$OUT/snp.filtered.vcf.gz" -Oz -o "$OUT/snp.pass.vcf.gz"
bcftools view -i 'FMT/DP[0]>15' "$OUT/indel.filtered.vcf.gz" -Oz -o "$OUT/indel.pass.vcf.gz"
'''
    filtering+='bcftools index -f -t "$OUT/snp.pass.vcf.gz"\nbcftools index -f -t "$OUT/indel.pass.vcf.gz"\n'
    filtering+='''bcftools concat -a "$OUT/snp.pass.vcf.gz" "$OUT/indel.pass.vcf.gz" -Ou | \\
  bcftools sort -Oz -o "$OUT/merged.vcf.gz"
bcftools index -f -t "$OUT/merged.vcf.gz"
bcftools view -r "$CHROMS"'''+('' if guarded else ' -f PASS')+''' "$OUT/merged.vcf.gz" -Oz -o "$FINAL"
bcftools index -f -t "$FINAL"'''
    if guarded:
        filtering=filtering.replace(' -Ou | '+chr(92)+'\n  bcftools sort -Oz',' -Oz')
    return [job('gatk-call','gatk',call,[ctx['BAM']],[ctx['RAW']],ctx),job('gatk-filter','scmut',filtering,[ctx['RAW']],[ctx['FINAL'],ctx['FINAL']+'.tbi'],ctx)]


def regular_jobs(tool, ctx, cfg, platform, level):
    if tool=='gatk':
        return gatk_jobs(ctx,platform)
    if tool=='freebayes':
        command='''mkdir -p "$OUT"
freebayes-parallel "$REGIONS" 16 -F 0.05 --min-coverage 10 -m 20 -q 20 -f "$REF" "$BAM" > "$RAW"'''
        return [job('freebayes-call','freebayes',command,[ctx['BAM'],ctx['REGIONS']],[ctx['RAW']],ctx),filtered_job(tool,ctx,cfg,platform,level)]
    if tool=='samtools':
        command='''mkdir -p "$OUT"
bcftools mpileup -Q 20 -q 20 -A -x -Ou --threads 16 -f "$REF" "$BAM" | \\
  bcftools call --threads 16 -mv -Ov -o "$RAW"'''
        return [job('samtools-call','samtools',command,[ctx['BAM']],[ctx['RAW']],ctx),filtered_job(tool,ctx,cfg,platform,level)]
    if tool=='varscan2':
        pileup='samtools mpileup -f "$REF" "$BAM"' if platform=='10x' else 'samtools mpileup -Q 20 -q 20 -f "$REF" "$BAM"'
        options=('--min-coverage 2 --p-value 0.95 --min-var-freq 0.01' if level=='pooled' else '--min-coverage 10 --p-value 0.95 --min-var-freq 0.05') if platform=='10x' else '--min-coverage 10 --p-value 0.05 --min-var-freq 0.05'
        command='mkdir -p "$OUT"\n'+pileup+' | java -jar "$VARSCAN" mpileup2snp '+options+' --output-vcf 1 --strand-filter 0 --min-avg-qual 20 > "$RAW"'
        return [job('varscan2-call','varscan2',command,[ctx['BAM'],ctx['VARSCAN']],[ctx['RAW']],ctx),filtered_job(tool,ctx,cfg,platform,level)]
    raise ValueError('Unsupported ordinary caller '+tool)


def stage_context(cfg, sample, root, unit, tool):
    out=root/tool/unit['id']
    resources=cfg['resources']
    return dict(OUT=str(out),BAM=unit.get('bam',''),REF=required(resources,'reference'),
        RAW=str(out/('raw.vcf.gz' if tool=='gatk' else 'raw.vcf')),FINAL=str(out/'filtered.vcf.gz'),
        CHROMS=CHROMS,REGIONS=required(resources,'freebayes_regions') if tool=='freebayes' else '',VARSCAN=required(resources,'varscan_jar') if tool=='varscan2' else '',
        SAMPLE=sample['id'],UNIT=unit['id'])


def cellsnp_jobs(cfg,sample,root,units,platform,level):
    resources=cfg['resources'];ref=required(resources,'reference');chrmap=required(resources,'chr_map')
    jobs=[];finals=[]
    if platform=='smartseq' or level=='single-cell':
        out=root/'cellsnp';raw=out/'cellSNP.cells.vcf.gz'
        rows=sample.get('cells',[])
        bamlist='\n'.join(u['bam'] for u in rows)+'\n'
        celllist='\n'.join(u['id'] for u in rows)+'\n'
        if platform=='smartseq':
            cmd='mkdir -p "$OUT"\nprintf %s '+quote(bamlist)+' > "$OUT/bams.txt"\nprintf %s '+quote(celllist)+' > "$OUT/cells.txt"\n'
            cmd+='''cellsnp-lite -S "$OUT/bams.txt" -i "$OUT/cells.txt" -O "$OUT" -p 16 \\
  --chrom 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M \\
  --cellTAG None --UMItag None --minCOUNT 10 --minMAF 0.01 --refseq "$REF" --gzip --genotype'''
            jobs.append(job('cellsnp-joint-call','cellsnp',cmd,[u['bam'] for u in rows],[str(raw)],dict(OUT=str(out),REF=ref)))
            for unit in rows:
                dest=out/unit['id']/'filtered.vcf.gz'
                values=dict(RAW=str(raw),FINAL=str(dest),UNIT=unit['id'],CHR_MAP=chrmap)
                cmd='mkdir -p '+quote(dest.parent)+'\nbcftools view -s "$UNIT" "$RAW" -Ou | bcftools view -i '+quote(FILTERS['cellsnp']['smartseq'])+' -Ou | bcftools annotate --rename-chrs "$CHR_MAP" -Oz -o "$FINAL"\nbcftools index -f -t "$FINAL"'
                jobs.append(job('cellsnp-filter-'+unit['id'],'scmut',cmd,[str(raw),chrmap],[str(dest)],values));finals.append((unit,dest))
        else:
            # The barcode-aware call follows the source pooled-candidate step.
            out=root/'cellsnp'/sample['id']
            values=dict(OUT=str(out),REF=ref,BAM=required(sample,'pooled_bam'),CELL_BAM=required(sample,'cellsnp_barcode_bam'),BARCODES=required(sample,'barcodes'),CHR_MAP=chrmap)
            cmd='''mkdir -p "$OUT/pooled" "$OUT/cells"
cellsnp-lite -s "$BAM" -O "$OUT/pooled" -p 16 --minMAF 0.1 --minCOUNT 20 --refseq "$REF" --cellTAG None --UMItag UB --gzip --genotype
cellsnp-lite -s "$CELL_BAM" -b "$BARCODES" -O "$OUT/cells" -R "$OUT/pooled/cellSNP.cells.vcf.gz" -p 16 --minMAF 0.1 --minCOUNT 20 --refseq "$REF" --gzip --genotype'''
            raw=out/'cells/cellSNP.cells.vcf.gz'
            jobs.append(job('cellsnp-barcode-call','cellsnp',cmd,[values['BAM'],values['CELL_BAM'],values['BARCODES']],[str(raw)],values))
            body=source_text('cellsnp_single_filter')
            body=re.sub(r'^samples = .*$', 'samples = '+repr([sample['id']]),body,flags=re.M)
            body=re.sub(r'^CELL_BASE = .*$', 'CELL_BASE = '+repr(str(root/'cellsnp-native')),body,flags=re.M)
            body=re.sub(r'^OUT_BASE = .*$', 'OUT_BASE = '+repr(str(root/'cellsnp-native/counts')),body,flags=re.M)
            link=root/'cellsnp-native'/(sample['id']+'_1a')
            linkcmd='mkdir -p '+quote(link.parent)+'\nln -s '+quote(out/'cells')+' '+quote(link)
            jobs.append(job('cellsnp-native-layout','scmut',linkcmd,[str(raw)],[],{}))
            jobs.append(python_job('cellsnp-singlecell-filter','scmut',body,[],[str(raw)],
                [str(root/'cellsnp-native/counts'/(sample['id']+'_cellSNP_filtered_cell_variant_count.tsv'))]))
        return jobs,finals
    for unit in units:
        out=root/'cellsnp'/unit['id'];raw=out/'cellSNP.cells.vcf.gz';final=out/'filtered.vcf.gz'
        values=dict(OUT=str(out),BAM=unit['bam'],REF=ref,CHR_MAP=chrmap,RAW=str(raw),FINAL=str(final))
        call='''mkdir -p "$OUT"
cellsnp-lite -s "$BAM" -O "$OUT" -p 16 --minMAF 0.1 --minCOUNT 20 --refseq "$REF" --cellTAG None --UMItag UB --gzip --genotype'''
        if level=='cell-type':call=call.replace('--minMAF 0.1 --minCOUNT 20','--minMAF 0.05 --minCOUNT 10')
        filt='bcftools annotate --rename-chrs "$CHR_MAP" "$RAW" -Ou | '
        if level=='cell-type':filt+='bcftools norm -f "$REF" -m -both -Ou | '
        filt+='bcftools view -i \'INFO/DP>10 && GT!="0/0"\' -Oz -o "$FINAL"\nbcftools index -f -t "$FINAL"'
        jobs += [job('cellsnp-call-'+unit['id'],'cellsnp',call,[unit['bam']],[str(raw)],values),job('cellsnp-filter-'+unit['id'],'scmut',filt,[str(raw),chrmap],[str(final)],values)]
        finals.append((unit,final))
    return jobs,finals


def strelka_jobs(cfg,sample,root,units,platform,level):
    resources=cfg['resources'];binary=required(resources,'strelka');ref=required(resources,'reference')
    jobs=[];finals=[]
    batches=[[unit] for unit in units] if platform=='smartseq' else [units]
    for batch in batches:
        tag=batch[0]['id'] if platform=='smartseq' else 'joint'
        out=root/'strelka2'/tag
        bams=[u['bam'] for u in batch] if platform=='smartseq' else required(sample,'strelka_joint_bams_'+('pooled' if level=='pooled' else 'celltype'))
        command='mkdir -p '+quote(out.parent)+'\n'+quote(binary)+' '+ ' '.join('--bam '+quote(b) for b in bams)
        if platform=='smartseq':command+=' --rna'
        command+=' --referenceFasta '+quote(ref)+' --runDir '+quote(out)+'\n'+quote(out/'runWorkflow.py')+' -m local -j 16'
        rawpaths=[]
        for unit in batch:
            if platform=='smartseq':raw=out/'results/variants/variants.vcf.gz'
            else:
                label=required(unit,'strelka_vcf_sample')
                raw=out/'results/variants'/('genome.'+label+'.vcf.gz')
            rawpaths.append(str(raw))
        jobs.append(job('strelka-configure-and-call-'+tag,'strelka2',command,bams,rawpaths))
        for unit,raw in zip(batch,rawpaths):
            ctx=stage_context(cfg,sample,root,unit,'strelka2');ctx['RAW']=raw
            jobs.append(job('strelka-filter-dir-'+unit['id'],'scmut','mkdir -p "$OUT"',[],[],ctx))
            jobs.append(filtered_job('strelka2',ctx,cfg,platform,level));finals.append((unit,Path(ctx['FINAL'])))
    return jobs,finals


def monovar_jobs(cfg,sample,root,units,platform,level):
    binary=required(cfg['resources'],'monovar');ref=required(cfg['resources'],'reference')
    jobs=[];finals=[]
    batches=[[u] for u in units] if platform=='smartseq' else [units]
    for batch in batches:
        tag=batch[0]['id'] if platform=='smartseq' else 'joint'
        out=root/'monovar'/tag;raw=out/'raw.vcf'
        bams=[u['bam'] for u in batch] if platform=='smartseq' else required(sample,'monovar_joint_bams')
        text='\n'.join(bams)+'\n'
        command='mkdir -p '+quote(out)+'\nprintf %s '+quote(text)+' > '+quote(out/'bams.txt')+'\n'
        command+='samtools mpileup -BQ0 -d 10000 -f '+quote(ref)+' -q 20 -b '+quote(out/'bams.txt')+' | '+quote(binary)+' -p 0.002 -a 0.2 -t 0.05 -m 16 -f '+quote(ref)+' -b '+quote(out/'bams.txt')+' -o '+quote(raw)
        jobs.append(job('monovar-call-'+tag,'monovar',command,bams,[str(raw)]))
        for unit in batch:
            ctx=stage_context(cfg,sample,root,unit,'monovar');ctx['RAW']=str(raw)
            final=Path(ctx['FINAL']);jobs.append(job('monovar-output-'+unit['id'],'scmut','mkdir -p "$OUT"',[],[],ctx))
            label=required(unit,'monovar_vcf_sample') if platform=='10x' else None
            jobs.append(filtered_job('monovar',ctx,cfg,platform,level,label));finals.append((unit,final))
    return jobs,finals


def native_jobs(tool,cfg,sample,root,level):
    # Full special-caller code is embedded so the wrapper does not replace the native method.
    resources=cfg['resources'];software=required(resources,'software_root');ref=required(resources,'reference')
    if tool=='scsnv':
        native_root=root/'scsnv';out=native_root/sample['id']
        body=source_text('scsnv_call')
        body=re.sub(r'^samples=\([^\n]*\)', 'samples=("$SCMUT_SAMPLE")',body,flags=re.M)
        body=body.replace('/home/tantingting/scmut',software)
        body=re.sub(r'^fastq_dir=.*$', 'fastq_dir="$SCMUT_FASTQ"',body,flags=re.M)
        body=re.sub(r'^out_root=.*$', 'out_root="$SCMUT_ROOT"',body,flags=re.M)
        body=re.sub(r'^\s*fastq_dir=.*$', '    fastq_dir="$SCMUT_FASTQ"',body,flags=re.M)
        body=body.replace('memusg -t ', '')
        body='mkdir -p "$SCMUT_ROOT/$SCMUT_SAMPLE"\n'+body
        values=dict(SCMUT_SAMPLE=sample['id'],SCMUT_FASTQ=required(sample,'fastq_dir'),SCMUT_ROOT=str(native_root))
        jobs=[job('scsnv-from-fastq','scsnv',body,[values['SCMUT_FASTQ']],[str(out/'snvs.csv')],values)]
        filtering=source_text('scsnv_filter').split('# Calculate the updated intersection.')[0]
        filtering=re.sub(r'^samples=\([^\n]*\)', 'samples=("$SCMUT_SAMPLE")',filtering,flags=re.M)
        filtering=filtering.replace('sample="${s}_cSCC"','sample="$SCMUT_SAMPLE"')
        filtering=re.sub(r'^(root|BASE)=.*$',lambda m:m[1]+'="$SCMUT_ROOT"',filtering,flags=re.M)
        final=out/'snv.dp10_af005.pos1.vcf.gz'
        jobs.append(job('scsnv-site-filter','scmut',filtering,[str(out/'snvs.csv')],[str(final)],values))
        if level=='single-cell':
            extract=source_text('scsnv_pairs')
            args=['--sample',sample['id'],'--indir',str(out),'--outdir',str(out/'cell_snv')]
            pairs=out/'cell_snv'/(sample['id']+'_scSNV_filtered_cell_mut.tsv.gz')
            jobs.append(python_job('scsnv-cell-pairs','scsnv',extract,args,[str(out/'snv.dp10_af005.filtered.csv')],[str(pairs)]))
            body=source_text('scsnv_valid')
            body=re.sub(r'^BASE_DIR = .*$', 'BASE_DIR = '+repr(str(native_root)),body,flags=re.M)
            body=re.sub(r'^samples = .*$', 'samples = '+repr([sample['id']]),body,flags=re.M)
            body=re.sub(r'^\s*reads_file = .*$', '    reads_file = '+repr(required(sample,'reads_csv')),body,flags=re.M)
            jobs.append(python_job('scsnv-valid-cells-min2','scsnv',body,[],[str(pairs),sample['reads_csv']],
                [str(out/'cell_snv'/(sample['id']+'_scSNV_cell_mut_validcell_min2c.tsv.gz'))]))
            return jobs,[]
        return jobs,[({'id':sample['id']},final)]
    if tool=='monopogen':
        native_root=root/'monopogen';out=native_root/sample['id']
        path=required(resources,'monopogen');panel=required(resources,'panel');regions=required(sample,'monopogen_regions')
        bam=required(sample,'pooled_bam');bamlist=required(sample,'monopogen_bamlist')
        cmd='mkdir -p '+quote(out)+'\n'
        cmd+='export LD_LIBRARY_PATH='+quote(str(Path(path).parent.parent/'apps'))+':${LD_LIBRARY_PATH:-}\n'
        app=required(resources,'monopogen_apps')
        cmd+='python '+quote(path)+' preProcess -b '+quote(bamlist)+' -o '+quote(out)+' -a '+quote(app)+' -t 16\n'
        cmd+='python '+quote(path)+' germline -a '+quote(app)+' -t 16 -s all -p '+quote(panel)+' -r '+quote(regions)+' -g '+quote(ref)+' -o '+quote(out)
        jobs=[job('monopogen-preprocess-germline','monopogen',cmd,[bam,bamlist,regions],[])]
        if level=='pooled':
            body=source_text('monopogen_filter').split('# Calculate the intersection.')[0]
            body=re.sub(r'^samples=\([^\n]*\)', 'samples=("$SCMUT_SAMPLE")',body,flags=re.M)
            body=re.sub(r'^MONO_BASE=.*$', 'MONO_BASE="$SCMUT_MONO"',body,flags=re.M)
            body=re.sub(r'^REF=.*$', 'REF="$SCMUT_REF"',body,flags=re.M)
            values=dict(SCMUT_SAMPLE=sample['id'],SCMUT_MONO=str(native_root),SCMUT_REF=ref)
            final=out/'germline'/(sample['id']+'.phased.depth_filtered.vcf.gz')
            jobs.append(job('monopogen-depth-filter','scmut',body,[],[str(final)],values))
            return jobs,[({'id':sample['id']},final)]
        barcode=required(sample,'barcode_reads')
        command=''
        for stage in ['featureInfo','cellScan','LDrefinement']:
            command+='python '+quote(path)+' somatic -a '+quote(app)+' -t 16 -s '+stage+' -k 1 -r '+quote(regions)+' -i '+quote(out)+' -l '+quote(barcode)+' -g '+quote(ref)+'\n'
        jobs.append(job('monopogen-somatic','monopogen',command,[barcode],[]))
        body=source_text('monopogen_single_filter')
        body=re.sub(r'^MONO_BASE <- .*$', 'MONO_BASE <- '+json.dumps(str(native_root)),body,flags=re.M)
        body=re.sub(r'^samples <- .*$', 'samples <- c('+json.dumps(sample['id'])+')',body,flags=re.M)
        command="Rscript - <<'SCMUT_R'\n"+body+"\nSCMUT_R"
        result=native_root/'all_samples_somatic_combine/All_samples_Monopogen_filtered_cell_snv_long.tsv'
        jobs.append(job('monopogen-cell-filter','plotting',command,[],[str(result)]))
        return jobs,[]
    raise ValueError(tool)


def scallele_jobs(cfg,sample,root,units,platform):
    ref=required(cfg['resources'],'reference');jobs=[];finals=[]
    out=root/'scallele'/sample['id'] if platform=='10x' else root/'scallele'
    for unit in units:
        prefix=out/(sample['id']+'_map_'+unit['id'] if platform=='10x' else unit['id'])
        cmd='mkdir -p '+quote(out)+'\nscAllele --AB 0.05 --min-map_quality 20 '
        if platform=='smartseq':cmd+='--DP 10 -n 16 '
        cmd+='-b '+quote(unit['bam'])+' -g '+quote(ref)+' -o '+quote(prefix)
        jobs.append(job('scallele-call-'+unit['id'],'scallele',cmd,[unit['bam']],[str(prefix)+'.vcf']))
    if platform=='10x':
        # Retain the original header repair, PASS/autosome filtering and recurrence calculation.
        body=source_text('scallele_single_filter')
        body=re.sub(r'^DIR=.*$', 'DIR="$SCMUT_DIR"',body,flags=re.M)
        body=re.sub(r'^REF=.*$', 'REF="$SCMUT_REF"',body,flags=re.M)
        body=re.sub(r'^samples=\([^\n]*\)', 'samples=("$SCMUT_SAMPLE")',body,flags=re.M)
        # The original filter expects per-sample raw VCFs below DIR/sample.
        values=dict(SCMUT_DIR=str(root/'scallele'),SCMUT_REF=ref,SCMUT_SAMPLE=sample['id'])
        body=body.replace('cell=${base#${s}_map_}', 'cell=${base#'+sample['id']+'_map_}')
        jobs.append(job('scallele-pass-filter','scmut',body,[],[str(root/'scallele/scAllele_all_samples_cell_snv_long.tsv')],values))
        rbody=source_text('scallele_recurrence')
        rbody=re.sub(r'^DIR <- .*$', 'DIR <- '+json.dumps(str(root/'scallele')),rbody,flags=re.M)
        jobs.append(job('scallele-min2-cells','plotting',"Rscript - <<'SCMUT_R'\n"+rbody+"\nSCMUT_R",[],[str(root/'scallele/scAllele_all_samples_cell_snv_long_min2cells.tsv')]))
    else:
        # Use the original Smart-seq2 header-repair/filter block with a private patient manifest.
        body=source_text('scallele_smart_filter')
        body=body.replace('/mnt/storage/tantingting/scmut_smart',str(root/'native'))
        body=body.replace('/home/tantingting/scmut/ref/bwa_ref/genome.fa',ref)
        body=re.sub(r'^MANIFEST=.*$', 'MANIFEST="$SCMUT_MANIFEST"',body,flags=re.M)
        manifest=root/'patient_manifest.tsv'
        cohort=required(sample,'cohort')
        target=root/'native/vcf'/(cohort+'_scRNA')/'scallele'
        setup='mkdir -p '+quote(target)+' '+quote(root/'native/vcf/1_TNBC_scRNA/scallele')+' '+quote(root/'native/vcf/2_BRCA_scRNA/scallele')+' '+quote(root/'native/vcf/filter_stat')+'\n'
        for unit in units:
            setup+='ln -s '+quote(out/(unit['id']+'.vcf'))+' '+quote(target/(unit['id']+'.vcf'))+'\n'
        text='scRNA_SRR\tSample\tWES_SRR\n'+''.join(u['id']+'\t'+sample['id']+'\t'+sample.get('wes_id','WES')+'\n' for u in units)
        setup+='printf %s '+quote(text)+' > '+quote(manifest)+'\n'
        outputs=[str(target/'PASS_norm_chr1_22'/(u['id']+'_PASS.vcf.gz')) for u in units]
        jobs.append(job('scallele-pass-filter','scmut',setup+body,[],outputs,dict(SCMUT_MANIFEST=str(manifest))))
        finals=[(u,Path(p)) for u,p in zip(units,outputs)]
    return jobs,finals


def scomatic_jobs(cfg,sample,root):
    ref=required(cfg['resources'],'reference');install=required(cfg['resources'],'scomatic')
    bam=required(sample,'pooled_bam');meta=required(sample,'cell_annotations');out=root/'scomatic'
    values=dict(OUT=str(out),REF=ref,BAM=bam,META=meta,SCOMATIC=install,SAMPLE=sample['id'])
    cmd='''mkdir -p "$OUT/counts" "$OUT/merged" "$OUT/calls"
python "$SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py" --bam "$BAM" --meta "$META" --id "$SAMPLE" --max_nM 5 --max_NH 1 --outdir "$OUT"
for bam in "$OUT"/*.bam; do
  cell_type=$(basename "$bam" | awk -F. '{print $(NF-1)}')
  mkdir -p "$OUT/counts/tmp_$cell_type"
  python "$SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py" --bam "$bam" --ref "$REF" --chrom all --out_folder "$OUT/counts" --min_bq 20 --tmp_dir "$OUT/counts/tmp_$cell_type" --nprocs 16 --min_mq 255 --min_af 0.05 --min_ac 2
done
python "$SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py" --tsv_folder "$OUT/counts" --outfile "$OUT/merged/$SAMPLE.tsv"
python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py" --infile "$OUT/merged/$SAMPLE.tsv" --outfile "$OUT/calls/$SAMPLE" --ref "$REF" --min_cell_types 1
python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py" --infile "$OUT/calls/$SAMPLE.calling.step1.tsv" --outfile "$OUT/calls/$SAMPLE" --editing "$SCOMATIC/RNAediting/AllEditingSites.hg38.txt" --pon "$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv"'''
    filtered=str(out/'calls'/required(sample,'scomatic_step2_filename'))
    jobs=[job('scomatic-call','scomatic',cmd,[bam,meta],[filtered],values)]
    # The caller's step2 table suffix depends on the source release and is explicit in the sample configuration.
    body=source_text('scomatic_filter')
    output=root/'scomatic/filtered'
    command='mkdir -p '+quote(output)+'\ncd '+quote(output)+'\npython - '+quote(filtered)+" <<'SCMUT_PY'\n"+body+'\nSCMUT_PY'
    jobs.append(job('scomatic-table-filter','scomatic',command,[],[str(output/'Epithelial.tsv'),str(output/'Myeloid.tsv'),str(output/'T.tsv')]))
    return jobs,[]


def evaluation_jobs(cfg,sample,root,tool,unit,final):
    truth=unit.get('truth_vcf',sample.get('truth_vcf'))
    bed=unit.get('callable_bed',sample.get('callable_bed'))
    if not truth or not bed:
        raise ValueError('Evaluation requires truth_vcf and callable_bed for '+unit['id'])
    sdf=required(cfg['resources'],'rtg_sdf');out=root/tool/unit['id']/'rtg'
    values=dict(TRUTH=truth,CALLS=str(final),BED=bed,SDF=sdf,RTG_OUT=str(out))
    cmd='''if [ -e "$RTG_OUT" ]; then echo "RTG output already exists: $RTG_OUT" >&2; exit 1; fi
rtg vcfeval -b "$TRUTH" -c "$CALLS" --bed-regions "$BED" -t "$SDF" -o "$RTG_OUT" --squash-ploidy --output-mode=annotate'''
    return job(tool+'-rtg-'+unit['id'],'rtg',cmd,[truth,bed,str(final),sdf],[str(out/'summary.txt')],values)


def plan(cfg,sample_id,selection,evaluate):
    sample=dict(required(cfg['samples'],sample_id));sample['id']=validate_id(sample_id)
    root=Path(required(cfg,'output_root')).expanduser().resolve()/PLATFORM/sample_id/LEVEL
    if LEVEL=='pooled':
        units=[dict(id=sample_id,bam=required(sample,'pooled_bam'),**{k:v for k,v in sample.items() if k.endswith('_vcf_sample')})]
    else:
        units=sample.get('cell_types' if LEVEL=='cell-type' else 'cells',[])
    for unit in units:validate_id(unit['id'])
    allowed=TOOLS[(PLATFORM,LEVEL)]
    chosen=allowed if selection=='all' else [allowed[int(t)-1] if t.isdigit() and 0<int(t)<=len(allowed) else t.lower() for t in selection.split(',')]
    if len(set(chosen))!=len(chosen) or any(t not in allowed for t in chosen):
        raise ValueError('Choose unique tools from '+', '.join(allowed))
    if not units and any(t not in ('scsnv','monopogen','scomatic') for t in chosen):
        raise ValueError('Configure cell_types or cells for the selected level')
    jobs=[];results=[]
    for tool in chosen:
        if tool in ['scsnv','monopogen']:new,finals=native_jobs(tool,cfg,sample,root,LEVEL)
        elif tool=='cellsnp':new,finals=cellsnp_jobs(cfg,sample,root,units,PLATFORM,LEVEL)
        elif tool=='strelka2':new,finals=strelka_jobs(cfg,sample,root,units,PLATFORM,LEVEL)
        elif tool=='monovar':new,finals=monovar_jobs(cfg,sample,root,units,PLATFORM,LEVEL)
        elif tool=='scallele':new,finals=scallele_jobs(cfg,sample,root,units,PLATFORM)
        elif tool=='scomatic':new,finals=scomatic_jobs(cfg,sample,root)
        else:
            new=[];finals=[]
            for unit in units:
                ctx=stage_context(cfg,sample,root,unit,tool)
                new.extend(regular_jobs(tool,ctx,cfg,PLATFORM,LEVEL));finals.append((unit,Path(ctx['FINAL'])))
        if evaluate and not finals:
            raise ValueError('--evaluate requires per-unit VCFs; '+tool+' produces native cell/site tables. Use the corresponding table-to-VCF evaluation code in Benchmark Evaluation.')
        jobs.extend(new)
        for unit,final in finals:
            results.append(dict(sample=sample_id,unit=unit['id'],tool=tool,vcf=str(final)))
            if evaluate:jobs.append(evaluation_jobs(cfg,sample,root,tool,unit,final))
    return jobs,root,results


def main():
    parser=argparse.ArgumentParser(description=PLATFORM+' '+LEVEL+' calling and filtering')
    parser.add_argument('--config',default='config/workflow.json')
    parser.add_argument('--sample')
    parser.add_argument('--tools',default='all',help='Comma-separated numbers or tool names')
    parser.add_argument('--list-tools',action='store_true')
    parser.add_argument('--evaluate',action='store_true')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    if args.list_tools:
        for i,t in enumerate(TOOLS[(PLATFORM,LEVEL)],1):print(str(i)+' '+t)
        return
    if not args.sample:parser.error('--sample is required')
    cfg=json.loads(Path(args.config).read_text())
    jobs,root,results=plan(cfg,args.sample,args.tools,args.evaluate)
    if jobs:jobs[0]['inputs'].append(required(cfg['resources'],'reference'))
    run_jobs(jobs,cfg,root,args.dry_run)
    if not args.dry_run:
        path=root/'variant_outputs.tsv';temporary=root/'variant_outputs.tsv.tmp'
        previous=[]
        if path.exists():
            with path.open() as stream:previous=list(csv.DictReader(stream,delimiter='\t'))
        updated={(r['sample'],r['unit'],r['tool']):r for r in previous}
        updated.update({(r['sample'],r['unit'],r['tool']):r for r in results})
        with temporary.open('w') as f:
            w=csv.DictWriter(f,fieldnames=['sample','unit','tool','vcf'],delimiter='\t');w.writeheader();w.writerows(updated.values())
        temporary.replace(path)
        print('Outputs: '+str(path))


if __name__=='__main__':
    try:main()
    except (ValueError,KeyError,RuntimeError,FileNotFoundError,subprocess.CalledProcessError,BlockingIOError) as error:
        print(str(error),file=sys.stderr)
        sys.exit(1)
