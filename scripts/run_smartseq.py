#!/usr/bin/env python3
"""Run source-derived calling and filtering stages for one sample."""
PLATFORM='smartseq'
LEVEL='single-cell'
FILTERS={'samtools': {'10x_pooled': 'FMT/GT!="0/0" && QUAL>20 && ( (INFO/DP4[2]+INFO/DP4[3]) / (INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3]) ) > 0.05 && ( (TYPE="snp" && INFO/DP>10) || (TYPE!="snp" && INFO/DP>15) )', '10x_celltype': 'FMT/GT!="0/0" && QUAL>20 && ( (INFO/DP4[2]+INFO/DP4[3]) / (INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3]) ) > 0.05 && ((TYPE="snp" && INFO/DP>10) || (TYPE!="snp" && INFO/DP>15))', 'smartseq': 'FMT/GT!="0/0" && QUAL>20 && (INFO/DP4[2]+INFO/DP4[3])/(INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3])>0.05 && ( (TYPE="snp" && INFO/DP>10) || (TYPE!="snp" && INFO/DP>15) )'}, 'freebayes': {'10x_pooled': 'FMT/GT!="0/0" && QUAL>20 && INFO/AF[0]>=0.05 && ( (TYPE="snp" && INFO/DP[0]>=10) || (TYPE!="snp" && INFO/DP[0]>=15) )', '10x_celltype': 'FMT/GT!="0/0" && QUAL>20 && INFO/AF[0]>=0.05 && ( (TYPE="snp" && INFO/DP[0]>10) || (TYPE!="snp" && INFO/DP[0]>15) )', 'smartseq': 'FMT/GT!="0/0" && QUAL>20 && INFO/AF[0]>0.05 && ( (TYPE="snp" && INFO/DP[0]>10) || (TYPE!="snp" && INFO/DP[0]>15) )'}, 'varscan2': {'10x_pooled': 'FMT/GT!="0/0" && FMT/DP>10 && (FMT/AD / (FMT/RD + FMT/AD)) > 0.05 && ( (TYPE="snp" && FMT/DP>10) || (TYPE!="snp" && FMT/DP>15) )', '10x_celltype': 'FMT/GT!="0/0" && (FMT/AD / (FMT/RD + FMT/AD)) > 0.05 && ( (TYPE="snp" && FMT/DP>10) || (TYPE!="snp" && FMT/DP>15) )', 'smartseq': 'FMT/GT!="0/0" && FMT/DP>10 && (FMT/AD / (FMT/RD + FMT/AD)) > 0.05'}, 'strelka2': {'10x_pooled': 'ALT!="." && FMT/GT!="0/0" && FMT/GT!="0|0" && QUAL>20 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP > 10) || (TYPE!="snp" && FMT/DPI > 15) )', '10x_celltype': 'ALT!="." && FMT/GT!="0/0" && FMT/GT!="0|0" && QUAL>20 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP > 10) || (TYPE!="snp" && FMT/DPI > 15) )', 'smartseq': 'ALT!="." && FMT/GT!="0/0" && FMT/GT!="0|0" && QUAL>20 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP > 10) || (TYPE!="snp" && FMT/DPI > 15) )'}, 'monovar': {'10x_pooled': 'FILTER="PASS" && FMT/GT!="0/0" && FMT/GT!="./." && (FMT/AD[0:0] + FMT/AD[0:1]) > 0 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP>10) || (TYPE!="snp" && FMT/DP>15) )', '10x_celltype': 'FILTER="PASS" && FMT/GT!="0/0" && FMT/GT!="./." && (FMT/AD[0:0] + FMT/AD[0:1]) > 0 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP>10) || (TYPE!="snp" && FMT/DP>15) )', 'smartseq': 'FMT/GT!="0/0" && FMT/GT!="./." && (FMT/AD[0:0] + FMT/AD[0:1]) > 0 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) > 0.05 && ( (TYPE="snp" && FMT/DP>10) || (TYPE!="snp" && FMT/DP>15) )'}, 'cellsnp': {'smartseq': 'FMT/GT!="0/0" && FMT/GT!="0|0" && FMT/GT!="." && FMT/GT!="./." && FMT/GT!=".|." && FMT/DP>10 && FMT/AD/FMT/DP>0.05'}}
SOURCES={}
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
