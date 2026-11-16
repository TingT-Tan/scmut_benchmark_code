#!/usr/bin/env python3
KIND='consensus'
NATIVE={}
NATIVE['196']=r'''
#!/bin/bash
combos=(
"Monopogen,Strelka2"
"Monopogen,Freebayes"
"Monopogen,GATK"
"Monopogen,scSNV"
"Monopogen,Varscan2"
"Strelka2,Freebayes"
"Strelka2,scSNV"
"Strelka2,Varscan2"
"Monopogen,Strelka2,Freebayes"
"Monopogen,Strelka2,GATK"
"Monopogen,Strelka2,scSNV"
"Monopogen,Strelka2,Varscan2"
"Monopogen,Freebayes,GATK"
"Monopogen,Freebayes,scSNV"
"Monopogen,Freebayes,Varscan2"
"Strelka2,Freebayes,GATK"
"Monopogen,Strelka2,Freebayes,GATK"
"Monopogen,Strelka2,Freebayes,scSNV"
"Monopogen,Strelka2,Freebayes,Varscan2"
) 


combos=(
"Strelka2,GATK"
"Freebayes,GATK"
"Freebayes,scSNV"
"Freebayes,Varscan2"
"Strelka2,Freebayes,scSNV"
"Strelka2,Freebayes,Varscan2"
"Monopogen,Strelka2,GATK,scSNV"
"Monopogen,Strelka2,GATK,Varscan2"
"Monopogen,Strelka2,scSNV,Varscan2"
"Monopogen,Freebayes,GATK,scSNV"
"Monopogen,Freebayes,GATK,Varscan2"
"Monopogen,Freebayes,scSNV,Varscan2"
"Strelka2,Freebayes,GATK,scSNV"
"Strelka2,Freebayes,GATK,Varscan2"
"Strelka2,Freebayes,scSNV,Varscan2"
)

POOL_BASE="/mnt/storage/tantingting/scmut/isec/pool_sites"
OUT_BASE="/mnt/storage/tantingting/scmut/isec/support_consensus_D1D2"
REF_SDF="/home/tantingting/scmut/ref/bwa_ref/genome.sdf"
mkdir -p "$OUT_BASE"

samples=("P4_cSCC" "P5_cSCC" "P6_cSCC" "P7_cSCC" "P9_cSCC" "P10_cSCC" "Acral1" "Acral2" "Acral3pre" "Acral3post" "Acral4" "Cutaneous1" "Cutaneous1Lym" "Cutaneous2")

get_truth_vcf() {
    local sample="$1"
    if [[ "$sample" == P*_cSCC ]]; then
        echo "/mnt/storage/tantingting/scmut/vcf/1_cSCC_WES/hap22/${sample}_hap_PASS3.vcf.gz"
    else
        case "$sample" in
            Acral1) echo "/mnt/storage/tantingting/scmut/vcf/2_melanoma_WES/hap/filter/SRR21849532_PASS.vcf.gz" ;;
            Acral2) echo "/mnt/storage/tantingting/scmut/vcf/2_melanoma_WES/hap/filter/SRR21849531_PASS.vcf.gz" ;;
            Acral3pre) echo "/mnt/storage/tantingting/scmut/vcf/2_melanoma_WES/hap/filter/SRR21849530_PASS.vcf.gz" ;;
            Acral3post) echo "/mnt/storage/tantingting/scmut/vcf/2_melanoma_WES/hap/filter/SRR21849529_PASS.vcf.gz" ;;
            Acral4) echo "/mnt/storage/tantingting/scmut/vcf/2_melanoma_WES/hap/filter/SRR21849528_PASS.vcf.gz" ;;
            Cutaneous1) echo "/mnt/storage/tantingting/scmut/vcf/2_melanoma_WES/hap/filter/SRR21849527_PASS.vcf.gz" ;;
            Cutaneous2) echo "/mnt/storage/tantingting/scmut/vcf/2_melanoma_WES/hap/filter/SRR21849526_PASS.vcf.gz" ;;
            Cutaneous1Lym) echo "/mnt/storage/tantingting/scmut/vcf/2_melanoma_WES/hap/filter/SRR21849525_PASS.vcf.gz" ;;
        esac
    fi
}

get_bed() {
    local sample="$1"
    if [[ "$sample" == P*_cSCC ]]; then
        local p=${sample%%_*}
        echo "/mnt/storage/tantingting/scmut/WES_gold/1_cSCC/region/whole/${p}_WES30_scRNA10_Shared_Genome.bed"
    else
        echo "/mnt/storage/tantingting/scmut/WES_gold/2_melanoma/region/whole/${sample}_Shared_Genome.bed"
    fi
}

make_vcf_from_keys() {
    local key_file="$1"
    local out_vcf="$2"
    {
        echo "##fileformat=VCFv4.2"
        for i in $(seq 1 22); do echo "##contig=<ID=chr${i}>"; done
        echo "##contig=<ID=chrX>"
        echo "##contig=<ID=chrY>"
        echo "##contig=<ID=chrM>"
        echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCONSENSUS"
        awk -F':' 'BEGIN{OFS="\t"}{print $1,$2,".",$3,$4,".","PASS",".","GT","0/1"}' "$key_file"
    } > "$out_vcf"
    bgzip -f "$out_vcf"
    tabix -f -p vcf "${out_vcf}.gz"
}

for sample in "${samples[@]}"; do
    echo "=============================="
    echo ">>> Processing $sample"
    echo "=============================="

    sample_out="${OUT_BASE}/${sample}"
    mkdir -p "$sample_out/tmp" "$sample_out/logs"

    truth_vcf=$(get_truth_vcf "$sample")
    bed_file=$(get_bed "$sample")

    [[ ! -f "$truth_vcf" ]] && { echo "[WARN] missing WES: $truth_vcf"; continue; }
    [[ ! -f "$bed_file" ]] && { echo "[WARN] missing BED: $bed_file"; continue; }

    for combo in "${combos[@]}"; do
        
        IFS=',' read -ra tools <<< "$combo"
        n_tools=${#tools[@]}
        combo_name=$(echo "$combo" | tr ',' '+')
        combo_tag=$(echo "$combo" | tr ',' '_')

        merged="${sample_out}/tmp/${sample}.${combo_tag}.merged.tsv"
        support_table="${sample_out}/${sample}.${combo_tag}.support_count.tsv"
        : > "$merged"

        miss=0
        for tool in "${tools[@]}"; do
            key_file="${POOL_BASE}/${sample}_${tool}/query.keys"
            if [[ ! -f "$key_file" ]]; then
                echo "[WARN] missing key: $sample $tool"
                miss=1
                break
            fi
            awk -v T="$tool" 'BEGIN{FS=":";OFS="\t"}{print $1,$2,$3,$4,T}' "$key_file" >> "$merged"
        done
        [[ "$miss" -eq 1 ]] && continue

        awk 'BEGIN{FS=OFS="\t"}{
            key=$1"\t"$2"\t"$3"\t"$4
            if(!(seen[key,$5]++)){
                count[key]++
                toolset[key]=(toolset[key]==""?$5:toolset[key]","$5)
            }
        }END{
            for(k in count) print k,count[k],toolset[k]
        }' "$merged" | sort -k1,1V -k2,2n > "$support_table"

        if [[ "$n_tools" -eq 2 ]]; then
            levels=("ge1" "eq2")
        elif [[ "$n_tools" -eq 3 ]]; then
            levels=("ge1" "ge2" "eq3")
        else
            levels=("ge1" "ge2" "ge3" "eq4")
        fi

        for level in "${levels[@]}"; do
            case "$level" in
                ge1) expr='$5>=1' ;;
                ge2) expr='$5>=2' ;;
                ge3) expr='$5>=3' ;;
                eq2) expr='$5==2' ;;
                eq3) expr='$5==3' ;;
                eq4) expr='$5==4' ;;
            esac

            key_txt="${sample_out}/${sample}.${combo_tag}.${level}.keys.txt"
            out_vcf="${sample_out}/${sample}.${combo_tag}.${level}.vcf"
            rtg_out="${sample_out}/rtg_${combo_tag}_${level}"

            awk 'BEGIN{FS=OFS="\t"} '"$expr"' {print $1":"$2":"$3":"$4}' "$support_table" > "$key_txt"
            n_sites=$(wc -l < "$key_txt")
            echo "[INFO] $sample | $combo_name | $level | $n_sites sites"

            [[ "$n_sites" -eq 0 ]] && continue

            make_vcf_from_keys "$key_txt" "$out_vcf"

            rm -rf "$rtg_out"

            rtg vcfeval -b "$truth_vcf" -c "${out_vcf}.gz" -t "$REF_SDF" --bed-regions "$bed_file" \
                -o "$rtg_out" --squash-ploidy --output-mode=annotate \
                > "${sample_out}/logs/${combo_tag}_${level}.log" 2>&1
        done
    done

    rm -rf "${sample_out}/tmp"
done

echo "All done."
'''
NATIVE['197']=r'''
BASE="/mnt/storage/tantingting/scmut/isec/support_consensus_D1D2"
OUT="${BASE}/rtg_summary_all.tsv"

echo -e "Sample\tCombo\tLevel\tTP\tFP\tFN\tPrecision\tSensitivity\tF_measure" > "$OUT"

for summary in ${BASE}/*/rtg_*/summary.txt; do
    [[ -f "$summary" ]] || continue
    sample=$(basename "$(dirname "$(dirname "$summary")")")
    name=$(basename "$(dirname "$summary")" | sed 's/^rtg_//')
    level=${name##*_}
    combo=${name%_*}
    line=$(grep -E '^[[:space:]]*None[[:space:]]+' "$summary" | tail -1)
    [[ -z "$line" ]] && continue
    TP=$(echo "$line" | awk '{print $3}')
    FP=$(echo "$line" | awk '{print $4}')
    FN=$(echo "$line" | awk '{print $5}')
    P=$(echo "$line" | awk '{print $6}')
    R=$(echo "$line" | awk '{print $7}')
    F1=$(echo "$line" | awk '{print $8}')
    echo -e "${sample}\t${combo}\t${level}\t${TP}\t${FP}\t${FN}\t${P}\t${R}\t${F1}" >> "$OUT"
done

echo "Done: $OUT"
'''
NATIVE['206']=r'''
#!/bin/bash
BASE=/mnt/storage/tantingting/scmut_smart
MANIFEST=${BASE}/bam/star_qc_unique_ge10k_scrna_sample_wes.tsv
POOL_BASE=/mnt/storage/tantingting/scmut_smart/vcf/site_D1D2
OUT_BASE=${BASE}/isec/support_consensus
REF_SDF=/home/tantingting/scmut/ref/bwa_ref/genome.sdf
RESULT=${OUT_BASE}/consensus_rtg_summary_rerun1.tsv
mkdir -p ${OUT_BASE}

combos=(
"monovar,freebayes"
"monovar,samtools"
"monovar,strelka2"
"monovar,varscan2"
"monovar,scallele"
"monovar,hap"
)

combos=(
"monovar,cellsnp"
"freebayes,samtools"
"freebayes,strelka2"
"freebayes,varscan2"
"samtools,strelka2"
"strelka2,varscan2"
)

combos=(
"monovar,freebayes,samtools"
"monovar,freebayes,strelka2"
"monovar,freebayes,varscan2"
"monovar,strelka2,varscan2"
"freebayes,samtools,strelka2"
"freebayes,samtools,varscan2"
"freebayes,strelka2,varscan2"
)

combos=(
"monovar,freebayes,samtools,strelka2"
"monovar,freebayes,samtools,varscan2"
"monovar,freebayes,strelka2,varscan2"
"monovar,samtools,strelka2,varscan2"
"freebayes,samtools,strelka2,varscan2"
)

echo -e "scRNA_SRR\tsample\tWES_SRR\tcombo\tlevel\ttruth_n\trecall" > ${RESULT}

make_vcf(){
	key_file=$1
	out_vcf=$2
	{
		echo "##fileformat=VCFv4.2"
		for i in $(seq 1 22); do echo "##contig=<ID=chr${i}>"; done
		echo "##contig=<ID=chrX>"
		echo "##contig=<ID=chrY>"
		echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
		echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCONSENSUS"
		awk -F':' 'BEGIN{OFS="\t"}{print $1,$2,".",$3,$4,".","PASS",".","GT","0/1"}' ${key_file}
	} > ${out_vcf}
	bgzip -f ${out_vcf}
	tabix -f -p vcf ${out_vcf}.gz
}

while read -r scRNA_SRR sample WES_SRR
do
	[[ "${scRNA_SRR}" == "scRNA_SRR" ]] && continue
	echo "===== ${scRNA_SRR} | ${sample} | ${WES_SRR} ====="

	if [[ ${sample} == PT* ]]; then
		DATASET=1_TNBC
		COHORT=TNBC
		truth_vcf=${BASE}/vcf/${DATASET}_WES/hap/filter/${WES_SRR}_PASS.vcf.gz
	elif [[ ${sample} == BC* ]]; then
		DATASET=2_BRCA
		COHORT=BRCA
		truth_vcf=${BASE}/vcf/${DATASET}_WES/hap/filter/${sample}_Tumor_${WES_SRR}_PASS.vcf.gz
	else
		echo "[WARN] unknown ${sample}"
		continue
	fi

	BED=${BASE}/WES_gold/smartseq_callable/${COHORT}/shared/${scRNA_SRR}_${WES_SRR}_WES30_scRNA1_shared.bed

	[[ ! -f ${truth_vcf} ]] && {
		echo "[WARN] missing WES ${truth_vcf}"
		continue
	}

	[[ ! -f ${BED} ]] && {
		echo "[WARN] missing BED ${BED}"
		continue
	}

	cell_out=${OUT_BASE}/${scRNA_SRR}
	mkdir -p ${cell_out}/tmp ${cell_out}/logs

	for combo in "${combos[@]}"
	do

	IFS=',' read -ra combo_tools <<< "${combo}"
	combo_tag=$(echo ${combo} | tr ',' '_')


	n=${#combo_tools[@]}

	if [[ ${n} -eq 2 ]]; then
	levels=("ge1" "eq2")
	elif [[ ${n} -eq 3 ]]; then
	levels=("ge1" "ge2" "eq3")
	elif [[ ${n} -eq 4 ]]; then
	levels=("ge1" "ge2" "ge3" "eq4")
	fi


	combo_complete=1
	for level in "${levels[@]}"
	do
		rtg_out=${cell_out}/rtg_${combo_tag}_${level}
		if [[ ! -s ${rtg_out}/summary.txt ]]; then
			combo_complete=0
			break
		fi
	done

	if [[ ${combo_complete} -eq 1 ]]; then
		echo "[SKIP COMBO] ${scRNA_SRR} ${combo} all levels completed"
		continue
	fi


	merged=${cell_out}/tmp/${combo_tag}.merged.tsv
	support=${cell_out}/${combo_tag}.support.tsv

	> ${merged}

	missing=0

	for tool in "${combo_tools[@]}"
	do
		key=${POOL_BASE}/${DATASET}_${scRNA_SRR}_${tool}.query.keys

		if [[ ! -f ${key} ]]; then
		echo "[WARN] missing ${key}"
		missing=1
		break
		fi

		awk -v T=${tool} '
		BEGIN{FS=":";OFS="\t"}
		{
		print $1,$2,$3,$4,T
		}' ${key} >> ${merged}
	done

	[[ ${missing} -eq 1 ]] && continue

	awk '
	BEGIN{FS=OFS="\t"}
	{
	k=$1"\t"$2"\t"$3"\t"$4
	if(!seen[k,$5]++){
	count[k]++
	tools[k]=(tools[k]==""?$5:tools[k]","$5)
	}
	}
	END{
	for(k in count)
	print k,count[k],tools[k]
	}
	' ${merged} | sort -k1,1V -k2,2n > ${support}

#####level
	for level in "${levels[@]}"
	do

	rtg_out=${cell_out}/rtg_${combo_tag}_${level}
	summary_file=${rtg_out}/summary.txt
	##################################################
	
	##################################################

	if [[ -s ${summary_file} ]]; then
		echo "[SKIP] ${scRNA_SRR} ${combo} ${level} already completed"
		continue
	fi
	echo "[RERUN] ${scRNA_SRR} ${combo} ${level}"

	case ${level} in
	ge1) expr='$5>=1';;
	ge2) expr='$5>=2';;
	ge3) expr='$5>=3';;
	eq2) expr='$5==2';;
	eq3) expr='$5==3';;
	eq4) expr='$5==4';;
	esac

	key_out=${cell_out}/tmp/${combo_tag}_${level}.keys
	out_vcf=${cell_out}/${combo_tag}_${level}.vcf

	awk 'BEGIN{FS=OFS="\t"}'"${expr}"'{print $1":"$2":"$3":"$4}' ${support} > ${key_out}

	n_sites=$(wc -l < ${key_out})

	if [[ ${n_sites} -eq 0 ]]; then
		echo "[EMPTY] ${scRNA_SRR} ${combo} ${level}: 0 consensus sites"
		continue
	fi

	make_vcf ${key_out} ${out_vcf}

	rm -rf ${rtg_out}

	rtg vcfeval -b ${truth_vcf} -c ${out_vcf}.gz --bed-regions ${BED} -t ${REF_SDF} -o ${rtg_out} \
	--squash-ploidy --output-mode=annotate > ${cell_out}/logs/${combo_tag}_${level}.log 2>&1

	if [[ -s ${rtg_out}/summary.txt ]]; then
		echo "[DONE] ${scRNA_SRR} ${combo} ${level}"
	else
		echo "[FAILED] ${scRNA_SRR} ${combo} ${level}"
	fi

done

done

#rm -rf ${cell_out}/tmp

done < ${MANIFEST}

echo "Done: ${RESULT}"
'''
NATIVE['207']=r'''
BASE=/mnt/storage/tantingting/scmut_smart/isec/support_consensus
OUT=${BASE}/consensus_rtg_all_long_0816.tsv

echo -e "Sample\tCombo\tLevel\ttruth_n\tTP\tFP\tFN\tPrecision\tSensitivity\tF_measure" > "$OUT"

for summary in ${BASE}/*/rtg_*/summary.txt; do

    [[ -f "$summary" ]] || continue

    
    rtg_dir=${summary%/summary.txt}

    
    rtg_name=${rtg_dir##*/}

    
    sample_dir=${rtg_dir%/*}
    sample=${sample_dir##*/}

    
    name=${rtg_name#rtg_}

    
    level=${name##*_}

    
    combo=${name%_*}

    
    awk -v sample="$sample" -v combo="$combo" -v level="$level" '
        /^[[:space:]]*None[[:space:]]+/ {
            TP=$3
            FP=$4
            FN=$5
            P=$6
            R=$7
            F1=$8
        }
        END {
            if (TP != "") {
                truth_n=TP+FN
                OFS="\t"
                print sample,combo,level,truth_n,TP,FP,FN,P,R,F1
            }
        }
        ' "$summary" >> "$OUT"

done

echo "Done: $OUT"
wc -l "$OUT"

# 93798 /mnt/storage/tantingting/scmut_smart/isec/support_consensus/consensus_rtg_all_long_0816.tsv
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


import argparse
import json
from pathlib import Path
import re
import shlex
import sys


def bind(text, key, value):
    pattern=r'^'+re.escape(key)+r'=[^\n]*$'
    text,n=re.subn(pattern,lambda m:key+'='+shlex.quote(str(value)),text,flags=re.M)
    if not n:raise ValueError('Source assignment not found: '+key)
    return text


def change_function(text, name, replacement):
    pattern=r'^'+name+r'\(\)\s*\{.*?^\}'
    text,n=re.subn(pattern,lambda m:replacement,text,flags=re.M|re.S)
    if n!=1:raise ValueError('Source function not found: '+name)
    return text


def simulation_plan(cfg,args):
    sample=cfg['samples'][args.sample];validate_id(args.sample)
    resources=cfg['resources'];settings=cfg.get('simulation',{})
    root=Path(cfg['output_root']).expanduser().resolve()/'simulation'/args.sample/args.mode
    if args.mode.startswith('smartseq'):
        if args.stage!='spike':raise ValueError('Use the original Smart-seq2 calling commands in Variant Calling on the simulated cell BAMs')
        key=args.mode.replace('-','_')
        region=required(settings,key+'_regions')
        cells=required(sample,'cells')
        if args.unit:
            cells=[u for u in cells if u['id']==args.unit]
            if not cells:raise ValueError('Unknown --unit in the selected sample')
        jobs=[]
        for cell in cells:
            validate_id(cell['id'])
            text=NATIVE[key]
            for name,value in [('REF',required(resources,'reference')),('OUTDIR',root/cell['id']),('REGION_DIR',region)]:text=bind(text,name,value)
            text=change_function(text,'get_bam','get_bam() {\n    printf "%s\\n" "$SCMUT_BAM"\n}')
            text=text.replace('samples=("$@")','samples=("$SCMUT_UNIT")')
            suffix='.final_sim_region.bed' if key=='smartseq_vaf' else '.final_sim_region_ge1.bed'
            inputs=[cell['bam'],str(Path(region)/(cell['id']+suffix)),resources['reference']]
            jobs.append(job('native-'+key+'-'+cell['id'],'simulation',text,inputs,[],dict(SCMUT_UNIT=cell['id'],SCMUT_BAM=cell['bam'])))
        return jobs,root
    if args.stage=='spike':
        i={'vaf':213,'dp':214,'dp-effective':215}[args.mode]
        text=NATIVE[str(i)]
        for key,value in [('REF',required(resources,'reference')),('OUTDIR',root)]:text=bind(text,key,value)
        if args.mode!='dp-effective':
            text=bind(text,'REGION_DIR',required(settings,'region_dir'))
            text=change_function(text,'get_bam','get_bam() {\n    printf "%s\\n" "$SCMUT_BAM"\n}')
            inputs=[required(sample,'pooled_bam'),required(settings,'region_dir')]
        else:
            text=bind(text,'OLDOUT',required(settings,'dp_cache_root'));inputs=[settings['dp_cache_root']]
        if args.mode!='dp-effective':
            suffix='.final_sim_region_dp40.bed' if args.mode=='vaf' else '.final_sim_region.bed'
            inputs.append(str(Path(settings['region_dir'])/(args.sample+suffix)))
        inputs.append(resources['reference'])
        text=text.replace('samples=("$@")','samples=("$SCMUT_SAMPLE")')
        return [job('native-'+args.mode+'-simulation','simulation',text,inputs,[],dict(SCMUT_SAMPLE=args.sample,SCMUT_BAM=sample.get('pooled_bam','')))],root
    if not args.bam:raise ValueError('--bam is required for --stage call')
    tools=args.tools.split(',')
    if len(set(tools))!=len(tools) or any(t not in ['gatk','freebayes','strelka2','varscan2','monopogen'] for t in tools):raise ValueError('Choose gatk,freebayes,strelka2,varscan2,monopogen')
    tag=Path(args.bam).name.removesuffix('.spikein.sorted.calmd.bam')
    jobs=[]
    for tool in tools:
        text=NATIVE['218' if tool=='gatk' else '217']
        if tool=='gatk':text=text[:text.index('\necho -e "Tag\\tSample')]
        else:text=text[:text.index('\n#for bam in ')]
        text=re.sub(r'^\s*(source .*conda.sh|conda activate .*|conda deactivate)\s*$', '',text,flags=re.M)
        bindings={'REF':required(resources,'reference'),'SIM_BASE':root,'ORIG_MANIFEST':required(settings,'original_vcf_manifest')}
        if tool!='gatk':bindings.update(REGION_FILE=resources.get('freebayes_regions',''),VARSCAN_JAR=resources.get('varscan_jar',''),STRELKA_BIN=resources.get('strelka',''),MONO_PATH=str(Path(resources.get('monopogen','')).parent.parent),PANEL_PATH=resources.get('panel',''),MONO_REGION=sample.get('monopogen_regions',''))
        for key,value in bindings.items():text=bind(text,key,value)
        fn='run_gatk_hap' if tool=='gatk' else 'run_'+tool
        text+='\n'+fn+' '+quote(tag)+' '+quote(args.sample)+' '+quote(args.bam)+'\n'
        jobs.append(job('native-simulation-call-'+tool,tool,text,[args.bam,settings['original_vcf_manifest']],[]))
    return jobs,root


def consensus_plan(cfg,args):
    sample=cfg['samples'][args.sample];validate_id(args.sample)
    root=Path(cfg['output_root']).expanduser().resolve()/'consensus'/args.platform/args.sample
    mapping={'gatk':'GATK','freebayes':'Freebayes','strelka2':'Strelka2','varscan2':'Varscan2','monopogen':'Monopogen','scsnv':'scSNV','cellsnp':'cellsnp','samtools':'samtools','monovar':'monovar','scallele':'scallele'}
    requested=args.tools.split(',')
    if len(set(requested))!=len(requested) or not 2<=len(requested)<=4:raise ValueError('Choose 2 to 4 unique tool names')
    names=[mapping[t] if args.platform=='10x' else ('hap' if t=='gatk' else t) for t in requested]
    jobs=[];keys=root/'keys';rows=[]
    with open(args.calls) as stream:
        rows=[r for r in csv.DictReader(stream,delimiter='\t') if r['sample']==args.sample and r['tool'] in requested]
    if not rows:raise ValueError('No matching rows in --calls')
    units={r['unit'] for r in rows}
    if args.platform=='10x' and units!={args.sample}:raise ValueError('10x consensus uses pooled VCFs; pass run_10x_pooled.py outputs')
    expected={(u,t) for u in units for t in requested}
    if {(r['unit'],r['tool']) for r in rows}!=expected or len(rows)!=len(expected):raise ValueError('Each unit must have exactly one VCF for every selected tool')
    cohort=sample.get('cohort','')
    cmd='mkdir -p '+quote(keys)+'\n'
    for row in rows:
        tool=mapping[row['tool']] if args.platform=='10x' else ('hap' if row['tool']=='gatk' else row['tool'])
        dest=keys/(args.sample+'_'+tool)/'query.keys' if args.platform=='10x' else keys/(cohort+'_'+row['unit']+'_'+tool+'.query.keys')
        cmd+='mkdir -p '+quote(dest.parent)+'\nbcftools query -f \'%CHROM:%POS:%REF:%ALT\\n\' '+quote(row['vcf'])+' > '+quote(dest)+'\n'
    jobs.append(job('consensus-input-keys','scmut',cmd,[r['vcf'] for r in rows],[]))
    text=NATIVE['196' if args.platform=='10x' else '206']
    text=re.sub(r'^combos=\(.*?^\)[^\n]*', 'combos=("'+','.join(names)+'")',text,flags=re.M|re.S)
    text=bind(text,'POOL_BASE',keys);text=bind(text,'OUT_BASE',root/'results');text=bind(text,'REF_SDF',required(cfg['resources'],'rtg_sdf'))
    inputs=[required(sample,'truth_vcf'),required(cfg['resources'],'rtg_sdf')]
    if args.platform=='10x':
        text=re.sub(r'^samples=\([^\n]*\)', 'samples=("'+args.sample+'")',text,flags=re.M)
        for fn,key in [('get_truth_vcf','truth_vcf'),('get_bed','callable_bed')]:
            text=change_function(text,fn,fn+'() {\n    printf "%s\\n" '+quote(required(sample,key))+'\n}')
        inputs.append(sample['callable_bed'])
    else:
        manifest=root/'manifest.tsv';beds=root/'callable'
        cells={u['id']:u for u in sample['cells']}
        setup='mkdir -p '+quote(beds)+'\n'
        lines='scRNA_SRR\tsample\tWES_SRR\n'
        for u in sorted(units):
            cell=cells[u];bed=required(cell,'callable_bed');inputs.append(bed)
            setup+='ln -sf '+quote(bed)+' '+quote(beds/(u+'.bed'))+'\n'
            lines+=u+'\t'+args.sample+'\t'+required(sample,'wes_id')+'\n'
        setup+='printf %s '+quote(lines)+' > '+quote(manifest)+'\n'
        jobs.append(job('consensus-manifest','scmut',setup,inputs,[]))
        text=bind(text,'MANIFEST',manifest)
        a=text.index('\tif [[ ${sample} == PT* ]]');b=text.index('\n\t[[ ! -f ${truth_vcf}',a)
        text=text[:a]+'\tDATASET='+quote(required(sample,'cohort'))+'\n\ttruth_vcf='+quote(sample['truth_vcf'])+'\n\tBED='+quote(beds)+'/${scRNA_SRR}.bed\n'+text[b:]
    jobs.append(job('native-consensus-and-rtg','rtg',text,inputs,[]))
    collect=NATIVE['197' if args.platform=='10x' else '207']
    collect=bind(collect,'BASE',root/'results')
    summary=root/'results'/('rtg_summary_all.tsv' if args.platform=='10x' else 'consensus_rtg_all_long_0816.tsv')
    jobs.append(job('native-consensus-summary','scmut',collect,[],[str(summary)]))
    return jobs,root


def main():
    parser=argparse.ArgumentParser(description=KIND+' workflow using the original Rmd procedures')
    parser.add_argument('--config',default='config/workflow.json');parser.add_argument('--sample',required=True)
    parser.add_argument('--dry-run',action='store_true')
    if KIND=='simulation':
        parser.add_argument('--mode',choices=['vaf','dp','dp-effective','smartseq-vaf','smartseq-dp'],default='vaf')
        parser.add_argument('--stage',choices=['spike','call'],default='spike');parser.add_argument('--bam');parser.add_argument('--unit')
        parser.add_argument('--tools',default='gatk,freebayes,strelka2,varscan2,monopogen')
    else:
        parser.add_argument('--platform',choices=['10x','smartseq'],required=True)
        parser.add_argument('--tools',required=True);parser.add_argument('--calls',required=True)
    args=parser.parse_args();cfg=json.loads(Path(args.config).read_text())
    jobs,root=simulation_plan(cfg,args) if KIND=='simulation' else consensus_plan(cfg,args)
    run_jobs(jobs,cfg,root,args.dry_run)


if __name__=='__main__':
    try:main()
    except (ValueError,KeyError,RuntimeError,FileNotFoundError,subprocess.CalledProcessError,BlockingIOError) as error:
        print(str(error),file=sys.stderr);sys.exit(1)
