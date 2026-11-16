#!/usr/bin/env python3
KIND='simulation'
NATIVE={}
NATIVE['213']=r'''
#!/bin/bash
set -euo pipefail

REF=/home/tantingting/scmut/ref/bwa_ref/genome.fa
REGION_DIR=/mnt/storage/tantingting/scmut/simu/region/final_sim_region_dp40
OUTDIR=/mnt/storage/tantingting/scmut/simu/direct_spike_vaf_dp40_final

N_SITE=10000
FLANK=0
MINDEPTH=40
MINMUTREADS=2
THREADS=16
SEED=1234
VAFS=(0.05 0.1 0.25 0.5 0.75)

mkdir -p ${OUTDIR}/{scripts,sites,window,bam,log,truth,summary}
mkdir -p ${OUTDIR}/cache

#samples=(P4_cSCC)
#samples=(
#P4_cSCC P5_cSCC P6_cSCC P7_cSCC P9_cSCC P10_cSCC
#Acral1 Acral2 Acral3pre Acral3post Acral4 Cutaneous1 Cutaneous1Lym Cutaneous2
#)

samples=("$@")

if [[ ${#samples[@]} -eq 0 ]]; then
    echo "[ERROR] No sample specified"
    exit 1
fi

get_bam() {
    sample=$1
    case "$sample" in
        P*_cSCC) echo "/mnt/storage/tantingting/scmut/bam/1_cSCC_scRNA/${sample}/outs/${sample}_map.bam" ;;
        Acral1) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Acral1_SRR20791206/outs/Acral1_SRR20791206_map.bam" ;;
        Acral2) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Acral2_SRR20791205/outs/Acral2_SRR20791205_map.bam" ;;
        Acral3pre) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Acral3pre_SRR20791203/outs/Acral3pre_SRR20791203_map.bam" ;;
        Acral3post) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Acral3post_SRR20791202/outs/Acral3post_SRR20791202_map.bam" ;;
        Acral4) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Acral4_SRR20791201/outs/Acral4_SRR20791201_map.bam" ;;
        Cutaneous1) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Cutaneous1_SRR20791198/outs/Cutaneous1_SRR20791198_map.bam" ;;
        Cutaneous1Lym) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Cutaneous1Lym_SRR20791196/outs/Cutaneous1Lym_SRR20791196_map.bam" ;;
        Cutaneous2) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Cutaneous2_SRR20791197/outs/Cutaneous2_SRR20791197_map.bam" ;;
        *) echo "NA" ;;
    esac
}

############################################

############################################

cat > ${OUTDIR}/scripts/make_master_sites.py << 'PY'
import sys, random, subprocess, bisect
import os
import pickle

region_bed, ref, n_site, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
random.seed(seed)
bases = ["A", "C", "G", "T"]
regions, cum = [], []
total = 0

with open(region_bed) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start, end = line.strip().split()[:3]
        start, end = int(start), int(end)
        if end <= start:
            continue
        total += end - start
        regions.append((chrom, start, end))
        cum.append(total)

if total == 0:
    sys.exit("No valid regions")

seen, tries, max_tries = set(), 0, n_site * 200

with open(out, "w") as w:
    while len(seen) < n_site and tries < max_tries:
        tries += 1
        i = bisect.bisect_left(cum, random.randint(1, total))
        chrom, start, end = regions[i]
        pos = random.randint(start, end - 1) + 1
        key = (chrom, pos)

        if key in seen:
            continue

        try:
            fa = subprocess.check_output(["samtools", "faidx", ref, f"{chrom}:{pos}-{pos}"], text=True).splitlines()
            ref_base = "".join(fa[1:]).upper()
        except Exception:
            continue

        if len(ref_base) != 1 or ref_base not in bases:
            continue

        alt = random.choice([b for b in bases if b != ref_base])
        seen.add(key)
        w.write(f"{chrom}\t{pos}\t{pos}\t{ref_base}\t{alt}\n")
PY

############################################




############################################

cat > ${OUTDIR}/scripts/direct_spike_bam.py << 'PY'
import sys, random, pysam, math
import bisect
from collections import defaultdict
import os
import pickle
from collections import defaultdict

bam_in, sites_file, ref_fa, bam_out, truth_out, min_depth, min_mutreads, seed, cache_file = sys.argv[1:10]
min_depth, min_mutreads, seed = int(min_depth), int(min_mutreads), int(seed)
random.seed(seed)

bam = pysam.AlignmentFile(bam_in, "rb")

sites = []
with open(sites_file) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        fs = line.strip().split()
        chrom, pos, end, vaf, ref_base, alt = fs[0], int(fs[1]), int(fs[2]), float(fs[3]), fs[4].upper(), fs[5].upper()
        sites.append((chrom, pos, vaf, ref_base, alt))
sites.sort(key=lambda x: (x[0], x[1])) 

# fragment_edit_map:
#   key = (chrom, pos, query_name)
#   value = alt

fragment_edit_map = {}
truth = []


############################################




############################################

##################################################

##################################################

if os.path.exists(cache_file):

    print(
        f"[INFO] load cache: {cache_file}",
        file=sys.stderr
    )

    with open(cache_file, "rb") as f:
        site_info = pickle.load(f)

    
    for chrom, pos, vaf, ref_base, alt in sites:

        if (chrom, pos) not in site_info:
            continue

        site_info[(chrom, pos)]["vaf"] = vaf

else:

    print(
        f"[INFO] build cache: {cache_file}",
        file=sys.stderr
    )
    site_info = {}

    
    for chrom, pos, vaf, ref_base, alt in sites:
        site_info[(chrom, pos)] = {
            "vaf": vaf,
            "ref": ref_base,
            "alt": alt,
            "cover": 0,
            "candidates": set()
        }
    chr_sites = defaultdict(set)

    for chrom, pos, vaf, ref_base, alt in sites:
        chr_sites[chrom].add(pos)

    print(
        f"[INFO] {len(site_info)} sites loaded",
        file=sys.stderr
    )

############################################

############################################

    for chrom in sorted(chr_sites.keys()):

        target_pos = sorted(chr_sites[chrom])

        
        start = min(target_pos) - 1
        end = max(target_pos)+1

        print(
            f"[INFO] scanning {chrom}",
            file=sys.stderr
        )

        try:
            iterator = bam.fetch(chrom, start, end)
        except ValueError:
            continue

        n_read = 0

        for r in iterator:

            if r.is_unmapped:
                continue

            if r.is_secondary:
                continue

            if r.is_supplementary:
                continue

            seq = r.query_sequence

            if seq is None:
                continue

            n_read += 1

        ####################################
        
        ####################################

            aligned = {
                rp + 1: qp
                for qp, rp in r.get_aligned_pairs(matches_only=False)
                if qp is not None and rp is not None
            }

            if not aligned:
                continue

        #read_start = min(aligned)
        #read_end = max(aligned)
        
            read_start = r.reference_start + 1
            read_end = r.reference_end
            if read_end is None:
                continue

        ####################################
        
        ####################################

            left = bisect.bisect_left(
                target_pos,
                read_start
            )

            right = bisect.bisect_right(
                target_pos,
                read_end
            )

            covered_sites = target_pos[left:right]

            if not covered_sites:
                continue

            for pos in covered_sites:

                qp = aligned.get(pos)

                if qp is None:
                    continue

                base = seq[qp].upper()

                if base not in "ACGT":
                    continue

                info = site_info[(chrom, pos)]

                info["cover"] += 1

                if base == info["ref"]:
                    info["candidates"].add(
                        r.query_name
                    )

        print(
            f"[INFO] {chrom}: {n_read} reads scanned",
            file=sys.stderr
        )
        
##################################################

##################################################

    with open(cache_file, "wb") as f:
        pickle.dump(site_info, f)
    print(
        f"[INFO] cache saved: {cache_file}",
        file=sys.stderr
    )

############################################

############################################

fragment_edit_map = {}
truth = []

for (chrom, pos), info in site_info.items():

    dp = info["cover"]

    if dp < min_depth:
        continue

    candidates = list(
        info["candidates"]
    )

    if len(candidates) == 0:
        continue

    target_alt_reads = max(
        min_mutreads,
        math.ceil(dp * info["vaf"])
    )

    n_choose = min(
        len(candidates),
        target_alt_reads
    )

    chosen_names = random.sample(
        candidates,
        n_choose
    )

    for name in chosen_names:
        fragment_edit_map[
            (chrom, pos, name)
        ] = info["alt"]

    truth.append(
        (
            chrom,
            pos,
            info["ref"],
            info["alt"],
            dp,
            len(candidates),
            n_choose,
            "NA"
        )
    )



print(
    f"[INFO] fragment_edit_map={len(fragment_edit_map)}",
    file=sys.stderr
)

############################################

#




#


#       (chrom, pos, query_name) -> alt
#

#       query_name -> chrom -> pos -> alt
#


############################################

bam.close()
bam = pysam.AlignmentFile(bam_in, "rb")
out = pysam.AlignmentFile(bam_out, "wb", template=bam)

edited_reads = 0
edited_bases = 0
site_edited_bases = defaultdict(int)

############################################

############################################

edit_by_qname = defaultdict(lambda: defaultdict(dict))

for chrom, pos, qname in fragment_edit_map:
    alt = fragment_edit_map[(chrom, pos, qname)]
    edit_by_qname[qname][chrom][pos] = alt

print(
    f"[INFO] edit_by_qname={len(edit_by_qname)}",
    file=sys.stderr
)

############################################

############################################

for r in bam:

    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.query_sequence is None:
        out.write(r)
        continue

    qname = r.query_name

    ########################################
    
    
    ########################################

    if qname not in edit_by_qname:
        out.write(r)
        continue

    chrom = r.reference_name

    if chrom is None:
        out.write(r)
        continue

    ########################################
    
    
    ########################################

    if chrom not in edit_by_qname[qname]:
        out.write(r)
        continue

    target_edits = edit_by_qname[qname][chrom]

    seq = list(r.query_sequence)
    quals = r.query_qualities
    changed = 0

    ########################################
    
    ########################################

    for qp, rp in r.get_aligned_pairs(matches_only=False):

        if qp is None or rp is None:
            continue

        pos = rp + 1

        if pos not in target_edits:
            continue

        alt = target_edits[pos]

        if qp < len(seq) and seq[qp].upper() in "ACGT" and seq[qp].upper() != alt:
            seq[qp] = alt
            changed += 1
            site_edited_bases[(chrom, pos)] += 1

    if changed > 0:
        r.query_sequence = "".join(seq)
        r.query_qualities = quals
        edited_reads += 1
        edited_bases += changed

    out.write(r)

out.close()
bam.close()

############################################



# real_vaf = mut_bases_edited / dp
############################################

with open(truth_out, "w") as w:
    w.write("chrom\tpos\tref\talt\tdp\tcandidate_fragments\tmut_fragments_planned\tmut_bases_edited\treal_vaf\n")
    for chrom, pos, ref_base, alt, dp, n_candidate, n_choose, _ in truth:
        mut_bases = site_edited_bases.get((chrom, pos), 0)
        real_vaf = mut_bases / dp if dp > 0 else 0
        w.write("\t".join(map(str, [chrom, pos, ref_base, alt, dp, n_candidate, n_choose, mut_bases, real_vaf])) + "\n")

sys.stderr.write(f"[INFO] sites_input={len(sites)} truth_sites={len(truth)} edited_reads={edited_reads} edited_bases={edited_bases}\n")
PY

manifest=${OUTDIR}/summary/prepare_manifest.tsv

if [[ ! -s ${manifest} ]]; then
    echo -e "Sample\tMode\tVAF\tMasterSites\tSites\tWindowBed\tOriginalSubsetBam\tSpikeinBam\tCalmdBam\tTruth" > ${manifest}
fi

############################################

############################################

for sample in "${samples[@]}"; do
    echo "========== ${sample} =========="

    bam=$(get_bam ${sample})
    region_bed=${REGION_DIR}/${sample}.final_sim_region_dp40.bed
    master_sites=${OUTDIR}/sites/${sample}.master_sites.txt
    sample_window=${OUTDIR}/window/${sample}.master.window.bed
    sample_subset=${OUTDIR}/bam/${sample}.master.original_subset.bam

    [[ "${bam}" != "NA" ]] || { echo "[SKIP] unknown sample: ${sample}"; continue; }
    [[ -f ${region_bed} ]] || { echo "[SKIP] no region: ${region_bed}"; continue; }
    [[ -f ${bam} ]] || { echo "[SKIP] no bam: ${bam}"; continue; }
    [[ -f ${bam}.bai ]] || samtools index ${bam}

    if [[ -s ${master_sites} ]]; then
        echo "[1] master sites exists, skip: ${master_sites}"
    else
        echo "[1] make master sites"
        python ${OUTDIR}/scripts/make_master_sites.py ${region_bed} ${REF} ${N_SITE} ${master_sites} ${SEED}
    fi

    n_master=$(wc -l < ${master_sites})
    echo "[INFO] master sites: ${n_master}"
    [[ ${n_master} -gt 0 ]] || { echo "[SKIP] no master sites: ${sample}"; continue; }

    if [[ -s ${sample_window} ]]; then
        echo "[2] master window exists, skip: ${sample_window}"
    else
        echo "[2] make master window"
        awk -v OFS="\t" -v flank=${FLANK} '{s=$2-1-flank; e=$3+flank; if(s<0)s=0; print $1,s,e}' ${master_sites} \
            | sort -k1,1 -k2,2n | bedtools merge > ${sample_window}
    fi

    if [[ -s ${sample_subset} && -s ${sample_subset}.bai ]]; then
        echo "[3] master subset BAM exists, skip: ${sample_subset}"
    else
        echo "[3] subset BAM once"
        samtools view -@ ${THREADS} -P -b -L ${sample_window} ${bam} > ${sample_subset}
        samtools index ${sample_subset}
    fi

    for vaf in "${VAFS[@]}"; do
        tag=${sample}.vaf${vaf}
        cache=${OUTDIR}/cache/${sample}.site_info.pkl
        sites=${OUTDIR}/sites/${tag}.sites.txt
        spikein=${OUTDIR}/bam/${tag}.spikein.bam
        spikein_sorted=${OUTDIR}/bam/${tag}.spikein.sorted.bam
        spikein_calmd=${OUTDIR}/bam/${tag}.spikein.sorted.calmd.bam
        truth=${OUTDIR}/truth/${tag}.truth.tsv

        echo "[4] ${sample} VAF=${vaf}"

        if [[ -s ${sites} ]]; then
            echo "[4] sites exists, skip: ${sites}"
        else
            awk -v OFS="\t" -v vaf=${vaf} '{print $1,$2,$3,vaf,$4,$5}' ${master_sites} > ${sites}
        fi

        if [[ -s ${spikein_calmd} && -s ${spikein_calmd}.bai && -s ${truth} ]]; then
            echo "[5] calmd spike-in BAM exists, skip: ${spikein_calmd}"
        else
            echo "[5] direct spike-in BAM"
            python ${OUTDIR}/scripts/direct_spike_bam.py ${sample_subset} ${sites} ${REF} ${spikein} \
                ${truth} ${MINDEPTH} ${MINMUTREADS} ${SEED} ${cache} > ${OUTDIR}/log/${tag}.log 2>&1

            [[ -s ${spikein} ]] || { echo "[ERROR] direct spike-in failed: check ${OUTDIR}/log/${tag}.log"; continue; }

            echo "[6] sort/index"
            samtools sort -@ ${THREADS} -o ${spikein_sorted} ${spikein}
            samtools index ${spikein_sorted}
            samtools quickcheck ${spikein_sorted}

            echo "[7] calmd/index"
            samtools calmd -@ ${THREADS} -b ${spikein_sorted} ${REF} > ${spikein_calmd}
            samtools index ${spikein_calmd}
            samtools quickcheck ${spikein_calmd}
        fi

        n_truth=$(($(wc -l < ${truth}) - 1))
        echo "[INFO] truth sites: ${n_truth}"

        if ! awk -v s="${sample}" -v v="${vaf}" '$1==s && $2=="VAF" && $3==v{found=1} END{exit found?0:1}' ${manifest}; then
            echo -e "${sample}\tVAF\t${vaf}\t${master_sites}\t${sites}\t${sample_window}\t${sample_subset}\t${spikein_sorted}\t${spikein_calmd}\t${truth}" >> ${manifest}
        fi
    done

    echo "[DONE] ${sample}"
done
'''
NATIVE['214']=r'''
#!/bin/bash
set -euo pipefail

############################################
# DP simulation design
#



#




#        DP >=10 / 25 / 50 / 100


#



############################################

REF=/home/tantingting/scmut/ref/bwa_ref/genome.fa
REGION_DIR=/mnt/storage/tantingting/scmut/simu/region/final_sim_region
OUTDIR=/mnt/storage/tantingting/scmut/simu/direct_spike_dp_gt10_vaf04_q20_final

N_CANDIDATE=50000
N_SITE=10000
FLANK=0
FIXED_VAF=0.4
DP_BINS=("10_25" "25_50" "50_100" "100plus")

MINMUTREADS=2
THREADS=16
SEED=1234
MIN_BQ=20

mkdir -p ${OUTDIR}/{scripts,sites,window,bam,log,truth,summary,cache}

samples=("$@")

if [[ ${#samples[@]} -eq 0 ]]; then
    echo "[ERROR] No sample specified"
    exit 1
fi

get_bam() {
    sample=$1
    case "$sample" in
        P*_cSCC) echo "/mnt/storage/tantingting/scmut/bam/1_cSCC_scRNA/${sample}/outs/${sample}_map.bam" ;;
        Acral1) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Acral1_SRR20791206/outs/Acral1_SRR20791206_map.bam" ;;
        Acral2) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Acral2_SRR20791205/outs/Acral2_SRR20791205_map.bam" ;;
        Acral3pre) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Acral3pre_SRR20791203/outs/Acral3pre_SRR20791203_map.bam" ;;
        Acral3post) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Acral3post_SRR20791202/outs/Acral3post_SRR20791202_map.bam" ;;
        Acral4) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Acral4_SRR20791201/outs/Acral4_SRR20791201_map.bam" ;;
        Cutaneous1) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Cutaneous1_SRR20791198/outs/Cutaneous1_SRR20791198_map.bam" ;;
        Cutaneous1Lym) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Cutaneous1Lym_SRR20791196/outs/Cutaneous1Lym_SRR20791196_map.bam" ;;
        Cutaneous2) echo "/mnt/storage/tantingting/scmut/bam/2_melanoma_scRNA/Cutaneous2_SRR20791197/outs/Cutaneous2_SRR20791197_map.bam" ;;
        *) echo "NA" ;;
    esac
}

############################################


############################################

cat > ${OUTDIR}/scripts/make_candidate_sites.py << 'PY'
import sys, random, subprocess, bisect

region_bed, ref, n_site, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
random.seed(seed)

bases = ["A", "C", "G", "T"]
regions, cum = [], []
total = 0

with open(region_bed) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start, end = line.strip().split()[:3]
        start, end = int(start), int(end)
        if end <= start:
            continue
        total += end - start
        regions.append((chrom, start, end))
        cum.append(total)

if total == 0:
    sys.exit("No valid regions")

seen = set()
tries = 0
max_tries = n_site * 300

with open(out, "w") as w:
    while len(seen) < n_site and tries < max_tries:
        tries += 1
        i = bisect.bisect_left(cum, random.randint(1, total))
        chrom, start, end = regions[i]

        
        pos = random.randint(start, end - 1) + 1
        key = (chrom, pos)

        if key in seen:
            continue

        try:
            fa = subprocess.check_output(
                ["samtools", "faidx", ref, f"{chrom}:{pos}-{pos}"],
                text=True
            ).splitlines()
            ref_base = "".join(fa[1:]).upper()
        except Exception:
            continue

        if len(ref_base) != 1 or ref_base not in bases:
            continue

        alt = random.choice([b for b in bases if b != ref_base])
        seen.add(key)
        w.write(f"{chrom}\t{pos}\t{pos}\t{ref_base}\t{alt}\n")

print(f"[INFO] candidate sites generated: {len(seen)}", file=sys.stderr)
PY

############################################


############################################

cat > ${OUTDIR}/scripts/build_site_cache.py << 'PY'
import sys, pysam, bisect, pickle
from collections import defaultdict

bam_in, candidate_sites, cache_file, min_bq = sys.argv[1:5]
min_bq = int(min_bq)

bam = pysam.AlignmentFile(bam_in, "rb")

site_info = {}
chr_sites = defaultdict(set)

with open(candidate_sites) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, pos, end, ref_base, alt = line.strip().split()[:5]
        pos = int(pos)
        ref_base = ref_base.upper()
        alt = alt.upper()

        site_info[(chrom, pos)] = {
            "ref": ref_base,
            "alt": alt,
            "cover": 0,
            "candidates": set()
        }
        chr_sites[chrom].add(pos)

print(f"[INFO] candidate sites loaded: {len(site_info)}", file=sys.stderr)

for chrom in sorted(chr_sites.keys()):
    target_pos = sorted(chr_sites[chrom])

    start = min(target_pos) - 1
    end = max(target_pos) + 1

    print(f"[INFO] scanning {chrom}", file=sys.stderr)

    try:
        iterator = bam.fetch(chrom, start, end)
    except ValueError:
        continue

    n_read = 0

    for r in iterator:
        #if r.is_unmapped or r.is_secondary or r.is_supplementary:
        #    continue
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail:
            continue
        if r.query_qualities is None:
            continue

        seq = r.query_sequence
        if seq is None:
            continue

        read_start = r.reference_start + 1
        read_end = r.reference_end

        if read_end is None:
            continue

        left = bisect.bisect_left(target_pos, read_start)
        right = bisect.bisect_right(target_pos, read_end)

        covered_sites = target_pos[left:right]
        if not covered_sites:
            continue

        n_read += 1

        aligned = {
            rp + 1: qp
            for qp, rp in r.get_aligned_pairs(matches_only=False)
            if qp is not None and rp is not None
        }

        if not aligned:
            continue

        for pos in covered_sites:
            qp = aligned.get(pos)

            if qp is None:
                continue
            if r.query_qualities[qp] < min_bq:
                continue
            base = seq[qp].upper()

            if base not in "ACGT":
                continue

            info = site_info[(chrom, pos)]
            info["cover"] += 1

            
            if base == info["ref"]:
                info["candidates"].add(r.query_name)

    print(f"[INFO] {chrom}: reads scanned={n_read}", file=sys.stderr)

bam.close()

with open(cache_file, "wb") as f:
    #pickle.dump(site_info, f)
    pickle.dump(site_info, f, protocol=pickle.HIGHEST_PROTOCOL)

n_ge10 = sum(1 for x in site_info.values() if x["cover"] >= 10)
n_ge25 = sum(1 for x in site_info.values() if x["cover"] >= 25)
n_ge50 = sum(1 for x in site_info.values() if x["cover"] >= 50)
n_ge100 = sum(1 for x in site_info.values() if x["cover"] >= 100)

print(f"[INFO] cache saved: {cache_file}", file=sys.stderr)
print(f"[INFO] DP>=10: {n_ge10}", file=sys.stderr)
print(f"[INFO] DP>=25: {n_ge25}", file=sys.stderr)
print(f"[INFO] DP>=50: {n_ge50}", file=sys.stderr)
print(f"[INFO] DP>=100: {n_ge100}", file=sys.stderr)
PY

############################################


############################################

cat > ${OUTDIR}/scripts/select_sites_by_dp.py << 'PY'
import sys, pickle, random, math

cache_file, dp_bin, n_site, vaf, out, seed = sys.argv[1:7]
n_site = int(n_site)
vaf = float(vaf)
seed = int(seed)

if dp_bin == "10_25":
    lo, hi = 10, 25
elif dp_bin == "25_50":
    lo, hi = 25, 50
elif dp_bin == "50_100":
    lo, hi = 50, 100
elif dp_bin == "100plus":
    lo, hi = 100, None
else:
    sys.exit(f"Unknown dp_bin: {dp_bin}")

random.seed(seed)

with open(cache_file, "rb") as f:
    site_info = pickle.load(f)

eligible = []

for (chrom, pos), info in site_info.items():
    dp = info["cover"]
    n_candidate = len(info["candidates"])

    in_bin = dp >= lo if hi is None else lo <= dp < hi
    need_alt = max(2, math.ceil(dp * vaf))

    if in_bin and n_candidate >= need_alt:
        eligible.append((chrom, pos, info["ref"], info["alt"], dp, n_candidate))

eligible.sort(key=lambda x: (x[0], x[1]))

print(f"[INFO] DPbin {dp_bin} eligible sites: {len(eligible)}", file=sys.stderr)

if len(eligible) == 0:
    sys.exit(f"No eligible sites for DP bin {dp_bin}")

if len(eligible) < n_site:
    print(f"[WARN] eligible sites < N_SITE: use all {len(eligible)} sites", file=sys.stderr)
    selected = eligible
else:
    selected = random.sample(eligible, n_site)
    selected.sort(key=lambda x: (x[0], x[1]))

with open(out, "w") as w:
    for chrom, pos, ref_base, alt, dp, n_candidate in selected:
        
        # chrom pos end vaf ref alt original_dp candidate_fragments
        w.write(f"{chrom}\t{pos}\t{pos}\t{vaf}\t{ref_base}\t{alt}\t{dp}\t{n_candidate}\n")

print(f"[INFO] selected sites: {len(selected)}", file=sys.stderr)
PY

############################################


############################################

cat > ${OUTDIR}/scripts/direct_spike_from_cache.py << 'PY'
import sys, random, pysam, math, pickle
from collections import defaultdict

bam_in, sites_file, bam_out, truth_out, min_mutreads, seed, cache_file, min_bq = sys.argv[1:9]
min_mutreads = int(min_mutreads)
seed = int(seed)
min_bq = int(min_bq)

random.seed(seed)

with open(cache_file, "rb") as f:
    site_info = pickle.load(f)

sites = []

with open(sites_file) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue

        fs = line.strip().split()
        chrom = fs[0]
        pos = int(fs[1])
        vaf = float(fs[3])
        ref_base = fs[4].upper()
        alt = fs[5].upper()

        sites.append((chrom, pos, vaf, ref_base, alt))

sites.sort(key=lambda x: (x[0], x[1]))

fragment_edit_map = {}
truth = []

############################################

############################################

for chrom, pos, vaf, ref_base, alt in sites:
    key = (chrom, pos)

    if key not in site_info:
        continue

    info = site_info[key]
    dp = info["cover"]
    candidates = list(info["candidates"])

    if len(candidates) == 0:
        continue

    target_alt_reads = max(min_mutreads, math.ceil(dp * vaf))
    n_choose = min(len(candidates), target_alt_reads)

    chosen_names = random.sample(candidates, n_choose)

    for name in chosen_names:
        fragment_edit_map[(chrom, pos, name)] = alt

    truth.append((chrom, pos, ref_base, alt, dp, len(candidates), n_choose))

print(f"[INFO] sites_input={len(sites)}", file=sys.stderr)
print(f"[INFO] truth_sites={len(truth)}", file=sys.stderr)
print(f"[INFO] fragment_edit_map={len(fragment_edit_map)}", file=sys.stderr)

############################################


############################################

edit_by_qname = defaultdict(lambda: defaultdict(dict))

for chrom, pos, qname in fragment_edit_map:
    edit_by_qname[qname][chrom][pos] = fragment_edit_map[(chrom, pos, qname)]

bam = pysam.AlignmentFile(bam_in, "rb")
out = pysam.AlignmentFile(bam_out, "wb", template=bam)

edited_reads = 0
edited_bases = 0
site_edited_bases = defaultdict(int)

for r in bam:
    #if r.is_unmapped or r.is_secondary or r.is_supplementary or r.query_sequence is None:
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.query_sequence is None:
        out.write(r)
        continue

    qname = r.query_name

    if qname not in edit_by_qname:
        out.write(r)
        continue

    chrom = r.reference_name

    if chrom is None or chrom not in edit_by_qname[qname]:
        out.write(r)
        continue

    target_edits = edit_by_qname[qname][chrom]

    seq = list(r.query_sequence)
    quals = r.query_qualities
    changed = 0

    for qp, rp in r.get_aligned_pairs(matches_only=False):
        if qp is None or rp is None:
            continue

        pos = rp + 1

        if pos not in target_edits:
            continue
        if quals is None or quals[qp] < min_bq:
            continue
        alt = target_edits[pos]

        if qp < len(seq) and seq[qp].upper() in "ACGT" and seq[qp].upper() != alt:
            seq[qp] = alt
            changed += 1
            site_edited_bases[(chrom, pos)] += 1

    if changed > 0:
        r.query_sequence = "".join(seq)
        r.query_qualities = quals
        edited_reads += 1
        edited_bases += changed

    out.write(r)

out.close()
bam.close()

############################################



# real_vaf = mut_bases_edited / dp
############################################

with open(truth_out, "w") as w:
    w.write("chrom\tpos\tref\talt\tdp\tcandidate_fragments\tmut_fragments_planned\tmut_bases_edited\treal_vaf\n")

    for chrom, pos, ref_base, alt, dp, n_candidate, n_choose in truth:
        mut_bases = site_edited_bases.get((chrom, pos), 0)
        real_vaf = mut_bases / dp if dp > 0 else 0

        w.write("\t".join(map(str, [
            chrom, pos, ref_base, alt,
            dp, n_candidate, n_choose,
            mut_bases, real_vaf
        ])) + "\n")

print(f"[INFO] edited_reads={edited_reads}", file=sys.stderr)
print(f"[INFO] edited_bases={edited_bases}", file=sys.stderr)
PY

manifest=${OUTDIR}/summary/prepare_manifest.tsv

if [[ ! -s ${manifest} ]]; then
    echo -e "Sample\tMode\tVAF\tDPBIN_Q20\tCandidateSites\tSelectedSites\tWindowBed\tOriginalSubsetBam\tSpikeinBam\tCalmdBam\tTruth" > ${manifest}
fi

############################################

############################################

for sample in "${samples[@]}"; do
    echo "========== ${sample} =========="

    bam=$(get_bam ${sample})
    region_bed=${REGION_DIR}/${sample}.final_sim_region.bed

    candidate_sites=${OUTDIR}/sites/${sample}.candidate_sites.txt
    candidate_window=${OUTDIR}/window/${sample}.candidate.window.bed
    sample_subset=${OUTDIR}/bam/${sample}.candidate.original_subset.bam
    cache=${OUTDIR}/cache/${sample}.site_info.pkl

    [[ "${bam}" != "NA" ]] || { echo "[SKIP] unknown sample: ${sample}"; continue; }
    [[ -f ${region_bed} ]] || { echo "[SKIP] no region: ${region_bed}"; continue; }
    [[ -f ${bam} ]] || { echo "[SKIP] no bam: ${bam}"; continue; }
    [[ -f ${bam}.bai ]] || samtools index ${bam}

    ########################################
    
    ########################################

    if [[ -s ${candidate_sites} ]]; then
        echo "[1] candidate sites exists, skip: ${candidate_sites}"
    else
        echo "[1] make candidate sites"
        python ${OUTDIR}/scripts/make_candidate_sites.py ${region_bed} ${REF} ${N_CANDIDATE} \
            ${candidate_sites} ${SEED} > ${OUTDIR}/log/${sample}.make_candidate.log 2>&1
    fi

    n_candidate=$(wc -l < ${candidate_sites})
    echo "[INFO] candidate sites: ${n_candidate}"

    [[ ${n_candidate} -gt 0 ]] || { echo "[SKIP] no candidate sites: ${sample}"; continue; }

    ########################################
    
    ########################################

    if [[ -s ${candidate_window} ]]; then
        echo "[2] candidate window exists, skip: ${candidate_window}"
    else
        echo "[2] make candidate window"
        awk -v OFS="\t" -v flank=${FLANK} '{s=$2-1-flank; e=$3+flank; if(s<0)s=0; print $1,s,e}' ${candidate_sites} \
            | sort -k1,1 -k2,2n | bedtools merge > ${candidate_window}
    fi

    if [[ -s ${sample_subset} && -s ${sample_subset}.bai ]]; then
        echo "[3] candidate subset BAM exists, skip: ${sample_subset}"
    else
        echo "[3] subset BAM once"
        samtools view -@ ${THREADS} -P -b -L ${candidate_window} ${bam} > ${sample_subset}
        samtools index ${sample_subset}
    fi

    ########################################
    
    
    ########################################

    if [[ -s ${cache} ]]; then
        echo "[4] site cache exists, skip: ${cache}"
    else
        echo "[4] build site cache"
        python ${OUTDIR}/scripts/build_site_cache.py ${sample_subset} ${candidate_sites} ${cache} ${MIN_BQ} \
            > ${OUTDIR}/log/${sample}.build_cache.log 2>&1
    fi

    ########################################
    
    ########################################
    for bin in "${DP_BINS[@]}"; do
        tag=${sample}.vaf${FIXED_VAF}.dpbin${bin}

        selected_sites=${OUTDIR}/sites/${tag}.selected_sites.txt
        spikein=${OUTDIR}/bam/${tag}.spikein.bam
        spikein_sorted=${OUTDIR}/bam/${tag}.spikein.sorted.bam
        spikein_calmd=${OUTDIR}/bam/${tag}.spikein.sorted.calmd.bam
        truth=${OUTDIR}/truth/${tag}.truth.tsv

        echo "[5] ${sample} VAF=${FIXED_VAF} DPbin=${bin}"

        if [[ -s ${selected_sites} ]]; then
            echo "[5] selected sites exists, skip: ${selected_sites}"
        else
            echo "[5] select sites by DP bin"
            python ${OUTDIR}/scripts/select_sites_by_dp.py ${cache} ${bin} ${N_SITE} ${FIXED_VAF} \
                ${selected_sites} ${SEED} > ${OUTDIR}/log/${tag}.select_sites.log 2>&1
        fi

        n_selected=$(wc -l < ${selected_sites})
        echo "[INFO] selected sites: ${n_selected}"

        [[ ${n_selected} -gt 0 ]] || { echo "[SKIP] no selected sites: ${tag}"; continue; }

        if [[ -s ${spikein_calmd} && -s ${spikein_calmd}.bai && -s ${truth} ]]; then
            echo "[6] calmd spike-in BAM exists, skip: ${spikein_calmd}"
        else
            echo "[6] direct spike-in from cache"
            python ${OUTDIR}/scripts/direct_spike_from_cache.py ${sample_subset} ${selected_sites} \
                ${spikein} ${truth} ${MINMUTREADS} ${SEED} ${cache} ${MIN_BQ} \
                > ${OUTDIR}/log/${tag}.spikein.log 2>&1

            [[ -s ${spikein} ]] || { echo "[ERROR] spike-in failed: check ${OUTDIR}/log/${tag}.spikein.log"; continue; }

            echo "[7] sort/index"
            samtools sort -@ ${THREADS} -o ${spikein_sorted} ${spikein}
            samtools index ${spikein_sorted}
            samtools quickcheck ${spikein_sorted}

            echo "[8] calmd/index"
            samtools calmd -@ ${THREADS} -b ${spikein_sorted} ${REF} > ${spikein_calmd}
            samtools index ${spikein_calmd}
            samtools quickcheck ${spikein_calmd}
        fi

        n_truth=$(($(wc -l < ${truth}) - 1))
        echo "[INFO] truth sites: ${n_truth}"

        if ! awk -v s="${sample}" -v v="${FIXED_VAF}" -v b="${bin}" '$1==s && $2=="DPBIN_Q20" && $3==v && $4==b{found=1} END{exit found?0:1}' ${manifest}; then
            echo -e "${sample}\tDPBIN_Q20\t${FIXED_VAF}\t${bin}\t${candidate_sites}\t${selected_sites}\t${candidate_window}\t${sample_subset}\t${spikein_sorted}\t${spikein_calmd}\t${truth}" >> ${manifest}
        fi
    done

    echo "[DONE] ${sample}"
done
'''
NATIVE['215']=r'''
#!/bin/bash

###########################################################
# DP effective-bin simulation v2
#



#








#      site_info.effective_q20.pkl

#      10-25, 25-50, 50-100, >=100


#




###########################################################

REF=/home/tantingting/scmut/ref/bwa_ref/genome.fa

OLDOUT=/mnt/storage/tantingting/scmut/simu/direct_spike_dp_gt10_vaf04_final
OUTDIR=/mnt/storage/tantingting/scmut/simu/direct_spike_dp_effbin_q20_vaf04_final

FIXED_VAF=0.4
N_SITE=10000
MINMUTREADS=2
THREADS=16
SEED=1234
MIN_BQ=20
MIN_MQ=0

DP_BINS=("10_25" "25_50" "50_100" "100plus")

mkdir -p ${OUTDIR}/{scripts,sites,bam,log,truth,summary,cache}

samples=("$@")
if [[ ${#samples[@]} -eq 0 ]]; then
    echo "[ERROR] No sample specified"
    echo "Usage: bash run_bam_simulation_final_dp_effbin_q20.sh Acral1"
    exit 1
fi

###########################################################
# Python 1: build_effective_cache.py
#



#


#

#   ref
#   alt
#   effective_dp

#   ref_count

#   candidate_fragments

#


###########################################################

cat > ${OUTDIR}/scripts/build_effective_cache.py << 'PY'
import sys, pysam, pickle
from collections import defaultdict

old_cache_file, bam_file, effective_cache_file, min_bq, min_mq = sys.argv[1:6]
min_bq = int(min_bq)
min_mq = int(min_mq)

with open(old_cache_file, "rb") as f:
    old = pickle.load(f)

sites_by_chr = defaultdict(list)
for (chrom, pos), info in old.items():
    sites_by_chr[chrom].append((pos, info["ref"].upper(), info["alt"].upper()))

bam = pysam.AlignmentFile(bam_file, "rb")
new = {}

for chrom in sorted(sites_by_chr):
    sites = sorted(sites_by_chr[chrom])
    pos_set = {x[0] for x in sites}

    start = min(pos_set) - 1
    end = max(pos_set)

    print(f"[INFO] scanning {chrom}", file=sys.stderr)

    for pos, ref, alt in sites:
        new[(chrom, pos)] = {
            "ref": ref,
            "alt": alt,
            "effective_dp": 0,
            "ref_count": 0,
            "candidate_fragments": set()
        }

    try:
        iterator = bam.fetch(chrom, start, end)
    except ValueError:
        continue

    for r in iterator:
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail:
            continue
        if r.mapping_quality < min_mq:
            continue
        if r.query_sequence is None or r.query_qualities is None:
            continue

        seq = r.query_sequence
        qual = r.query_qualities

        # qp: query position; rp: reference position, 0-based
        
        
        for qp, rp in r.get_aligned_pairs(matches_only=False):
            if qp is None or rp is None:
                continue

            pos = rp + 1
            if pos not in pos_set:
                continue

            if qual[qp] < min_bq:
                continue

            base = seq[qp].upper()
            if base not in "ACGT":
                continue

            info = new[(chrom, pos)]
            info["effective_dp"] += 1

            
            if base == info["ref"]:
                info["ref_count"] += 1
                info["candidate_fragments"].add(r.query_name)

bam.close()

with open(effective_cache_file, "wb") as f:
    pickle.dump(new, f, protocol=pickle.HIGHEST_PROTOCOL)

print(f"[INFO] saved effective cache: {effective_cache_file}", file=sys.stderr)

for lo, hi, label in [(10,25,"10_25"), (25,50,"25_50"), (50,100,"50_100")]:
    n = sum(1 for x in new.values() if lo <= x["effective_dp"] < hi)
    print(f"[INFO] DPbin {label}: {n}", file=sys.stderr)

n = sum(1 for x in new.values() if x["effective_dp"] >= 100)
print(f"[INFO] DPbin 100plus: {n}", file=sys.stderr)
PY

###########################################################
# Python 2: select_sites_by_effbin.py
#


#


#   2. candidate_fragments >= ceil(effective_dp * VAF)
#


###########################################################

cat > ${OUTDIR}/scripts/select_sites_by_effbin.py << 'PY'
import sys, pickle, random, math

cache_file, dp_bin, n_site, vaf, min_mutreads, out, seed = sys.argv[1:8]
n_site = int(n_site)
vaf = float(vaf)
min_mutreads = int(min_mutreads)
seed = int(seed)
random.seed(seed)

if dp_bin == "10_25":
    lo, hi = 10, 25
elif dp_bin == "25_50":
    lo, hi = 25, 50
elif dp_bin == "50_100":
    lo, hi = 50, 100
elif dp_bin == "100plus":
    lo, hi = 100, None
else:
    sys.exit(f"Unknown dp_bin: {dp_bin}")

with open(cache_file, "rb") as f:
    cache = pickle.load(f)

eligible = []

for (chrom, pos), info in cache.items():
    dp = info["effective_dp"]
    ref_count = info["ref_count"]
    n_candidate = len(info["candidate_fragments"])
    need_alt = max(min_mutreads, math.ceil(dp * vaf))

    in_bin = dp >= lo if hi is None else lo <= dp < hi

    if in_bin and n_candidate >= need_alt:
        eligible.append((chrom, pos, info["ref"], info["alt"], dp, ref_count, n_candidate, need_alt))

eligible.sort(key=lambda x: (x[0], x[1]))

print(f"[INFO] {dp_bin} eligible={len(eligible)}", file=sys.stderr)

if len(eligible) == 0:
    sys.exit(f"No eligible sites for {dp_bin}")

if len(eligible) > n_site:
    selected = random.sample(eligible, n_site)
    selected.sort(key=lambda x: (x[0], x[1]))
else:
    selected = eligible
    print(f"[WARN] eligible < N_SITE, use all: {len(selected)}", file=sys.stderr)

with open(out, "w") as w:
    for chrom, pos, ref, alt, dp, ref_count, n_candidate, need_alt in selected:
        # chrom pos end vaf ref alt effective_dp_q20 ref_count_q20 candidate_fragments_q20 need_alt dp_bin
        w.write(f"{chrom}\t{pos}\t{pos}\t{vaf}\t{ref}\t{alt}\t{dp}\t{ref_count}\t{n_candidate}\t{need_alt}\t{dp_bin}\n")

print(f"[INFO] selected={len(selected)}", file=sys.stderr)
PY

###########################################################
# Python 3: direct_spike_from_effective_cache.py
#


#






#

#   effective_dp_q20
#   ref_count_q20
#   candidate_fragments_q20
#   mut_fragments_planned
#   mut_bases_edited
#   real_vaf = mut_bases_edited / effective_dp_q20
###########################################################

cat > ${OUTDIR}/scripts/direct_spike_from_effective_cache.py << 'PY'
import sys, random, pysam, pickle
from collections import defaultdict

bam_in, sites_file, bam_out, truth_out, seed, cache_file, min_bq = sys.argv[1:8]
seed = int(seed)
min_bq = int(min_bq)
random.seed(seed)

with open(cache_file, "rb") as f:
    cache = pickle.load(f)

sites = []
with open(sites_file) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        fs = line.strip().split()
        chrom = fs[0]
        pos = int(fs[1])
        vaf = float(fs[3])
        ref = fs[4].upper()
        alt = fs[5].upper()
        dp = int(fs[6])
        ref_count = int(fs[7])
        n_candidate = int(fs[8])
        need_alt = int(fs[9])
        dp_bin = fs[10]
        sites.append((chrom, pos, vaf, ref, alt, dp, ref_count, n_candidate, need_alt, dp_bin))

fragment_edit_map = {}
truth = []

for chrom, pos, vaf, ref, alt, dp, ref_count, n_candidate, need_alt, dp_bin in sites:
    info = cache.get((chrom, pos))
    if info is None:
        continue

    candidates = list(info["candidate_fragments"])
    if len(candidates) < need_alt:
        continue

    chosen = random.sample(candidates, need_alt)

    for qname in chosen:
        fragment_edit_map[(chrom, pos, qname)] = alt

    truth.append((chrom, pos, ref, alt, dp, ref_count, len(candidates), need_alt, dp_bin))

edit_by_qname = defaultdict(lambda: defaultdict(dict))
for chrom, pos, qname in fragment_edit_map:
    edit_by_qname[qname][chrom][pos] = fragment_edit_map[(chrom, pos, qname)]

bam = pysam.AlignmentFile(bam_in, "rb")
out = pysam.AlignmentFile(bam_out, "wb", template=bam)

site_edited_bases = defaultdict(int)
edited_reads = 0
edited_bases = 0

for r in bam:
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.query_sequence is None:
        out.write(r)
        continue

    qname = r.query_name
    if qname not in edit_by_qname:
        out.write(r)
        continue

    chrom = r.reference_name
    if chrom is None or chrom not in edit_by_qname[qname]:
        out.write(r)
        continue

    if r.query_qualities is None:
        out.write(r)
        continue

    target = edit_by_qname[qname][chrom]
    seq = list(r.query_sequence)
    quals = r.query_qualities
    changed = 0

    for qp, rp in r.get_aligned_pairs(matches_only=False):
        if qp is None or rp is None:
            continue

        pos = rp + 1
        if pos not in target:
            continue

        
        if quals[qp] < min_bq:
            continue

        alt = target[pos]
        if qp < len(seq) and seq[qp].upper() in "ACGT" and seq[qp].upper() != alt:
            seq[qp] = alt
            changed += 1
            site_edited_bases[(chrom, pos)] += 1

    if changed > 0:
        r.query_sequence = "".join(seq)
        r.query_qualities = quals
        edited_reads += 1
        edited_bases += changed

    out.write(r)

out.close()
bam.close()

with open(truth_out, "w") as w:
    w.write("chrom\tpos\tref\talt\teffective_dp_q20\tref_count_q20\tcandidate_fragments_q20\tmut_fragments_planned\tmut_bases_edited\treal_vaf\tdp_bin\n")
    for chrom, pos, ref, alt, dp, ref_count, n_candidate, need_alt, dp_bin in truth:
        mut_bases = site_edited_bases.get((chrom, pos), 0)
        real_vaf = mut_bases / dp if dp > 0 else 0
        w.write("\t".join(map(str, [
            chrom, pos, ref, alt,
            dp, ref_count, n_candidate,
            need_alt, mut_bases, real_vaf, dp_bin
        ])) + "\n")

print(f"[INFO] truth_sites={len(truth)}", file=sys.stderr)
print(f"[INFO] edited_reads={edited_reads}", file=sys.stderr)
print(f"[INFO] edited_bases={edited_bases}", file=sys.stderr)
PY

###########################################################
# Manifest
###########################################################

manifest=${OUTDIR}/summary/prepare_manifest.tsv
if [[ ! -s ${manifest} ]]; then
    echo -e "Sample\tMode\tVAF\tDP_bin\tOldCandidateSites\tOldSubsetBam\tEffectiveCache\tSelectedSites\tSpikeinBam\tCalmdBam\tTruth" > ${manifest}
fi

###########################################################
# Main loop
###########################################################

for sample in "${samples[@]}"; do
    echo "========== ${sample} =========="

    old_candidate_sites=${OLDOUT}/sites/${sample}.candidate_sites.txt
    old_subset_bam=${OLDOUT}/bam/${sample}.candidate.original_subset.bam
    old_cache=${OLDOUT}/cache/${sample}.site_info.pkl
    eff_cache=${OUTDIR}/cache/${sample}.site_info.effective_q${MIN_BQ}.pkl

    [[ -s ${old_candidate_sites} ]] || { echo "[SKIP] no old candidate sites: ${old_candidate_sites}"; continue; }
    [[ -s ${old_subset_bam} ]] || { echo "[SKIP] no old subset bam: ${old_subset_bam}"; continue; }
    [[ -s ${old_subset_bam}.bai ]] || samtools index ${old_subset_bam}
    [[ -s ${old_cache} ]] || { echo "[SKIP] no old cache: ${old_cache}"; continue; }

    #######################################################
    # Step 1
    
    
    #######################################################

    if [[ -s ${eff_cache} ]]; then
        echo "[1] effective cache exists, skip: ${eff_cache}"
    else
        echo "[1] build effective cache from old 50000 sites"
        python ${OUTDIR}/scripts/build_effective_cache.py ${old_cache} ${old_subset_bam} ${eff_cache} \
            ${MIN_BQ} ${MIN_MQ} > ${OUTDIR}/log/${sample}.build_effective_cache.log 2>&1
    fi

    #######################################################
    # Step 2-5
    
    
    #######################################################

    for bin in "${DP_BINS[@]}"; do
        tag=${sample}.vaf${FIXED_VAF}.dpbin${bin}.q${MIN_BQ}

        selected_sites=${OUTDIR}/sites/${tag}.selected_sites.txt
        spikein=${OUTDIR}/bam/${tag}.spikein.bam
        spikein_sorted=${OUTDIR}/bam/${tag}.spikein.sorted.bam
        spikein_calmd=${OUTDIR}/bam/${tag}.spikein.sorted.calmd.bam
        truth=${OUTDIR}/truth/${tag}.truth.tsv

        echo "[2] ${sample} DPbin=${bin}"

        if [[ -s ${selected_sites} ]]; then
            echo "[2] selected sites exists, skip: ${selected_sites}"
        else
            echo "[2] select sites by effective DP bin"
            python ${OUTDIR}/scripts/select_sites_by_effbin.py ${eff_cache} ${bin} ${N_SITE} ${FIXED_VAF} \
                ${MINMUTREADS} ${selected_sites} \
                $((SEED + ${#bin})) \
                > ${OUTDIR}/log/${tag}.select_sites.log 2>&1
        fi

        n_selected=$(wc -l < ${selected_sites})
        echo "[INFO] selected sites: ${n_selected}"

        [[ ${n_selected} -gt 0 ]] || { echo "[SKIP] no selected sites: ${tag}"; continue; }

        if [[ -s ${spikein_calmd} && -s ${spikein_calmd}.bai && -s ${truth} ]]; then
            echo "[3] calmd spike-in BAM exists, skip: ${spikein_calmd}"
        else
            echo "[3] direct spike-in"
            python ${OUTDIR}/scripts/direct_spike_from_effective_cache.py ${old_subset_bam} \
                ${selected_sites} ${spikein} ${truth} \
                $((SEED + ${#bin})) \
                ${eff_cache} ${MIN_BQ} > ${OUTDIR}/log/${tag}.spikein.log 2>&1

            [[ -s ${spikein} ]] || { echo "[ERROR] spike-in failed: ${tag}"; continue; }

            echo "[4] sort/index"
            samtools sort -@ ${THREADS} -o ${spikein_sorted} ${spikein}
            samtools index ${spikein_sorted}
            samtools quickcheck ${spikein_sorted}

            echo "[5] calmd/index"
            samtools calmd -@ ${THREADS} -b ${spikein_sorted} ${REF} > ${spikein_calmd}
            samtools index ${spikein_calmd}
            samtools quickcheck ${spikein_calmd}
        fi

        n_truth=$(($(wc -l < ${truth}) - 1))
        echo "[INFO] truth sites: ${n_truth}"

        if ! awk -v s="${sample}" -v v="${FIXED_VAF}" -v b="${bin}" '$1==s && $2=="DPBIN_Q20" && $3==v && $4==b{found=1} END{exit found?0:1}' ${manifest}; then
            echo -e "${sample}\tDPBIN_Q20\t${FIXED_VAF}\t${bin}\t${old_candidate_sites}\t${old_subset_bam}\t${eff_cache}\t${selected_sites}\t${spikein_sorted}\t${spikein_calmd}\t${truth}" >> ${manifest}
        fi
    done

    echo "[DONE] ${sample}"
done
'''
NATIVE['217']=r'''
#!/bin/bash
export PATH=/home/tantingting/scmut/tools:$PATH
REF=/home/tantingting/scmut/ref/bwa_ref/genome.fa

REGION_FILE=/home/tantingting/scmut/ref/bwa_ref/10w_regions_1.txt
VARSCAN_JAR=/home/tantingting/scmut/tools/VarScan.v2.4.6.jar
STRELKA_BIN=/home/tantingting/scmut/tools/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py

MONO_PATH=/home/tantingting/scmut/tools/Monopogen
PANEL_PATH=/home/tantingting/scmut/ref/1000g_panel/
MONO_REGION=/mnt/storage/tantingting/scmut/vcf/2_melanoma_scRNA/monopo/region.lst
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MONO_PATH}/apps

SIM_BASE=/mnt/storage/tantingting/scmut/simu/direct_spike_vaf_dp40_final
BAM_DIR=${SIM_BASE}/bam
OUT_BASE=${SIM_BASE}/caller
ORIG_MANIFEST=/mnt/storage/tantingting/scmut/isec/vcf_manifest_pool_D1D2.tsv

THREADS=16
CHR_LIST=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22

mkdir -p ${OUT_BASE}/{freebayes,strelka2,varscan2,monopogen,newonly,tmp,log,summary}

get_sample() {
    tag=$1
    echo "$tag" | sed -E 's/\.vaf[0-9.]+.*//; s/\.dp[0-9]+.*//'
}

get_orig_vcf() {
    sample=$1
    tool=$2
    awk -v s="$sample" -v t="$tool" '
    BEGIN{FS="[ \t]+"}
    NR>1 && $1==s && tolower($2)==tolower(t) {print $3; exit}
    ' "$ORIG_MANIFEST"
}

subtract_original() {
    tag=$1
    sample=$2
    tool=$3
    sim_vcf=$4

    out=${OUT_BASE}/newonly/${tag}.${tool}.newonly.vcf.gz
    tmp=${OUT_BASE}/tmp/${tag}.${tool}
    mkdir -p "$tmp"

    orig_vcf=$(get_orig_vcf "$sample" "$tool" || true)

    sim_norm=${tmp}/sim.norm.vcf.gz
    orig_norm=${tmp}/orig.norm.vcf.gz

    bcftools view -r "$CHR_LIST" "$sim_vcf" | bcftools norm -f "$REF" -m -both | \
    bcftools view -v snps -Oz -o "$sim_norm"
    tabix -f -p vcf "$sim_norm"

    if [[ -z "${orig_vcf}" || ! -s "${orig_vcf}" ]]; then
        
		echo "[WARN] original VCF not found: ${sample} ${tool}, keep simulated VCF only"
        cp "$sim_norm" "$out"
        cp "${sim_norm}.tbi" "${out}.tbi"
    else
        
		bcftools view -r "$CHR_LIST" "$orig_vcf" | bcftools norm -f "$REF" -m -both | \
        bcftools view -v snps -Oz -o "$orig_norm"
        tabix -f -p vcf "$orig_norm"

        bcftools isec -c none -n=1 -w1 "$sim_norm" "$orig_norm" -Oz -o "$out"
        tabix -f -p vcf "$out"
    fi

    echo -e "${tag}\t${sample}\t${tool}\t${out}\t$(bcftools view -H "$out" | wc -l)" >> ${OUT_BASE}/summary/newonly_counts.tsv
    echo "[NEWONLY] ${tool} ${tag}: $(bcftools view -H "$out" | wc -l)"
}

run_freebayes() {
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate freebayes

    tag=$1
    sample=$2
    bam=$3

    out=${OUT_BASE}/freebayes
    raw=${out}/${tag}.freebayes.raw.vcf.gz
    pass=${out}/${tag}.freebayes.pass.vcf.gz

    if [[ ! -s "$raw" || ! -s "${raw}.tbi" ]]; then
        echo "[FreeBayes CALL] $tag"
        freebayes-parallel ${REGION_FILE} ${THREADS} -F 0.05 --min-coverage 10 -m 20 -q 20 \
          -f ${REF} ${bam} | bgzip -c > ${raw}
        tabix -f -p vcf "$raw"
    else
        echo "[FreeBayes SKIP CALL] $tag"
    fi

    if [[ ! -s "$pass" || ! -s "${pass}.tbi" ]]; then
        echo "[FreeBayes FILTER] $tag"
        bcftools norm -f "$REF" -m -both "$raw" | \
        bcftools view -r "$CHR_LIST" -v snps -i 'FMT/GT!="0/0" && QUAL>20' -Oz -o "$pass"
        tabix -f -p vcf "$pass"
    fi

    subtract_original "$tag" "$sample" "Freebayes" "$pass"
}

run_strelka2() {
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate py2
	tag=$1
    sample=$2
    bam=$3

    out=${OUT_BASE}/strelka2
    run_dir=${out}/${tag}
    raw=${run_dir}/results/variants/variants.vcf.gz
    pass=${out}/${tag}.strelka2.pass.vcf.gz

    if [[ ! -s "$pass" || ! -s "${pass}.tbi" ]]; then
        echo "[Strelka2 CALL] $tag"
        rm -rf "$run_dir"

        "$STRELKA_BIN" --bam "$bam" --rna --referenceFasta "$REF" --runDir "$run_dir"

        python ${run_dir}/runWorkflow.py -m local -j ${THREADS} > ${OUT_BASE}/log/${tag}.strelka2.log 2>&1

        echo "[Strelka2 FILTER] $tag"
        bcftools view -r "$CHR_LIST" -i 'ALT!="." && FILTER="PASS"' "$raw" | \
        bcftools norm -f "$REF" -m -both | \
        bcftools view -v snps -i 'FMT/GT!="0/0" && FMT/GT!="0|0" && QUAL>20' -Oz -o "$pass"
        tabix -f -p vcf "$pass"
    else
        echo "[Strelka2 SKIP] $tag"
    fi

    subtract_original "$tag" "$sample" "Strelka2" "$pass"
}

run_varscan2() {
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate scmut
	tag=$1
    sample=$2
    bam=$3

    out=${OUT_BASE}/varscan2
    raw=${out}/${tag}.varscan2.vcf
    rawgz=${raw}.gz
    pass=${OUT_BASE}/varscan2/${tag}.varscan2.pass.vcf.gz

    if [[ ! -s "$rawgz" || ! -s "${rawgz}.tbi" ]]; then
        echo "[VarScan2 CALL] $tag"
        samtools mpileup -Q 20 -q 20 -f ${REF} ${bam} | java -jar ${VARSCAN_JAR} mpileup2snp \
          --min-coverage 10 --output-vcf 1 --strand-filter 0 --p-value 0.05 --min-avg-qual 20 > ${raw}

        bgzip -f "$raw"
        tabix -f -p vcf "$rawgz"
    else
        echo "[VarScan2 SKIP CALL] $tag"
    fi

    if [[ ! -s "$pass" || ! -s "${pass}.tbi" ]]; then
        echo "[VarScan2 FILTER] $tag"
        bcftools view -r "$CHR_LIST" "$rawgz" | bcftools norm -f "$REF" -m -both | \
        bcftools view -v snps -i 'FMT/GT!="0/0"' -Oz -o "$pass"
        tabix -f -p vcf "$pass"
    fi

    subtract_original "$tag" "$sample" "VarScan2" "$pass"
}

run_monopogen() {
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate Monopogen
	tag=$1
    sample=$2
    bam=$3

    out=${OUT_BASE}/monopogen/${tag}
    bam_lst=${out}/bam.lst
    pass=${OUT_BASE}/monopogen/${tag}.monopogen.pass.vcf.gz
    mkdir -p "$out" ${OUT_BASE}/monopogen

    if [[ ! -s "$pass" || ! -s "${pass}.tbi" ]]; then
        echo "[Monopogen CALL] $tag"
        echo -e "${tag},${bam}" > "$bam_lst"

        python ${MONO_PATH}/src/Monopogen.py preProcess -b "$bam_lst" -o "$out" -a ${MONO_PATH}/apps \
            -t ${THREADS}

        python ${MONO_PATH}/src/Monopogen.py germline -a ${MONO_PATH}/apps -t ${THREADS} -s all \
            -p ${PANEL_PATH} -r ${MONO_REGION} -g ${REF} -o "$out"

        gl_all=${out}/${tag}.gl.all.vcf.gz
        phased_all=${out}/${tag}.phased.all.vcf.gz
        gl_pass=${out}/${tag}.gl.depth_pass.vcf.gz
        keep=${out}/${tag}.keep.keys.txt

        bcftools concat -a -Oz -o "$gl_all" \
            $(find "${out}/germline" -maxdepth 1 -name "chr*.gl.vcf.gz" | sort -V)
        tabix -f -p vcf "$gl_all"

        bcftools concat -a -Oz -o "$phased_all" \
            $(find "${out}/germline" -maxdepth 1 -name "chr*.phased.vcf.gz" | sort -V)
        tabix -f -p vcf "$phased_all"

        bcftools view -v snps "$gl_all" -Oz -o "$gl_pass"
        tabix -f -p vcf "$gl_pass"

        bcftools query -f '%CHROM:%POS:%REF:%ALT\n' "$gl_pass" | sort -u > "$keep"

        bcftools view -T <(awk -F: 'BEGIN{OFS="\t"}{print $1,$2,$2,$3,$4}' "$keep") "$phased_all" | \
        bcftools view -r "$CHR_LIST" -v snps -i 'GT!="0/0" && GT!="0|0" && GT!="./." && GT!=".|."' \
            -Oz -o "$pass"
        tabix -f -p vcf "$pass"

        rm -f "$gl_all" "$gl_all.tbi" "$phased_all" "$phased_all.tbi" "$gl_pass" "$gl_pass.tbi"
    else
        echo "[Monopogen SKIP] $tag"
    fi

    subtract_original "$tag" "$sample" "Monopogen" "$pass"
}

echo -e "Tag\tSample\tTool\tNewOnlyVCF\tN_newonly" > ${OUT_BASE}/summary/newonly_counts.tsv

#for bam in ${BAM_DIR}/*.spikein.sorted.calmd.bam; do
#    [[ -s "$bam" ]] || continue
#
#    tag=$(basename "$bam" .spikein.sorted.calmd.bam)
#    sample=$(get_sample "$tag")
#
#    echo
#    echo "======================================"
#    echo "BAM:    $bam"
#    echo "TAG:    $tag"
#    echo "SAMPLE: $sample"
#    echo "======================================"
#
#    run_freebayes "$tag" "$sample" "$bam"
#    run_strelka2 "$tag" "$sample" "$bam"
#    run_varscan2 "$tag" "$sample" "$bam"
#    run_monopogen "$tag" "$sample" "$bam"
#
#    echo "[DONE] $tag"
#done
#
#echo
#echo "ALL DONE"
#echo "Summary: ${OUT_BASE}/summary/newonly_counts.tsv"

############################################
# Test one BAM only
############################################

bam=${BAM_DIR}/Acral1.vaf0.05.spikein.sorted.calmd.bam

if [[ ! -s "$bam" ]]; then
    echo "[ERROR] test BAM not found: $bam"
    exit 1
fi

tag=$(basename "$bam" .spikein.sorted.calmd.bam)
sample=$(get_sample "$tag")

echo
echo "======================================"
echo "[TEST ONE BAM]"
echo "BAM:    $bam"
echo "TAG:    $tag"
echo "SAMPLE: $sample"
echo "======================================"

run_freebayes "$tag" "$sample" "$bam"
run_strelka2 "$tag" "$sample" "$bam"
run_varscan2 "$tag" "$sample" "$bam"
run_monopogen "$tag" "$sample" "$bam"

echo
echo "[TEST DONE] $tag"
echo "Summary: ${OUT_BASE}/summary/newonly_counts.tsv"
'''
NATIVE['218']=r'''
#!/bin/bash
export PATH=/home/tantingting/scmut/tools:$PATH

REF=/home/tantingting/scmut/ref/bwa_ref/genome.fa
SIM_BASE=/mnt/storage/tantingting/scmut/simu/direct_spike_vaf_dp40_final
BAM_DIR=${SIM_BASE}/bam
OUT_BASE=${SIM_BASE}/caller
ORIG_MANIFEST=/mnt/storage/tantingting/scmut/isec/vcf_manifest_pool_D1D2.tsv

THREADS=16
JAVA_OPTS="-Xmx24g"
CHR_LIST=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22

mkdir -p ${OUT_BASE}/{hap,newonly,log,summary,tmp}

get_orig_vcf() {
    sample=$1
    tool=$2
    awk -v s="$sample" -v t="$tool" '
    BEGIN{FS="[ \t]+"}
    NR>1 && $1==s && tolower($2)==tolower(t) {print $3; exit}
    ' "$ORIG_MANIFEST"
}

subtract_original() {
    tag=$1
    sample=$2
    tool=$3
    sim_vcf=$4

    out=${OUT_BASE}/newonly/${tag}.${tool}.newonly.vcf.gz
    orig_vcf=$(get_orig_vcf "$sample" "$tool" || true)

    if [[ -z "$orig_vcf" || ! -s "$orig_vcf" ]]; then
        echo "[WARN] original VCF not found: ${sample} ${tool}, keep simulated VCF only"
        cp "$sim_vcf" "$out"
        cp "${sim_vcf}.tbi" "${out}.tbi"
    else
        bcftools isec -c none -n=1 -w1 "$sim_vcf" "$orig_vcf" -Oz -o "$out"
        tabix -f -p vcf "$out"
    fi

    n=$(bcftools view -H "$out" | wc -l)
    echo -e "${tag}\t${sample}\t${tool}\t${out}\t${n}" >> ${OUT_BASE}/summary/newonly_counts.gatk_hap.tsv
    echo "[NEWONLY] ${tool} ${tag}: ${n}"
}

run_gatk_hap() {
    tag=$1
    sample=$2
    bam=$3

    out=${OUT_BASE}/gatk_hap/${tag}
    mkdir -p "$out"

    split_bam=${out}/${tag}.split.bam
    raw_vcf=${out}/${tag}.gatk.raw.vcf.gz
    pass_vcf=${out}/${tag}.gatk.pass.vcf.gz
    norm_vcf=${out}/${tag}.gatk.norm.vcf.gz
    snp_vcf=${out}/${tag}.gatk.snp.vcf.gz
    filt_vcf=${out}/${tag}.gatk.snp.filtered.vcf.gz

    if [[ ! -s "$split_bam" || ! -s "${split_bam}.bai" ]]; then
        echo "[GATK SplitNCigarReads] $tag"
        gatk --java-options "$JAVA_OPTS" SplitNCigarReads -R "$REF" -I "$bam" \
            -O "$split_bam" --create-output-bam-index true
    else
        echo "[GATK SKIP SPLIT] $tag"
    fi

    if [[ ! -s "$raw_vcf" || ! -s "${raw_vcf}.tbi" ]]; then
        echo "[GATK HaplotypeCaller] $tag"
        gatk --java-options "$JAVA_OPTS" HaplotypeCaller -R "$REF" -I "$split_bam" -O "$raw_vcf" \
            --native-pair-hmm-threads "$THREADS" --standard-min-confidence-threshold-for-calling 30.0
    else
        echo "[GATK SKIP CALL] $tag"
    fi

    if [[ ! -s "$pass_vcf" || ! -s "${pass_vcf}.tbi" ]]; then
        echo "[GATK FILTER SNP] $tag"

        bcftools view -r "$CHR_LIST" "$raw_vcf" -Oz -o ${out}/${tag}.chr.vcf.gz
        tabix -f -p vcf ${out}/${tag}.chr.vcf.gz

        bcftools norm -f "$REF" -m -both ${out}/${tag}.chr.vcf.gz -Oz -o "$norm_vcf"
        tabix -f -p vcf "$norm_vcf"

        gatk SelectVariants -R "$REF" -V "$norm_vcf" -select-type SNP -O "$snp_vcf"

        gatk VariantFiltration -R "$REF" -V "$snp_vcf" -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "(vc.hasAttribute('MQRankSum') && MQRankSum < -12.5)" --filter-name "MQRS12.5" \
            -filter "(vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)" --filter-name "RPRS8" \
            -O "$filt_vcf"

        bcftools view -f PASS "$filt_vcf" -Oz -o "$pass_vcf"
        tabix -f -p vcf "$pass_vcf"
    else
        echo "[GATK SKIP FILTER] $tag"
    fi

    subtract_original "$tag" "$sample" "GATK" "$pass_vcf"
}

echo -e "Tag\tSample\tTool\tNewOnlyVCF\tN_newonly" > ${OUT_BASE}/summary/newonly_counts.gatk_hap.tsv

BATCH=${1:-1}

if [[ "$BATCH" == "1" ]]; then
    samples=(Acral1 Acral2 Acral3pre Acral3post)
elif [[ "$BATCH" == "2" ]]; then
    samples=(Acral4 Cutaneous1 Cutaneous1Lym Cutaneous2)
elif [[ "$BATCH" == "3" ]]; then
    samples=(P4_cSCC P5_cSCC P6_cSCC P7_cSCC P9_cSCC P10_cSCC)
else
    echo "[ERROR] BATCH must be 1, 2, or 3"
    exit 1
fi

for sample in "${samples[@]}"; do
    for bam in ${BAM_DIR}/${sample}*.spikein.sorted.calmd.bam; do
        [[ -s "$bam" ]] || continue

        tag=$(basename "$bam" .spikein.sorted.calmd.bam)

        echo
        echo "======================================"
        echo "BATCH:  $BATCH"
        echo "BAM:    $bam"
        echo "TAG:    $tag"
        echo "SAMPLE: $sample"
        echo "======================================"

        run_gatk_hap "$tag" "$sample" "$bam"

        echo "[DONE] $tag"
    done
done

echo
echo "ALL DONE batch ${BATCH}"
echo "Summary: ${OUT_BASE}/summary/newonly_counts.gatk_hap.tsv"
'''
NATIVE['smartseq_vaf']=r'''
#!/bin/bash
set -euo pipefail

REF=/home/tantingting/scmut/ref/bwa_ref/genome.fa
REGION_DIR=/mnt/storage/tantingting/scmut_smart/simu/region/vaf_sim_region
OUTDIR=/mnt/storage/tantingting/scmut_smart/simu/direct_spike_vaf_dp40_final
MANIFEST=/mnt/storage/tantingting/scmut_smart/bam/star_qc_unique_ge10k_scrna_sample_wes.tsv

N_SITE=3000
FLANK=0
MINDEPTH=40
MINMUTREADS=2
THREADS=16
SEED=1234
VAFS=(0.05 0.1 0.25 0.5 0.75)

mkdir -p ${OUTDIR}/{scripts,sites,window,bam,log,truth,summary}
mkdir -p ${OUTDIR}/cache


samples=("$@")

if [[ ${#samples[@]} -eq 0 ]]; then
    echo "[ERROR] No sample specified"
    exit 1
fi

get_bam(){

    srr=$1

    sample=$(awk -v s="$srr" '$1==s{print $2}' ${MANIFEST})

    if [[ $sample == PT* ]]; then
        echo "/mnt/storage/tantingting/scmut_smart/bam/1_TNBC_scRNA/md_bam/${srr}_sortmd.bam"
    elif [[ $sample == BC* ]]; then
        echo "/mnt/storage/tantingting/scmut_smart/bam/2_BRCA_scRNA/md_bam/${srr}_sortmd.bam"
    else
        echo "NA"
    fi
}

############################################

############################################

cat > ${OUTDIR}/scripts/make_master_sites.py << 'PY'
import sys, random, subprocess, bisect
import os
import pickle

region_bed, ref, n_site, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
random.seed(seed)
bases = ["A", "C", "G", "T"]
regions, cum = [], []
total = 0

with open(region_bed) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, start, end = line.strip().split()[:3]
        start, end = int(start), int(end)
        if end <= start:
            continue
        total += end - start
        regions.append((chrom, start, end))
        cum.append(total)

if total == 0:
    sys.exit("No valid regions")

seen, tries, max_tries = set(), 0, n_site * 200

with open(out, "w") as w:
    while len(seen) < n_site and tries < max_tries:
        tries += 1
        i = bisect.bisect_left(cum, random.randint(1, total))
        chrom, start, end = regions[i]
        pos = random.randint(start, end - 1) + 1
        key = (chrom, pos)

        if key in seen:
            continue

        try:
            fa = subprocess.check_output(["samtools", "faidx", ref, f"{chrom}:{pos}-{pos}"], text=True).splitlines()
            ref_base = "".join(fa[1:]).upper()
        except Exception:
            continue

        if len(ref_base) != 1 or ref_base not in bases:
            continue

        alt = random.choice([b for b in bases if b != ref_base])
        seen.add(key)
        w.write(f"{chrom}\t{pos}\t{pos}\t{ref_base}\t{alt}\n")
PY

############################################




############################################

cat > ${OUTDIR}/scripts/direct_spike_bam.py << 'PY'
import sys, random, pysam, math
import bisect
from collections import defaultdict
import os
import pickle
from collections import defaultdict

bam_in, sites_file, ref_fa, bam_out, truth_out, min_depth, min_mutreads, seed, cache_file = sys.argv[1:10]
min_depth, min_mutreads, seed = int(min_depth), int(min_mutreads), int(seed)
random.seed(seed)

bam = pysam.AlignmentFile(bam_in, "rb")

sites = []
with open(sites_file) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        fs = line.strip().split()
        chrom, pos, end, vaf, ref_base, alt = fs[0], int(fs[1]), int(fs[2]), float(fs[3]), fs[4].upper(), fs[5].upper()
        sites.append((chrom, pos, vaf, ref_base, alt))
sites.sort(key=lambda x: (x[0], x[1]))

# fragment_edit_map:
#   key = (chrom, pos, query_name)
#   value = alt

fragment_edit_map = {}
truth = []


############################################




############################################

##################################################

##################################################

if os.path.exists(cache_file):

    print(
        f"[INFO] load cache: {cache_file}",
        file=sys.stderr
    )

    with open(cache_file, "rb") as f:
        site_info = pickle.load(f)


    for chrom, pos, vaf, ref_base, alt in sites:

        if (chrom, pos) not in site_info:
            continue

        site_info[(chrom, pos)]["vaf"] = vaf

else:

    print(
        f"[INFO] build cache: {cache_file}",
        file=sys.stderr
    )
    site_info = {}


    for chrom, pos, vaf, ref_base, alt in sites:
        site_info[(chrom, pos)] = {
            "vaf": vaf,
            "ref": ref_base,
            "alt": alt,
            "cover": 0,
            "candidates": set()
        }
    chr_sites = defaultdict(set)

    for chrom, pos, vaf, ref_base, alt in sites:
        chr_sites[chrom].add(pos)

    print(
        f"[INFO] {len(site_info)} sites loaded",
        file=sys.stderr
    )

############################################

############################################

    for chrom in sorted(chr_sites.keys()):

        target_pos = sorted(chr_sites[chrom])


        start = min(target_pos) - 1
        end = max(target_pos)+1

        print(
            f"[INFO] scanning {chrom}",
            file=sys.stderr
        )

        try:
            iterator = bam.fetch(chrom, start, end)
        except ValueError:
            continue

        n_read = 0

        for r in iterator:

            if r.is_unmapped:
                continue

            if r.is_secondary:
                continue

            if r.is_supplementary:
                continue

            seq = r.query_sequence

            if seq is None:
                continue

            n_read += 1

        ####################################

        ####################################

            aligned = {
                rp + 1: qp
                for qp, rp in r.get_aligned_pairs(matches_only=False)
                if qp is not None and rp is not None
            }

            if not aligned:
                continue

        #read_start = min(aligned)
        #read_end = max(aligned)

            read_start = r.reference_start + 1
            read_end = r.reference_end
            if read_end is None:
                continue

        ####################################

        ####################################

            left = bisect.bisect_left(
                target_pos,
                read_start
            )

            right = bisect.bisect_right(
                target_pos,
                read_end
            )

            covered_sites = target_pos[left:right]

            if not covered_sites:
                continue

            for pos in covered_sites:

                qp = aligned.get(pos)

                if qp is None:
                    continue

                base = seq[qp].upper()

                if base not in "ACGT":
                    continue

                info = site_info[(chrom, pos)]

                info["cover"] += 1

                if base == info["ref"]:
                    info["candidates"].add(
                        r.query_name
                    )

        print(
            f"[INFO] {chrom}: {n_read} reads scanned",
            file=sys.stderr
        )
        
##################################################

##################################################

    with open(cache_file, "wb") as f:
        pickle.dump(site_info, f)
    print(
        f"[INFO] cache saved: {cache_file}",
        file=sys.stderr
    )

############################################

############################################

fragment_edit_map = {}
truth = []

for (chrom, pos), info in site_info.items():

    dp = info["cover"]

    if dp < min_depth:
        continue

    candidates = list(
        info["candidates"]
    )

    if len(candidates) == 0:
        continue

    target_alt_reads = max(
        min_mutreads,
        math.ceil(dp * info["vaf"])
    )

    n_choose = min(
        len(candidates),
        target_alt_reads
    )

    chosen_names = random.sample(
        candidates,
        n_choose
    )

    for name in chosen_names:
        fragment_edit_map[
            (chrom, pos, name)
        ] = info["alt"]

    truth.append(
        (
            chrom,
            pos,
            info["ref"],
            info["alt"],
            dp,
            len(candidates),
            n_choose,
            "NA"
        )
    )



print(
    f"[INFO] fragment_edit_map={len(fragment_edit_map)}",
    file=sys.stderr
)

############################################

#




#


#       (chrom, pos, query_name) -> alt
#

#       query_name -> chrom -> pos -> alt
#


############################################

bam.close()
bam = pysam.AlignmentFile(bam_in, "rb")
out = pysam.AlignmentFile(bam_out, "wb", template=bam)

edited_reads = 0
edited_bases = 0
site_edited_bases = defaultdict(int)

############################################

############################################

edit_by_qname = defaultdict(lambda: defaultdict(dict))

for chrom, pos, qname in fragment_edit_map:
    alt = fragment_edit_map[(chrom, pos, qname)]
    edit_by_qname[qname][chrom][pos] = alt

print(
    f"[INFO] edit_by_qname={len(edit_by_qname)}",
    file=sys.stderr
)

############################################

############################################

for r in bam:

    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.query_sequence is None:
        out.write(r)
        continue

    qname = r.query_name

    ########################################


    ########################################

    if qname not in edit_by_qname:
        out.write(r)
        continue

    chrom = r.reference_name

    if chrom is None:
        out.write(r)
        continue

    ########################################


    ########################################

    if chrom not in edit_by_qname[qname]:
        out.write(r)
        continue

    target_edits = edit_by_qname[qname][chrom]

    seq = list(r.query_sequence)
    quals = r.query_qualities
    changed = 0

    ########################################

    ########################################

    for qp, rp in r.get_aligned_pairs(matches_only=False):

        if qp is None or rp is None:
            continue

        pos = rp + 1

        if pos not in target_edits:
            continue

        alt = target_edits[pos]

        if qp < len(seq) and seq[qp].upper() in "ACGT" and seq[qp].upper() != alt:
            seq[qp] = alt
            changed += 1
            site_edited_bases[(chrom, pos)] += 1

    if changed > 0:
        r.query_sequence = "".join(seq)
        r.query_qualities = quals
        edited_reads += 1
        edited_bases += changed

    out.write(r)

out.close()
bam.close()

############################################



# real_vaf = mut_bases_edited / dp
############################################

with open(truth_out, "w") as w:
    w.write("chrom\tpos\tref\talt\tdp\tcandidate_fragments\tmut_fragments_planned\tmut_bases_edited\treal_vaf\n")
    for chrom, pos, ref_base, alt, dp, n_candidate, n_choose, _ in truth:
        mut_bases = site_edited_bases.get((chrom, pos), 0)
        real_vaf = mut_bases / dp if dp > 0 else 0
        w.write("\t".join(map(str, [chrom, pos, ref_base, alt, dp, n_candidate, n_choose, mut_bases, real_vaf])) + "\n")

sys.stderr.write(f"[INFO] sites_input={len(sites)} truth_sites={len(truth)} edited_reads={edited_reads} edited_bases={edited_bases}\n")
PY

manifest=${OUTDIR}/summary/prepare_manifest.tsv

if [[ ! -s ${manifest} ]]; then
    echo -e "Sample\tMode\tVAF\tMasterSites\tSites\tWindowBed\tOriginalSubsetBam\tSpikeinBam\tCalmdBam\tTruth" > ${manifest}
fi

############################################

############################################

for sample in "${samples[@]}"; do
    echo "========== ${sample} =========="

    bam=$(get_bam ${sample})
    region_bed=${REGION_DIR}/${sample}.final_sim_region.bed
    master_sites=${OUTDIR}/sites/${sample}.master_sites.txt
    sample_window=${OUTDIR}/window/${sample}.master.window.bed
    sample_subset=${OUTDIR}/bam/${sample}.master.original_subset.bam

    [[ "${bam}" != "NA" ]] || { echo "[SKIP] unknown sample: ${sample}"; continue; }
    [[ -f ${region_bed} ]] || { echo "[SKIP] no region: ${region_bed}"; continue; }
    [[ -f ${bam} ]] || { echo "[SKIP] no bam: ${bam}"; continue; }
    [[ -f ${bam}.bai ]] || samtools index ${bam}

    if [[ -s ${master_sites} ]]; then
        echo "[1] master sites exists, skip: ${master_sites}"
    else
        echo "[1] make master sites"
        python ${OUTDIR}/scripts/make_master_sites.py ${region_bed} ${REF} ${N_SITE} ${master_sites} ${SEED}
    fi

    n_master=$(wc -l < ${master_sites})
    echo "[INFO] master sites: ${n_master}"
    [[ ${n_master} -gt 0 ]] || { echo "[SKIP] no master sites: ${sample}"; continue; }

    if [[ -s ${sample_window} ]]; then
        echo "[2] master window exists, skip: ${sample_window}"
    else
        echo "[2] make master window"
        awk -v OFS="\t" -v flank=${FLANK} '{s=$2-1-flank; e=$3+flank; if(s<0)s=0; print $1,s,e}' ${master_sites} \
            | sort -k1,1 -k2,2n \
            | bedtools merge \
            > ${sample_window}
    fi

    if [[ -s ${sample_subset} && -s ${sample_subset}.bai ]]; then
        echo "[3] master subset BAM exists, skip: ${sample_subset}"
    else
        echo "[3] subset BAM once"
        samtools view -@ ${THREADS} -P -b -L ${sample_window} ${bam} > ${sample_subset}
        samtools index ${sample_subset}
    fi

    for vaf in "${VAFS[@]}"; do
        tag=${sample}.vaf${vaf}
        cache=${OUTDIR}/cache/${sample}.site_info.pkl
        sites=${OUTDIR}/sites/${tag}.sites.txt
        spikein=${OUTDIR}/bam/${tag}.spikein.bam
        spikein_sorted=${OUTDIR}/bam/${tag}.spikein.sorted.bam
        spikein_calmd=${OUTDIR}/bam/${tag}.spikein.sorted.calmd.bam
        truth=${OUTDIR}/truth/${tag}.truth.tsv

        echo "[4] ${sample} VAF=${vaf}"

        if [[ -s ${sites} ]]; then
            echo "[4] sites exists, skip: ${sites}"
        else
            awk -v OFS="\t" -v vaf=${vaf} '{print $1,$2,$3,vaf,$4,$5}' ${master_sites} > ${sites}
        fi

        if [[ -s ${spikein_calmd} && -s ${spikein_calmd}.bai && -s ${truth} ]]; then
            echo "[5] calmd spike-in BAM exists, skip: ${spikein_calmd}"
        else
            echo "[5] direct spike-in BAM"
            python ${OUTDIR}/scripts/direct_spike_bam.py \
                ${sample_subset} \
                ${sites} \
                ${REF} \
                ${spikein} \
                ${truth} \
                ${MINDEPTH} \
                ${MINMUTREADS} \
                ${SEED} \
                ${cache} \
                > ${OUTDIR}/log/${tag}.log 2>&1

            [[ -s ${spikein} ]] || { echo "[ERROR] direct spike-in failed: check ${OUTDIR}/log/${tag}.log"; continue; }

            echo "[6] sort/index"
            samtools sort -@ ${THREADS} -o ${spikein_sorted} ${spikein}
            samtools index ${spikein_sorted}
            samtools quickcheck ${spikein_sorted}

            echo "[7] calmd/index"
            samtools calmd -@ ${THREADS} -b ${spikein_sorted} ${REF} > ${spikein_calmd}
            samtools index ${spikein_calmd}
            samtools quickcheck ${spikein_calmd}
        fi

        n_truth=$(($(wc -l < ${truth}) - 1))
        echo "[INFO] truth sites: ${n_truth}"

        if ! awk -v s="${sample}" -v v="${vaf}" '$1==s && $2=="VAF" && $3==v{found=1} END{exit found?0:1}' ${manifest}; then
            echo -e "${sample}\tVAF\t${vaf}\t${master_sites}\t${sites}\t${sample_window}\t${sample_subset}\t${spikein_sorted}\t${spikein_calmd}\t${truth}" >> ${manifest}
        fi
    done

    echo "[DONE] ${sample}"
done
'''
NATIVE['smartseq_dp']=r'''
#!/bin/bash

# Smart-seq2 DP effective-bin simulation
#


#    Q>=20
#    A/C/G/T base
#    remove duplicate/secondary/supplementary/qcfail

#    DP5-10
#    DP10-25
#    DP25-50
#    DP>=50



REF=/home/tantingting/scmut/ref/bwa_ref/genome.fa
OUTDIR=/mnt/storage/tantingting/scmut_smart/simu/direct_spike_dp_effbin_q20_vaf04_final
MANIFEST=/mnt/storage/tantingting/scmut_smart/bam/star_qc_unique_ge10k_scrna_sample_wes.tsv
REGION_DIR=/mnt/storage/tantingting/scmut_smart/simu/region/dp_sim_region

N_CANDIDATE=50000
FLANK=0
FIXED_VAF=0.4
N_SITE=3000
MINMUTREADS=2
THREADS=16
SEED=1234
MIN_BQ=20
MIN_MQ=0

DP_BINS=("5_10" "10_25" "25_50" "50plus")
#DP_BINS=("5_10" "10_25")

mkdir -p ${OUTDIR}/{scripts,sites,window,bam,cache,truth,log,summary}

samples=("$@")
if [[ ${#samples[@]} -eq 0 ]]; then
    echo "[ERROR] No sample specified"
    exit 1
fi

get_bam(){

    srr=$1

    sample=$(awk -v s="$srr" '$1==s{print $2}' ${MANIFEST})

    if [[ $sample == PT* ]]; then
        echo "/mnt/storage/tantingting/scmut_smart/bam/1_TNBC_scRNA/md_bam/${srr}_sortmd.bam"
    elif [[ $sample == BC* ]]; then
        echo "/mnt/storage/tantingting/scmut_smart/bam/2_BRCA_scRNA/md_bam/${srr}_sortmd.bam"
    else
        echo "NA"
    fi
}

###########################################################
# Step 1:

#

#   final_sim_region_ge1.bed
#




#

#   candidate_sites.txt
#



###########################################################

cat > ${OUTDIR}/scripts/make_candidate_sites.py << 'PY'
import sys, random, subprocess, bisect

region_bed, ref, n_site, out, seed = (
    sys.argv[1],
    sys.argv[2],
    int(sys.argv[3]),
    sys.argv[4],
    int(sys.argv[5])
)

random.seed(seed)
bases=["A","C","G","T"]

regions=[]
cum=[]
total=0

with open(region_bed) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue

        chrom,start,end=line.strip().split()[:3]
        start=int(start)
        end=int(end)

        if end<=start:
            continue

        total += end-start
        regions.append((chrom,start,end))
        cum.append(total)

seen=set()

with open(out,"w") as w:
    tries=0
    max_tries=n_site*300
    while len(seen)<n_site and tries<max_tries:
        tries+=1
        idx=bisect.bisect_left(
            cum,
            random.randint(1,total)
        )

        chrom,start,end=regions[idx]
        pos=random.randint(start,end-1)+1
        if (chrom,pos) in seen:
            continue
        try:
            fa=subprocess.check_output(
                [
                    "samtools",
                    "faidx",
                    ref,
                    f"{chrom}:{pos}-{pos}"
                ],
                text=True
            ).splitlines()
            ref_base="".join(fa[1:]).upper()
        except:
            continue
        if ref_base not in bases:
            continue
        alt=random.choice(
            [x for x in bases if x!=ref_base]
        )
        seen.add((chrom,pos))
        w.write(
            f"{chrom}\t{pos}\t{pos}\t{ref_base}\t{alt}\n"
        )

print(
    f"[INFO] candidate sites generated: {len(seen)}",
    file=sys.stderr
)
PY


###########################################################

#

#






#   6. Mapping quality >= MIN_MQ
#   7. Base quality >= MIN_BQ

#

#   ref
#   alt



###########################################################


cat > ${OUTDIR}/scripts/build_effective_cache.py << 'PY'
import sys
import pysam
import pickle
from collections import defaultdict

bam_file=sys.argv[1]
candidate_file=sys.argv[2]
out_cache=sys.argv[3]

MIN_BQ=int(sys.argv[4])
MIN_MQ=int(sys.argv[5])

sites={}
chr_sites=defaultdict(set)

with open(candidate_file) as f:
    for line in f:
        chrom,pos,end,ref,alt=line.strip().split()
        pos=int(pos)
        sites[(chrom,pos)]={
            "ref":ref,
            "alt":alt,
            "effective_dp":0,
            "ref_count":0,
            "candidate_fragments":set()
        }
        chr_sites[chrom].add(pos)

bam=pysam.AlignmentFile(
    bam_file,
    "rb"
)
for chrom in chr_sites:
    positions=sorted(chr_sites[chrom])
    pos_set=set(positions)
    print(
        f"[INFO] scanning {chrom}",
        file=sys.stderr
    )
    try:
        reads=bam.fetch(
            chrom,
            min(positions)-1,
            max(positions)
        )
    except:
        continue
    for r in reads:
        if (
            r.is_unmapped
            or r.is_secondary
            or r.is_supplementary
            or r.is_duplicate
            or r.is_qcfail
        ):
            continue
        if r.mapping_quality < MIN_MQ:
            continue
        if (
            r.query_sequence is None
            or r.query_qualities is None
        ):
            continue
        seq=r.query_sequence
        qual=r.query_qualities

        for qp,rp in r.get_aligned_pairs(
            matches_only=False
        ):
            if qp is None or rp is None:
                continue
            pos=rp+1
            if pos not in pos_set:
                continue
            if qual[qp] < MIN_BQ:
                continue
            base=seq[qp].upper()
            if base not in "ACGT":
                continue
            info=sites[(chrom,pos)]
            info["effective_dp"] += 1
            if base == info["ref"]:
                info["ref_count"] += 1
                info["candidate_fragments"].add(
                    r.query_name
                )
bam.close()
with open(out_cache,"wb") as f:

    pickle.dump(
        sites,
        f,
        protocol=pickle.HIGHEST_PROTOCOL
    )
print(
    f"[INFO] saved {out_cache}",
    file=sys.stderr
)
for label,func in [
    ("DP1plus",
     lambda x:x["effective_dp"]>=1),
    ("DP5plus",
     lambda x:x["effective_dp"]>=5),
    ("DP10plus",
     lambda x:x["effective_dp"]>=10),
    ("DP25plus",
     lambda x:x["effective_dp"]>=25),
    ("DP50plus",
     lambda x:x["effective_dp"]>=50),
]:
    print(
        label,
        sum(
            func(x)
            for x in sites.values()
        ),
        file=sys.stderr
    )
PY

###########################################################
# Python 2: select_sites_by_effbin.py
#


#


#   2. candidate_fragments >= ceil(effective_dp * VAF)
#


###########################################################

cat > ${OUTDIR}/scripts/select_sites_by_effbin.py << 'PY'
import sys, pickle, random, math

cache_file, dp_bin, n_site, vaf, min_mutreads, out, seed = sys.argv[1:8]
n_site = int(n_site)
vaf = float(vaf)
min_mutreads = int(min_mutreads)
seed = int(seed)
random.seed(seed)

if dp_bin == "5_10":
    lo, hi = 5,10
elif dp_bin == "10_25":
    lo, hi = 10,25
elif dp_bin == "25_50":
    lo, hi = 25,50
elif dp_bin == "50plus":
    lo, hi = 50,None
else:
    sys.exit(f"Unknown dp_bin: {dp_bin}")

with open(cache_file, "rb") as f:
    cache = pickle.load(f)

eligible = []

for (chrom, pos), info in cache.items():
    dp = info["effective_dp"]
    ref_count = info["ref_count"]
    n_candidate = len(info["candidate_fragments"])
    need_alt = max(min_mutreads, math.ceil(dp * vaf))

    in_bin = dp >= lo if hi is None else lo <= dp < hi

    if in_bin and n_candidate >= need_alt:
        eligible.append((chrom, pos, info["ref"], info["alt"], dp, ref_count, n_candidate, need_alt))

eligible.sort(key=lambda x: (x[0], x[1]))

print(f"[INFO] {dp_bin} eligible={len(eligible)}", file=sys.stderr)

if len(eligible) == 0:
    sys.exit(f"No eligible sites for {dp_bin}")

if len(eligible) > n_site:
    selected = random.sample(eligible, n_site)
    selected.sort(key=lambda x: (x[0], x[1]))
else:
    selected = eligible
    print(f"[WARN] eligible < N_SITE, use all: {len(selected)}", file=sys.stderr)

with open(out, "w") as w:
    for chrom, pos, ref, alt, dp, ref_count, n_candidate, need_alt in selected:
        # chrom pos end vaf ref alt effective_dp_q20 ref_count_q20 candidate_fragments_q20 need_alt dp_bin
        w.write(f"{chrom}\t{pos}\t{pos}\t{vaf}\t{ref}\t{alt}\t{dp}\t{ref_count}\t{n_candidate}\t{need_alt}\t{dp_bin}\n")

print(f"[INFO] selected={len(selected)}", file=sys.stderr)
PY

###########################################################
# Python 3: direct_spike_from_effective_cache.py
#


#






#

#   effective_dp_q20
#   ref_count_q20
#   candidate_fragments_q20
#   mut_fragments_planned
#   mut_bases_edited
#   real_vaf = mut_bases_edited / effective_dp_q20
###########################################################

cat > ${OUTDIR}/scripts/direct_spike_from_effective_cache.py << 'PY'
import sys, random, pysam, pickle
from collections import defaultdict

bam_in, sites_file, bam_out, truth_out, seed, cache_file, min_bq = sys.argv[1:8]
seed = int(seed)
min_bq = int(min_bq)
random.seed(seed)

with open(cache_file, "rb") as f:
    cache = pickle.load(f)

sites = []
with open(sites_file) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        fs = line.strip().split()
        chrom = fs[0]
        pos = int(fs[1])
        vaf = float(fs[3])
        ref = fs[4].upper()
        alt = fs[5].upper()
        dp = int(fs[6])
        ref_count = int(fs[7])
        n_candidate = int(fs[8])
        need_alt = int(fs[9])
        dp_bin = fs[10]
        sites.append((chrom, pos, vaf, ref, alt, dp, ref_count, n_candidate, need_alt, dp_bin))

fragment_edit_map = {}
truth = []

for chrom, pos, vaf, ref, alt, dp, ref_count, n_candidate, need_alt, dp_bin in sites:
    info = cache.get((chrom, pos))
    if info is None:
        continue

    candidates = list(info["candidate_fragments"])
    if len(candidates) < need_alt:
        continue

    chosen = random.sample(candidates, need_alt)

    for qname in chosen:
        fragment_edit_map[(chrom, pos, qname)] = alt

    truth.append((chrom, pos, ref, alt, dp, ref_count, len(candidates), need_alt, dp_bin))

edit_by_qname = defaultdict(lambda: defaultdict(dict))
for chrom, pos, qname in fragment_edit_map:
    edit_by_qname[qname][chrom][pos] = fragment_edit_map[(chrom, pos, qname)]

bam = pysam.AlignmentFile(bam_in, "rb")
out = pysam.AlignmentFile(bam_out, "wb", template=bam)

site_edited_bases = defaultdict(int)
edited_reads = 0
edited_bases = 0

for r in bam:
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.query_sequence is None:
        out.write(r)
        continue

    qname = r.query_name
    if qname not in edit_by_qname:
        out.write(r)
        continue

    chrom = r.reference_name
    if chrom is None or chrom not in edit_by_qname[qname]:
        out.write(r)
        continue

    if r.query_qualities is None:
        out.write(r)
        continue

    target = edit_by_qname[qname][chrom]
    seq = list(r.query_sequence)
    quals = r.query_qualities
    changed = 0

    for qp, rp in r.get_aligned_pairs(matches_only=False):
        if qp is None or rp is None:
            continue

        pos = rp + 1
        if pos not in target:
            continue


        if quals[qp] < min_bq:
            continue

        alt = target[pos]
        if qp < len(seq) and seq[qp].upper() in "ACGT" and seq[qp].upper() != alt:
            seq[qp] = alt
            changed += 1
            site_edited_bases[(chrom, pos)] += 1

    if changed > 0:
        r.query_sequence = "".join(seq)
        r.query_qualities = quals
        edited_reads += 1
        edited_bases += changed

    out.write(r)

out.close()
bam.close()

with open(truth_out, "w") as w:
    w.write("chrom\tpos\tref\talt\teffective_dp_q20\tref_count_q20\tcandidate_fragments_q20\tmut_fragments_planned\tmut_bases_edited\treal_vaf\tdp_bin\n")
    for chrom, pos, ref, alt, dp, ref_count, n_candidate, need_alt, dp_bin in truth:
        mut_bases = site_edited_bases.get((chrom, pos), 0)
        real_vaf = mut_bases / dp if dp > 0 else 0
        w.write("\t".join(map(str, [
            chrom, pos, ref, alt,
            dp, ref_count, n_candidate,
            need_alt, mut_bases, real_vaf, dp_bin
        ])) + "\n")

print(f"[INFO] truth_sites={len(truth)}", file=sys.stderr)
print(f"[INFO] edited_reads={edited_reads}", file=sys.stderr)
print(f"[INFO] edited_bases={edited_bases}", file=sys.stderr)
PY

###########################################################
# Manifest
###########################################################

manifest=${OUTDIR}/summary/prepare_manifest.tsv

if [[ ! -s ${manifest} ]]; then

echo -e \
"Sample\tMode\tVAF\tDP_bin\tCandidateSites\tSubsetBam\tEffectiveCache\tSelectedSites\tSpikeinBam\tCalmdBam\tTruth" \
> ${manifest}

fi

###########################################################
# Main loop
###########################################################

for sample in "${samples[@]}"; do
    echo "========== ${sample} =========="

    bam=$(get_bam ${sample})
    region_bed=${REGION_DIR}/${sample}.final_sim_region_ge1.bed
    candidate_sites=${OUTDIR}/sites/${sample}.candidate_sites.txt
    sample_subset=${OUTDIR}/bam/${sample}.candidate.original_subset.bam
    eff_cache=${OUTDIR}/cache/${sample}.site_info.effective_q${MIN_BQ}.pkl
	candidate_window=${OUTDIR}/window/${sample}.candidate.window.bed


    [[ -s ${bam} ]] || {
        echo "[SKIP] BAM missing: ${bam}"
        continue
    }


    [[ -s ${region_bed} ]] || {
        echo "[SKIP] region missing: ${region_bed}"
        continue
    }
#######################################################
# Step 0
# final region -> 15000 candidate sites
#######################################################

if [[ -s ${candidate_sites} ]]; then
    echo "[0] candidate exists"

else

    echo "[0] generate candidate sites"
    python ${OUTDIR}/scripts/make_candidate_sites.py \
        ${region_bed} \
        ${REF} \
        ${N_CANDIDATE} \
        ${candidate_sites} \
        ${SEED}

fi

#######################################################
# Step 0.5
# subset BAM
#######################################################

if [[ -s ${candidate_window} ]]; then
    echo "[2] candidate window exists, skip: ${candidate_window}"
else
    echo "[2] make candidate window"
    awk -v OFS="\t" -v flank=${FLANK} '{s=$2-1-flank; e=$3+flank; if(s<0)s=0; print $1,s,e}' ${candidate_sites} \
        | sort -k1,1 -k2,2n \
        | bedtools merge \
        > ${candidate_window}
fi

if [[ -s ${sample_subset} && -s ${sample_subset}.bai ]]; then
    echo "[3] candidate subset BAM exists, skip: ${sample_subset}"
else
    echo "[3] subset BAM once"
    samtools view -@ ${THREADS} -P -b -L ${candidate_window} ${bam} > ${sample_subset}
    samtools index ${sample_subset}
fi

    #######################################################
    # Step 1
    #######################################################

    if [[ -s ${eff_cache} ]]; then
        echo "[1] effective cache exists, skip: ${eff_cache}"
    else
        echo "[1] build effective cache from 50000 sites"
    python ${OUTDIR}/scripts/build_effective_cache.py \
        ${sample_subset} \
        ${candidate_sites} \
        ${eff_cache} \
        ${MIN_BQ} \
        ${MIN_MQ} \
        > ${OUTDIR}/log/${sample}.build_effective_cache.log 2>&1
        fi

    #######################################################
    # Step 2-5


    #######################################################

    for bin in "${DP_BINS[@]}"; do
        tag=${sample}.vaf${FIXED_VAF}.dpbin${bin}.q${MIN_BQ}

		selected_sites=${OUTDIR}/sites/${tag}.selected_sites.txt
        spikein=${OUTDIR}/bam/${tag}.spikein.bam
        spikein_sorted=${OUTDIR}/bam/${tag}.spikein.sorted.bam
        spikein_calmd=${OUTDIR}/bam/${tag}.spikein.sorted.calmd.bam
        truth=${OUTDIR}/truth/${tag}.truth.tsv

        echo "[2] ${sample} DPbin=${bin}"

        if [[ -s ${selected_sites} ]]; then
            echo "[2] selected sites exists, skip: ${selected_sites}"
        else
            echo "[2] select sites by effective DP bin"
            python ${OUTDIR}/scripts/select_sites_by_effbin.py \
                ${eff_cache} \
                ${bin} \
                ${N_SITE} \
                ${FIXED_VAF} \
                ${MINMUTREADS} \
                ${selected_sites} \
                $((SEED + ${#bin})) \
                > ${OUTDIR}/log/${tag}.select_sites.log 2>&1
        fi

        n_selected=$(wc -l < ${selected_sites})
        echo "[INFO] selected sites: ${n_selected}"

        [[ ${n_selected} -gt 0 ]] || { echo "[SKIP] no selected sites: ${tag}"; continue; }

        if [[ -s ${spikein_calmd} && -s ${spikein_calmd}.bai && -s ${truth} ]]; then
            echo "[3] calmd spike-in BAM exists, skip: ${spikein_calmd}"
        else
            echo "[3] direct spike-in"
            python ${OUTDIR}/scripts/direct_spike_from_effective_cache.py \
                ${sample_subset} \
                ${selected_sites} \
                ${spikein} \
                ${truth} \
                $((SEED + ${#bin})) \
                ${eff_cache} \
                ${MIN_BQ} \
                > ${OUTDIR}/log/${tag}.spikein.log 2>&1

            [[ -s ${spikein} ]] || { echo "[ERROR] spike-in failed: ${tag}"; continue; }

            echo "[4] sort/index"
            samtools sort -@ ${THREADS} -o ${spikein_sorted} ${spikein}
            samtools index ${spikein_sorted}
            samtools quickcheck ${spikein_sorted}

            echo "[5] calmd/index"
            samtools calmd -@ ${THREADS} -b ${spikein_sorted} ${REF} > ${spikein_calmd}
            samtools index ${spikein_calmd}
            samtools quickcheck ${spikein_calmd}
        fi

        n_truth=$(($(wc -l < ${truth}) - 1))
        echo "[INFO] truth sites: ${n_truth}"

        if ! awk \
			-v s="${sample}" \
			-v v="${FIXED_VAF}" \
			-v b="${bin}" \
			'$1==s && $2=="DPBIN_Q20" && $3==v && $4==b{found=1} END{exit found?0:1}' \
			${manifest}
		then
			echo -e "${sample}\tDPBIN_Q20\t${FIXED_VAF}\t${bin}\t${candidate_sites}\t${subset_bam}\t${eff_cache}\t${selected_sites}\t${spikein_sorted}\t${spikein_calmd}\t${truth}" \
			>> ${manifest}

fi
    done

    echo "[DONE] ${sample}"
done
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
