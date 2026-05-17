#!/usr/bin/env bash
#
# find_orthologs.sh
# -----------------
# Reciprocal Best Hits (RBH) ortholog finder for two species, given CDS FASTAs.
#
# Pipeline:
#   1. Translate CDS -> protein (frame 1, stop codons trimmed)
#   2. Build BLAST protein databases for each species
#   3. blastp A -> B and blastp B -> A
#   4. Keep the top hit per query (by bitscore)
#   5. Intersect to produce reciprocal best hits = putative orthologs
#
# Dependencies (must be on PATH):
#   - NCBI BLAST+ (makeblastdb, blastp)
#   - EMBOSS transeq             (alternative: seqkit translate)
#
# Usage:
#   ./find_orthologs.sh -a speciesA.cds.fa -b speciesB.cds.fa [-o outdir] [-e 1e-5] [-t 4]

set -euo pipefail

# ---------- defaults ----------
OUTDIR="ortholog_results"
EVALUE="1e-5"
THREADS=4

usage() {
    cat <<EOF
Usage: $0 -a <speciesA.cds.fa> -b <speciesB.cds.fa> [options]

Required:
  -a FILE   CDS FASTA for species A (nucleotide)
  -b FILE   CDS FASTA for species B (nucleotide)

Options:
  -o DIR    Output directory          (default: $OUTDIR)
  -e NUM    E-value threshold         (default: $EVALUE)
  -t INT    BLAST threads             (default: $THREADS)
  -h        Show this help
EOF
    exit 1
}

# ---------- parse args ----------
while getopts "a:b:o:e:t:h" opt; do
    case $opt in
        a) CDS_A="$OPTARG" ;;
        b) CDS_B="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "${CDS_A:-}" || -z "${CDS_B:-}" ]] && usage
[[ -f "$CDS_A" ]] || { echo "ERROR: $CDS_A not found"; exit 1; }
[[ -f "$CDS_B" ]] || { echo "ERROR: $CDS_B not found"; exit 1; }

# ---------- check deps ----------
for tool in makeblastdb blastp transeq; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not found on PATH"; exit 1; }
done

# Use absolute paths before changing directory
CDS_A=$(readlink -f "$CDS_A")
CDS_B=$(readlink -f "$CDS_B")

mkdir -p "$OUTDIR"
cd "$OUTDIR"

NAME_A=$(basename "$CDS_A" | sed -E 's/\.(fa|fasta|fna|cds)(\.gz)?$//I')
NAME_B=$(basename "$CDS_B" | sed -E 's/\.(fa|fasta|fna|cds)(\.gz)?$//I')

echo ">>> [1/5] Translating CDS to protein..."
transeq -sequence "$CDS_A" -outseq "${NAME_A}.pep.fa" -frame 1 -trim -clean >/dev/null 2>&1
transeq -sequence "$CDS_B" -outseq "${NAME_B}.pep.fa" -frame 1 -trim -clean >/dev/null 2>&1
# transeq appends "_1" to every header; strip it so IDs match the original CDS
sed -i -E 's/^(>[^ ]+)_1\b/\1/' "${NAME_A}.pep.fa"
sed -i -E 's/^(>[^ ]+)_1\b/\1/' "${NAME_B}.pep.fa"

echo ">>> [2/5] Building BLAST databases..."
makeblastdb -in "${NAME_A}.pep.fa" -dbtype prot -out "db_${NAME_A}" >/dev/null
makeblastdb -in "${NAME_B}.pep.fa" -dbtype prot -out "db_${NAME_B}" >/dev/null

FMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

echo ">>> [3/5] blastp ${NAME_A} -> ${NAME_B} ..."
blastp -query "${NAME_A}.pep.fa" -db "db_${NAME_B}" \
       -outfmt "$FMT" -evalue "$EVALUE" \
       -num_threads "$THREADS" \
       -out "${NAME_A}_vs_${NAME_B}.tsv"

echo ">>> [4/5] blastp ${NAME_B} -> ${NAME_A} ..."
blastp -query "${NAME_B}.pep.fa" -db "db_${NAME_A}" \
       -outfmt "$FMT" -evalue "$EVALUE" \
       -num_threads "$THREADS" \
       -out "${NAME_B}_vs_${NAME_A}.tsv"

echo ">>> [5/5] Finding reciprocal best hits..."
# Best hit per query = highest bitscore (column 12); break ties by lowest evalue (col 11).
# Sort by query asc, bitscore desc, evalue asc, then keep first row per query.
sort -k1,1 -k12,12gr -k11,11g "${NAME_A}_vs_${NAME_B}.tsv" \
    | awk -F'\t' '!seen[$1]++ {print $1"\t"$2"\t"$3"\t"$11"\t"$12}' \
    > "best_${NAME_A}_to_${NAME_B}.tsv"

sort -k1,1 -k12,12gr -k11,11g "${NAME_B}_vs_${NAME_A}.tsv" \
    | awk -F'\t' '!seen[$1]++ {print $1"\t"$2"\t"$3"\t"$11"\t"$12}' \
    > "best_${NAME_B}_to_${NAME_A}.tsv"

# Reciprocal: A's best hit in B is X, and X's best hit in A is back to A.
awk -F'\t' 'BEGIN{
    OFS="\t"
    print "queryA","bestB","pident_A2B","evalue_A2B","bitscore_A2B","pident_B2A","evalue_B2A","bitscore_B2A"
}
NR==FNR { best_b2a[$1]=$2; pid[$1]=$3; ev[$1]=$4; bs[$1]=$5; next }
($2 in best_b2a) && (best_b2a[$2] == $1) {
    print $1, $2, $3, $4, $5, pid[$2], ev[$2], bs[$2]
}' "best_${NAME_B}_to_${NAME_A}.tsv" "best_${NAME_A}_to_${NAME_B}.tsv" \
    > "RBH_orthologs.tsv"

N=$(( $(wc -l < RBH_orthologs.tsv) - 1 ))
echo ""
echo ">>> Done."
echo ">>> Reciprocal best hit ortholog pairs: $N"
echo ">>> Output: $(pwd)/RBH_orthologs.tsv"
