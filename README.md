# Finding Orthologs Between Two Species with Reciprocal Best Hits

A walkthrough for building a simple, reproducible ortholog-finding pipeline in bash. Given the coding sequences (CDS) of two species, this tutorial shows you how to identify pairs of genes that are likely **orthologs** — genes in different species that descend from a single gene in their last common ancestor.

The full script lives in [`find_orthologs.sh`](./find_orthologs.sh). This README explains what it does, why each step matters, and how to use and extend it.

---

## Why orthologs?

When you compare genes across species, the relationships fall into two main categories:

- **Orthologs** are genes that diverged because the *species* diverged. They usually keep the same function — the human and mouse versions of a gene are typically orthologs.
- **Paralogs** are genes that diverged because the *gene* duplicated within a genome. They often drift toward new functions.

If you want to answer questions like *"what is the mouse equivalent of this human gene?"* or *"is this enzyme conserved in my non-model organism?"*, you are asking about orthology. Getting the right pairs matters: a sloppy match between paralogs can lead you to compare genes that no longer do the same thing.

## The Reciprocal Best Hits (RBH) method

RBH is the workhorse method for one-to-one ortholog detection. The idea is simple:

> If gene **A₁** in species A finds gene **B₅** as its best match in species B, *and* gene **B₅** in species B finds **A₁** as its best match in species A, then A₁ and B₅ are reciprocal best hits — putative orthologs.

```
   Species A                  Species B
   ─────────                  ─────────
      A1  ──── best hit ────►  B5
                                │
      A1  ◄─── best hit ────────┘

      → A1 and B5 are reciprocal best hits ✓
```

It's fast, it's easy to reason about, and it gives high-precision pairs. It misses some orthologs (especially after lineage-specific duplications), but for a first pass between two genomes it is hard to beat.

## What you'll need

Two pieces of software, both standard in bioinformatics:

| Tool | Purpose | Install |
|---|---|---|
| **NCBI BLAST+** | sequence search | `sudo apt install ncbi-blast+` *(Debian/Ubuntu)* or `brew install blast` *(macOS)* |
| **EMBOSS** | provides `transeq` for translation | `sudo apt install emboss` *(Debian/Ubuntu)* or `brew install brewsci/bio/emboss` *(macOS)* |

You will also need two **CDS FASTA files**, one per species — nucleotide sequences of protein-coding genes, with one entry per gene. You can usually download these from Ensembl, NCBI RefSeq, or the genome project for your organism. Headers should be unique, e.g.:

```
>ENSG00000139618
ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA...
>ENSG00000141510
ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCA...
```

## Quick start

```bash
git clone https://github.com/<you>/<repo>.git
cd <repo>
chmod +x find_orthologs.sh

./find_orthologs.sh \
    -a human.cds.fa \
    -b mouse.cds.fa \
    -o human_mouse_orthologs \
    -t 8
```

That's it. Twenty to thirty minutes later (depending on genome size and thread count), you'll have a tab-separated file of ortholog pairs in `human_mouse_orthologs/RBH_orthologs.tsv`.

## How it works, step by step

The script runs five stages. Knowing what each does will help you debug and customize.

### 1. Translate CDS to protein

```bash
transeq -sequence speciesA.cds.fa -outseq speciesA.pep.fa -frame 1 -trim
```

We translate nucleotide CDS into amino acid sequences using EMBOSS `transeq`. Why translate? Because **protein-level comparison is dramatically more sensitive than nucleotide comparison** when species are even moderately distant. The genetic code is redundant — many codons map to the same amino acid — so two true orthologs can diverge at the DNA level while their proteins remain nearly identical. blastp sees through that noise.

`transeq` appends `_1` to every header (its frame indicator); the script strips that so the IDs stay consistent with the original CDS.

### 2. Build BLAST databases

```bash
makeblastdb -in speciesA.pep.fa -dbtype prot -out db_speciesA
```

BLAST needs an indexed database to search against. We build one for each species, in protein mode (`-dbtype prot`).

### 3. BLAST each species against the other

```bash
blastp -query speciesA.pep.fa -db db_speciesB ...
blastp -query speciesB.pep.fa -db db_speciesA ...
```

Two `blastp` runs, in opposite directions. Output is BLAST's tabular format 6, which includes bitscore — the metric we'll use to pick best hits.

A subtle point worth knowing: BLAST has a `-max_target_seqs 1` option that *looks* like it should give you the top hit, but it doesn't reliably do so. Internally, BLAST applies that limit *during* the search heuristic rather than after final ranking, and the kept hit may not be the best one. So we ask BLAST for all hits and pick the best one ourselves.

### 4. Pick the best hit per query

```bash
sort -k1,1 -k12,12gr -k11,11g blast.tsv | awk '!seen[$1]++' > best_hits.tsv
```

Sort by query ID, then by bitscore descending (`-k12,12gr`), then by e-value ascending as a tiebreaker. The `awk '!seen[$1]++'` idiom keeps only the first row for each query — which, after that sort, is the best hit.

### 5. Intersect for reciprocal best hits

```awk
NR==FNR { best_b2a[$1]=$2; next }
($2 in best_b2a) && (best_b2a[$2] == $1) { print }
```

Read all of B→A's best hits into a hash, then walk through A→B's best hits and keep only the rows whose target points back to the original query. The rows that survive are your reciprocal best hits.

## Understanding the output

`RBH_orthologs.tsv` has one row per ortholog pair:

| Column | Meaning |
|---|---|
| `queryA` | Gene ID in species A |
| `bestB` | Its ortholog in species B |
| `pident_A2B` | % identical residues in the A→B alignment |
| `evalue_A2B` | E-value of the A→B hit |
| `bitscore_A2B` | Bitscore of the A→B hit |
| `pident_B2A` etc. | Same three statistics for the reverse direction |

Higher bitscore and lower e-value mean stronger evidence. For most downstream uses, you can filter on **percent identity** (e.g. keep pairs with `pident ≥ 30%`) and **alignment coverage**. To add a coverage filter, edit the BLAST output format string in the script to include `qlen` and `slen`, then compute `length / min(qlen, slen)` in the final awk step.

## Tuning the run

The script exposes the three knobs you're most likely to touch:

```bash
./find_orthologs.sh -a A.fa -b B.fa \
    -e 1e-10 \     # stricter e-value
    -t 16 \        # more threads
    -o my_outdir   # custom output directory
```

For closely related species (say, human/mouse) an e-value of `1e-10` or even `1e-20` is fine. For more distant pairs (vertebrate/insect, plant/animal) loosen it to `1e-5` or `1e-3` so you don't miss real but weakly-conserved orthologs.

## Limitations and when to use something else

RBH is a great first pass, but it has known blind spots:

**It assumes one-to-one orthology.** After a gene duplication in one lineage, two paralogs in species A may both descend from a single gene in species B. RBH will pick only one of them and ignore the other real ortholog. If your biology involves whole-genome duplications (teleost fish, many plants), you will miss real signal.

**It is sensitive to incomplete annotations.** If a gene is missing from one CDS file, its true ortholog in the other species will be paired with whatever is "next best" — or end up with no pair at all.

**It does not build gene trees.** For phylogenetically rigorous orthology, dedicated tools like [OrthoFinder](https://github.com/davidemms/OrthoFinder), [ProteinOrtho](https://gitlab.com/paulklemm_PHD/proteinortho), or [OMA](https://omabrowser.org/) construct gene families and resolve one-to-many and many-to-many relationships. They are slower but more complete.

A reasonable workflow: use this RBH script for a fast first pass and sanity check, then move to a tree-based tool if you need to handle gene families properly.

## Citing

If you use this in a publication, please cite the underlying tools:

- **BLAST+**: Camacho et al. (2009) *BMC Bioinformatics* 10:421
- **EMBOSS**: Rice, Longden & Bleasby (2000) *Trends in Genetics* 16(6):276–277
