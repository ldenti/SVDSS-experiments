# PingPong - SVs evaluation

## Prerequisites
Install tools:
```
conda install minimap2 pbsv sniffles cutesv svim svim-asm dipcall samtools bcftools pysam biopython
pip3 install truvari

git clone https://github.com/philres/ngmlr.git
cd ngmlr
mkdir build ; cd build
cmake ..
make

# pbmm2 will be loaded by snakemake from envs/pbmm2.yaml due to conflicts
```

Dowwnload (hg19) reference and annotations:
```
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
# Extract chromosomes
samtools faidx hs37d5.fa 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y > hs37d5.chroms.fa
# Map non-ACGT characters to N:
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' hs37d5.chroms.fa

# Get hg19 PAR region from dipcall repository:
wget https://raw.githubusercontent.com/lh3/dipcall/master/data/hs37d5.PAR.bed

# Get trf from pbsv repository:
wget https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_hs37d5.trf.bed

# Get GIAB VCF and tiers:
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier2_v0.6.bed
```

### config.yaml
`config.yaml` is used to tweak snakemake execution:
* `name`: custom name for the run
* `fa`: reference genome (i.e., _hs37d5.chroms.fa_)
* `fq`: HiFi sample
* `trf`: trf annotation
* `par`: PAR regions
* `hap1`: first haplotype
* `hap2`: second haplotype
* `tier1`: GIAB tier 1
* `tier2`: GIAB tier 2
* `out`: output directory (everything will go here)
* `threads`: number of threads for each rule that requires more threads

## HG007
Download [corrected HiFi reads](https://console.cloud.google.com/storage/browser/_details/brain-genomics-public/research/deepconsensus/publication/deepconsensus_predictions/hg007_15kb/three_smrt_cells/HG007_230654_115437_2fl_DC_hifi_reads.fastq) and haplotypes built from those:
* [HG007.asm.hap1.p_ctg.fa](https://console.cloud.google.com/storage/browser/brain-genomics-public/research/deepconsensus/publication/analysis/genome_assembly/hg007_15kb/two_smrt_cells/dc/HG007.asm.hap1.p_ctg.fa)
* [HG007.asm.hap2.p_ctg.fa](https://console.cloud.google.com/storage/browser/brain-genomics-public/research/deepconsensus/publication/analysis/genome_assembly/hg007_15kb/two_smrt_cells/dc/HG007.asm.hap2.p_ctg.fa)

Change `config.yaml` accordingly and run:
```
# Run mappers and callers
snakemake --use-conda [-n] -p -j 16
# Run dipcall and truvari
snakemake -s Snakefile.analysis [-n] -p -j 8
```

## HG002


## CHM13




