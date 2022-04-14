# SVDSS experiments

## Prerequisites
Install tools:
```
conda install minimap2 pbsv=2.6.2 sniffles=1.0.12 cutesv=1.0.11 svim=1.4.2 dipcall samtools bcftools pysam biopython numpy pandas seaborn
pip3 install truvari

# We need ngmlr v0.2.8 (master) due to a bug in previous releases (not yet in conda at the time of the experiments)
git clone https://github.com/philres/ngmlr.git
cd ngmlr
mkdir build ; cd build
cmake ..
make

# pbmm2 and debreak will be loaded by snakemake from envs/ folder due to conflicts

# clone and install SVDSS from https://github.com/Parsoa/SVDSS
```

Download (hg38) reference and annotations:
```
# Reference
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > hg38.chroms.fa
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' hg38.chroms.fa

# GIAB tiers
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
sed -e 's/^/chr/' HG002_SVs_Tier1_v0.6.bed > Tier1_v0.6.hg19.bed
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
liftOver Tier1_v0.6.hg19.bed hg19ToHg38.over.chain Tier1_v0.6.bed Tier1_v0.6.unlifted.bed
for c in $(seq 1 22) X Y ; do grep -P "^chr${c}\t" Tier1_v0.6.hg38.bed ; done | sort -k1,1 -k2,2n > Tier1_v0.6.hg38.noalt.bed
bedtools complement -i Tier1_v0.6.hg38.noalt.bed -g hg38.chroms.genome > Tier23_v0.6.hg38.noalt.bed

# Dipcall PAR
wget https://raw.githubusercontent.com/lh3/dipcall/master/data/hs38.PAR.bed

# pbsv TRF regions
wget https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
```

### config.yaml
The file `config.yaml` is used to tweak snakemake execution:
* **name**: custom name for the run
* **fa**: reference genome (i.e., _hs37d5.chroms.fa_)
* **fq**: HiFi sample
* **trf**: trf annotation
* **par**: PAR regions
* **hap1**: first haplotype
* **hap2**: second haplotype
* **tier1**: GIAB tier 1
* **out**: output directory (everything will go here)
* **pingpong_bin**: path to PingPong binary
* **ngmlr_bin**: path to ngmlr (v0.2.8) binary
* **threads**: number of threads for each rule that requires more threads
* **aligners**: list of aligners to use
* **callers**: list of callers to use

### Truth VCFs
For users convenience, we provide the three VCFs computed with dipcall against hg38 and used as groundtruth in our analysis (see `truthsets` folder).

## HG007
Download [corrected HiFi reads](https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/deepconsensus_predictions/hg002_15kb/two_smrt_cells/HG002_15kb_222723_002822_2fl_DC_hifi_reads.fastq) and corresponding haplotypes:
* [HG007.asm.hap1.p_ctg.fa](https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/analysis/genome_assembly/hg007_15kb/two_smrt_cells/dc/HG007.asm.hap1.p_ctg.fa)
* [HG007.asm.hap2.p_ctg.fa](https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/analysis/genome_assembly/hg007_15kb/two_smrt_cells/dc/HG007.asm.hap2.p_ctg.fa)

Change `config.yaml` accordingly and run:
```
# Mappers + callers
snakemake --use-conda [-n] -p -j 16
# Dipcall and truvari
snakemake -s Snakefile.dipvari [-n] -p -j 8
```

#### Coverage titration for 5x and 10x
```
samtools view -b -s 0.3 /path/to/out/dir/pbmm2.bam > pbmm2.5x.bam
samtools fastq pbmm2.5x.bam > pbmm2.5x.fq
# change config.yaml
snakemake --use-conda -p -j 16
snakemake -s Snakefile.dipvari -p -j 8


samtools view -b -s 0.6 /path/to/out/dir/pbmm2.bam > pbmm2.10x.bam
samtools fastq pbmm2.10x.bam > pbmm2.10x.fq
# change config.yaml
snakemake --use-conda -p -j 16
snakemake -s Snakefile.dipvari -p -j 8
```

To plot results, assuming `${OUT-*}` to be the output folders set in `config.yaml`, run:
```
python3 plot_coverage.py ${OUT-5x}/pbmm2/results.csv ${OUT-10x}/pbmm2/results.csv ${OUT}/pbmm2/results.csv
```

#### Other plots

Assuming `${OUT}` to be the output folder set in `config.yaml`:
```
# VCF length distribution
python3 plot_vcflen.py ${OUT}

# Upset plot for Ext. Tier 2 region
python3 plot_upset.py ${OUT}/pbmm2/truvari/ tier23conf > hg007.tier23conf.ppexclusive.vcf

# Exclusive calls from SVDSS heterozygosity
python3 check_exclusive.py get hg007.tier23conf.ppexclusive.vcf ${OUT}/pbmm2/truvari/dipcall.tier23conf.vcf > hg007.tier23conf.ppexclusive.dist
python3 check_exclusive.py plot hg007.tier23conf.ppexclusive.dist
```

## HG002
Download [corrected HiFi reads](https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/deepconsensus_predictions/hg002_15kb/two_smrt_cells/HG002_15kb_222723_002822_2fl_DC_hifi_reads.fastq) and corresponding haplotypes:
* [HG002.asm.hap1.p_ctg.fa](https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/analysis/genome_assembly/hg002_15kb/two_smrt_cells/dc/HG002.asm.hap1.p_ctg.fa)
* [HG002.asm.hap2.p_ctg.fa](https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/analysis/genome_assembly/hg002_15kb/two_smrt_cells/dc/HG002.asm.hap2.p_ctg.fa)

Update `config.yaml` accordingly and run `snakemake`.

#### CMRG analysis
```
# Download the CMRG callset
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi

# Run truvari against the CMRG callset (assuming ${HG2-OUT} to be the output folder for the HG002 analysis)
for t in pp cutesv pbsv svim sniffles debreak ; do truvari bench -b HG002_GRCh38_CMRG_SV_v1.00.vcf.gz -c ${HG2-OUT}/pbmm2/${t}.vcf.gz -o ${t} -f ~/data/hg38-refanno/hg38.chroms.fa -r 1000 -p 0.00 -s 20 -S 20 ; done

# CMRG venn and supervenn
python3 medical_analysis.py .
```

## CHM13
Download [HiFi reads](https://github.com/marbl/CHM13#hifi-data) and [corresponding haplotypes](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz).

Update `config.yaml` accordingly and run `snakemake`.
