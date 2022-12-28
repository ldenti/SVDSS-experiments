# SVDSS experiments

## Prerequisites
### Tools
```bash
conda install minimap2 pbsv=2.6.2 sniffles=1.0.12 cutesv=1.0.11 svim=1.4.2 dipcall samtools bcftools pysam biopython numpy pandas seaborn
pip3 install truvari

# We need ngmlr v0.2.8 (master) due to a bug in previous releases (not yet in conda at the time of the experiments)
git clone https://github.com/philres/ngmlr.git
cd ngmlr
git checkout a2a31fb6a63547be29c5868a7747e0c6f6e9e41f
mkdir build ; cd build
cmake ..
make

# pbmm2 and debreak will be loaded by snakemake from envs/ folder due to conflicts

# clone and install SVDSS from https://github.com/Parsoa/SVDSS
```

### Data
In our experimental evaluation we used hg38.

#### Reference and annotations
```bash
# Reference
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > hg38.chroms.fa
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' hg38.chroms.fa
samtools faidx hg38.chroms.fa
cut -f1,2 hg38.chroms.fa.fai | sort -k 1 > hg38.chroms.fa.genome

# Dipcall PAR
wget https://raw.githubusercontent.com/lh3/dipcall/master/data/hs38.PAR.bed

# pbsv TRF regions
wget https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
```

#### Tiers
In our paper we used the Tier 1 from GIAB and then we created an Extended Tier 2, that is everything outside Tier 1. Tier 1 is downloaded from GIAB FTP server and lifted to hg38. Extended Tier 2 (here called Tier23) is computed by complementing Tier 1.
```bash
# Get Tier 1 and lift to hg38
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
sed -e 's/^/chr/' HG002_SVs_Tier1_v0.6.bed > Tier1_v0.6.hg19.bed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz
./liftOver Tier1_v0.6.hg19.bed hg19ToHg38.over.chain Tier1_v0.6.hg38.bed Tier1_v0.6.unlifted.bed

# Clean and sort bed to be able to complement it
for c in $(seq 1 22) X Y ; do grep -P "^chr${c}\t" Tier1_v0.6.hg38.bed ; done | sort -k1,1 -k2,2n > Tier1_v0.6.hg38.noalt.bed
bedtools complement -i Tier1_v0.6.hg38.noalt.bed -g hg38.chroms.genome > Tier23_v0.6.hg38.noalt.bed
```

#### Truth VCFs
For users convenience, we provide the three VCFs computed with dipcall against hg38 and used as groundtruth in our analysis (see `truthsets` folder).

#### config.yaml
The experiments are provided as [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows. The file `config.yaml` is used to tweak snakemake execution:
* **name**: custom name for the run
* **fa**: reference genome (i.e., _hg38.chroms.fa_)
* **fq**: HiFi sample
* **trf**: trf annotation
* **par**: PAR regions
* **hap1**: first haplotype
* **hap2**: second haplotype
* **tier1**: GIAB Tier 1
* **tier23**: Extended Tier 2
* **out**: output directory (everything will go here)
* **pingpong_bin**: path to SVDSS binary
* **ngmlr_bin**: path to ngmlr (v0.2.8) binary
* **threads**: number of threads for each rule that requires more threads
* **aligners**: list of aligners to use
* **callers**: list of callers to use

## HG007
Download [corrected HiFi reads](https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/deepconsensus_predictions/hg007_15kb/three_smrt_cells/HG007_230654_115437_2fl_DC_hifi_reads.fastq) and corresponding haplotypes:
* [HG007.asm.hap1.p_ctg.fa](https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/analysis/genome_assembly/hg007_15kb/two_smrt_cells/dc/HG007.asm.hap1.p_ctg.fa)
* [HG007.asm.hap2.p_ctg.fa](https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/analysis/genome_assembly/hg007_15kb/two_smrt_cells/dc/HG007.asm.hap2.p_ctg.fa)

Change `config.yaml` accordingly and run:
```bash
# Mappers + callers
snakemake --use-conda [-n] -p -j 16
# Dipcall and truvari
snakemake -s Snakefile.dipvari [-n] -p -j 8
```

#### Coverage titration for 5x and 10x
```bash
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
```bash
python3 plot_coverage.py ${OUT-5x}/pbmm2/results.csv ${OUT-10x}/pbmm2/results.csv ${OUT}/pbmm2/results.csv
```

#### Other plots
Assuming `${OUT}` to be the output folder set in `config.yaml`:
```bash
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
To compare the callsets to the CMRG callset:
```bash
# Download the CMRG callset
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi

# Run truvari against the CMRG callset (assuming ${HG2-OUT} to be the output folder for the HG002 analysis)
for t in pp cutesv pbsv svim sniffles debreak ; do truvari bench -b HG002_GRCh38_CMRG_SV_v1.00.vcf.gz -c ${HG2-OUT}/pbmm2/${t}.vcf.gz -o ${t} -f ~/data/hg38-refanno/hg38.chroms.fa -r 1000 -p 0.00 -s 20 -S 20 ; done

# Plot the CMRG venn and supervenn
python3 medical_analysis.py . # . is the folder with the folders created by truvari in the previous cycle
```

#### GIAB vs dipcall analysis
To compare the tools against the GIAB callset, download the hg19 data and run the pipeline.

Download (hg19) reference and annotations:
```bash
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
# Extract chromosomes
samtools faidx hs37d5.fa 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y > hs37d5.chroms.fa
# Map non-ACGT characters to N:
sed -i '/^[^>]/ y/BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrsuvwxyz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/' hs37d5.chroms.fa

# Get hg19 PAR region from dipcall repository:
wget https://raw.githubusercontent.com/lh3/dipcall/master/data/hs37d5.PAR.bed
```
Download GIAB Tier 1 and complement it:
```bash
# Get GIAB VCF and tiers:
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed
# Here we just complement the Tier 1 to get the Extended Tier 2
bedtools complement -i HG002_SVs_Tier1_v0.6.bed -g hs37d5.chroms.fa.fai > HG002_SVs_Tier23_v0.6.bed
```

Then, update `config.yaml` accordingly (see `config-hg19.yaml`, i.e., everything should point to hg19 and a new output directory) and:
```bash
# run everything on the older reference release
snakemake --use-conda -p -j 16
# compare the new results with dipcall (against hg19)
snakemake -s Snakefile.dipvari -p -j 8
# compare the new results with the GIAB callset, GIAB vs dipcall, and dipcall vs GIAB
snakemake -s Snakefile.giabvari -p -j 8

# Heterozygosity plot
python3 plot_hetebars.py {dipcall.clean.vcf.gz} {GIAB.clean.vcf.gz} ${OUT}/pbmm2/
```

## CHM13
Download [HiFi reads](https://github.com/marbl/CHM13#hifi-data) and [corresponding haplotypes](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz).

Update `config.yaml` accordingly and run `snakemake`.
