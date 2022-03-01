configfile: "config.yaml"

from os.path import join as pjoin

REF = config["fa"]
SAMPLE = config["fq"]
TRF = config["trf"]

RUN = config["name"]

ODIR = config["out"] # os.getcwd()
THREADS = config["threads"]

PP_BIN = config["pingpong_bin"]
NGMLR_BIN = config["ngmlr_bin"]

aligners = []
for aligner in config["aligners"]:
    aligners.append(aligner)
callers = []
for caller in config["callers"]:
    callers.append(caller)



rule run:
    input:
        expand(pjoin(ODIR, "{aligner}", "{caller}.vcf.gz"),
               aligner = aligners,
               caller = callers),


################
### ALIGNERS ###
################

rule minimap2:
    input:
        fa = REF,
        fq = SAMPLE
    output:
        bam = pjoin(ODIR, "minimap2.md.bam")
    params:
        sam = pjoin(ODIR, "minimap2.sam"),
        name = RUN
    threads: THREADS
    benchmark: pjoin(ODIR, "benchmark", "minimap2.txt")
    shell:
        """
        minimap2 -ax map-hifi --MD --eqx -Y -R \'@RG\\tID:{params.name}\' -t {threads} {input.fa} {input.fq} -o {params.sam}
        samtools view -u {params.sam} | samtools sort -T {output.bam}.sort-tmp > {output.bam}
        samtools index {output.bam}
        """

rule pbmm2:
    input:
        fa = REF,
        fq = SAMPLE
    output:
        bam = pjoin(ODIR, "pbmm2.md.bam")
    params:
        bam = pjoin(ODIR, "pbmm2.bam"),
        name = RUN
    threads: THREADS
    benchmark: pjoin(ODIR, "benchmark", "pbmm2.txt")
    conda: "envs/pbmm2.yaml"
    shell:
        """
        pbmm2 align -j {threads} --preset CCS --sort --rg \'@RG\\tID:{params.name}\' --sample {params.name} {input.fa} {input.fq} {params.bam}
        samtools index {params.bam}
        samtools calmd -b {params.bam} {input.fa} > {output.bam}
        samtools index {output.bam}
        """

rule ngmlr:
    input:
        fa = REF,
        fq = SAMPLE
    output:
        bam = pjoin(ODIR, "ngmlr.md.bam")
    params:
        sam = pjoin(ODIR, "ngmlr.sam"),
        name = RUN
    threads: THREADS
    benchmark: pjoin(ODIR, "benchmark", "ngmlr.txt")
    shell:
        """
        {NGMLR_BIN} -t {threads} --rg-id {params.name} -r {input.fa} -q {input.fq} -o {params.sam}
        samtools view -uS {params.sam} | samtools sort -T {output.bam}.sort-tmp > {output.bam}
        samtools index {output.bam}
        """

# This is just for SVDSS
rule primaligns:
    input:
        bam = pjoin(ODIR, "{aligner}.md.bam")
    output:
        bam = pjoin(ODIR, "{aligner}.md.primary.bam")
    threads: 1
    shell:
        """
        samtools view -b -F 2304 {input.bam} > {output.bam}
        samtools index {output.bam}
        """

###############
### CALLERS ###
###############

rule pbsv_discover:
    input:
        bed = TRF,
        bam = pjoin(ODIR, "{aligner}.md.bam")
    output:
        svsig = pjoin(ODIR, "{aligner}", "pbsv", "variations.svsig.gz")
    threads: 1
    benchmark: pjoin(ODIR, "benchmark", "{aligner}", "pbsv-1.txt")
    shell:
        """
	    pbsv discover --tandem-repeats {input.bed} {input.bam} {output.svsig}
        """

rule pbsv_call:
    input:
        fa = REF,
        svsig = pjoin(ODIR, "{aligner}", "pbsv", "variations.svsig.gz")
    output:
        vcf = pjoin(ODIR, "{aligner}", "pbsv", "variations.vcf")
    threads: THREADS
    benchmark: pjoin(ODIR, "benchmark", "{aligner}", "pbsv-2.txt")
    shell:
        """
        pbsv call -j {threads} --ccs -t INS,DEL {input.fa} {input.svsig} {output.vcf}
        """
# -m 20

rule pbsv_post:
    input:
        vcf = pjoin(ODIR, "{aligner}", "pbsv", "variations.vcf")
    output:
        vcf = pjoin(ODIR, "{aligner}", "pbsv.vcf")
    shell:
        """
        mv {input.vcf} {output.vcf}
        """


rule cutesv:
    input:
        fa = REF,
        bam = pjoin(ODIR, "{aligner}.md.bam")
    output:
        vcf = pjoin(ODIR, "{aligner}", "cutesv", "variations.vcf")
    params:
        wdir = pjoin(ODIR, "{aligner}", "cutesv")
    threads: THREADS
    benchmark: pjoin(ODIR, "benchmark", "{aligner}", "cutesv.txt")
    shell:
        """
        mkdir -p {params.wdir}
        cuteSV -t {threads} -s 2 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 {input.bam} {input.fa} {output.vcf} {params.wdir}
        """
# -l 30

rule cutesv_post:
    input:
        vcf = pjoin(ODIR, "{aligner}", "cutesv", "variations.vcf")
    output:
        vcf = pjoin(ODIR, "{aligner}", "cutesv.vcf")
    shell:
        """
        grep -v 'INV\|DUP\|BND' {input.vcf} > {output.vcf}
        sed -i '4 a ##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand ? (line manually added)">' {output.vcf}
        """

rule svim:
    input:
        fa = REF,
        bam = pjoin(ODIR, "{aligner}.md.bam")
    output:
        odir = directory(pjoin(ODIR, "{aligner}", "svim"))
    threads: 1
    benchmark: pjoin(ODIR, "benchmark", "{aligner}", "svim.txt")
    shell:
        """
        svim alignment --min_sv_size 30 --cluster_max_distance 1.4 --interspersed_duplications_as_insertions --tandem_duplications_as_insertions {output.odir} {input.bam} {input.fa}
        """

rule svim_post:
    input:
        idir = pjoin(ODIR, "{aligner}", "svim")
    output:
        vcf = pjoin(ODIR, "{aligner}", "svim.vcf")
    threads: 1
    shell:
        """
        cat {input.idir}/variants.vcf \
        | sed 's/INS:NOVEL/INS/g' \
        | sed 's/DUP:INT/INS/g' \
        | sed 's/DUP:TANDEM/INS/g' \
        | awk '{{ if($1 ~ /^#/) {{ print $0 }} else {{ if($5=="<DEL>" || $5=="<INS>") {{ print $0 }} }} }}' \
        | grep -v 'SUPPORT=1;' \
        | sed 's/q5/PASS/g' > {output.vcf}
        """
# | grep -v 'SUPPORT=1;\|SUPPORT=2;\|SUPPORT=3;' \


# rule sniffles2:
#     input:
#         fa = REF,
#         bam = pjoin(ODIR, "{aligner}.md.bam"),
#         bed = TRF
#     output:
#         vcf = pjoin(ODIR, "{aligner}", "sniffles2", "variations.vcf")
#     threads: THREADS
#     benchmark: pjoin(ODIR, "benchmark", "{aligner}", "sniffles2.txt")
#     conda: "envs/sniffles2.yaml"
#     shell:
#         """
#         sniffles --input {input.bam} --vcf {output.vcf} --reference {input.fa} --tandem-repeats {input.bed} --threads {threads} --minsupport 2 --minsvlen 30
#         """

# rule sniffles2_post:
#     input:
#         vcf = pjoin(ODIR, "{aligner}", "sniffles2", "variations.vcf")
#     output:
#         vcf = pjoin(ODIR, "{aligner}", "sniffles2.vcf")
#     shell:
#         """
#         cat <(cat {input.vcf} | grep "^#") <(cat {input.vcf} | grep -vE "^#" | grep 'DUP\|INS\|DEL' | sed 's/DUP/INS/g' | sort -k1,1 -k2,2g) > {output.vcf}
#         """

rule debreak:
    input:
        fa = REF,
        bam = pjoin(ODIR, "{aligner}.md.bam"),
        bed = TRF
    output:
        vcf = pjoin(ODIR, "{aligner}", "debreak", "debreak.vcf")
    params:
        odir = pjoin(ODIR, "{aligner}", "debreak"),
        aligner = "{aligner}"
    threads: THREADS
    benchmark: pjoin(ODIR, "benchmark", "{aligner}", "debreak.txt")
    conda: "envs/debreak.yaml"
    shell:
        """
        debreak --bam {input.bam} --outpath {params.odir} --rescue_large_ins --poa --ref {input.fa} --min_size 30 --min_support 2 --aligner params.aligner --thread {threads}
        """

rule debreak_post:
    input:
        debreak = pjoin(ODIR, "{aligner}", "debreak", "debreak.vcf"),
        pp = pjoin(ODIR, "{aligner}", "pp.vcf")
    output:
        vcf = pjoin(ODIR, "{aligner}", "debreak.vcf")
    threads: 1
    shell:
        """
        bash clean_debreak.sh {input.debreak} {input.pp} > {output.vcf}
        """


#############
### SVDSS ###
#############
rule pp_index:
    input:
        fa = REF
    output:
        fmd = REF + ".fmd"
    threads: 1
    benchmark: pjoin(ODIR, "benchmark", "pp-index.txt")
    shell:
        """
        {PP_BIN} index --fastq {input.fa} --index {output.fmd}
        """

rule pp_smooth:
    input:
        fa = REF,
        bam = pjoin(ODIR, "{aligner}.md.primary.bam")
    output:
        bam = pjoin(ODIR, "{aligner}", "pp", "smoothed.selective.bam")
    params:
        wd = pjoin(ODIR, "{aligner}", "pp")
    threads: THREADS
    benchmark: pjoin(ODIR, "benchmark", "{aligner}", "pp-2-smooth.txt")
    shell:
        """
        {PP_BIN} smooth --reference {input.fa} --bam {input.bam} --workdir {params.wd} --threads {threads}
        """

rule pp_sortsmoothed:
    input:
        bam = pjoin(ODIR, "{aligner}", "pp", "smoothed.selective.bam")
    output:
        bam = pjoin(ODIR, "{aligner}", "pp", "smoothed.selective.sorted.bam")
    threads: 1
    benchmark: pjoin(ODIR, "benchmark", "{aligner}", "pp-3-sort.txt")
    shell:
        """
        samtools sort -T {output.bam}.sort-tmp {input.bam} > {output.bam}
        samtools index {output.bam}
        """

rule pp_search:
    input:
        fmd = REF + ".fmd",
        bam = pjoin(ODIR, "{aligner}", "pp", "smoothed.selective.sorted.bam")
    output:
        sfs = pjoin(ODIR, "{aligner}", "pp", "solution_batch_0.assembled.sfs")
    params:
        wd = pjoin(ODIR, "{aligner}", "pp")
    threads: THREADS
    benchmark: pjoin(ODIR, "benchmark", "{aligner}", "pp-4-search.txt")
    shell:
        """
        {PP_BIN} search --index {input.fmd} --bam {input.bam} --threads {threads} --workdir {params.wd} --assemble
        """

rule pp_call:
    input:
        fa = REF,
        bam = pjoin(ODIR, "{aligner}", "pp", "smoothed.selective.sorted.bam"),
        sfs = pjoin(ODIR, "{aligner}", "pp", "solution_batch_0.assembled.sfs")
    output:
        vcf = pjoin(ODIR, "{aligner}", "pp.vcf")
    params:
        wd = pjoin(ODIR, "{aligner}", "pp")
    threads: THREADS
    benchmark: pjoin(ODIR, "benchmark", "{aligner}", "pp-5-call.txt")
    shell:
        """
        n=$(ls {params.wd}/solution_batch_*.assembled.sfs | wc -l)
        {PP_BIN} call --reference {input.fa} --bam {input.bam} --threads {threads} --workdir {params.wd} --batches ${{n}} --min-cluster-weight 2
        mv {params.wd}/svs_poa.vcf {output.vcf}
        """


############
### POST ###
############

rule vcf2gz:
    input:
        "{fname}.vcf"
    output:
        "{fname}.vcf.gz"
    params:
        tmp_prefix = "{fname}.bcftools-sort-tmp"
    threads: 1
    shell:
        """
        bcftools sort -T {params.tmp_prefix} -Oz {input} > {output}
        tabix -p vcf {output}
        """
