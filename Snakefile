SAMPLES = ["SRR2584857_1"]
GENOME = ["ecoli-rel606"]

rule make_bams:
    input:
        expand("outputs/{sample}.x.{genome}.bam.sorted",
               sample=SAMPLES, genome=GENOME)

rule map_reads:
    input:
        reads="{reads}.fastq.gz",
        ref="{genome}.fa.gz"
    output: "outputs/{reads}.x.{genome}.sam"
    shell: """
        minimap2 -ax sr {input.ref} {input.reads} > {output}
    """

rule sam_to_bam:
    input: "outputs/{reads}.x.{genome}.sam",
    output: "outputs/{reads}.x.{genome}.bam",
    shell: """
        samtools view -b -F 4 {input} > {output}
     """

rule sort_bam:
    input: "outputs/{reads}.x.{genome}.bam"
    output: "outputs/{reads}.x.{genome}.bam.sorted"
    shell: """
        samtools sort {input} > {output}
    """

rule uncompress_genome:
    input: "{genome}.fa.gz"
    output: "outputs/{genome}.fa"
    shell: """
        "gunzip -k {input} > {output}
    """

rule call_variants:
    input:
        ref="outputs/{genome}.fa",
        bam="outputs/{reads}.x.{genome}.bam.sorted",
    output:
        pileup="outputs/{reads}.x.{genome}.pileup",
        bcf="outputs/{reads}.x.{genome}.bcf",
        vcf="outputs/{reads}.x.{genome}.vcf",
    shell: """
        bcftools mpileup -Ou -f {output.ref} {input.bam} > {output.pileup}
        bcftools call -mv -Ob {output.pileup} -o {output.bcf}
        bcftools view {output.bcf} > {output.vcf}
    """