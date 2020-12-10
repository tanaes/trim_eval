from os.path import join

configfile: "config.yaml"

indir = config['indir']
outdir = config['outdir']

samples = config['samples']
depths = config['depths']
trims = config['trims']

rule sample:
    input:
        fwd = join(indir, "{sample}.R1.fastq.gz"),
        rev = join(indir, "{sample}.R2.fastq.gz")
    output:
        fwd = temp(join(outdir,
                        'trimmed',
                        "{sample}.{depth}.R1.fastq.gz")),
        rev = temp(join(outdir,
                        'trimmed',
                        "{sample}.{depth}.R2.fastq.gz"))
    conda:
        "env.yaml"
    shell:
        """
        seqtk sample -s100 {input.fwd} {wildcards.depth} > {output.fwd}
        seqtk sample -s100 {input.rev} {wildcards.depth} > {output.rev}
        """


rule trim:
    input:
        fwd = rules.sample.output.fwd,
        rev = rules.sample.output.rev
    output: 
        fwd = join(outdir,
                   'trimmed',
                   "{sample}.{depth}.{trim}.R1.fastq.gz"),
        rev = join(outdir,
                   'trimmed',
                   "{sample}.{depth}.{trim}.R2.fastq.gz")
    conda:
        "env.yaml"
    shell:
        """
        seqtk trimfq -e {wildcards.trim} {input.fwd} > {output.fwd}
        seqtk trimfq -e {wildcards.trim} {input.rev} > {output.rev}
        """


rule assemble:
    input:
        fwd = rules.trim.output.fwd,
        rev = rules.trim.output.rev
    output:
        assembly = join(outdir,
                        'assembled',
                        "{sample}.{depth}.{trim}",
                        'final.contigs.fa')
    conda:
        "env.yaml"
    resources:
        mem_mb=50000
    threads:
        4
    log:
        join(outdir, "logs", "assemble.{sample}.{depth}.{trim}.log")
    shell:
        """
        outdir=$(dirname "{output.assembly}")

        megahit -1 {input.fwd} -2 {input.rev} -o $outdir \
        2> {log} 1>&2
        """


rule quast:
    input:
        rules.assemble.output.assembly
    output:
        report = join(outdir,
                      'quast',
                      "{sample}.{depth}.{trim}",
                      'report.tsv')
    conda:
        "env.yaml"
    log:
        join(outdir, "logs", "assemble.{sample}.{depth}.{trim}.log")
    shell:
        """
        outdir=$(dirname "{output.report}")

        quast.py -t {threads} -o $outdir \
        {input} 2> {log} 1>&2
        """


rule multiqc:
    input:
        expand(rules.quast.output.report,
               sample=samples,
               depth=depths,
               trim=trims)
    output:
        join(outdir, "assembled/multiqc.html")
    wrapper:
        "0.49.0/bio/multiqc"
