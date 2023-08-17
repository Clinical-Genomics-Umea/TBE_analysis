from pathlib import Path
import re
import pyfastx
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# COMMON FUNCTIONS

def get_samples(sample_dir: str) -> list[str]:
    return [
        x.stem for x in Path(SAMPLE_DIR).iterdir() 
        if x.is_dir() 
        and not x.stem.startswith(".")
    ]

def get_fastx(wildcards):
    in_folder = Path(f"{SAMPLE_DIR}/{wildcards.sample}")
    return [file for file in Path(in_folder).iterdir() if re.search(r"fq|fastq|fa|fasta|fna", file.name)]


def fastx_file_to_df(fastx_file: str) -> pd.DataFrame:
    fastx = pyfastx.Fastx(fastx_file)
    reads = list(zip(*[[x[0], x[1]] for x in fastx]))

    df = (
        pd.DataFrame({"name": reads[0], "sequence": reads[1]})
        .assign(read_len=lambda x: x.sequence.str.len())
        .sort_values("read_len", ascending=False)
    )
    
    return df

def revcomp(dna_seq):
    return dna_seq[::-1].translate(str.maketrans("ATGC", "TACG"))


def run_blastn(contigs: str, db: str, temp_file: str) -> pd.DataFrame:
    # In order for blastn to find the files
    os.environ["BLASTDB"] = config["BLASTDB_ENVIRON_VARIABLE"]
    df = pd.read_csv(contigs)
    if df.shape[0] == 0:
        return df
    
    matches = []
    
    for contig in df.itertuples():
        with open(temp_file, "w+") as f:
            print(f">{contig.name}", file=f)
            print(contig.sequence, file=f)
        command = [
            "blastn", "-query", temp_file, "-db", db, "-max_target_seqs", "1",
            "-outfmt", "6 stitle sacc pident slen sstart send sstrand"
        ]
        
        match = subprocess.check_output(command, universal_newlines=True).strip()
        matches.append(match)
        
    df = df.assign(matches=matches).loc[lambda x: x.matches != ""]
    if df.shape[0] == 0:
        return df
    
    df[["match_name", "accession", "percent_identity", "subject_len", "align_start", "align_end", "subject_strand"]] = (
       df.matches.str.split("\t", expand=True).loc[:, 0:6]
    )
    df = ( 
        df
        .assign(subject_strand=lambda x: [y[0] for y in x.subject_strand.str.split("\n")])
        .assign(sequence=lambda x: [y.sequence if y.subject_strand == "plus" else revcomp(y.sequence) for y in x.itertuples()])
        .assign(align_start_temp=lambda x: x.align_start)
        .assign(align_start=lambda x: [y.align_start if y.subject_strand == "plus" else y.align_end for y in x.itertuples()])
        .assign(align_end=lambda x: [y.align_end if y.subject_strand == "plus" else y.align_start_temp for y in x.itertuples()])
        .drop(columns="align_start_temp")
        .sort_values("align_start")
    )
    return df

# IO
configfile: "config.yaml"

RESULTS = config["RESULTS"]
SAMPLE_DIR = config["SAMPLES"]
SAMPLES = get_samples(SAMPLE_DIR)


# -- RULES --- #
rule all:
    input:
        expand(f"{RESULTS}/{{sample}}/BLASTN/{{sample}}_contigs_blastn.csv", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/COVERAGE/{{sample}}_coverage_plot.svg", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/CONSENSUS/{{sample}}_consensus.fa", sample=SAMPLES),


rule merge:
    input:
        fastq_files=get_fastx
    output:
        merged=f"{RESULTS}/{{sample}}/MERGED_FASTQ/{{sample}}_merged.fastq.gz"
    shell:
        """
        cat {input.fastq_files} > {output.merged}
        """

rule porechop:
    input:
        fastq=rules.merge.output.merged
    output:
        trimmed=f"{RESULTS}/{{sample}}/PORECHOP/{{sample}}_trimmed.fastq.gz",
    log:
        f"{RESULTS}/{{sample}}/logs/porechop.log",
    shell:
        """
        porechop \
        -i {input.fastq} \
        -o {output.trimmed} \
        > {log} 2>&1
        """

rule flye:
    input:
        fastq=rules.porechop.output.trimmed
    params:
        outdir=f"{RESULTS}/{{sample}}/FLYE",
        genome_size=config["GENOME_SIZE"]
    output:
        assembly=f"{RESULTS}/{{sample}}/FLYE/assembly.fasta",
    threads:
        config["THREADS"]
    shell:
        """
        flye \
        --meta \
        --threads {threads} \
        --genome-size {params.genome_size} \
        --out-dir {params.outdir} \
        --nano-raw {input.fastq}
        """
        #--asm-coverage 50 \

rule wrangle_flye:
    input:
        contigs=rules.flye.output.assembly
    output:
        csv=f"{RESULTS}/{{sample}}/FLYE/{{sample}}_contigs.csv"
    run:
        df = fastx_file_to_df(input.contigs)
        df = df.assign(sample_id=wildcards.sample)
        df.to_csv(output.csv, index=False)
        
rule blastn:
    input:
        contigs=rules.wrangle_flye.output.csv
    output:
        blast=f"{RESULTS}/{{sample}}/BLASTN/{{sample}}_contigs_blastn.csv"
    params:
        temp_file=f"{RESULTS}/{{sample}}/BLASTN/temp.blastn",
        db=config["BLASTN_DB"],
    run:
        df = run_blastn(contigs=input.contigs, db=params.db, temp_file=params.temp_file)
        df.to_csv(output.blast, index=False)
        
        os.remove(params.temp_file)



rule minimap2:
    input:
        fastq=rules.porechop.output.trimmed
    output:
        bam=f"{RESULTS}/{{sample}}/MINIMAP2/{{sample}}.bam",
    params:
        reference=config["LANGAT_REF"],
    threads:
        config["THREADS"]
    shell:
        """
        minimap2 -t {threads} -a -x map-ont {params.reference} {input.fastq} | \
        samtools view -h -F 4 -F 256 -F 2048 | samtools sort -o {output.bam}
        """

rule mpileup:
    input:
        bam=rules.minimap2.output.bam
    output:
        mpileup=f"{RESULTS}/{{sample}}/COVERAGE/{{sample}}_mpileup.tsv",
    shell:
        """
        samtools mpileup -o {output.mpileup} -a {input.bam}
        """

rule consensus:
    input:
        mpileup=rules.mpileup.output.mpileup
    output:
        consensus=f"{RESULTS}/{{sample}}/CONSENSUS/{{sample}}_consensus.fa",
    shell:
        """
        cat {input.mpileup} | ivar consensus -t .7 -m 10 -p {wildcards.sample}_consensus
        mv {wildcards.sample}_consensus.fa {output.consensus}
        rm {wildcards.sample}_consensus.qual.txt
        """


rule plot_coverage:
    input:
        mpileup=rules.mpileup.output.mpileup,
    params:
        treshold=config["COVERAGE_THRESHOLD"],
    output:
        coverage_plot_svg=f"{RESULTS}/{{sample}}/COVERAGE/{{sample}}_coverage_plot.svg",
        coverage_plot_png=f"{RESULTS}/{{sample}}/COVERAGE/{{sample}}_coverage_plot.png",
    run:
        COVERAGE_THRESHOLD = params.treshold
        
        def parse_mpileup(mpileup: str) -> pd.DataFrame:
            return pd.read_csv(
                mpileup,
                sep="\t",
                header=None,
                names=["id", "pos", "base", "coverage", "x", "y"],
            )

        def coverage_plot(mpileup: str, rolling_average: int = 10):

            mpileup_df = parse_mpileup(mpileup).assign(
                depth=lambda x: x.coverage.rolling(rolling_average).mean()
            )

            mean_coverage = mpileup_df.coverage.mean()
            coverage = (
                sum(1 if x > COVERAGE_THRESHOLD else 0 for x in mpileup_df.coverage)
                / mpileup_df.shape[0]
                * 100
            )

            sns.set_theme()
            coverage_plot = plt.figure(figsize=(15, 8))
            sns.lineplot(data=mpileup_df, x="pos", y="depth")
            zero = plt.axhline(y=0, color="red")
            zero.set_label("Zero")
            mean = plt.axhline(y=mean_coverage, color="green")
            mean.set_label(f"Mean coverage: {mean_coverage: .1f}X")
            plt.legend(loc="upper right")
            plt.title(f"Percent bases with coverage above {COVERAGE_THRESHOLD}X: {coverage: .1f}%")
            plt.suptitle(f"Ref: {mpileup_df.iloc[0].id} | Sample: {wildcards.sample}")
            #plt.text(0.05, 0.85, f"Sequence identity: {identity: .1f}", bbox={'facecolor': 'oldlace', 'alpha': 0.8, 'pad': 8})
            plt.close()
            return coverage_plot

        plot = coverage_plot(input.mpileup)
        plot.savefig(output.coverage_plot_svg)
        plot.savefig(output.coverage_plot_png)


# If viralFlye
rule viralFlye:
    input:
        fastq=rules.porechop.output.trimmed,
        flye_dir=rules.flye.params.outdir,
    params:
        outdir=f"{RESULTS}/{{sample}}/VIRALFLYE",
        hmm=config["HMM"],
    output:
        # TODO
        assembly=f"{RESULTS}/{{sample}}/VIRAFLY/TODO.fasta",
    threads:
        config["THREADS"]
    shell:
        """
        viralFlye.py \
        --dir {input.flye_dir} \
        --hmm {params.HMM} \
        --reads {input.fastq} \
        --outdir {params.outdir} 
        """
