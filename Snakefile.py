from pathlib import Path
import re
import pyfastx
import pandas as pd
import numpy as np

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
            "-outfmt", "6 stitle sacc pident slen"
        ]
        
        match = subprocess.check_output(command, universal_newlines=True).strip()
        matches.append(match)
        
    df = df.assign(matches=matches).loc[lambda x: x.matches != ""]
    if df.shape[0] == 0:
        return df
    
    df[["match_name", "accession", "percent_identity", "sequence_len"]] = (
        df.matches.str.split("\t", expand=True).loc[:, :3]
    )
    df = df.assign(sequence_len=lambda x: [y[0] for y in x.sequence_len.str.split("\n")])
    
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
        expand(f"{RESULTS}/{{sample}}/MINIMAP2/{{sample}}.bam", sample=SAMPLES),


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
    output:
        assembly=f"{RESULTS}/{{sample}}/FLYE/assembly.fasta",
    threads:
        config["THREADS"]
    shell:
        """
        flye \
        --meta \
        --threads {threads} \
        --out-dir {params.outdir} \
        --nano-raw {input.fastq}
        """

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
        reference=config["TBE_REF"]
    threads:
        config["THREADS"]
    shell:
        """
        minimap2 -t {threads} -a -x map-ont {params.reference} {input.fastq} | \
        samtools view -h -F 4 -F 256 | samtools sort -o {output.bam}
        """


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
