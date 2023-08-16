from pathlib import Path
import re

# COMMON FUNCTIONS

def get_samples(sample_dir: str) -> list[str]:
    return [
        x.stem for x in Path(SAMPLE_DIR).iterdir() 
        if x.is_dir() 
        and not x.stem.startswith(".")
    ]

def get_fastx(wildcards):
    in_folder = Path(f"{SAMPLE_DIR}/{wildcards.sample}")
    files = [file for file in Path(in_folder).iterdir() if re.search(r"fq|fastq|fa|fasta|fna", file.name)]
    print(files)
    return files


# IO
configfile: "config.yaml"

RESULTS = config["RESULTS"]
SAMPLE_DIR = config["SAMPLES"]
SAMPLES = get_samples(SAMPLE_DIR)


# RULES
rule all:
    input:
        expand(f"{RESULTS}/{{sample}}/MERGED_FASTQ/{{sample}}_merged.fastq.gz", sample=SAMPLES),


# --- CONCAT FASTQ ---

rule concatenate_fastq:
    input:
        fastq_files=get_fastx
    output:
        merged=f"{RESULTS}/{{sample}}/merged_fastq/{{sample}}_merged.fastq.gz"
    shell:
        """
        cat {input.fastq_files} > {output.merged}
        """

# --- PORECHOP ---

rule porechop:
    input:
        fastq=get_fastx
    output:
        chopped=f"{RESULTS}/{{sample}}/{{sample}}_raw_porechop.fastq.gz",
    log:
        f"{RESULTS}/{{sample}}/logs/log_porechops.txt",
    message:
        """
        Running porechop on {input.fastq}
        """
    shell:
        """
        porechop \
        -i {input.fastq} \
        -o {output.chopped} \
        | tee {log}
        """
