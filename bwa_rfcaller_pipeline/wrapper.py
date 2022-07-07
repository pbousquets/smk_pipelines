#! /usr/bin/env python3
import click
import subprocess
import yaml
from os import remove, mkdir, chdir
from os.path import dirname, exists
from pathlib import Path


def run_validations(fastqs, comparison, fasta, dbsnp, targets, pon, ploidy, threads, other_threads, cores, memory, max_memory):
    if fastqs:
        assert exists(fastqs), f"MissingFileError: FastQs metadata file does not exist: {fastqs} \n"
    assert exists(comparison) and ("rfcaller_vcf" in targets or "all" in targets), f"MissingFileError: RFcaller is expected to run, but comparison file does not exist: {comparison} \n"
    assert exists(fasta), f"MissingFileError: Fasta file does not exist: {fasta} \n"
    assert exists(dbsnp), f"dbSNP file does not exist: {dbsnp} \n"
    assert pon.lower() in ["hg38", "hg19"] or exists(pon), f"PoNError: The pon variable isn't any of hg38 or hg19 and doesn't not exist either: {pon} \n"
    assert ploidy.lower() in ["grch38", "grch19"] or exists(ploidy), f"PloidyError: The ploidy variable isn't any of GRCh38 or GRCh37 and doesn't not exist either: {ploidy} \n"
    assert (threads + other_threads <= cores), "CoresError: The BWA threads and other_threads sum cannot exceed the max total cores \n"
    assert (max_memory > memory), "MemoryError: The max_memory variable ({max_memory}) must be greater that the memory variable ({memory})"

    ploidy_exists = True if exists(ploidy) else False
    pon_exists = True if exists(pon) else False

    return [ploidy_exists, pon_exists]

@click.command(context_settings={'show_default': True})
@click.option('--fastqs', help='\b\nInput with fastqs metadata. Fields: PATIENT, SAMPLE, EXPERIMENT, F1, F2')
@click.option('--comparison', required = True, help='\b\nTwo-column file with experiment names for somatic variant calling. If skipping alignment, add the paths to bam files in the additional fields. Fields: TUMOR_EXPERIMENT, NORMAL_EXPERIMENT [TUMOR_BAM, NORMAL_BAM]')
@click.option('--fasta', required = True,  help='Path to reference genome fasta')
@click.option('--dbsnp', required = True,  help='Path to dbsnp')
@click.option('--outdir', default = ".", help='\b\nOutput directory. If missing, it will be created')

@click.option('--platform', default = "ILLUMINA-NOVASEQ-6000", help='\b\nName of the sequencing platform')
@click.option('--center', default = "MACROGEN", help='Name of the platform center')
@click.option('--pon',  default= "hg38", help='\b\nPanel of normals. hg38 and hg19 are built-in, or path to custom pon', type = str)
@click.option('--ploidy',  default= "GRCh38", help='\b\nBcftools ploidy file (GRCh38 or GRCh37)', type = str)

@click.option('--threads', default=24, help='Number of threads per job in BWA',  type = int)
@click.option('--other_threads', default=5, help='\b\nNumber of threads for minor steps (merge, sort, ...)', type = int)
@click.option('--cores', default=60, help='\b\nMax number of cores provided to snakemake', type = int)
@click.option('--memory', default=40, help='\b\nMax RAM memory allowed to sort (Gb) per FASTQ pair', type = int)
@click.option('--max_memory', default=100, help='Max total RAM memory allowed (Gb)', type = int)

@click.option('--verbose', default = False, help="Increase verbosity",is_flag=True)
@click.option('--clean/--no-clean', is_flag=True, default = True, help="\b\nWhether to remove the config file and logs or not")
@click.option('--dryrun', is_flag=True, default = False, help="\b\nMake validations only and launch snakemake in dry-run mode")
@click.argument('targets', nargs=-1, type = click.Choice(["all", "rg_bams", "discordants_split", "bqsr_bams", "rfcaller_vcf"], case_sensitive=False))

def run(fastqs, comparison, fasta, dbsnp, targets, platform, center, pon, ploidy, threads, other_threads, cores, memory, max_memory, verbose, clean, dryrun, outdir):
    """Simple wrapper for launching BWA-mem2/RFcaller with snakemake."""

    ploidy_exists, pon_exists = run_validations(fastqs, comparison, fasta, dbsnp, targets, pon, ploidy, threads, other_threads, cores, memory, max_memory)

    targets = " ".join(["rg_bams", "discordants_split", "bqsr_bams", "rfcaller_vcf"]) if "all" in targets else " ".join(targets)

    click.echo(f"Running BWA/RFCaller with {cores} cores and {max_memory}G RAM.")

    config = {
        "comparison": str(Path(comparison).absolute()),
        "fasta": str(Path(fasta).absolute()),
        "dbsnp": str(Path(dbsnp).absolute()),
        "PL": platform,
        "CN": center,
        "pon": str(Path(pon).absolute()) if pon_exists else pon,
        "ploidy": str(Path(ploidy).absolute()) if ploidy_exists else ploidy,
        "threads": threads,
        "other_threads": other_threads,
        "memory": memory,
        "workdir": str(Path(outdir).absolute())
    }

    if fastqs:
        config["fastqs"] = str(Path(fastqs).absolute())

    snakefile_dir = str(Path(dirname(__file__)).absolute())

    f = open('config.yaml', 'w+')
    yaml.dump(config, f)
    cmd = f"snakemake --snakefile {snakefile_dir}/Snakefile --cores {cores} --resources mem_gb={max_memory} --configfile config.yaml "

    if not verbose:
        cmd += "--quiet "

    if dryrun:
        cmd += "-np "

    cmd += f"{targets}" 
    print(cmd)

    with open("run.log","w+") as log:
        try:
            if verbose:
                subprocess.check_call(cmd, shell = True)
            else:
                subprocess.check_call(cmd, stdout=log, shell = True)
                
        except subprocess.CalledProcessError as e:
            print(f'Error in snakemake invocation: {e}', file=log)
            if verbose:
                print(f'Error in snakemake invocation: {e}')
            return e.returncode

    if clean:
        remove("config.yaml")
        remove("run.log")

if __name__ == '__main__':
    run()


