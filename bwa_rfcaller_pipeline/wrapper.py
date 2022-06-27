import click
import subprocess
from pyrsistent import s
import yaml
from os import remove, mkdir, chdir
from os.path import dirname, exists
from pathlib import Path


def run_validations(input, comparison, fasta, dbsnp, targets, pon, ploidy, threads, other_threads, cores, memory, max_memory):
    assert exists(input), f"MissingFileError: Input file does not exist: {input} \n"
    assert exists(comparison) and ("rfcaller_vcf" in targets or "all" in targets), f"MissingFileError: RFcaller is expected to run, but comparison file does not exist: {comparison} \n"
    assert exists(fasta), f"MissingFileError: Fasta file does not exist: {fasta} \n"
    assert exists(dsnp), f"dbSNP file does not exist: {dbsnp} \n"
    assert pon.lower() in ["hg38", "hg19"] or exists(pon), f"PoNError: The pon variable isn't any of hg38 or hg19 and doesn't not exist either: {pon} \n"
    assert ploidy.lower() in ["grch38", "grch19"] or exists(ploidy), f"PloidyError: The ploidy variable isn't any of GRCh38 or GRCh37 and doesn't not exist either: {ploidy} \n"
    assert (threads + other_threads <= cores), "CoresError: The BWA threads and other_threads sum cannot exceed the max total cores \n"
    assert (max_memory > memory), "MemoryError: The max_memory variable ({max_memory}) must be greater that the memory variable ({memory})"

    ploidy_exists = True if exists(ploidy) else False
    pon_exists = True if exists(pon) else False

    return [ploidy_exists, pon_exists]

@click.command()
@click.option('--input', required = True, help='\b\nInput with fastqs metadata. Fields: PATIENT, SAMPLE, EXPERIMENT, F1, F2')
@click.option('--comparison', required = True, help='\b\nTwo-column file with experiment names for somatic variant calling. Fields: TUMOR_EXPERIMENT, NORMAL_EXPERIMENT')
@click.option('--fasta', required = True,  help='Path to reference genome fasta')
@click.option('--dbsnp', required = True,  help='Path to dbsnp')
@click.argument('targets', nargs=-1, type = click.Choice(["all", "rg_bams", "discordants_split", "bqsr_bams", "rfcaller_vcf"], case_sensitive=False))
@click.option('--PL', default = "ILLUMINA-NOVASEQ-6000", help='Name of the sequencing platform')
@click.option('--CN', default = "MACROGEN", help='Name of the platform center')
@click.option('--PL', default = "SM1", help='Sample ID')
@click.option('--LB', default = "LB1", help='Library ID')
@click.option('--pon',  default= "hg38", help='\b\nPanel of normals. hg38 and hg19 are built-in, or path to custom pon', type = str)
@click.option('--ploidy',  default= "GRCh38", help='Bcftools ploidy file (GRCh38 or GRCh37)', type = str)
@click.option('--threads', default=24, help='Number of threads per job in BWA',  type = int)
@click.option('--other_threads', default=5, help='\b\nNumber of threads for minor steps (merge, sort, ...)', type = int)
@click.option('--cores', default=60, help='Max number of cores provided to snakemake', type = int)
@click.option('--memory', default=40, help='\b\nMax RAM memory allowed to sort (Gb) per FASTQ pair', type = int)
@click.option('--max_memory', default=100, help='Max total RAM memory allowed (Gb)', type = int)
@click.option('--verbose', default = False, help="Increase verbosity",is_flag=True)
@click.option('--clean/--no-clean', is_flag=True, default = True, help="\b\nWhether to remove the config file and logs or not")
@click.option('--outdir', default = ".", help='\b\nOutput directory. If missing, it will be created')

def run(input, comparison, fasta, dbsnp, targets, PL, CN, SM, LB, pon, ploidy, threads, other_threads, cores, memory, max_memory, verbose, clean, outdir):
    """Simple wrapper for launching BWA-mem2/RFcaller with snakemake."""

    ploidy_exists, pon_exists = run_validations(input, comparison, fasta, dbsnp, targets, pon, ploidy, threads, other_threads, cores, memory, max_memory)

    targets = " ".join(["rg_bams", "discordants_split", "bqsr_bams", "rfcaller_vcf"]) if "all" in targets else " ".join(targets)

    click.echo(f"Running BWA/RFCaller with {cores} cores and {max_memory}G RAM.")

    config = {
        "input": str(Path(input).absolute()),
        "comparison": str(Path(comparison).absolute()),
        "fasta": str(Path(fasta).absolute()),
        "dbsnp": str(Path(dbsnp).absolute()),
        "PL": PL,
        "SM": SM,
        "LB": LB,
        "CN": CN,
        "pon": str(Path(pon).absolute()) if pon_exists else pon,
        "ploidy": str(Path(ploidy).absolute()) if ploidy_exists else ploidy,
        "threads": threads,
        "other_threads": other_threads,
        "memory": memory
    }
    
    snakefile_dir = str(Path(dirname(__file__)).absolute())

    if outdir != ".":
        if not exists(outdir):
            mkdir(outdir)
        chdir(outdir)

    f = open('config.yaml', 'w+')
    yaml.dump(config, f)
    cmd = f"snakemake --snakefile {snakefile_dir}/Snakefile --cores {cores} --resources mem_gb = {max_memory} --configfile config.yaml {targets}"

    if not verbose:
        cmd += "--quiet"
    print(cmd)
    with open("run.log","w+") as log:
        try:
            subprocess.check_call(cmd, stdout=log, shell = True)
        except subprocess.CalledProcessError as e:
            print(f'Error in snakemake invocation: {e}', file=log)
            return e.returncode

    if clean:
        remove("config.yaml")
        remove("run.log")

if __name__ == '__main__':
    run()

