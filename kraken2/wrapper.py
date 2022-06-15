#! /usr/bin/env python3
import click
import subprocess
import yaml
from os import remove, mkdir, chdir
from os.path import dirname, exists
from pathlib import Path

@click.command()
@click.option('--input', required = True, help='Input with bam list')
@click.option('--db', required = True, help='Path to Kracken database. Alternatively, type "default" if using the viral db in our docker version.')
@click.option('--outdir', default = ".", help='Output directory. If missing, it will be created')
@click.option('--extra_contigs', default = "none", help='Path to text-file with extra contigs to analyze. Leave none for only unmapped ("*")')
@click.option('--taxa', default = "S", type = click.Choice(["D","P","C","O","F","G","S","S1"]), help='Taxa for abundance estimation with bracken. Default is S(pecies).')
@click.option('--read_length',  default=150, help='Read length', type = int)
@click.option('--qual',  default=20, help='Kraken minimum base quality argument', type = int)
@click.option('--confidence', default=0.6, help='Kraken confidence argument', type = float)
@click.option('--threads', default=4, help='Number of threads per job (4 recommended)', type = int)
@click.option('--kraken_threads', default=10, help='Number of threads for kraken', type = int)
@click.option('--cores', default=16, help='Max number of cores provided to snakemake', type = int)
@click.option('--verbose', default = False, help="Increase verbosity",is_flag=True)
@click.option('--clean/--no-clean', is_flag=True, default = True, help="Whether remove the config file nad logs or not")

def run(input, db, extra_contigs, read_length, threads, verbose, outdir, qual, confidence, taxa, cores, kraken_threads, clean):
    """Simple wrapper for launching Kraken/Bracken with snakemake."""
    assert 0 <= confidence <= 1, "confidence parameter must be within [0-1]" 
    assert (threads <= cores) & (kraken_threads <= cores), "The number of threads per job cannot exceed the max total cores"
    assert exists(input), f"Input file does not exist: {input}"
    assert exists(extra_contigs) or extra_contigs == "none", f"Extra contigs file does not exist: {extra_contigs}"
    assert exists(db) or db == "default", f"Extra contigs file does not exist: {db}"

    click.echo(f"Running Kraken with {cores} cores")

    config = {
        "input": str(Path(input).absolute()),
        "db": str(Path(db).absolute()) if db != "default" else db,
        "contigs": str(Path(extra_contigs).absolute()) if extra_contigs != "none" else extra_contigs,
        "taxa": taxa,
        "read_length": read_length,
        "qual": qual,
        "confidence":confidence,
        "threads": threads,
        "kraken_threads": kraken_threads
    }
    snakefile_dir = str(Path(dirname(__file__)).absolute())

    if outdir != ".":
        if not exists(outdir):
            mkdir(outdir)
        chdir(outdir)

    f = open('config.yaml', 'w+')
    yaml.dump(config, f)
    cmd = f"snakemake --snakefile {snakefile_dir}/snakefile --cores {cores} --configfile config.yaml "

    if not verbose:
        cmd += "--quiet"

    with open("kraken.log","wb") as log:
        try:
            subprocess.check_call(cmd, stdout=log, shell = True)
        except subprocess.CalledProcessError as e:
            print(f'Error in snakemake invocation: {e}', file=log)
            return e.returncode

    if clean:
        remove("config.yaml")
        remove("kraken.log")
if __name__ == '__main__':
    run()
 