#!/opt/conda/envs/inStrain/bin/python3
import asyncio
import aiofiles
import os
import click
import pathlib

async def add_prefix_to_fasta(input_fasta: str, output_fasta: str, prefix: str):
    """Asynchronously reads a FASTA file and adds a prefix to headers."""
    async with aiofiles.open(input_fasta, mode='r') as infile, aiofiles.open(output_fasta, mode='w') as outfile:
        async for line in infile:
            if line.startswith(">"):
                await outfile.write(f">{prefix}{line[1:]}")  # Modify header
            else:
                await outfile.write(line)  # Keep sequence unchanged

async def process_bins(bin_files, output_dir, prefix_template):
    """Processes multiple FASTA bins asynchronously."""
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists
    tasks = []

    for bin_file in bin_files:
        prefix = prefix_template.format(binname=pathlib.Path(bin_file).name.replace(".","_")) 
        output_fasta = os.path.join(output_dir, f"{os.path.basename(bin_file)}")
        tasks.append(add_prefix_to_fasta(bin_file, output_fasta, prefix))

    await asyncio.gather(*tasks)  # Run all tasks concurrently
    click.echo(f"Processed {len(bin_files)} bins. Saved in '{output_dir}'.")

@click.command()
@click.argument("bin_files", type=click.Path(exists=True), nargs=-1)
@click.option("--output-dir", "-o", type=click.Path(), default="output_bins", help="Directory to save modified FASTA files.")
@click.option("--prefix", "-p", default="{binname}_", help="Prefix template for FASTA headers.")
def cli(bin_files, output_dir, prefix):
    """CLI tool to add prefixes to FASTA headers in multiple bins."""
    if not bin_files:
        click.echo("Error: No FASTA files provided.", err=True)
        return
    asyncio.run(process_bins(bin_files, output_dir, prefix))

if __name__ == "__main__":
    cli()