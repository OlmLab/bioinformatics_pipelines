#!/opt/conda/envs/roadmap6/bin/python3
import re
from pathlib import Path   
import argparse

TEMPLATE="""#Default config file for eukdetect. Copy and edit for analysis

#Directory where EukDetect output should be written
output_dir: <OUTPUT_DIRECTORY_PATH>
  
#Indicate whether reads are paired (true) or single (false)
paired_end: <IS_PAIRED_END> 

#filename excluding sample name. no need to edit if paired_end = false
fwd_suffix: <FWD_SUFFIX> 

#filename excludign sample name. no need to edit if paired_end = false
rev_suffix: <REV_SUFFIX>

#file name excluding sample name. no need to edit if paired_end = true 
se_suffix: <SE_SUFFIX> 

#length of your reads. pre-trimming reads not recommended
readlen: <READ_LENGTH>

#full path to directory with raw fastq files
fq_dir: <FASTQ_DIRECTORY_PATH>

#full path to folder with eukdetect database files
database_dir: <DATABASE_DIRECTORY_PATH>

#name of database. Default is original genomes only database name
database_prefix: "ncbi_eukprot_met_arch_markers.fna"

#full path to eukdetect installation folder
eukdetect_dir: <EUKDETECT_DIRECTORY_PATH>

#list sample names here. fastqs must correspond to {samplename}{se_suffix} for SE reads or {samplename}{fwd_suffix} and {samplename}{rev_suffix} for PE
#each sample name should be preceded by 2 spaces and followed by a colon character
samples:
<SAMPLENAMES>
"""


class Config:
    def __init__(
        self,
        sample_name: str,
        read1: Path,
        read2: Path = None,
        output_dir: Path = Path("./eukdetect_output"),
        db_path: Path = Path("./eukdb"),
        read_length: int|None = None,
        eukdetect_dir: Path = Path("."),
    ):
        self.sample_name = sample_name
        self.read1 = read1
        self.read2 = read2
        self.output_dir = output_dir
        self.db_path = db_path
        self.eukdetect_dir = eukdetect_dir
        if read2:
            self.paired_end = "true"
        else:
            self.paired_end = "false"
        self.file_type= "compressed" if self.read1.suffix in [".gz", ".bz2", ".xz"] else "uncompressed"

        self.read_length = read_length if read_length is not None else self._infer_read_length()
        if self.paired_end == "true":
            self.fwd_suffix = self.read1.name.replace(sample_name, "")
            self.rev_suffix = self.read2.name.replace(sample_name, "")
            self.se_suffix = None
        else:
            self.fwd_suffix = None
            self.rev_suffix = None
            self.se_suffix = self.read1.name.replace(sample_name, "")
        
    
    def _infer_read_length(self) -> int:
        if self.file_type == "compressed":
            import gzip
            open_func = gzip.open
        else:
            open_func = open
        with open_func(self.read1, 'rt') as f:
            first_read = f.readline().strip()  # Skip header
            seq = f.readline().strip()  # Read sequence line
            return len(seq)
    
    def generate_config(self) -> str:
        fq_dir = self.read1.parent.resolve()
        samplenames = f"  {self.sample_name}:"
        config_content = TEMPLATE
        config_content = re.sub(r"<OUTPUT_DIRECTORY_PATH>", str(self.output_dir.resolve()), config_content)
        config_content = re.sub(r"<IS_PAIRED_END>", self.paired_end, config_content)
        config_content = re.sub(r"<FWD_SUFFIX>", str(self.fwd_suffix) if self.fwd_suffix else "", config_content)
        config_content = re.sub(r"<REV_SUFFIX>", str(self.rev_suffix) if self.rev_suffix else "", config_content)
        config_content = re.sub(r"<SE_SUFFIX>", str(self.se_suffix) if self.se_suffix else "", config_content)
        config_content = re.sub(r"<READ_LENGTH>", str(self.read_length), config_content)
        config_content = re.sub(r"<FASTQ_DIRECTORY_PATH>", str(fq_dir), config_content)
        config_content = re.sub(r"<DATABASE_DIRECTORY_PATH>", str(self.db_path.resolve()), config_content)
        config_content = re.sub(r"<EUKDETECT_DIRECTORY_PATH>", str(self.eukdetect_dir.resolve()), config_content)
        config_content = re.sub(r"<SAMPLENAMES>", samplenames, config_content)
        return config_content

def write_config_to_file(config: Config, filepath: Path):
    config_content = config.generate_config()
    with open(filepath, 'w') as f:
        f.write(config_content)

        
    
        
def parse_args():
    parser = argparse.ArgumentParser(description="Generate EukDetect config file.")
    parser.add_argument("--sample-name", required=True, help="Sample name.")
    parser.add_argument("--read1", required=True, type=Path, help="Path to read 1 fastq file.")
    parser.add_argument("--read2", type=Path, help="Path to read 2 fastq file (if paired-end).")
    parser.add_argument("--output-dir", type=Path, default=Path("./eukdetect_output"), help="Output directory for EukDetect results.")
    parser.add_argument("--db-path", type=Path, default=Path("./eukdb"), help="Path to EukDetect database directory.")
    parser.add_argument("--read-length", type=int, help="Length of reads. If not provided, inferred from read1.")
    parser.add_argument("--eukdetect-dir", type=Path, default=Path("."), help="Path to EukDetect installation directory.")
    parser.add_argument("--config-out", required=True, type=Path, help="Output path for generated config file.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    config = Config(
        sample_name=args.sample_name,
        read1=args.read1,
        read2=args.read2,
        output_dir=args.output_dir,
        db_path=args.db_path,
        read_length=args.read_length,
        eukdetect_dir=args.eukdetect_dir,
    )
    write_config_to_file(config, args.config_out)
    
    
        
    