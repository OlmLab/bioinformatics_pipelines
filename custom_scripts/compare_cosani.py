import click
import inStrain.SNVprofile as SNVprofile
import pandas as pd
import inStrain
import numpy as np
import json


def get_cos_ani(f1:np.ndarray,f2:np.ndarray)->float:
    """
    Calculate the cosine ANI (Average Nucleotide Identity) between two samples.
    
    Args:
        f1: A numpy array of the frequency tensor of sample 1. The shape of the frequency tensor is (4, n), where n is the number of nucleotide positions.
        f2: A numpy array of the frequency tensor of sample 2.  The shape of the frequency tensor is (4, n), where n is the number of nucleotide positions.
    
    Returns:
        A float value of the cosine ANI.
    """

    
    
    return np.nan_to_num(np.sum(np.multiply(f1,f2))/(np.sqrt(np.sum(np.square(f1)))*np.sqrt(np.sum(np.square(f2))))*100)

def get_popani(f1:np.ndarray,f2:np.ndarray)->float:
    """
    Calculate the population ANI (Average Nucleotide Identity) between two samples.
    
    Args:
        f1: A numpy array of the frequency tensor of sample 1. The shape of the frequency tensor is (4, n), where n is the number of nucleotide positions.
        f2: A numpy array of the frequency tensor of sample 2.  The shape of the frequency tensor is (4, n), where n is the number of nucleotide positions.
        
    Returns:
        A float value of the population ANI.
    """
    if f1.shape!=f2.shape:
        raise ValueError("shape mismatch")
    np.any(np.multiply(f1,f2),axis=0)
    return np.sum(np.any(np.multiply(f1,f2),axis=0))/f1.shape[1]*100

def get_conani(f1:np.ndarray,f2:np.ndarray)->float:
    """
    Calculate the consensus ANI (Average Nucleotide Identity) between two samples.
    
    Args:
        f1: A numpy array of the frequency tensor of sample 1. The shape of the frequency tensor is (4, n), where n is the number of nucleotide positions.
        f2: A numpy array of the frequency tensor of sample 2.  The shape of the frequency tensor is (4, n), where n is the number of nucleotide positions.
        
    Returns:
        A float value of the consensus ANI.
    """
    if np.all(np.logical_or(f1==0,f2==0)):
        return np.zeros(f1.shape[1])
    if f1.shape!=f2.shape:
        raise ValueError("shape mismatch")
    return np.sum(np.argmax(f1,axis=0)==np.argmax(f2,axis=0))/f1.shape[1]*100


def compare(
    profile_1:SNVprofile.SNVprofile,
    profile_2:SNVprofile.SNVprofile,
    stb_file: str,
    cov_threshold: float = 10,
    breadth_threshold: float = 0.5
    
            )->dict:
    """
    Compare two SNV profiles and calculate the cosine ANI.
    
    Args:
        profile_1: The first SNV profile.
        profile_2: The second SNV profile. 
        stb_file: The path to the STB file that contains the binning information.  
        cov_threshold: The coverage threshold to call a genome as present.
        breadth_threshold: The breadth threshold to call a genome as present.
    Returns:
        None
    """

    sample_1_genomes=set(profile_1.get("genome_level_info")["genome"].to_list())
    sample_2_genomes=set(profile_2.get("genome_level_info")["genome"].to_list())
    union_genomes=sample_1_genomes.union(sample_2_genomes)
    covT_1 = profile_1.get("covT")
    covT_2 = profile_2.get("covT")
    sc_stat_1 = profile_1.get_nonredundant_scaffold_table()
    sc_stat_1.set_index("scaffold", inplace=True)
    sc_stat_2 = profile_2.get_nonredundant_scaffold_table()
    sc_stat_2.set_index("scaffold", inplace=True)
    b2s=pd.read_table(stb_file,header=None).groupby(1)[0].agg(list).reset_index()
    b2s=dict(zip(b2s[1],b2s[0]))
    output={
        profile_1.location.split("/")[-1]:{"genomes":list(sample_1_genomes)},
        profile_2.location.split("/")[-1]:{"genomes":list(sample_2_genomes)},
        "shared_genomes":{},
    }
    for bin in union_genomes:
        if bin not in output[profile_1.location.split("/")[-1]]["genomes"] or bin not in output[profile_2.location.split("/")[-1]]["genomes"]:
            continue
        scratch={"cosani":[],"popani":[],"conani":[], "len":[]}
        for scaffold in b2s[bin]:
            if not (covT_1.get(scaffold) and covT_2.get(scaffold)):
                continue
            sc1=inStrain.profile.profile_utilities.mm_counts_to_counts_shrunk(covT_1.get(scaffold))
            sc2=inStrain.profile.profile_utilities.mm_counts_to_counts_shrunk(covT_2.get(scaffold))
            common_index=sc1.index.intersection(sc2.index)
            if len(common_index) == 0:
                continue 
            sc1=sc1[common_index]
            sc2=sc2[common_index]
            vect1=np.zeros((4,len(common_index)),dtype=np.float64)
            vect1[0,:]=1
            vect2=np.zeros((4,len(common_index)),dtype=np.float64)
            vect2[0,:]=1
            raw_snp_table_1=profile_1.get("raw_snp_table")
            sc1_snps=raw_snp_table_1[raw_snp_table_1["scaffold"]==scaffold]
            sc1_snpprof=sc1_snps[["A","T","C","G"]].to_numpy().T
            sc1_snp_index=sc1.index.get_indexer(sc1_snps["position"])
            vect1[:,sc1_snp_index]=sc1_snpprof
            ### Normalize the vectors
            # Normalize the frequency tensors to sum to 1
            vect1= vect1 / np.sum(vect1, axis=0, keepdims=True)
            raw_snp_table_2=profile_2.get("raw_snp_table")
            sc2_snps=raw_snp_table_2[raw_snp_table_2["scaffold"]==scaffold]
            sc2_snpprof=sc2_snps[["A","T","C","G"]].to_numpy().T
            sc2_snp_index=sc2.index.get_indexer(sc2_snps["position"])
            vect2[:,sc2_snp_index]=sc2_snpprof
            # Normalize the frequency tensors to sum to 1
            vect2= vect2 / np.sum(vect2, axis=0, keepdims=True)
            scratch["cosani"].append(get_cos_ani(vect1,vect2))
            scratch["popani"].append(get_popani(vect1,vect2))
            scratch["conani"].append(get_conani(vect1,vect2))
            scratch["len"].append(len(common_index))
        output["shared_genomes"][bin]={"cosani":np.sum([scratch["cosani"][i]*scratch["len"][i] for i in range(len(scratch["len"]))])/np.sum(scratch["len"]),
                                        "popani":np.sum([scratch["popani"][i]*scratch["len"][i] for i in range(len(scratch["len"]))])/np.sum(scratch["len"]),
                                        "conani":np.sum([scratch["conani"][i]*scratch["len"][i] for i in range(len(scratch["len"]))])/np.sum(scratch["len"]),
                                        }
    
    
    return output
        
def compare_stats(compare_profile:dict,
                      strain_pop_treshold:float,
                      strain_cos_treshold:float,
                      strain_con_treshold:float,
                      )->dict:
    """
    Get the comparison statistics from the comparison profile.
    
    Args:
        compare_profile: The comparison profile from compare function.
        strain_treshold: The strain threshold.
    Returns:
        None
    """
    shared_genomes={
        "popani":list(),
        "conani":list(),
        "cosani":list(),
    }
    for genome in compare_profile["shared_genomes"]:
        if compare_profile["shared_genomes"][genome]["cosani"]>strain_cos_treshold:
            shared_genomes["cosani"].append(genome)
        if compare_profile["shared_genomes"][genome]["popani"]>strain_pop_treshold:
            shared_genomes["popani"].append(genome)
        if compare_profile["shared_genomes"][genome]["conani"]>strain_con_treshold:
            shared_genomes["conani"].append(genome)
            
        
    output={
        profile:{"all_genomes":len(compare_profile[profile]["genomes"]),
                 "shared_genomes":{"popani":len(shared_genomes["popani"]),
                                    "conani":len(shared_genomes["conani"]),
                                    "cosani":len(shared_genomes["cosani"])},
                 "sharing_rate":{"popani":len(shared_genomes["popani"])/len(compare_profile[profile]["genomes"]),
                                 "conani":len(shared_genomes["conani"])/len(compare_profile[profile]["genomes"]),
                                 "cosani":len(shared_genomes["cosani"])/len(compare_profile[profile]["genomes"])},
                 }
        for profile in compare_profile if profile != "shared_genomes"
    }
        
        
    return output
    
@click.group()
def cli():
    """Compare SNV profiles and analyze strain sharing statistics."""
    pass

@cli.command('compare')
@click.option('--profile_1', type=click.Path(exists=True), required=True, help='Path to the first SNV profile.')
@click.option('--profile_2', type=click.Path(exists=True), required=True, help='Path to the second SNV profile.')
@click.option('--stb_file', type=click.Path(exists=True), required=True, help='Path to the STB file that contains the binning information.')
@click.option('--output_file', type=click.Path(), required=True, help='Path to save the comparison results as a JSON file.')
def compare_command(profile_1: str, profile_2: str, stb_file: str, output_file: str):
    """Compare two SNV profiles and calculate ANI metrics."""
    profile_1 = SNVprofile.SNVprofile(profile_1)
    profile_2 = SNVprofile.SNVprofile(profile_2)
    
    res = compare(
        profile_1=profile_1,
        profile_2=profile_2,
        stb_file=stb_file
    )
    
    with open(output_file, "w") as f:
        json.dump(res, f, indent=4)
    
    click.echo(f"Comparison results saved to {output_file}")

@cli.command('stats')
@click.option('--compare_profile', type=click.Path(exists=True), required=True, help='Path to the comparison profile JSON file.')
@click.option('--strain_pop_treshold', type=float, required=False, help='The population ANI threshold to call a genome as shared.', default=99)
@click.option('--strain_cos_treshold', type=float, required=False, help='The cosine ANI threshold to call a genome as shared.',default=99)
@click.option('--strain_con_treshold', type=float, required=False, help='The consensus ANI threshold to call a genome as shared.',default=99)
@click.option('--output_file', type=click.Path(), help='Path to save the stats results (defaults to input filename with _stats suffix).')
def stats_command(compare_profile: str, strain_pop_treshold: float, strain_cos_treshold: float, strain_con_treshold: float, output_file: str = None):
    """Generate statistics from a comparison profile."""
    with open(compare_profile, "r") as f:
        profile_data = json.load(f)
    
    res = compare_stats(
        compare_profile=profile_data,
        strain_pop_treshold=strain_pop_treshold,
        strain_cos_treshold=strain_cos_treshold,
        strain_con_treshold=strain_con_treshold
    )
    
    if output_file is None:
        output_file = compare_profile.replace(".json", "_stats.json")
    
    with open(output_file, "w") as f:
        json.dump(res, f, indent=4)
    
    click.echo(f"Statistics results saved to {output_file}")

if __name__ == "__main__":
    cli()