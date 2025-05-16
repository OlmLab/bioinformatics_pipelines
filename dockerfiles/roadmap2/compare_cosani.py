#!/opt/conda/envs/inStrain/bin/python3
import click
import inStrain.SNVprofile as SNVprofile
import pandas as pd
import inStrain
import numpy as np
import json
from multiprocessing import Pool
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
    null_model:np.ndarray = None,
    min_freq: float = 0.05,
    n_procs: int = 4,
    cov_threshold: float = 0,
    breadth_threshold: float = 0,   
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
    profile_1_genomes=profile_1.get("genome_level_info")
    profile_1_genomes=profile_1_genomes[(profile_1_genomes["coverage"]>cov_threshold)&(profile_1_genomes["breadth"]>breadth_threshold)]
    sample_1_genomes=set(profile_1_genomes["genome"].to_list())
    del profile_1_genomes
    profile_2_genomes=profile_2.get("genome_level_info")
    profile_2_genomes=profile_2_genomes[(profile_2_genomes["coverage"]>cov_threshold)&(profile_2_genomes["breadth"]>breadth_threshold)]
    sample_2_genomes=set(profile_2_genomes["genome"].to_list())
    del profile_2_genomes
    union_genomes=sample_1_genomes.union(sample_2_genomes)
    covT_1 = profile_1.get("covT")
    covT_2 = profile_2.get("covT")
    b2s=pd.read_table(stb_file,header=None).groupby(1)[0].agg(list).reset_index()
    b2s=dict(zip(b2s[1],b2s[0]))
    output={
        profile_1.location.split("/")[-1]:{"genomes":list(sample_1_genomes)},
        profile_2.location.split("/")[-1]:{"genomes":list(sample_2_genomes)},
        "shared_genomes":{},
    }
    raw_snp_table_1=profile_1.get("raw_snp_table")
    raw_snp_table_2=profile_2.get("raw_snp_table")
    args={
        "scaffold":[],
        "bin_":[],
        "scaff_covT_1":[],
        "scaff_covT_2":[],
        "raw_snp_table_1":[],
        "raw_snp_table_2":[],
        "min_freq":[],
        "null_model":[],
    }
    for bin_ in union_genomes:
        for scaffold in b2s[bin_]:
            if not (covT_1.get(scaffold) and covT_2.get(scaffold)):
                continue
            args["scaffold"].append(scaffold)
            args["bin_"].append(bin_)
            args["scaff_covT_1"].append(covT_1.get(scaffold))
            args["scaff_covT_2"].append(covT_2.get(scaffold))
            args["raw_snp_table_1"].append(raw_snp_table_1)
            args["raw_snp_table_2"].append(raw_snp_table_2)
            args["min_freq"].append(min_freq)
            args["null_model"].append(null_model)
    

    with Pool(n_procs) as pool:
        res=pool.starmap(
            compare_scaffold,
            zip(
            args["scaffold"],
            args["bin_"],
            args["scaff_covT_1"],
            args["scaff_covT_2"],
            args["raw_snp_table_1"],
            args["raw_snp_table_2"],
            args["min_freq"],
            args["null_model"],
            )
        )
    output={
        "shared_genomes":{},
        profile_1.location.split("/")[-1]:{"genomes":list(sample_1_genomes)},
        profile_2.location.split("/")[-1]:{"genomes":list(sample_2_genomes)},
    }
    df=pd.DataFrame(res,columns=["bin","scaffold","common_index","cosani","popani","conani"])
    df["weghted_sum_cosani"]=df["cosani"]*df["common_index"]
    df["weghted_sum_popani"]=df["popani"]*df["common_index"]
    df["weghted_sum_conani"]=df["conani"]*df["common_index"]
    df=df.groupby("bin").agg({"weghted_sum_cosani":"sum","weghted_sum_popani":"sum","weghted_sum_conani":"sum","common_index":"sum"}).reset_index()
    df["cosani"]=df["weghted_sum_cosani"]/df["common_index"]
    df["popani"]=df["weghted_sum_popani"]/df["common_index"]
    df["conani"]=df["weghted_sum_conani"]/df["common_index"]
    df=df[["bin","cosani","popani","conani"]]
    output["shared_genomes"]=df.set_index("bin").to_dict(orient="index")
    return output
def compare_scaffold(scaffold:str,
                     bin_:str,
                        scaff_covT_1,
                        scaff_covT_2,
                        raw_snp_table_1:pd.DataFrame,
                        raw_snp_table_2:pd.DataFrame,
                        min_freq:float = 0.01,
                        null_model:np.ndarray = None,
                        )->tuple:
    
    common_index=inStrain.readComparer.calc_mm2overlap_v3(scaff_covT_1,scaff_covT_2)[0][0]
    if len(common_index) == 0:
        return bin_,scaffold,0,0,0,0 
    common_index_list=list(common_index)
    common_index_map=dict(zip(common_index_list,range(len(common_index_list))))
    vect1=np.zeros((4,len(common_index)),dtype=np.float64)
    vect1[0,:]=1
    sc1_snps=raw_snp_table_1[(raw_snp_table_1["scaffold"]==scaffold) & (raw_snp_table_1["position"].isin(common_index))]
    sc2_snps=raw_snp_table_2[(raw_snp_table_2["scaffold"]==scaffold) & (raw_snp_table_2["position"].isin(common_index))]
    allpos=pd.concat([sc1_snps,sc2_snps]).drop_duplicates(subset=["position"], keep="first")[["position","ref_base"]]
    allpos=pd.concat([allpos["position"],pd.get_dummies(allpos['ref_base'],dtype=int).reindex(columns=['A', 'T', 'C', 'G'], fill_value=0)],axis=1)
    vect1[:,allpos["position"].map(common_index_map).to_numpy()]=allpos[["A","T","C","G"]].to_numpy().T
    vect2=vect1.copy() ### Up to here, the vect1 and vect2 are the same
    tres_1=get_treshold_matrix(sc1_snps[["A","T","C","G"]].to_numpy(),null_model,min_freq=min_freq).T
    sc1_snpprof=sc1_snps[["A","T","C","G"]].to_numpy().T
    sc1_snpprof=sc1_snpprof/np.sum(sc1_snpprof, axis=0, keepdims=True)
    sc1_snpprof[sc1_snpprof<tres_1]=0
    sc1_snpprof=sc1_snpprof/np.sum(sc1_snpprof, axis=0, keepdims=True)
    vect1[:, sc1_snps["position"].map(common_index_map).to_numpy()]=sc1_snpprof
    tres_2=get_treshold_matrix(sc2_snps[["A","T","C","G"]].to_numpy(),null_model,min_freq=min_freq).T
    sc2_snpprof=sc2_snps[["A","T","C","G"]].to_numpy().T
    sc2_snpprof=sc2_snpprof/np.sum(sc2_snpprof, axis=0, keepdims=True)
    sc2_snpprof[sc2_snpprof<tres_2]=0
    sc2_snpprof=sc2_snpprof/np.sum(sc2_snpprof, axis=0, keepdims=True)
    vect2[:,sc2_snps["position"].map(common_index_map).to_numpy()]=sc2_snpprof
    ans=bin_,scaffold,len(common_index), get_cos_ani(vect1,vect2),get_popani(vect1,vect2),get_conani(vect1,vect2)
    return ans


def get_treshold_matrix(
        counts_matrix:np.ndarray,
        null_model:np.ndarray,
        min_freq:float = 0.01,
                        )->np.ndarray:
    """
    Get the treshold matrix from the counts matrix.
    
    Args:
        counts_matrix: The counts matrix.
        null_model: The null model.
        min_freq: The minimum frequency to call a base present.
    Returns:
        A numpy array of the treshold matrix.
    """
    coverage=np.sum(counts_matrix,axis=1)
    tres=null_model[coverage,:]
    row_indices = np.repeat(np.arange(counts_matrix.shape[0]), counts_matrix.shape[1])
    flat_counts = counts_matrix.flatten()
    flat_counts = np.clip(flat_counts, 0, null_model.shape[1]-1)
    result = np.zeros_like(counts_matrix, dtype=float)
    result.flat = tres[row_indices, flat_counts]
    tres = result
    tres[tres<min_freq]=min_freq
    return tres
    
      
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
@click.option('--cpus', type=int, default=8, help='Number of CPU cores to use for parallel processing.')
@click.option('--min_freq', type=float, default=0.05, help='Minimum frequency to call a base present.')
@click.option('--cov_threshold', type=float, default=0, help='Coverage threshold to call a genome as present.')
@click.option('--breadth_threshold', type=float, default=0, help='Breadth threshold to call a genome as present.')
@click.option('--null_model', type=click.Path(exists=True),default="/data/NullModel.txt" ,required=False, help='Path to the null model file.')
def compare_command(profile_1: str, profile_2: str, stb_file: str, output_file: str, cpus: int, min_freq: float, cov_threshold: float, breadth_threshold: float, null_model: str):
    """Compare two SNV profiles and calculate ANI metrics."""
    profile_1 = SNVprofile.SNVprofile(profile_1)
    profile_2 = SNVprofile.SNVprofile(profile_2)
    null_model = pd.read_table(null_model,delimiter="\t")
    null_model=np.vstack([np.zeros((1,19)),null_model.to_numpy()])
    null_model[:,0]=0
    res = compare(
        profile_1=profile_1,
        profile_2=profile_2,
        stb_file=stb_file,
        null_model=null_model,  
        n_procs=cpus,
        min_freq=min_freq,
        cov_threshold=cov_threshold,
        breadth_threshold=breadth_threshold,
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
