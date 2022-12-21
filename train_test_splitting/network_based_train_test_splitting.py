import pandas as pd
import numpy as np
import networkx as nx
from operator import itemgetter
from datetime import datetime
import matplotlib.pyplot as plt

print(datetime.now())


#--generating heterogenous networks for train-test splitting--
    
#-Louvain method-
import community  #Package name is community but refer to python-louvain on pypi
                   # -> to install via conda: conda install -c conda-forge python-louvain
                   # -> to install via pip: pip install python-louvain

def component(G):
    size = [len(i) for i in sorted(nx.connected_components(G),key=len,reverse=True)]
    members = [i for i in sorted(nx.connected_components(G),key=len,reverse=True)]
    return members, size


def partition(family, init_dtp, G, partition_params):
    #split the input graph G into partitions using Louvain community detection method     
    parts = community.best_partition(G, random_state=42)
    values = [parts.get(node) for node in G.nodes()]
    
    #obtain node clusters of initial partitions
    df = pd.DataFrame({"cluster_id":["C{}".format(str(i)) for i in values],"cluster_members":[node for node in G.nodes()]})
    df = df.groupby(by="cluster_id")["cluster_members"].apply(list).reset_index()
    df.insert(loc=1,column="cluster_size",value=[len(i) for i in df["cluster_members"]])    
    df = df.sort_values(by="cluster_size",ascending=False).reset_index(drop=True)
    #df.to_csv(fr"fully_dissimilar_split\{family}\{family}_node_clusters.tsv", sep="\t", index=None) 

    #df = pd.read_csv(fr"fully_dissimilar_split\{family}\{family}_node_clusters.tsv", sep="\t", converters={"cluster_members":lambda x:x.strip("[]").replace("'","").split(", ")})

    #group node clusters with a reasonable size and merge their members to obtain acceptable ratios of train and test samples
    group1 = [node for ind in range(partition_params[0], round(len(df)*partition_params[1]),partition_params[2]) for node in df["cluster_members"][ind]]
    group2 = [node for node in G.nodes() if node not in group1]    
    
    #obtain subgraphs of grouped nodes based on edges in the input graph G
    group1_G = nx.subgraph(G,group1).copy()
    group2_G = nx.subgraph(G,group2).copy()
    
    #obtain the list of edges to be removed to disconnect similar train and test samples
    gr1_gr2_mergeG = nx.compose_all([group1_G,group2_G])
    diffG = nx.difference(G,gr1_gr2_mergeG)
    edge_list = list(diffG.edges())

    #select compounds and datapoints to be removed
    rm_cmp = []
    rm_dtp_c =[]
    rm_dtp_p = []
    for edges in edge_list: 
        if (edges[0].startswith("CHEMBL")) and (edges[1].startswith("CHEMBL")):
                rm_cmp.extend([edges[0],edges[1]])
        
        elif (edges[0].startswith("CHEMBL")) and (not edges[1].startswith("CHEMBL")):
            rm_dtp_c.append(edges[0])
            rm_dtp_p.append(edges[1])
        elif (not edges[0].startswith("CHEMBL")) and (edges[1].startswith("CHEMBL")):
            rm_dtp_c.append(edges[1])
            rm_dtp_p.append(edges[0])
    
    #if train and test samples are connected via bioactivity datapoints, keep their nodes and remove only bioactivity edges between these nodes              
    rm_dtp = pd.DataFrame({"compound_id":rm_dtp_c,"target_id":rm_dtp_p})    
    new_dtp = pd.concat([init_dtp,rm_dtp])
    new_dtp.drop_duplicates(subset=["compound_id","target_id"],keep=False,inplace=True)

    #if train and test samples are connected via two similar compounds, remove both compound nodes 
    new_dtp2 = new_dtp.loc[~new_dtp["compound_id"].isin(rm_cmp)]

    return new_dtp2, rm_cmp


def fully_dissimilar_split(family, partition_params, tr_ts):
    #obtain compound-compound similarity graph
    cmp = pd.read_csv(fr"initial_files/{family}/{family}_compound_tanimoto_sim_0.5thr.tsv", sep="\t")
    G_cmp = nx.from_pandas_edgelist(cmp,source="compound1",target="compound2",edge_attr="similarity")
    
    #obtain protein-protein similarity graph   
    prot = pd.read_csv(fr"initial_files/{family}/{family}_protein-pairwise-sims.tsv", sep="\t")
    G_prot = nx.from_pandas_edgelist(prot,source="target1",target="target2")

    #obtain protein-compound bioactivity graph
    dtp = pd.read_csv(fr"initial_files/{family}/{family}_dataset.tsv", sep="\t")
    G_dtp = nx.from_pandas_edgelist(dtp,source="compound_id",target="target_id")
    
    #merge compound-compound similarity, protein-protein similarity, and protein-compound bioactivity graphs
    G_merged = nx.compose_all([G_cmp,G_prot,G_dtp])
    comp_merged = component(G_merged)
           
    #apply community detection for the largest component of the merged graph, which is the first one 
    subG_merged_comp1 = nx.subgraph(G_merged,list(comp_merged[0][0])).copy() 
    subG_merged_comp1_partition = partition(family, dtp, subG_merged_comp1, partition_params)
    new_dtp = subG_merged_comp1_partition[0]
    
    #update merged graph after dividing the largest component into sub-components
    newG_dtp = nx.from_pandas_edgelist(new_dtp,source="compound_id",target="target_id")
    newG_cmp = nx.subgraph(G_cmp,list(set(G_cmp.nodes())-set(subG_merged_comp1_partition[1]))).copy()

    newG_merged = nx.compose_all([newG_cmp,G_prot,newG_dtp])
    
    #obtain components of newly merged graph
    comp_newG_merged = component(newG_merged)

    #select train/test nodes
    if tr_ts == "default":
        train_final_nodes = comp_newG_merged[0][0] 
        test_final_nodes = [node for node_list in comp_newG_merged[0][1:] for node in node_list]
    else:  
        train_final_nodes = list(comp_newG_merged[0][0]) + [node for node_list in comp_newG_merged[0][tr_ts:] for node in node_list] 
        test_final_nodes = [node for node_list in comp_newG_merged[0][1:tr_ts] for node in node_list]

    #obtain final bioactivity files based on selected train/test nodes
    for nodes,dataset in zip([train_final_nodes,test_final_nodes],["train","test"]): 
        subG_member = nx.subgraph(newG_merged,nodes).copy()
        member_dtp_cmps = []
        member_dtp_prots = []
        for edg in subG_member.edges():
            if (edg[0].startswith("CHEMBL")) and (not edg[1].startswith("CHEMBL")):
                member_dtp_cmps.append(edg[0])
                member_dtp_prots.append(edg[1])
            elif (not edg[0].startswith("CHEMBL")) and (edg[1].startswith("CHEMBL")):
                member_dtp_cmps.append(edg[1])
                member_dtp_prots.append(edg[0]) 
        df_member = pd.DataFrame({"compound_id":member_dtp_cmps,"target_id":member_dtp_prots})
        df_member_dtps = dtp.merge(df_member,on=["compound_id","target_id"])
        
        #this part is only for hydrolases, some of datapoints were removed to balance train-test bioactivity distributions
        if family == "hydrolases":
            if dataset == "train":
                df_member_dtps = df_member_dtps.sample(frac=1, random_state=42)  
                df_train_thr = df_member_dtps.loc[(df_member_dtps["pchembl_value"]>4.2)&(df_member_dtps["pchembl_value"]<5.0)]
                df_member_dtps = pd.concat([df_member_dtps,df_train_thr.iloc[:2400]])  #df_shuffle
                df_member_dtps.drop_duplicates(keep=False,inplace=True)
    
            elif dataset == "test":
                df_member_dtps = df_member_dtps.sample(frac=1, random_state=42)  
                df_test_thr1 = df_member_dtps.loc[(df_member_dtps["pchembl_value"]>6.2)&(df_member_dtps["pchembl_value"]<6.6)]
                df_test_thr2 = df_member_dtps.loc[(df_member_dtps["pchembl_value"]>7.6)&(df_member_dtps["pchembl_value"]<10)]            
                df_member_dtps = pd.concat([df_member_dtps,df_test_thr1.iloc[:300],df_test_thr2.iloc[:300]])  #df_shuffle
                df_member_dtps.drop_duplicates(keep=False,inplace=True)

        df_member_dtps.to_csv(fr"final_files/fully_dissimilar_split/{family}_{dataset}.tsv", sep="\t", index=None)

    print(datetime.now())
    

# fds_epigenetic_regulators = fully_dissimilar_split("epigenetic-regulators", [6,0.2,2], "default")
# fds_hydrolases = fully_dissimilar_split("hydrolases", [0,0.38,1], -40)
# fds_ion_channels = fully_dissimilar_split("ion-channels", [4,0.14,2], "default")
# fds_membrane_receptors = fully_dissimilar_split("membrane-receptors", [0,0.44,1], "default")
# fds_other_enzymes = fully_dissimilar_split("other-enzymes", [5,0.24,2], 43)
# fds_oxidoreductases = fully_dissimilar_split("oxidoreductases", [0,0.44,1], "default")
# fds_proteases = fully_dissimilar_split("proteases", [0,0.39,1], "default")
# fds_transcription_factors = fully_dissimilar_split("transcription-factors", [4,0.1,2], "default")
# fds_transferases = fully_dissimilar_split("transferases", [0,0.36,1], "default")
# fds_transporters = fully_dissimilar_split("transporters", [10,0.28,2], -17)




"To obtain dissimilar_compound and random_splits, you need to obtain fully_dissimilar_splits first!"        


def dissimilar_compound_split(family, params):
    #initial bioactivity data is based on train/test samples of fully_dissimilar_split sets
    tr = pd.read_csv(fr"../datasets/large_scale/fully_dissimilar_split/{family}_train.tsv", sep="\t")
    ts = pd.read_csv(fr"../datasets/large_scale/fully_dissimilar_split/{family}_test.tsv", sep="\t")
    tr_ts = pd.concat([tr,ts])

    cmp = pd.read_csv(fr"initial_files/{family}/{family}_compound_tanimoto_sim_0.5thr.tsv", sep="\t")
    cmp_flt = cmp.loc[(cmp["compound1"].isin(tr_ts["compound_id"].unique()))&(cmp["compound2"].isin(tr_ts["compound_id"].unique()))]
    
    #obtain components of compounds from compound-compound similarity graph
    G_cmp = nx.from_pandas_edgelist(cmp_flt,source="compound1",target="compound2",edge_attr="similarity")
    comp_cmp = component(G_cmp)
    
    #select test compounds
    new_ts_cmp = [cmp for ind in range(params[0],len(comp_cmp[1]), params[1]) for cmp in comp_cmp[0][ind]]
    
    #obtain final bioactivity files based on selected test compounds
    new_ts = tr_ts.loc[tr_ts["compound_id"].isin(new_ts_cmp)]
    new_tr = tr_ts.loc[~(tr_ts["compound_id"].isin(new_ts_cmp))]
        
    new_tr.to_csv(fr"final_files/dissimilar_compound_split/{family}_train.tsv", sep="\t", index=None)
    new_ts.to_csv(fr"final_files/dissimilar_compound_split/{family}_test.tsv", sep="\t", index=None)
    
    print(datetime.now())


# dcs_epigenetic_regulators = dissimilar_compound_split("epigenetic-regulators", [4,5])
# dcs_hydrolases = dissimilar_compound_split("hydrolases", [7,7])
# dcs_ion_channels = dissimilar_compound_split("ion-channels", [6,7])
# dcs_membrane_receptors = dissimilar_compound_split("membrane-receptors", [9,4])
# dcs_other_enzymes = dissimilar_compound_split("other-enzymes", [2,4])
# dcs_oxidoreductases = dissimilar_compound_split("oxidoreductases", [4,9])
# dcs_proteases = dissimilar_compound_split("proteases", [3,5])
# dcs_transcription_factors = dissimilar_compound_split("transcription-factors", [9,9])
# dcs_transferases = dissimilar_compound_split("transferases", [6,5])
# dcs_transporters = dissimilar_compound_split("transporters", [8,5])

    




def random_split(family, rand_st):
    #initial bioactivity data is based on train/test samples of fully_dissimilar_split sets
    tr = pd.read_csv(fr"../datasets/large_scale/fully_dissimilar_split/{family}_train.tsv", sep="\t")
    ts = pd.read_csv(fr"../datasets/large_scale/fully_dissimilar_split/{family}_test.tsv", sep="\t")
    tr_ts = pd.concat([tr,ts])
    print(len(tr),len(ts),len(tr_ts),len(tr)+len(ts))
    
    #obtain final bioactivity files using random splitting 
    new_ts = tr_ts.sample(n=len(ts), random_state=rand_st)    
    new_tr = pd.concat([tr_ts,new_ts]).drop_duplicates(keep=False)
        
    new_tr.to_csv(fr"final_files/random_split/{family}_train.tsv", sep="\t", index=None)
    new_ts.to_csv(fr"final_files/random_split/{family}_test.tsv", sep="\t", index=None)

    print(datetime.now())


# rs_epigenetic_regulators = random_split("epigenetic-regulators", 1)
# rs_hydrolases = random_split("hydrolases", 5)
# rs_ion_channels = random_split("ion-channels", 1)
# rs_membrane_receptors = random_split("membrane-receptors", 2)
# rs_other_enzymes = random_split("other-enzymes", 2)
# rs_oxidoreductases = random_split("oxidoreductases", 1)
# rs_proteases = random_split("proteases", 1)
# rs_transcription_factors = random_split("transcription-factors", 3)
# rs_transferases = random_split("transferases", 1)
# rs_transporters = random_split("transporters", 5)

 





















