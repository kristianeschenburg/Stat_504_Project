#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 11:39:13 2018

@author: kristianeschenburg
"""

import numpy as np
import os
import pandas as pd

"""
CD into the directory where the ./gene_expression_matrix/ directort is.
"""

col_file = './gene_expression_matrix_2016-03-03/columns-samples.csv'
col_names = pd.read_csv(col_file)

print col_names.head()

subj_id = list(set(col_names['donor_id']))
regions = set(col_names['structure_acronym'])

subj_region_map = {k: {}.fromkeys(regions) for k in subj_id}

for index,row in col_names.iterrows():
    
    subject = row['donor_id']
    reg = row['structure_acronym']
    rna_id = row['rnaseq_profile_id']
    
    subj_region_map[subject][reg] = str(rna_id)
    
row_file = './gene_expression_matrix_2016-03-03/rows-genes.csv'
row_names = pd.read_csv(row_file)
gene_file = './gene_expression_matrix_2016-03-03/fpkm_table_normalized.csv'
gene_data = pd.read_csv(gene_file)

fwm = region_table(gene_data,subj_region_map,'FWM')
hip = region_table(gene_data,subj_region_map,'HIP')
pcx = region_table(gene_data,subj_region_map,'PCx')
tcx = region_table(gene_data,subj_region_map,'TCx')

"""

fwm.to_csv('./gene_expression_matrix_2016-03-03/GeneData_FWM.csv)
hip.to_csv('./gene_expression_matrix_2016-03-03/GeneData_HIP.csv)
pcx.to_csv('./gene_expression_matrix_2016-03-03/GeneData_PCX.csv)
tcx.to_csv('./gene_expression_matrix_2016-03-03/GeneData_TCX.csv)
"""

def region_table(gene_data,s2r_map,region):
    
    """
    Given table of gene expression data, subset table by region.
    """
    
    # get unique Gene ID
    gene_ids = list(gene_data['gene_id \\ rnaseq_profile_id'])
    
    dataArray = np.empty(shape=(len(s2r_map.keys()),len(gene_ids)))
    
    # gene subject ID
    for num,subj in enumerate(s2r_map.keys()):

        # gene subject-region RNA ID
        rna_id = s2r_map[subj][region]
        
        # if subject had this region sampled, get data
        if rna_id:
            rna_data = gene_data[rna_id]
        # otherwise create vector of all zeros
        else:
            rna_data = np.zeros(shape=len(gene_ids),)
            rna_data.fill(np.nan)
        
        dataArray[num,:] = rna_data

    G = pd.DataFrame(data=dataArray,columns=gene_ids)
    G.index = s2r_map.keys()
    
    return G