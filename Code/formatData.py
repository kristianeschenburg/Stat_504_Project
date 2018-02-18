#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 11:39:13 2018

@author: kristianeschenburg
"""

import numpy as np
import pandas as pd

def region_table(gene_data,s2r_map,region):
    
    """
    Given a table of (Gene IDs x RNA IDs), subset gene data to only RNA IDs
    corresponding to a specific region.
    
    Parameters:
    - - - - -
        gene_data : rna expression data of dims (Gene IDs x RNA IDs)
        s2r_map : subj_id to RNA ID maps
        region : unique region acronym name
    """
    
    # get unique Gene IDs
    gene_ids = list(gene_data['gene_id \\ rnaseq_profile_id'])
    
    # initialize empty data array
    # not as intuitive, but faster
    dataArray = np.empty(shape=(len(s2r_map.keys()),len(gene_ids)))
    
    # loop over subject IDs
    for num,subj in enumerate(s2r_map.keys()):

        # get unique subject-specific RNA ID for region
        rna_id = s2r_map[subj][region]
        
        # if subject had this region sampled, get expression data
        if rna_id:
            rna_data = gene_data[rna_id]
        # otherwise create vector of all NaN
        else:
            rna_data = np.zeros(shape=len(gene_ids),)
            rna_data.fill(np.nan)
        
        # fill array with expression data
        dataArray[num,:] = rna_data

    # initialize / fill data frame
    G = pd.DataFrame(data=dataArray,columns=gene_ids)
    # allow indexing of rows by subject ID
    # ex. G.loc[subject_ID]
    G.index = s2r_map.keys()
    
    return G


if __name__ == '__main__':
    
    """
    Assumes you have a directory structure like this:
        
        YourMainDir/
            Code/
                formatData.py
            Data/
                gene_expression_matrix/
                    columns-samples.csv
                    rows-genes.csv
                    fpkm_table_normalized.csv
    """
    
    dataDir = '../Data/gene_expression_matrix/'
    
    # Load file mapping RNA-IDs to subject donor IDs
    col_file = ''.join([dataDir,'columns-samples.csv'])
    col_names = pd.read_csv(col_file)
    
    # We're interested in rnaseq_profile_id, donor_id, and structure acronym
    print col_names.head()
    
    # Get lists of donor IDs and unique region names
    subj_id = list(set(col_names['donor_id']))
    regions = set(col_names['structure_acronym'])
    
    # Map donor ID to RNA IDs
    # {subject_ID: {Structure_Name: RNA_ID_#}}
    subj_region_map = {k: {}.fromkeys(regions) for k in subj_id}
    
    for index,row in col_names.iterrows():
        
        subject = row['donor_id']
        reg = row['structure_acronym']
        rna_id = row['rnaseq_profile_id']
        
        subj_region_map[subject][reg] = str(rna_id)
        
    # Load file mapping gene IDs to gene description
    row_file = ''.join([dataDir,'rows-genes.csv'])
    row_names = pd.read_csv(row_file)
    
    # Load file of (Gene ID x RNA ID) expression values
    gene_file = ''.join([dataDir,'fpkm_table_normalized.csv'])
    gene_data = pd.read_csv(gene_file)
    
    # fo each region, generate unique data frame
    fwm = region_table(gene_data,subj_region_map,'FWM')
    hip = region_table(gene_data,subj_region_map,'HIP')
    pcx = region_table(gene_data,subj_region_map,'PCx')
    tcx = region_table(gene_data,subj_region_map,'TCx')
    
    # for each region, save data frame to csv
    fwm.to_csv(''.join([dataDir,'GeneData_FWM.csv']))
    hip.to_csv(''.join([dataDir,'GeneData_HIP.csv']))
    pcx.to_csv(''.join([dataDir,'GeneData_PCX.csv']))
    tcx.to_csv(''.join([dataDir,'GeneData_TCX.csv']))
