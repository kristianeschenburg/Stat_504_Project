# Stats 504 Project
Code and data for Statistics 504 Project

We'll be making inference on the relationship of various clinical covariates with Alzheimer's disease.  The data is provided by the Allen Brain Institute, based on their "Aging, Dementia, and TBI" study, which is a sub-project of the broader "Adult Changes in Thought" (ACT) study.  The ACT study is organized by the University of Washington and the Kaiser Permanente Health Research Institute for the longitudinal analysis of the relationship between brain aging and incident dementia in the Seattle metropolitan area.

Included in this Allen Institute dataset are study participant metadata, quantitative histology and protein measurements of neuropathology, and RNA-sequencing (seq) analyses of the hippocampus and neocortex.

Preprocessing by the Allen Institute included between-sample normalization of RNA-seq counts within a given region (hippocampus, parietal cortex, etc.) -- both the unnormalized and normalized data are available, but we use the normalized data for our analysis.  Due to the large file size of the RNA-seq data, we do not include the table here -- but it can be downloaded from here: http://aging.brain-map.org/api/v2/well_known_file_download/502999992 (note, this will download the .zip file directly).  More information on the the Aging, TBI, and Dementia Study can be found here: http://aging.brain-map.org/overview/home.

Detailed documentation on the actual data used in this analysis can be found here: http://aging.brain-map.org/download/index
