# Validation-and-quality-assessment-of-macromolecular-structures-using-complex-network-analysis

Protein model analysis 
The R script (File name: ProtModAna.r) was used to analyze protein models. The script reads PDB file, converts 3D protein model into graph and calculates mean node degree, Zscore, poor and long subgraphs.

Smooth mean and standard deviation was calculated using data source (raw_data.csv), with a window length [0.8N, 1.2N], where N is the length of protein (File name: SmoothData.csv).  

USAGE
1) Test run

Save "ProtModAna.r" and "SmoothData.csv" in to the same directory and run R script "ProtModAna.r".

2) Custom analyze

User can change the path to the "SmoothData.csv" data file. Note, this file is necessary to calculate Zscore.
To validate protein model change the line where R script reads pdb file (see line 19 and 20 in "ProtModAna.r")

