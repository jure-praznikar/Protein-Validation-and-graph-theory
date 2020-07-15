Validation and quality assessment of macromolecular structures using complex network analysis

Validation of three-dimensional structures is at the core of structural determination methods. The local validation criteria, such as deviations from ideal bond length and bonding angles, Ramachandran plot outliers and clashing contacts, are a standard part of structure analysis before structure deposition, whereas the global and regional packing may not yet have been addressed. In the last two decades, three-dimensional models of macromolecules such as proteins have been successfully described by a network of nodes and edges. Amino acid residues as nodes and close contact between the residues as edges have been used to explore basic network properties, to study protein folding and stability and to predict catalytic sites. Using complex network analysis, we introduced common network parameters to distinguish between correct and incorrect three-dimensional protein structures. The analysis showed that correct structures have a higher average node degree, higher graph energy, and lower shortest path length than their incorrect counterparts. Thus, correct protein models are more densely intra-connected, and in turn, the transfer of information between nodes/amino acids is more efficient. Moreover, protein graph spectra were used to investigate model bias in protein structure.


For details see paper:
https://www.nature.com/articles/s41598-019-38658-9

Protein model analysis 
The R script (File name: ProtModAna.r) was used to analyze protein models. The script reads PDB file, converts 3D protein model into graph and calculates mean node degree, Zscore, poor and long subgraphs.

Smooth mean and standard deviation was calculated using data source (raw_data.csv), with a window length [0.8N, 1.2N], where N is the length of protein (File name: SmoothData.csv).  

USAGE
1) Test run

Save "ProtModAna.r" and "SmoothData.csv" in to the same directory and run R script "ProtModAna.r".

2) Custom analyze

User can change the path to the "SmoothData.csv" data file. Note, this file is necessary to calculate Zscore.
To validate protein model change the line where R script reads pdb file (see line 19 and 20 in "ProtModAna.r")

