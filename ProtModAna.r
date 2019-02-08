#    This is free script: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    This script is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this script.  If not, see <http://www.gnu.org/licenses/>.

#R script for network analysis of protein models
#load lib
library(igraph)
library(bio3d)
#read PDB file
#pdb <- read.pdb("177l.pdb")
pdb <- read.pdb("179l.pdb")
inds <- atom.select(pdb, elety = "CA") # CA atoms
my.atoms<-pdb$atom[inds$atom,c("x", "y", "z")] # CA atom coordinates
numatoms<-dim(my.atoms)[1]
#adj matrix
adj <- matrix(nrow = numatoms, ncol = numatoms,0) # adj matrix
adjpoor<- matrix(nrow = numatoms, ncol = numatoms,0) # matrix-poor CA distances
adjlong<- matrix(nrow = numatoms, ncol = numatoms,0) # matrix-long CA distances
# define function,  Euclidean distance of two vectors
euc.dist <- function(x1,x2) sqrt(sum((x1-x2)^2)) 
#set threshold's
threshold<-7.0 #Upper threshold=7.0
thresholdLONG<-3.9 # threshold/long CA distances
thresholdPOORup<-3.7 # Upper threshold-poor CA distances
thresholdPOORdown<-3.0 # Down threshold-poor CA distances
#init
ca.poor<-0
ca.long<-0
#
cat("\n")
start_time <- Sys.time()
cat(" .... Calculation in progress ....\n")
# Create adj matrix - only upper triangle (no diagonal elements)
adj <- matrix(0,nrow = numatoms, ncol = numatoms)
DD<-dist(my.atoms)
adj[lower.tri(adj, diag=FALSE)] <- DD #write distances to adj 
adj<-adj+t(adj) # add upper triangle, is symmetric
adjpoor<-adj # write distances to adjpoor
templong<-adj[1:numatoms-1,2:numatoms] #use original adj matrix
#A  d  A  A      d  A  A
#A  A  d  A      A  d  A
#A  A  A  d      A  A  d
#A  A  A  A
#only "d" elements represent sequential CA atoms
#in smaller matrix are sequential CA atoms ("d") represented in diagonal
#adjlong<-adj # write distances to adjlong
#set values to zero or one if below threshold
adj[adj>threshold]<-0 #above threshold=0
adj[adj>1.0E-5]<-1 #else=1
#
adjpoor[adjpoor > thresholdPOORup | adjpoor < thresholdPOORdown]<-0
adjpoor[adjpoor > 1.0E-5]<-1
ca.poor<-length(which(adjpoor==1))/2 #divide by two, it is symmetric
#long
templong[upper.tri(templong,diag=FALSE) | lower.tri(templong,diag=FALSE)]<-0
diag(templong)[which(diag(templong)<3.9)]<-0
diag(templong)[which(diag(templong)>1.0E-5)]<-1
adjlong[1:numatoms-1,2:numatoms]<-templong
adjlong<-adjlong+t(adjlong)# is symmetric
ca.long<-length(which(adjlong==1))/2 #divide by two, it is symmetric
#
end_time <- Sys.time()
cat("Time (sec) to construct adj matrix: ", round(end_time-start_time,5),"\n")
#
#make graph object from adj matrix
g=graph.adjacency(adj,mode="undirected",weighted=NULL) #graph
gp=graph.adjacency(adjpoor,mode="undirected",weighted=NULL) # poor CA graph
gl=graph.adjacency(adjlong,mode="undirected",weighted=NULL) # long CA graph
#Global parameters
Mnd<-mean(degree(g)) # Mean node degree
mydata<-read.csv("./SmoothData.csv",header=TRUE) # read data from file
index<-which(mydata[,1]==vcount(g))
degreeDB<-mydata[index,2] # mean degree from Data Base, read second column
stdDB<-mydata[index,3] # stdev from Data Base, read third column
Zscore<-(Mnd-degreeDB)/stdDB # calculate Z-score
#Print on Screen
cat("\n")
cat(" ----- Global parameters -----\n")
cat("Number of nodes (CA atoms):",format(vcount(g),digits=5),"\n")
cat("Calculated Node degree:",format(Mnd,,digits=3),"\n")
cat("Node degree from data base:",format(degreeDB,,digits=3),"\n")
cat("Standard deviation from data base:",format(stdDB,,digits=3),"\n")
cat("Z-score:",format(Zscore,digits=2),"\n")
#Local parameters
clup <- components(gp) # poor connected components
clul <- components(gl) # long connected components
percp<-ca.poor/(vcount(g)-1)*100 # calculate percentage
percl<-ca.long/(vcount(g)-1)*100
cat("\n")
cat("--- Count poor and long CA distances ---\n")
cat("Poor CA edges (total/percentage): ",ca.poor,"/",format(percp,digits=3),"%\n")
cat("Long CA edges (total/percentage): ",ca.long,"/",format(percl,digits=3),"%\n")
# POOR -  Print on Screen
cat("\n")
cat("--- Poor subgraphs ( Number of nodes >= 4 ) ---\n")
cat(" Max. subgraph has ",max(clup$csize)," nodes.\n")
cat(" List of subgraphs:\n")
k<-which(clup$csize >= 4)
if (length(k)>0){
   for (i in 1:length(k)) {
     ak<-which(clup$membership == k[i])
     cat(i,")",ak,"\n") 
   }
   cat("! Note, first residue in input file has index=1 !\n")
} else {
  cat("! List of POOR subgraphs is empty !\n")
}
# LONG -  Print on Screen
cat("\n")
cat("--- Long subgraphs ( Number of nodes >= 4 ) ---\n")
cat(" Max. subgraph has ",max(clul$csize)," nodes.\n")
cat(" List of subgraphs:\n")
k<-which(clul$csize >= 4)
if (length(k)>0){
   for (i in 1:length(k)) {
     ak<-which(clul$membership == k[i])
     cat(i,")",ak,"\n") 
   }
   cat("! Note, first residue in input file has index=1 !\n")
} else {
  cat("! List of LONG subgraphs is empty !\n")
}
cat("\n")
cat(" .... The End ....\n")

