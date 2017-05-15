#Copyright 2016 Lingjiao Chen, Arun Kumar, Jeffrey Naughton, and Jignesh M. Patel, Version 0.8.
#All rights reserved.
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

##############################################################
print('Kmeans Clustering on dataset Expedia')

##############################################################
# Load in Packages
print('Loading Packages...')
library("Matrix");
library("matrixStats");
library(data.table);
MyDir = getwd();
SrcDir = paste(getwd(),'/src',sep="");
setwd(SrcDir)
file.sources = list.files(pattern="*.r",recursive=TRUE)
sapply(file.sources,source,.GlobalEnv)
setwd(MyDir)
# Parameter Setting
options(digits = 15)
options(scipen=20)
print('Loading Packages is done')

##############################################################
# Load in Data
print("Loading data...");
Sdir = './Data/Expedia/MLSSparse.txt';
S =readMM(file=Sdir)+0;
dS=ncol(S);

FK1dir = './Data/Expedia/MLFK1.csv';
R1dir = './Data/Expedia/MLR1Sparse.txt'
JSet1 = as.matrix(read.table(FK1dir,header=TRUE)); 
nS = nrow(JSet1);
nR1 = max(JSet1);
FK1 = sparseMatrix(i=c(1:nS),j=JSet1,x=1,dims=c(nS,nR1));
R1S = readMM(file=R1dir)+0;
dR1 = ncol(R1S);


FK2dir = './Data/Expedia/MLFK2.csv';
R2dir = './Data/Expedia/MLR2Sparse.txt';
JSet2 = as.matrix(read.table(FK2dir,header=TRUE)); 
nS = nrow(JSet2);
nR2 = max(JSet2);
FK2 = sparseMatrix(i=c(1:nS),j=JSet2,x=1,dims=c(nS,nR2));
R2S = readMM(file=R2dir)+0;
dR2 = ncol(R2S);

Ydir = './Data/Expedia/MLY.csv';
Ytest = as.matrix(( fread(Ydir)));
Y = (Ytest)+0;
print("Loading data is done");

#############################################################
# Creating Materialized data table and normalized matrix
print("Creating Materialized data table and normalized matrix");
T = cbind(S,FK1%*%R1S,FK2%*%R2S);
TNM = NormalMatrix(EntTable=list(S),AttTables=list(R1S,R2S),KFKDs=list(FK1,FK2),Sparse=TRUE);
print("Creating Materialized data table and normalized matrix is done");

#############################################################
# Parameter 
Center_Number = 5;
Max_Iter = 20;
k_center = as.matrix(t(T[1:Center_Number,1:(dS+dR1+dR2)])); # (dS+dR) x Center number
Result_eps = 1e-6;

#############################################################
# Materilaized Version
print("Materialized Version:");
Material_Time1 = Sys.time();
Material_Result = KmeansClustering(T,Max_Iter,Center_Number,k_center,nS);
Material_RunTime = Sys.time()-Material_Time1;
print("Materialized Version is done\n");
print("Materialization Runtime:");
print(Material_RunTime);

#############################################################
# Morpheus 
print("Morpheus:\n");
Morpheus_Time1 = Sys.time();
Morpheus_Result = KmeansClustering(TNM,Max_Iter,Center_Number,k_center,nS);
Morpheus_RunTime = Sys.time()-Morpheus_Time1;
print("Morpheus is done\n");
print("Morpheus Runtime: ");
print(Morpheus_RunTime);

if(norm(as.matrix(Material_Result[[1]] - Morpheus_Result[[1]]),'F')<Result_eps*(1+norm(as.matrix(Material_Result[[1]]),'F')) )
{
	print("Pass verification (Both approches return the same ML model).");
}else
{
	print("Fail to pass verification (Both approches return different models).");
}
