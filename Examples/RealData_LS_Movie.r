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
print('Logistic Regression on dataset Movie')

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
##############################################################
print("Loading data...");
S = matrix(0,0,0);

FK1dir = './Data/MovieLens1M/MLFK1.csv';
R1dir = './Data/MovieLens1M/MLR1Sparse.txt';
JSet1 = as.matrix(read.table(FK1dir,header=TRUE)); 
nS = nrow(JSet1);
nR1 = max(JSet1);
FK1 = sparseMatrix(i=c(1:nS),j=JSet1,x=1,dims=c(nS,nR1));
R1S = readMM(file=R1dir)+0;


FK2dir = './Data/MovieLens1M/MLFK2.csv';
R2dir = './Data/MovieLens1M/MLR2Sparse.txt';
JSet2 = as.matrix(read.table(FK2dir,header=TRUE)); 
nS = nrow(JSet2);
nR2 = max(JSet2);
FK2 = sparseMatrix(i=c(1:nS),j=JSet2,x=1,dims=c(nS,nR2));
R2S = readMM(file=R2dir)+0;

Ydir = './Data/MovieLens1M/MLY.csv';
Ytest = as.matrix(( fread(Ydir)));
Y = Ytest+0;
print("Loading data is done");

#############################################################
# Creating Materialized data table and normalized matrix
print("Creating Materialized data table and normalized matrix");
T = cbind(FK1%*%R1S,FK2%*%R2S);
TNM = NormalMatrix(EntTable=list(S),AttTables=list(R1S,R2S),KFKDs=list(FK1,FK2),Sparse=TRUE);
print("Creating Materialized data table and normalized matrix is done");

#############################################################
# Parameter 
winit = Matrix(rnorm(ncol(T),2),ncol(T),1,sparse=TRUE);
gamma0 = 0.000001;
Max_Iter = 1;
Result_eps = 1e-6;

#############################################################
# Materilaized Version
print("Materialized Version:");
Material_Time1 = Sys.time();
Material_Result = LinearRegression(T,Max_Iter, winit, gamma0, Y);
Material_RunTime = Sys.time()-Material_Time1;
print("Materialized Version is done\n");
print("Materialization Runtime:");
print(Material_RunTime);

#############################################################
# Morpheus 
print("Morpheus:\n");
Morpheus_Time1 = Sys.time();
Morpheus_Result = LinearRegression(TNM,Max_Iter, winit, gamma0, Y);
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
