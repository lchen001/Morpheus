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
print('Logistic Regression on dataset Yelp')

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
print("Reading in data...");
S = matrix(0,0,0);
dS = nrow(S);

FK1dir = './Data/Yelp/MLFK1.csv';
R1dir = './Data/Yelp/MLR1Sparse.txt';
JSet1 = as.matrix(read.table(FK1dir,header=TRUE)); 
nS = nrow(JSet1);
nR1 = max(JSet1);
FK1 = sparseMatrix(i=c(1:nS),j=JSet1,x=1,dims=c(nS,nR1));
R1S = readMM(file=R1dir)+0;
dR1 = ncol(R1S);

FK2dir = './Data/Yelp/MLFK2.csv';
R2dir = './Data/Yelp/MLR2Sparse.txt';

JSet2 = as.matrix(read.table(FK2dir,header=TRUE)); 
nS = nrow(JSet2);
nR2 = max(JSet2);
FK2 = sparseMatrix(i=c(1:nS),j=JSet2,x=1,dims=c(nS,nR2));
R2S = readMM(file=R2dir)+0;
dR2 = ncol(R2S);

Ydir = './Data/Yelp/MLY.csv';
Ytest = as.matrix(( fread(Ydir)));
Y = (Ytest > 2.5)+0;
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
Max_Iter = 20;
Result_eps = 1e-6;

#############################################################
# Materilaized Version
print("Materialized Version:");
Material_Time1 = Sys.time();
Material_Result = LogisticRegression(T,Max_Iter, winit, gamma0, Y);
Material_RunTime = Sys.time()-Material_Time1;
print("Materialized Version is done\n");
print("Materialization Runtime:");
print(Material_RunTime);

#############################################################
# Morpheus 
print("Morpheus:\n");
Morpheus_Time1 = Sys.time();
Morpheus_Result = LogisticRegression(TNM,Max_Iter, winit, gamma0, Y);
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










library("Matrix");
library("matrixStats");
library(data.table);
source('./NMPackageTest/NormalMatrix.r')
# Parameter Setting
options(digits = 15)
options(scipen=20)
IterSet=c(1,5,10,20);
print("LR for real data!")

##############################################################
print("Reading in data...");
S = matrix(0,0,0);

FK1dir = './RealTest/Yelp/MLFK1.csv';
JSet1 = as.matrix(read.table(FK1dir,header=TRUE)); 
nS = nrow(JSet1);
nR1 = max(JSet1);
FK1 = sparseMatrix(i=c(1:nS),j=JSet1,x=1,dims=c(nS,nR1));
R1S = readMM(file='./RealTest/Yelp/MLR1Sparse.txt')+0;


FK2dir = './RealTest/Yelp/MLFK2.csv';
JSet2 = as.matrix(read.table(FK2dir,header=TRUE)); 
nS = nrow(JSet2);
nR2 = max(JSet2);
FK2 = sparseMatrix(i=c(1:nS),j=JSet2,x=1,dims=c(nS,nR2));
R2S = readMM(file='./RealTest/Yelp/MLR2Sparse.txt')+0;

Ydir = './RealTest/Yelp/MLY.csv';
Ytest = as.matrix(( fread(Ydir)));
Y = (Ytest > 2.5)+0;
x=Y;
print("Finish Reading!");
#############################################################
T = cbind(FK1%*%R1S,FK2%*%R2S);
TNM = NormalMatrix(EntTable=list(S),AttTables=list(R1S,R2S),KFKDs=list(FK1,FK2),Sparse=TRUE);
print("Finish Creating T and TNM!");
#############################################################

winit = Matrix(rnorm(ncol(T),2),ncol(T),1,sparse=TRUE);

gamma0 = 0.000001;
for( k_MaxIter in 1:length(IterSet))
{
	#Material
	print('LR Material No Disk');
	# Warm up the cache
	gc();
	gc();
	gc();
	gc();
	gc();
	w = winit;
	time1 = Sys.time(); 
	for( k in 1: 2 )
	{
		w = w - gamma0 * (t(T) %*% (x/(1+exp(T%*%w))) );
	}
	time2 = Sys.time();
	M4 =  difftime(time2, time1, units="secs");
	print(M4)
	print('Finish Warm-up');
						
	# Material LR
	gc();
	gc();
	gc();
	gc();
	gc();
	w = winit;
	time1 = Sys.time(); 
	for( k in 1: IterSet[k_MaxIter] )
	{
		w = w - gamma0 * (t(T) %*% (x/(1+exp(T%*%w))) );
	}
	time2 = Sys.time();
	M4 =  difftime(time2, time1, units="secs");
	print(M4)
	# Write the time interval 1
	timeF=as.numeric(M4);
	LogName=paste('./DataLog_RD/LR/Yelp/Log_LR_M_Yelp_iter_',IterSet[k_MaxIter],'_1.txt',sep="");#Log_LR_M_S_10000000x20_R_500000x80_x_100x1_3.txt
	write.table(as.matrix(timeF),LogName, sep=",", row.names = FALSE,  col.names=FALSE);
		
				gc();
				gc();
				gc();
				gc();
				gc();
				w = winit;
				time1 = Sys.time(); 
				for( k in 1: IterSet[k_MaxIter] )
				{
					w = w - gamma0 * (t(T) %*% (x/(1+exp(T%*%w))) );
				}
				time2 = Sys.time();
				M4 =  difftime(time2, time1, units="secs");
				print(M4)
				# Write the time interval 1
				timeF=as.numeric(M4);
					LogName=paste('./DataLog_RD/LR/Yelp/Log_LR_M_Yelp_iter_',IterSet[k_MaxIter],'_2.txt',sep="");#Log_LR_M_S_10000000x20_R_500000x80_x_100x1_3.txt
					write.table(as.matrix(timeF),LogName, sep=",", row.names = FALSE,  col.names=FALSE);
		
				gc();
				gc();
				gc();
				gc();
				gc();
				w = winit;
				time1 = Sys.time(); 
				for( k in 1: IterSet[k_MaxIter] )
				{
					w = w - gamma0 * (t(T) %*% (x/(1+exp(T%*%w))) );
				}
				time2 = Sys.time();
				M4 =  difftime(time2, time1, units="secs");
				print(M4)
				# Write the time interval 1
				timeF=as.numeric(M4);
				LogName=paste('./DataLog_RD/LR/Yelp/Log_LR_M_Yelp_iter_',IterSet[k_MaxIter],'_3.txt',sep="");#Log_LR_M_S_10000000x20_R_500000x80_x_100x1_3.txt
				write.table(as.matrix(timeF),LogName, sep=",", row.names = FALSE,  col.names=FALSE);
				
				
				




				#Factor
				print('LR Factor No Disk');
				# Warm up the cache
				gc();
				gc();
				gc();
				gc();
				gc();
				w=winit;
				time1 = Sys.time(); 
				for( k in 1: 2 )
				{
					print(k);
					w = w - gamma0 * (t(TNM) %*% (x/(1+exp(TNM%*%w))) );
				}
				time2 = Sys.time();
				M4 =  difftime(time2, time1, units="secs");
				print(M4)
				print('Finish Warm-up');
				
				# Factor LR
				gc();
				gc();
				gc();
				gc();
				gc();
				w=winit;				
				time1 = Sys.time(); 
				for( k in 1: IterSet[k_MaxIter] )
				{
					w = w - gamma0 * (t(TNM) %*% (x/(1+exp(TNM%*%w))) );
				}
				time2 = Sys.time();
				M4 =  difftime(time2, time1, units="secs");
				print(M4)
				# Write the time interval 1
				timeF=as.numeric(M4);
				LogName=paste('./DataLog_RD/LR/Yelp/Log_LR_F_Yelp_iter_',IterSet[k_MaxIter],'_1.txt',sep="");#Log_LR_M_S_10000000x20_R_500000x80_x_100x1_3.txt
				write.table(as.matrix(timeF),LogName, sep=",", row.names = FALSE,  col.names=FALSE);
		
				gc();
				gc();
				gc();
				gc();
				gc();
				w=winit;				
				time1 = Sys.time(); 
				for( k in 1: IterSet[k_MaxIter] )
				{
					w = w - gamma0 * (t(TNM) %*% (x/(1+exp(TNM%*%w))) );
				}
				time2 = Sys.time();
				M4 =  difftime(time2, time1, units="secs");
				print(M4)
				# Write the time interval 1
				timeF=as.numeric(M4);
				LogName=paste('./DataLog_RD/LR/Yelp/Log_LR_F_Yelp_iter_',IterSet[k_MaxIter],'_2.txt',sep="");#Log_LR_M_S_10000000x20_R_500000x80_x_100x1_3.txt
				write.table(as.matrix(timeF),LogName, sep=",", row.names = FALSE,  col.names=FALSE);
		

				gc();
				gc();
				gc();
				gc();
				gc();
				w=winit;				
				time1 = Sys.time(); 
				for( k in 1: IterSet[k_MaxIter] )
				{
					w = w - gamma0 * (t(TNM) %*% (x/(1+exp(TNM%*%w))) );
				}
				time2 = Sys.time();
				M4 =  difftime(time2, time1, units="secs");
				print(M4)
				# Write the time interval 1
				timeF=as.numeric(M4);
					LogName=paste('./DataLog_RD/LR/Yelp/Log_LR_F_Yelp_iter_',IterSet[k_MaxIter],'_3.txt',sep="");#Log_LR_M_S_10000000x20_R_500000x80_x_100x1_3.txt
					write.table(as.matrix(timeF),LogName, sep=",", row.names = FALSE,  col.names=FALSE);
				
}	
