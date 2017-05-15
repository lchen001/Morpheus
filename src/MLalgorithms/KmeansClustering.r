#Copyright 2016 Lingjiao Chen, Version 0.8.
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

KmeansClustering <- function(Data,Max_Iter,Center_Number,k_center,nS)
{
# k_center is the k centers, and YA is the assignment
	All1 = matrix(1, nS, 1);
	All1_k = matrix(1,1,Center_Number);
	All1_C = t(matrix(1,1,nrow(k_center)));	
	T2 = rowSums(Data^{2}) %*% All1_k;
	T22 = Data * 2;
	for( k in 1: Max_Iter )
	{
		Dist = T2 - as.matrix(T22 %*% k_center ) +  All1 %*% colSums(k_center ^2);
		YA = (Matrix(Dist == (rowMins((Dist)) %*% All1_k),sparse=TRUE))+0;
		k_center = as.matrix(  ( t(Data) %*% YA ) / ( (All1_C) %*% colSums(YA) )  );
	}
	return(list(k_center ,YA));
}