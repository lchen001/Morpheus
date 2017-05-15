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

library(Matrix);
NormalMatrix <- setClass (
	"NormalMatrix",
	slots = c (
		EntTable = "list",
		AttTables = "list",
		KFKDs = "list",
		Trans = "logical",
		Sparse = "logical"
		),
	prototype=prototype(Trans=FALSE,Sparse=FALSE)	
	)

#############################################
# Transpose Rule
't.NormalMatrix' <- function( NM1)
{
	NM1@Trans=!NM1@Trans;
	return(NM1);
}
#############################################

#############################################
# Element-wise Scalar Operators
setMethod("+", c("numeric","NormalMatrix"), function(e1,e2) {
	Rnum = length(e2@AttTables);
	Sempty = nrow(e2@EntTable[[1]])*ncol(e2@EntTable[[1]]);
	if(Sempty!=0)
	{
		e2@EntTable[[1]]=e2@EntTable[[1]]+e1;
	}
	for(i in 1:Rnum)
	{
		e2@AttTables[[i]] = e2@AttTables[[i]]+e1;
	}	
	return(e2);
})

setMethod("+", c("NormalMatrix","numeric"), function(e1,e2) {
	Rnum = length(e1@AttTables);
	Sempty = nrow(e1@EntTable[[1]])*ncol(e1@EntTable[[1]]);
	if(Sempty!=0)
	{
		e1@EntTable[[1]]=e1@EntTable[[1]]+e2;
	}
	for(i in 1:Rnum)
	{
		e1@AttTables[[i]] = e1@AttTables[[i]]+e2;
	}	
	return(e1);
})

setMethod("-", c("numeric","NormalMatrix"), function(e1,e2) {
	Rnum = length(e2@AttTables);
	Sempty = nrow(e2@EntTable[[1]])*ncol(e2@EntTable[[1]]);
	if(Sempty!=0)
	{
		e2@EntTable[[1]]=e2@EntTable[[1]]-e1;
	}
	for(i in 1:Rnum)
	{
		e2@AttTables[[i]] = e2@AttTables[[i]]-e1;
	}	
	return(e2);
})

setMethod("-", c("NormalMatrix","numeric"), function(e1,e2) {
	Rnum = length(e1@AttTables);
	Sempty = nrow(e1@EntTable[[1]])*ncol(e1@EntTable[[1]]);
	if(Sempty!=0)
	{
		e1@EntTable[[1]]=e1@EntTable[[1]]-e2;
	}
	for(i in 1:Rnum)
	{
		e1@AttTables[[i]] = e1@AttTables[[i]]-e2;
	}	
	return(e1);
})

setMethod("*", c("numeric","NormalMatrix"), function(e1,e2) {
	Rnum = length(e2@AttTables);
	Sempty = nrow(e2@EntTable[[1]])*ncol(e2@EntTable[[1]]);
	if(Sempty!=0)
	{
		e2@EntTable[[1]]=e2@EntTable[[1]]*e1;
	}
	for(i in 1:Rnum)
	{
		e2@AttTables[[i]] = e2@AttTables[[i]]*e1;
	}	
	return(e2);
})

setMethod("*", c("NormalMatrix","numeric"), function(e1,e2) {
	Rnum = length(e1@AttTables);
	Sempty = nrow(e1@EntTable[[1]])*ncol(e1@EntTable[[1]]);
	if(Sempty!=0)
	{
		e1@EntTable[[1]]=e1@EntTable[[1]]*e2;
	}
	for(i in 1:Rnum)
	{
		e1@AttTables[[i]] = e1@AttTables[[i]]*e2;
	}	
	return(e1);
})

setMethod("/", c("numeric","NormalMatrix"), function(e1,e2) {
	Rnum = length(e2@AttTables);
	Sempty = nrow(e2@EntTable[[1]])*ncol(e2@EntTable[[1]]);
	if(Sempty!=0)
	{
		e2@EntTable[[1]]=e2@EntTable[[1]]/e1;
	}
	for(i in 1:Rnum)
	{
		e2@AttTables[[i]] = e2@AttTables[[i]]/e1;
	}	
	return(e2);
})

setMethod("/", c("NormalMatrix","numeric"), function(e1,e2) {
	Rnum = length(e1@AttTables);
	Sempty = nrow(e1@EntTable[[1]])*ncol(e1@EntTable[[1]]);
	if(Sempty!=0)
	{
		e1@EntTable[[1]]=e1@EntTable[[1]]/e2;
	}
	for(i in 1:Rnum)
	{
		e1@AttTables[[i]] = e1@AttTables[[i]]/e2;
	}	
	return(e1);
})

setMethod("^", c("numeric","NormalMatrix"), function(e1,e2) {
	Rnum = length(e2@AttTables);
	Sempty = nrow(e2@EntTable[[1]])*ncol(e2@EntTable[[1]]);
	if(Sempty!=0)
	{
		e2@EntTable[[1]]=e2@EntTable[[1]]^e1;
	}
	for(i in 1:Rnum)
	{
		e2@AttTables[[i]] = e2@AttTables[[i]]^e1;
	}	
	return(e2);
})

setMethod("^", c("NormalMatrix","numeric"), function(e1,e2) {
	Rnum = length(e1@AttTables);
	Sempty = nrow(e1@EntTable[[1]])*ncol(e1@EntTable[[1]]);
	if(Sempty!=0)
	{
		e1@EntTable[[1]]=e1@EntTable[[1]]^e2;
	}
	for(i in 1:Rnum)
	{
		e1@AttTables[[i]] = e1@AttTables[[i]]^e2;
	}	
	return(e1);
})

#############################################

#############################################
#Aggregation
setMethod("sum", c("NormalMatrix"), function(x) {
	Sempty = nrow(x@EntTable[[1]])*ncol(x@EntTable[[1]]);
	y = ifelse(Sempty!=0,sum(x@EntTable[[1]]),0 );
	Rnum = length(x@AttTables);
	for(i in 1:Rnum)
	{
		y = y + colSums(x@KFKDs[[i]])%*%rowSums(x@AttTables[[i]]);
	}	
	return(y);
})

setMethod("rowSums", c("NormalMatrix"), function(x) {
	if(x@Trans==TRUE)
	{
		return(t(NormalMatrixcolSums(x)));
	}
	else
	{
		return(NormalMatrixrowSums(x));
	}
})

setMethod("colSums", c("NormalMatrix"), function(x) {
	if(x@Trans==TRUE)
	{
		return(t(NormalMatrixrowSums(x)));
	}
	else
	{
		return(NormalMatrixcolSums(x));
	}
})

NormalMatrixrowSums <- function(x){
	Sempty = nrow(x@EntTable[[1]])*ncol(x@EntTable[[1]]);
	i = 1;
	if(Sempty!=0)
	{
		y = rowSums(x@EntTable[[1]])+as.numeric(x@KFKDs[[i]]%*%rowSums(x@AttTables[[i]]));;
	}
	else
	{
		y = as.numeric(x@KFKDs[[i]]%*%rowSums(x@AttTables[[i]]));
	}
	Rnum = length(x@AttTables);
	i = 2;
	while(i <= Rnum)
	{
		y = y + as.numeric(x@KFKDs[[i]]%*%rowSums(x@AttTables[[i]]));
		i = i + 1;
	}	
	return(y);
}

NormalMatrixcolSums <- function(x){
	Sempty = nrow(x@EntTable[[1]])*ncol(x@EntTable[[1]]);
	Rnum = length(x@AttTables);
	if(Sempty!=0)
	{
		y = c(colSums(x@EntTable[[1]]), colSums(x@KFKDs[[1]])%*%x@AttTables[[1]]);
	}
	else
	{
		y = colSums(x@KFKDs[[1]])%*%x@AttTables[[1]];
	}	
	i = 2;
	while(i <= Rnum)
	{
		y = c(y, colSums(x@KFKDs[[i]])%*%x@AttTables[[i]] );
		i = i + 1;
	}	
	return(y);
}


#############################################

#############################################
# Multiplication

# Left Matrix Multiplication
NormalMatrixLMM <- function(x,y)
{
	Sempty = nrow(x@EntTable[[1]])*ncol(x@EntTable[[1]]);
	Rnum = length(x@AttTables);
	if(x@Sparse==TRUE) # EntTable is empty
	{
		if(Sempty==0) # EntTable is empty
		{
			Out1 = x@KFKDs[[1]] %*% (x@AttTables[[1]] %*% y[1:(ncol(x@AttTables[[1]]) ),] );
		}
		else
		{
			Out1 = x@EntTable[[1]]%*% y[1:(ncol(x@EntTable[[1]]) ),] + x@KFKDs[[1]] %*% (x@AttTables[[1]] %*% y[(1+ncol(x@EntTable[[1]])):(ncol(x@AttTables[[1]])+ncol(x@EntTable[[1]]) ),] );
		}		
		StartCol = ncol(x@EntTable[[1]])+ncol(x@AttTables[[1]]);
		i = 2;
		while(i <= Rnum) # multi-table joins
		{
			Out1 = Out1 + (x@KFKDs[[i]]) %*% (x@AttTables[[i]] %*% y[(1+StartCol):(StartCol+ncol(x@AttTables[[i]]) ),] );	
			StartCol=StartCol+(ncol(x@AttTables[[i]]) );
			i = i + 1;
		}
		return (Out1);
	}
	else
	{
		if(Sempty==0) # EntTable is empty
		{
			Out1 = as.matrix(x@KFKDs[[1]] %*% (x@AttTables[[1]] %*% y[1:(ncol(x@AttTables[[1]]) ),] ));
		}
		else
		{
			Out1 = x@EntTable[[1]]%*% y[1:(ncol(x@EntTable[[1]]) ),] + as.matrix(x@KFKDs[[1]] %*% (x@AttTables[[1]] %*% y[(1+ncol(x@EntTable[[1]])):(ncol(x@AttTables[[1]])+ncol(x@EntTable[[1]]) ),] ));
		}		
		StartCol = ncol(x@EntTable[[1]])+ncol(x@AttTables[[1]]);
		i = 2;			
		while(i <= Rnum) # multi-table joins
		{
			Out1 = Out1 + as.matrix((x@KFKDs[[i]]) %*% (x@AttTables[[i]] %*% y[(1+StartCol):(StartCol+ncol(x@AttTables[[i]]) ),] ));	
			StartCol=StartCol+(ncol(x@AttTables[[i]]) );
			i = i + 1;
		}
		return (Out1);

	}
}

# Right Matrix Multiplication
NormalMatrixRMM <- function(x,y)
{
	Sempty = nrow(y@EntTable[[1]])*ncol(y@EntTable[[1]]);
	Rnum = length(y@AttTables);
	if(y@Sparse==TRUE)# Sparse 
	{
		if(Sempty==0) # EntTable is empty
		{
			Out1 = (x %*% y@KFKDs[[1]])%*%y@AttTables[[1]];
		}
		else
		{
			Out1 = cbind(x%*% y@EntTable[[1]],(x %*% y@KFKDs[[1]])%*%y@AttTables[[1]]);	
		}	
		i = 2;	
		while( i <= Rnum) # multi-table joins
		{
			Out1 = cbind(Out1,(x %*% y@KFKDs[[i]])%*%y@AttTables[[i]]);
			i = i + 1;	
		}		
		return (Out1);
	}
	else
	{
		if(Sempty==0) # EntTable is empty
		{
			Out1 = as.matrix(x %*% y@KFKDs[[1]])%*%y@AttTables[[1]];
		}
		else
		{
			Out1 = cbind(x%*% y@EntTable[[1]],as.matrix(x %*% y@KFKDs[[1]])%*%y@AttTables[[1]]);	
		}		
		i = 2;	
		while( i <= Rnum) # multi-table joins
		{
			Out1 = cbind(Out1,as.matrix(x %*% y@KFKDs[[i]])%*%y@AttTables[[i]]);
			i = i + 1;	
		}		
		return (Out1);
	}
}

# API for LMM and RMM with different matrix types.

#LMM API.
setMethod("%*%", c("NormalMatrix","dgCMatrix"), function(x,y) 
{
	if(x@Trans==FALSE)
	{
		return( NormalMatrixLMM(x,y) );
	}
	else
	{
		return( t(NormalMatrixRMM(t(y),x)) );
	}
})

setMethod("%*%", c("NormalMatrix","matrix"), function(x,y) 
{
	if(x@Trans==FALSE)
	{
		return( NormalMatrixLMM(x,y) );
	}
	else
	{
		return( t(NormalMatrixRMM(t(y),x)) );
	}
})

setMethod("%*%", c("NormalMatrix","dgeMatrix"), function(x,y) 
{
	if(x@Trans==FALSE)
	{
		return( NormalMatrixLMM(x,y) );
	}
	else
	{
		return( t(NormalMatrixRMM(t(y),x)) );
	}
})

setMethod("%*%", c("NormalMatrix","lgCMatrix"), function(x,y) 
{
	if(x@Trans==FALSE)
	{
		return( NormalMatrixLMM(x,y) );
	}
	else
	{
		return( t(NormalMatrixRMM(t(y),x)) );
	}
})

#RMM API.
setMethod("%*%", c("dgCMatrix","NormalMatrix"), function(x,y) {
	if(y@Trans==FALSE)
	{
		return( NormalMatrixRMM(x,y) );
	}
	else
	{
		return( t(NormalMatrixLMM(y,t(x))) );
	}
})

setMethod("%*%", c("dgeMatrix","NormalMatrix"), function(x,y) {
	if(y@Trans==FALSE)
	{
		return( NormalMatrixRMM(x,y) );
	}
	else
	{
		return( t(NormalMatrixLMM(y,t(x))) );
	}
})

setMethod("%*%", c("lgCMatrix","NormalMatrix"), function(x,y) {
	if(y@Trans==FALSE)
	{
		return( NormalMatrixRMM(x,y) );
	}
	else
	{
		return( t(NormalMatrixLMM(y,t(x))) );
	}
})

setMethod("%*%", c("matrix","NormalMatrix"), function(x,y) {
	if(y@Trans==FALSE)
	{
		return( NormalMatrixRMM(x,y) );
	}
	else
	{
		return( t(NormalMatrixLMM(y,t(x))) );
	}
})

#############################################


#############################################
#Derived Multiplication Operators.
#Include Crossprod and double multiplication.
#############################################

#############################################
# Cross Product
NormalMatrixCrossprod <- function(x)
{
	Sempty = nrow(x@EntTable[[1]])*ncol(x@EntTable[[1]]);
	Rnum = length(x@AttTables);
	if(x@Sparse==TRUE)# Sparse 
	{
		if(Sempty==0) # EntTable is empty
		{
			Out1 = ( crossprod( Diagonal( x=colSums(x@KFKDs[[1]])^{1/2} ) %*% x@AttTables[[1]]) ) ;
			i = 2;	
			while( i <= Rnum) # multi-table joins
			{
				# Construct Y_21
				j = 1;
				Y21 = crossprod( x@AttTables[[i]], crossprod(crossprod(x@KFKDs[[j]], x@KFKDs[[i]] ), x@AttTables[[j]]) ) ;
				j = 2;
				while( j <=(i-1))
				{
					j = j + 1;
					Y21 = cbind(Y21, crossprod( x@AttTables[[i]], crossprod(crossprod(x@KFKDs[[j]], x@KFKDs[[i]] ), x@AttTables[[j]]) )) ;
				}
				Out1 = (rbind( cbind(Out1,t(Y21) ) , cbind( Y21, (crossprod( Diagonal( x=colSums(x@KFKDs[[i]])^{1/2} ) %*% x@AttTables[[i]]) ) )) );
				i = i + 1;
			}
			return( Out1 );
		}
		else
		{
			Y21 = crossprod(x@AttTables[[1]], ( crossprod(x@KFKDs[[1]], x@EntTable[[1]])  ) );
			Out1 = (rbind( cbind(crossprod(x@EntTable[[1]]),t(Y21) ) , cbind( Y21, (crossprod( Diagonal( x=colSums(x@KFKDs[[1]])^{1/2} ) %*% x@AttTables[[1]]) ) )) );
			i = 2;	
			while( i <= Rnum) # multi-table joins
			{
				# Construct Y_21
				Y21 = crossprod(x@AttTables[[i]], ( crossprod(x@KFKDs[[i]], x@EntTable[[1]])  ) );
				for( j in 1:(i-1))
				{
					Y21 = cbind(Y21, crossprod( x@AttTables[[i]], crossprod(crossprod(x@KFKDs[[j]], x@KFKDs[[i]] ), x@AttTables[[j]]) )) ;
				}
				Out1 = (rbind( cbind(Out1,t(Y21) ) , cbind( Y21, (crossprod( Diagonal( x=colSums(x@KFKDs[[i]])^{1/2} ) %*% x@AttTables[[i]]) ) )) );
				i = i + 1;
			}
			return( Out1 );
		}		
	}
	else
	{
		if(Sempty==0) # EntTable is empty
		{
			Out1 = as.matrix( crossprod( Diagonal( x=colSums(x@KFKDs[[1]])^{1/2} ) %*% x@AttTables[[1]]) ) ;
			i = 2;	
			while( i <= Rnum) # multi-table joins
			{
				# Construct Y_21
				j = 1;
				Y21 = crossprod( x@AttTables[[i]], as.matrix(crossprod(crossprod(x@KFKDs[[j]], x@KFKDs[[i]] ), x@AttTables[[j]]) )) ;
				j = 2;
				while( j <=(i-1))
				{
					j = j + 1;
					Y21 = cbind(Y21, crossprod( x@AttTables[[i]], as.matrix(crossprod(crossprod(x@KFKDs[[j]], x@KFKDs[[i]] ), x@AttTables[[j]]) )) ) ;
				}
				Out1 = (rbind( cbind(Out1,t(Y21) ) , cbind( Y21, as.matrix(crossprod( Diagonal( x=colSums(x@KFKDs[[i]])^{1/2} ) %*% x@AttTables[[i]]) ) )) );
				i = i + 1;			
			}
			return( Out1 );
		}
		else
		{
			Y21 = crossprod(x@AttTables[[1]], as.matrix( crossprod(x@KFKDs[[1]], x@EntTable[[1]])  ) );
			Out1 = (rbind( cbind(crossprod(x@EntTable[[1]]),t(Y21) ) , cbind( Y21, as.matrix(crossprod( Diagonal( x=colSums(x@KFKDs[[1]])^{1/2} ) %*% x@AttTables[[1]]) ) )) );
			i = 2;	
			while( i <= Rnum) # multi-table joins
			{
				# Construct Y_21
				Y21 = crossprod(x@AttTables[[i]], ( crossprod(x@KFKDs[[i]], x@EntTable[[1]])  ) );
				for( j in 1:(i-1))
				{
					Y21 = cbind(Y21, crossprod( x@AttTables[[i]], as.matrix(crossprod(crossprod(x@KFKDs[[j]], x@KFKDs[[i]] ), x@AttTables[[j]]) ))) ;
				}
				Out1 = (rbind( cbind(Out1,t(Y21) ) , cbind( Y21, as.matrix(crossprod( Diagonal( x=colSums(x@KFKDs[[i]])^{1/2} ) %*% x@AttTables[[i]]) ) )) );
				i = i + 1;			
			}
			return( Out1 );
		}	
	}
}

# Transpose Cross Product
NormalMatrixCrossprodTrans <- function(x)
{
	Sempty = nrow(x@EntTable[[1]])*ncol(x@EntTable[[1]]);
	Rnum = length(x@AttTables);
	if(x@Sparse==TRUE)# Sparse 
	{
		if(Sempty==0) # EntTable is empty
		{
			Out1 = tcrossprod(x@KFKDs[[1]] %*% tcrossprod(x@AttTables[[1]]), x@KFKDs[[1]]);
		}
		else
		{
			Out1 = tcrossprod(x@EntTable[[1]]) + tcrossprod(x@KFKDs[[1]] %*% tcrossprod(x@AttTables[[1]]), x@KFKDs[[1]]);
		}
		i = 2;	
		while( i <= Rnum) # multi-table joins
		{
			Out1 = Out1 + tcrossprod(x@KFKDs[[i]] %*% tcrossprod(x@AttTables[[i]]), x@KFKDs[[i]]);
			i = i + 1;
		}
		return( Out1 );
	}
	else
	{
		if(Sempty==0) # EntTable is empty
		{
			Out1 = as.matrix(tcrossprod(x@KFKDs[[1]] %*% tcrossprod(x@AttTables[[1]]), x@KFKDs[[1]]));
		}
		else
		{
			Out1 = tcrossprod(x@EntTable[[1]]) + as.matrix(tcrossprod(x@KFKDs[[1]] %*% tcrossprod(x@AttTables[[1]]), x@KFKDs[[1]]));
		}
		i = 2;	
		while( i <= Rnum) # multi-table joins
		{
			Out1 = Out1 + as.matrix(tcrossprod(x@KFKDs[[i]] %*% tcrossprod(x@AttTables[[i]]), x@KFKDs[[i]]));
			i = i + 1;
		}
		return( Out1 );
	}	
}

# Crossprod API
setMethod("crossprod", c("NormalMatrix"), function(x) {
	if(x@Trans==FALSE)
	{
		return( NormalMatrixCrossprod(x) );
	}
	else
	{
		return( (NormalMatrixCrossprodTrans(x)) );
	}
})

setMethod("tcrossprod", c("NormalMatrix"), function(x) {
	if(x@Trans==FALSE)
	{
		return( NormalMatrixCrossprodTrans(x) );
	}
	else
	{
		return( (NormalMatrixCrossprod(x)) );
	}
})
#############################################


#############################################
#Double Multiplication

#############################################
# Cross Product Alike, i.e., t(A)%*% B 
# Warning: Only test the case when A = B.
NormalMatrixCrossprod2M <- function(x,y)
{
	Sempty1 = nrow(x@EntTable[[1]])*ncol(x@EntTable[[1]]);
	Rnum1 = length(x@AttTables);
	Sempty2 = nrow(y@EntTable[[1]])*ncol(y@EntTable[[1]]);
	Rnum2 = length(y@AttTables);
	if(x@Sparse==TRUE)# Sparse 
	{
		# First compute attribute tables
		i = 1;
		j = 1;
		Out1 = t(x@AttTables[[i]]) %*%( t(x@KFKDs[[i]]) %*% y@KFKDs[[j]] ) %*% (y@AttTables[[j]]);
		j = j+1;
		while( 2<=j && j<= Rnum2)# Enter this loop only if B is multi-table.
		{
			Out1 = cbind(Out1,t(x@AttTables[[i]]) %*%( t(x@KFKDs[[i]]) %*% y@KFKDs[[j]] ) %*% (y@AttTables[[j]]));
			j = j+1;
		}
		i = i+1;
		while( 2<=i && i<=Rnum1 ) # Enter this loop only if A is multi-table.
		{
			j = 1;
			Temp = t(x@AttTables[[i]]) %*%( t(x@KFKDs[[i]]) %*% y@KFKDs[[j]] ) %*% (y@AttTables[[j]]);
			j = j+1;
			while( 2<=j && j<= Rnum2)# Enter this loop only if B is multi-table.
			{
				Temp = cbind(Temp,t(x@AttTables[[i]]) %*%( t(x@KFKDs[[i]]) %*% y@KFKDs[[j]] ) %*% (y@AttTables[[j]]));
				j = j+1;
			}
			Out1 = rbind(Out1,Temp);
			i=i+1;
		}
		# Finish attribute tables.
		
		# Consider Entity table.
		if(Sempty1 == 0 && Sempty2 != 0 ) # No Entity table in A.
		{
			j = 1;
			Temp = t(x@AttTables[[j]]) %*%( t(x@KFKDs[[j]]) %*% y@EntTable[[1]] ) ;
			j = j+1;
			while( 2<=j && j<= Rnum1)# Enter this loop only if A is multi-table.
			{
				Temp = rbind(Temp,t(x@AttTables[[j]]) %*%( t(x@KFKDs[[j]]) %*% y@EntTable[[1]] ));
				j = j+1;
			}			
			Out1 = cbind(Temp,Out1);
		}
		if(Sempty1 != 0 && Sempty2 == 0 ) # No Entity table in B.
		{
			j = 1;
			Temp = (t(x@EntTable[[1]]) %*% y@KFKDs[[j]]) %*% y@AttTables[[j]]  ;
			j = j+1;
			while( 2<=j && j<= Rnum2)# Enter this loop only if B is multi-table.
			{
				Temp = cbind(Temp,(t(x@EntTable[[1]]) %*% y@KFKDs[[j]]) %*% y@AttTables[[j]] );
				j = j+1;
			}			
			Out1 = rbind(Temp,Out1);
		}
		if(Sempty1 != 0 && Sempty2 != 0 ) # Both with Entity table.
		{
			# EntTable in A.
			j = 1;
			Temp = t(x@AttTables[[j]]) %*%( t(x@KFKDs[[j]]) %*% y@EntTable[[1]] ) ;
			j = j+1;
			while( 2<=j && j<= Rnum1)# Enter this loop only if A is multi-table.
			{
				Temp = rbind(Temp,t(x@AttTables[[j]]) %*%( t(x@KFKDs[[j]]) %*% y@EntTable[[1]] ) );
				j = j+1;
			}			
			Out1 = cbind(Temp,Out1);
			# EntTable in B.
			j = 1;
			Temp = cbind( (t(x@EntTable[[1]]) %*% y@EntTable[[1]]) , (t(x@EntTable[[1]]) %*% y@KFKDs[[j]]) %*% y@AttTables[[j]] ) ;
			j = j+1;
			while( 2<=j && j<= Rnum2)# Enter this loop only if B is multi-table.
			{
				Temp = cbind(Temp,(t(x@EntTable[[1]]) %*% y@KFKDs[[j]]) %*% y@AttTables[[j]] );
				j = j+1;
			}			
			Out1 = rbind(Temp,Out1);
		}
			
			return(Out1);	
	}
	else
	{
		# First compute attribute tables
		i = 1;
		j = 1;
		Out1 = t(x@AttTables[[i]]) %*%as.matrix(( t(x@KFKDs[[i]]) %*% y@KFKDs[[j]] ) %*% (y@AttTables[[j]]));
		j = j+1;
		while( 2<=j && j<= Rnum2)# Enter this loop only if B is multi-table.
		{
			Out1 = cbind(Out1,t(x@AttTables[[i]]) %*%as.matrix(( t(x@KFKDs[[i]]) %*% y@KFKDs[[j]] ) %*% (y@AttTables[[j]])));
			j = j+1;
		}
		i = i+1;
		while( 2<=i && i<=Rnum1 ) # Enter this loop only if A is multi-table.
		{
			j = 1;
			Temp = t(x@AttTables[[i]]) %*%as.matrix(( t(x@KFKDs[[i]]) %*% y@KFKDs[[j]] ) %*% (y@AttTables[[j]]));
			j = j+1;
			while( 2<=j && j<= Rnum2)# Enter this loop only if B is multi-table.
			{
				Temp = cbind(Temp,t(x@AttTables[[i]]) %*%as.matrix(( t(x@KFKDs[[i]]) %*% y@KFKDs[[j]] ) %*% (y@AttTables[[j]])));
				j = j+1;
			}
			Out1 = rbind(Out1,Temp);
			i=i+1;
		}
		# Finish attribute tables.
		
		# Consider Entity table.
		if(Sempty1 == 0 && Sempty2 != 0 ) # No Entity table in A.
		{
			j = 1;
			Temp = t(x@AttTables[[j]]) %*%as.matrix( t(x@KFKDs[[j]]) %*% y@EntTable[[1]] ) ;
			j = j+1;
			while( 2<=j && j<= Rnum1)# Enter this loop only if A is multi-table.
			{
				Temp = rbind(Temp,t(x@AttTables[[j]]) %*%as.matrix( t(x@KFKDs[[j]]) %*% y@EntTable[[1]] ) );
				j = j+1;
			}			
			Out1 = cbind(Temp,Out1);
		}
		if(Sempty1 != 0 && Sempty2 == 0 ) # No Entity table in B.
		{
			j = 1;
			Temp = as.matrix(t(x@EntTable[[1]]) %*% y@KFKDs[[j]]) %*% y@AttTables[[j]]  ;
			j = j+1;
			while( 2<=j && j<= Rnum2)# Enter this loop only if B is multi-table.
			{
				Temp = cbind(Temp,as.matrix(t(x@EntTable[[1]]) %*% y@KFKDs[[j]]) %*% y@AttTables[[j]] );
				j = j+1;
			}			
			Out1 = rbind(Temp,Out1);
		}
		if(Sempty1 != 0 && Sempty2 != 0 ) # Both with Entity table.
		{
			# EntTable in A.
			j = 1;
			Temp = t(x@AttTables[[j]]) %*%as.matrix( t(x@KFKDs[[j]]) %*% y@EntTable[[1]] ) ;
			j = j+1;
			while( 2<=j && j<= Rnum1)# Enter this loop only if A is multi-table.
			{
				Temp = rbind(Temp,t(x@AttTables[[j]]) %*%as.matrix( t(x@KFKDs[[j]]) %*% y@EntTable[[1]] ) );
				j = j+1;
			}			
			Out1 = cbind(Temp,Out1);
			# EntTable in B.
			j = 1;
			Temp = cbind( (t(x@EntTable[[1]]) %*% y@EntTable[[1]]) , as.matrix(t(x@EntTable[[1]]) %*% y@KFKDs[[j]]) %*% y@AttTables[[j]] ) ;
			j = j+1;
			while( 2<=j && j<= Rnum2)# Enter this loop only if B is multi-table.
			{
				Temp = cbind(Temp,as.matrix(t(x@EntTable[[1]]) %*% y@KFKDs[[j]]) %*% y@AttTables[[j]] );
				j = j+1;
			}			
			Out1 = rbind(Temp,Out1);
		}
		return( Out1 );
	}	
}

#Note: Double Multiplication is not fully optimized.
#############################################




#############################################
# M:N joins
NormalMatrixJoin <- function(S,R,SJA,RJA,IsSparse=FALSE){
	# Input : S, R, SJA, RJA
	nS = nrow(S);
	dS = ncol(S);
	nR = nrow(R);
	dR = ncol(R);

	SJAS = sparseMatrix(c(1:nS),SJA );
	RJAS = sparseMatrix(c(1:nR),RJA );
	if( ncol(SJAS) > ncol(RJAS) )
	{
		SJAS = SJAS[,1:ncol(RJAS)];
	}
	if( ncol(SJAS) < ncol(RJAS) )
	{
		RJAS = RJAS[,1:ncol(SJAS)];
	}
	JM = SJAS %*%t(RJAS);
	nT = sum(colSums((JM)));
	if(nT!=0)
	{
		JMS = summary(JM);
		SJM = sparseMatrix(c(1:nT),JMS[,1],dims=c(nT,nS),x=1);
		RJM = sparseMatrix(c(1:nT),JMS[,2],dims=c(nT,nR),x=1);
	}
	if(nT==0)
	{
		SJM = Matrix(0,nT,nS);
		RJM = Matrix(0,nT,nR);
	}
	S0 = matrix(0,0,0);
	TNM = NormalMatrix(EntTable=list(S0),AttTables=list(S,R),KFKDs=list(SJM,RJM),Sparse=IsSparse);
	return (TNM);
}
#############################################
#Handle conflicts with regular matrix package
setMethod(Ops, c("NormalMatrix","NormalMatrix"), function(e1,e2) {
	return ();
})