#!/bin/bash  
# clear

# Logistic Regression
R -f  Examples/RealData_LR_Expedia.r
R -f  Examples/RealData_LR_Movie.r
R -f  Examples/RealData_LR_Yelp.r

# Linear Regression
R -f  Examples/RealData_LS_Expedia.r
R -f  Examples/RealData_LS_Movie.r
R -f  Examples/RealData_LS_Yelp.r

# Kmeans Clustering
R -f  Examples/RealData_KmeansCLA_Expedia.r
R -f  Examples/RealData_KmeansCLA_Movie.r
R -f  Examples/RealData_KmeansCLA_Yelp.r
 