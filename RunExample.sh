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
 
