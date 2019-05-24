%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ben Engelhard, Princeton University (2019).

This package is provided free without any warranty; you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. If this code is used, please cite: B Engelhard et al. Specialized coding of sensory, motor, and cognitive variables in VTA dopamine neurons. Nature, 2019.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Encoding Model

This package preprocesses joint behavioral and neuronal data and then processes it through an encoding model in order to obtain a quantitative measure of the contributions of the behavioral variables to the activity of single neurons, as described in the following paper: B Engelhard et al. Specialized coding of sensory, motor, and cognitive variables in midbrain dopamine neurons. Nature, 2019. Please see the paper for details on the encoding model.

Function list:
make_predictor_matrix_generalcase.m
process_encoding_model.m
find_non_empty_cells.m
get_CV_R2.m

Data files list:
spline_basis30_int.mat

Instructions:
First, the data has to be formatted to be used in the make_predictor_matrix_generalcase function. See the function header for specific details. After the predictor matrix is generated, it can be directly processed with the process_encoding_model function. In order to obtain a measure of significance for the relationship between behavioral variables and the neural activity, the obtained F-statistic for each behavioral variable should then be compared to a distribution of F-statistics obtained from a (reasonably large) number of shuffled data instantiations (i.e. the neural activity shuffled but with the same matrix of predictors). 

Notes: 
Currently, event variables are convolved with a basis set composed of 7 splines and 30 timepoints. If you wish to change this, you may use the following package to generate a different spline basis set:
Ramsay JO (2014). fdaM: Functional Data Analysis. MATLAB package, URL http://www. psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/.
The generated spline basis set should be named spline_basis and saved in a matfile named 'spline_basis30_int.mat'.


