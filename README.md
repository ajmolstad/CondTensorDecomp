# CondTensorDecomp
Simulation code for "Conditional probability tensor decompositions for multivariate categorical response regression" 

To run the simulations displayed in Figure 5, execute the Run_R3.sh file. This will initiate 2400 jobs, each of which corresponds to a single replicate. For each replicate, the .sh will first run the R3.R script, which generates data from the appropriate data generating model and then sources the Main_Random.R script. The Main_Random.R script is where model fitting and performance evaluation is done. The Random suffix refers to the method of initialization for the B's at each candidate tuning parameter value. The only function needed is in the Functions directory: this must be sourced before running the Main_Random.R script. 

Please direct any questions to amolstad@ufl.edu. 
