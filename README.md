# CondTensorDecomp
Simulation code for "Conditional probability tensor decompositions for multivariate categorical response regression" 

To run the simulations displayed in Figure 5, execute the Run_R3.sh file. This will initiate 2400 jobs, each of which is a single replicate. This will first call the R3_Random.R script, which generates data from the appropriate data generating model and then sources the Main_Random.R script. The Random suffix refers to the method of initialization for the B's. The only function needed is in the Functions directory. 
