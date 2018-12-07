title: 
Multivariate multiple test procedures based on nonparametric copula estimation

authors: 
André Neumann(*), Taras Bodnar, Dietmar Pfeifer, and Thorsten Dickhaus
(*)Author responsible for the code (email: neumann@uni-bremen.de).


Section 2.2: 
Load the file "comparison_of_copula_approximation_methods.R". Then run the main 
function "BC_Comparison()". The results are also provided in the file 
"2017-05-31-BC_Comparison.RData", and can be plotted (and saved to file)
via the function "plot_results".

Section 4:
Load the file "simulation_study.R". Then run the main function 
"BC_SimulationStudy()". The computations can take very long (multiple days on a 
single computer). Consequently, the results are also provided in the files 
"2018-01-12-BC_SimulationStudy.RData" and "2018-01-12-BC_SimulationStudy.xlsx".

Section 5:
No source code is provided. The data can be loaded from the file 
"BC_Application_Data.RData".


list of configurations (sessionInfo()):
R version 3.4.3 (2017-11-30)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252    
    LC_MONETARY=German_Germany.1252
[4] LC_NUMERIC=C                    LC_TIME=German_Germany.1252    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base 

other attached packages:
 [1] copula_0.999-18   doRNG_1.6.6       rngtools_1.2.4    pkgmaker_0.22     
     registry_0.5      doParallel_1.0.11
 [7] iterators_1.0.9   foreach_1.4.4     mvtnorm_1.0-7     formatR_1.5      

loaded via a namespace (and not attached):
 [1] magrittr_1.5      ADGofTest_0.3     xtable_1.8-2      pspline_1.0-18    
     lattice_0.20-35   pcaPP_1.9-73     
 [7] stringr_1.2.0     tools_3.4.3       grid_3.4.3        digest_0.6.15     
     numDeriv_2016.8-1 Matrix_1.2-12    
[13] codetools_0.2-15  stringi_1.1.6     compiler_3.4.3    gsl_1.9-10.3      
     stats4_3.4.3      stabledist_0.7-1