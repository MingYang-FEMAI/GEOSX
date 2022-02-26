

============== =============================================================== ======= ========================================================================================= 
Name           Type                                                            Default Description                                                                               
============== =============================================================== ======= ========================================================================================= 
method         geosx_LinearSolverParameters_Multiscale_Coarsening_Metis_Method kway    METIS partitioning method, one of: ``kway``, ``recursive``                                
minCommonNodes integer                                                         3       Minimum number of nodes shared between two cells when constructing the connectivity graph 
ufactor        integer                                                         30      METIS ufactor parameter, affects partitioning balance/edgecut tradeoff                    
============== =============================================================== ======= ========================================================================================= 


