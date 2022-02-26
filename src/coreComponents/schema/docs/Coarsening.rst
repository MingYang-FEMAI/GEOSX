

============= ================================================================ ======= ============================================================================================================ 
Name          Type                                                             Default Description                                                                                                  
============= ================================================================ ======= ============================================================================================================ 
minGlobalDof  globalIndex                                                      0       limit of coarsening across all ranks (i.e. trim the grid hierarchy globally)                                 
minLocalDof   localIndex                                                       0       Limit of coarsening on current rank (i.e. keep a local coarsening ratio of 1 once this problem size reached) 
partitionType geosx_LinearSolverParameters_Multiscale_Coarsening_PartitionType metis   Partition type for generating coarse aggregates. Available options are: ``metis``, ``rib``, ``cartesian``    
ratio         real64_array                                                     {0}     Coarsening ratio (number of fine cells per coarse cell)                                                      
Metis         node                                                             unique  :ref:`XML_Metis`                                                                                             
============= ================================================================ ======= ============================================================================================================ 


