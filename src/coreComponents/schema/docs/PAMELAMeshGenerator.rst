

================= ============ ======== ============================================================ 
Name              Type         Default  Description                                                  
================= ============ ======== ============================================================ 
fieldNamesInGEOSX string_array {}       Name of the fields within GEOSX                              
fieldsToImport    string_array {}       Fields to be imported from the external mesh file            
file              path         required path to the mesh file                                        
logLevel          integer      0        Log level                                                    
name              string       required A name is required for any non-unique nodes                  
scale             R1Tensor     {1,1,1}  Scale the coordinates of the vertices by given scale factors 
translate         R1Tensor     {0,0,0}  Translate the coordinates of the vertices by a given vector  
================= ============ ======== ============================================================ 


