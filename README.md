The folder named "Functions" contain the functions used in some of the other scripts. 
It is therefore important that when running the scripts that this folder is added to the path where the script is being run. Otherwise an error will occur.
Scripts with "NewmarkStatic" in their name are functions using the Newton-Raphson algorithm.
Script with "V_physical", "V_StandardBasis" or "V_TaylorBasis" are functions using the CDM and either the full physical coordinates,
the standard basis or the Taylor basis in the calculations.

The folder named "Two-bar simple truss" contains the script used to find the results regarding the truss structure.
Scripts named "Non_linear_CDM_constantc" or "Non_linear_Newmark_constantC" are the scripts for analysing the two-bar simple truss structure.
with either the CDM or the Newmark algorithm. "simple_1_dof_nnlinear_static is used to inspect the static aspect of the Two-bar structure.

The folder named "Nonlinear cable" contains the script used to find the results for the nonlinear cable structure. 
The script with "standardBasis" in their name are used to simmulate the responses using the standard basis 
and the ones with "TaylorBasis" in their name are used to simmulate the responses using the Taylor basis. 
If "physical" is in the script name it is used to calculate the physical response. 
The physical response is also calculated in some of the script with "StandardBasis" or "TaylorBasis". 
Files with "Stab" are files used to analyse the stability of the structure.

The folder named "Bridge section" contains the scripts used to calculate the results regarding the bridge section structure. 
For the scripts with "error" in their name it is important that the Folder "Data" is added to the path for which the script is run.
