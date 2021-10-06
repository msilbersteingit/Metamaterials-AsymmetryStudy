# Metamaterials-AsymmetryStudy
Codes used for the paper covering a study of symmetry and asymmetry in lattice-based metamaterials

The file primary_analysis_code.m contains all the code used in the body of the paper, while secondary_analysis_code.m additionally contains functions used to generate plots in the supplementary information of the paper.  The files design_sets_3x3.mat and design_sets_5x5.mat include all of the lattice designs (in the form of connectivity arrays) used in our study.  Using this code package requires the installation of (1) MATLAB and (2) ANSYS.  In order to use either primary_analysis_code.m or secondary_analysis_code.m (and the critical functions within), all of the MATLAB subfunctions and text files found in this repository must be included in the same directory on your machine.

Included in primary_analysis_code.m is the implementation of the finite-element beam model (written in ANSYS APDL) used to study each lattice design.  In order to use this code, the filepaths on lines 7,28,39,50,58,66,74,82, and 297 of Truss2D_NxN.txt must be edited to match the locations of these files on your machine.  Additionally, line 2861 of primary_analysis_code.m must be edited to reflect the filepath of ANSYS on your machine (the location should be similar to the filepath shown), and line 2873 must be edited to reflect the ANSYS license type you are using (aa_r in the existing code indicates a Research license).  In secondary_analysis_code.m these two lines are 4275 and 4287, respectively.


