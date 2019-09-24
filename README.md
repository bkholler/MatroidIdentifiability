# MatroidIdentifiability
The files in this repository are the supplementary files for "Identifiability in Phylogenetics Using Algebraic Matroids". 
Any reference to an algorithm or result is a reference to that paper. 
Below is a brief description of all of the files in this repository. 
To use any of these files one should download the whole repository and keep all of the files in one folder. 
All of the actual functions and algorithms are in PhylogeneticMatroids.m.

1. PhylogeneticMatroids.m
    - Contains all of the algorithms and methods described in "Identifiability in Phylogenetics Using Algebraic Matroids"
    - Also contains some basic methods to compute parameterizations of phylogenetic models and their jacobians
    
2. unrooted4, unrooted5, unrooted6
    - Each of these files is just a list of the unrooted n leaf trees (for n = 4, 5, 6 respectively) given by their splits


3. fourLeafOrbs, fiveLeafOrbs, sixLeafOrbs
    - A list of pairs of pairs of integers in the form {{i_1, i_2}, {i_3, i_4}} 
    where i_j corresponds to the i_jth tree in unrooted4, unrooted5, or unrooted6 respectively
    - The pairs of pairs of trees are orbit representatives for the orbits of the natural symmetric group action 
    on pairs of pairs trees (permuting the leaves of each tree) 
    - Essentially these files just encode the pairs of pairs of trees symmetry

4. certsCFN
    - Contains a list of the certificates for each pair of pairs of trees in the form {i, cert} 
    where i specifies the position of the pair of pairs in sixLeafOrbs
    - 

  
  
5. CFN_6Leaf_Mixtures.nb
    - This notebook file verifies that the sets in certsCFN are in fact certificates as described in Algorithm 3.3
    - Running the ParallelDo loop in full may take a long time. We ran it in batches of about 1000 on normal computers. 
    
6. K3P_4LeafMixtures.nb
    - This notebook is similar to CFN_6LeafMixtures but for 4 leaf K3P mixtures
    - The only other difference is that the certificates are included in the notebook instead of kept in a separate file

7. K2P_Networks.nb
    - This file contains the proof of Lemmas 5.7, 5.8 and 5.9 for the K2P model
    
8. K3P_Networks.nb
    - This file contains the proof of Lemmas 5.7, 5.8 and 5.9 for the K3P model
  
