# ChromeBat
An metahueristic approach to solving the 3D Genome Reconstruction Problem

# Usage
Enter either command

>python ChromeBat.py contact_matrix

>python ChromeBat.py contact_matrix parameter_file

where contact_matrix is text file representing a square contact matrix derived from a HiC experiement.
If no parameter_file is passed the algorithm will use default parameters.


# Parameters
All parameters must be specified in a text file like in the parameters.txt example.
The parameters include
+ ouput_file: This is the name both of the outputted files will have.
+ alpha: This is the conversion factor used to convert the contact matrix to a distance matrix. If no alpha value is found a search across alphas in [0,1] will be performed
+ num_bats: How many bats the algorithm will simulate.
+ generations: How many iterations the algorithm will perform.
+ min_freq: The minimuim frequency that a bat can have. Low frequency means a bat will explore more then exploit.
+ max_freq: The maximum frequency a bat will be simulated with. Bats with high frequency expoit more then explore.
+ volume: A bat's volume determines how willing it is to accept new solutions. A loud bat will accept solutions with high probability.
+ pulse: When a bat pulses it teleports to the current best known solution. High pulse means it teleports with high probability. 
+ perturbation: This determines the size of the random walk a bat takes after pulsing
+ structs: This is how many structures in addition to the alpha search structures the algorithm should generate. If an alpha search is performed only the structures generated with the optimal alpha will be written to files.


# Output
Chromebat.py will produce two files per structure generated with desired alpha whose names are specified by the output_file parameter in the parameter file.

These files are
  + output_fileX.pdb : contains the model that may be visualized using PyMol
  + output_fileX.log : contains the input file name, Spearman's and Pearson's Correlation Coeffiecents and the Root Mean Squared Error

  Where X indicates the a structures number. For example if structs=2 and an alpha search is performed 3 structures with the optimal alpha will be generated, and 3 pairs of .pdb,.log files will be generated.


