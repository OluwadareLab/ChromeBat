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
+ alpha: This is the conversion factor used.
+ num_bats: How many bats the algorithm will simulate.
+ generations: How many iterations the algorithm will perform.
+ upper_bound: The upper bound of the search space.
+ lower_bound: The lower bound of the search space.
+ min_freq: The minimuim frequency that a bat can have. Low frequency means a bat will explore more then exploit.
+ max_freq: The maximum frequency a bat will be simulated with. Bats with high frequency expoit more then explore.
+ volume: A bat's volume determines how willing it is to accept new solutions. A loud bat will accept solutions with high probability.
+ pulse: When a bat pulses it teleports to the current best known solution. High pulse means it teleports with high probability. 


# Output

Chromebat.py will produce two files whose names are specified by the output_file parameter in the parameter file.

These files are
  + output_file.pdb : contains the model that may be visualized using PyMol
  + output_file.log : contains the input file name, Spearman's and Pearson's Correlation Coeffiecents and the Roor Mean Squared Error

# Comments

The approach is based on the bat algorithm "A New Metaheuristic Bat-Inspired Algorithm, in: Nature Inspired Cooperative Strategies for Optimization".
Additionally, I used https://github.com/buma/BatAlgorithm as reference and starting point.
The algorithm features numerous hyperparameters whose effects I describe in the parameters.txt.
In general, generations and num_bats are the most important parameters.
I found the best results with 10,000 generations and 30 bats, which took around 5 minutes on my computer.
