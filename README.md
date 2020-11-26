# ChromeBat
An metahueristic approach to solving the 3D Genome Reconstruction Problem

# Usage
Enter either command

>python ChromeBat.py contact_matrix

>python ChromeBat.py contact_matrix parameter_file

where contact_matrix is text file representing a square contact matrix derived from a HiC experiement.
If no parameter_file is passed the algorithm will use default parameters.
The parameter_file should be formatted like the parameters.txt.

# Comments

The approach is based on the bat algorithm "A New Metaheuristic Bat-Inspired Algorithm, in: Nature Inspired Cooperative Strategies for Optimization".
Additionally, I used https://github.com/buma/BatAlgorithm as reference and starting point.
The algorithm features numerous hyperparameters whose effects I describe in the parameters.txt.
In general, generations and num_bats are the most important parameters.
I found the best results with 10,000 generations and 30 bats, which took around 5 minutes on my computer.
