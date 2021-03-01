------------------------------------------------------------------------------------------------------------------------------------
# ChromeBat: A Bio-Inspired Approach to Genome Reconstruction
------------------------------------------------------------------------------------------------------------------------------------
**OluwadareLab,**
**University of Colorado, Colorado Springs**

----------------------------------------------------------------------
**Developed by:** <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Brandon Collins<br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Department of Computer Science <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;University of Colorado, Colorado Springs <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Email: bcollin3@uccs.edu <br /><br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Philip Brown, PhD<br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Department of Computer Science <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;University of Colorado, Colorado Springs <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Email: pbrown2@uccs.edu 
		 <br /><br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Oluwatosin Oluwadare, PhD <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Department of Computer Science <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;University of Colorado, Colorado Springs <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Email: ooluwada@uccs.edu


----------------------------------------------------------------------

**1.	Content of folders:**
----------------------------------------------------------------------
* src: Source code and parameter file for ChromeBat
* results: output structures on GM06990 and GM12878

**3.	Input matrix file format:**
-----------------------------------------------------------

Square Matrix Input format: The square matrix is a whitespace seperated N by N intra-chromosomal contact matrix derived from Hi-C data, where N is the number of equal-sized regions of a chromosome.


**4.	Dependencies Installation:**
-----------------------------------------------------------

&nbsp;&nbsp;&nbsp;&nbsp;python - 3.8.5 <br />
&nbsp;&nbsp;&nbsp;&nbsp;scipy - 1.5.2 <br />
&nbsp;&nbsp;&nbsp;&nbsp;numpy - 1.19.2 <br />
&nbsp;&nbsp;&nbsp;&nbsp;scikit-learn - 0.23.2 <br />


**5.	Usage:**
-----------------------------------------------------------
>python ChromeBat.py contact_matrix parameter_file

where contact_matrix is text file representing a whitespace delimited square contact matrix derived from a HiC experiement.




**6.	Parameters:**
-----------------------------------------------------------
All parameters must be specified in a text file like in the parameters_heavy.txt example.
Default values of parameters_heavy.txt are given in [].
However, during the parameter search phase of the algoirthm the heavy parameter file will run 30 processes concurrently.
Because this is impractical to run on most local machines we also provide a parameters_light.txt file that will run only 6 processes concurrently.

Searched Parameters:
+ alpha: [0.1,0.3,0.5,0.7,0.9,1.0] This is the conversion factor used to convert the contact matrix to a distance matrix. If no alpha value is found a search across alphas in [0,1] will be performed
+ perturbation: [0.002,0.004,0.006,0.008,0.01] This determines the size of the random walk a bat takes after pulsing

These parameters may have comma delimited values in the parameter file, doing so will result in a search over all combinations of possible purturbation and alpha values.

Normal Parameters:
+ ouput_file: [bat] This is the name both of the outputted files will have.
+ num_bats: [10] How many bats the algorithm will simulate.
+ generations: [10] How many iterations the algorithm will perform.
+ min_freq: [0] The minimuim frequency that a bat can have. Low frequency means a bat will explore more then exploit.
+ max_freq: [0.1] The maximum frequency a bat will be simulated with. Bats with high frequency expoit more then explore.
+ volume: [0.9] A bat's volume determines how willing it is to accept new solutions. A loud bat will accept solutions with high probability.
+ pulse: [0.9] When a bat pulses it teleports to the current best known solution. High pulse means it teleports with high probability. 
+ structs: [10] This is how many structures in addition to the alpha search structures the algorithm should generate. If an alpha search is performed only the structures generated with the optimal alpha will be written to files.

These parameters may only have 1 value.

**7.	Output:**
-----------------------------------------------------------
Chromebat.py will produce two files per structure generated whose names are specified by the output_file parameter in the parameter file.

These files are
  + output_fileX.pdb : contains the model that may be visualized using PyMol
  + output_fileX.log : contains the input file name, Spearman's and Pearson's Correlation Coeffiecents and the Root Mean Squared Error

  Where X indicates the a structures number. For example if structs=2 and an alpha/perturbation search is performed 3 structures with the optimal alpha/perturbation will be generated, and 3 pairs of .pdb,.log files will be generated.


