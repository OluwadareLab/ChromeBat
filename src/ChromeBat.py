#ChromeBat
#Written By Brandon Collins for CS5435
#11/24/2020


import argparse
import numpy as np
import random
import re
import scipy.stats
import sklearn.metrics
import time
from multiprocessing import Pool
from multiprocessing import cpu_count
from itertools import repeat

#regular expressions
comment_re=re.compile("#.*")
numeric_re=re.compile('\d+(\.\d+)?')

#accepts a canidate solution(xyz coordinates of all loci) can be have shape (3n,) or (n,3)
#returns the distance matrix between all of said loci
def sol2distVec(sol):
    #Ensuring xyz_col to have shape (n,3)
    if len(np.shape(sol))==1:
        xyz_col=np.array(np.split(sol,len(sol)/3))
    else:
        xyz_col=sol
    n=len(xyz_col)
    #we construct rows_m,cols_m to have shape (n,n,3)
    #
    rows=np.reshape(xyz_col,(n,1,3))
    rows_m=np.tile(rows,(1,n,1))
    cols=np.reshape(xyz_col,(1,n,3))
    cols_m=np.tile(cols,(n,1,1))
    return np.sum((cols_m-rows_m)**2,2)**0.5

#accepts contact matrix and scaling parameter alpha
#returns the distance matrix
def contactToDistance(contact,alpha=.5):
	with np.errstate(divide="ignore"):
		matrix=1/(contact**alpha)
		#np.fill_diagonal(matrix,0)
	return matrix



#this is the optimization is performed
#initilize a bat object and then call fly()
#parameters are described on github
class Bat:
	params=set(["matrix","structs","alpha","func","num_bats","upper_bound","lower_bound","min_freq","max_freq","volume","generations","pulse","perturbation"])
	def __init__(self,matrix,num_bats=20,lower_bound=0,upper_bound=3,min_freq=0,max_freq=2,volume=0.5,generations=100,pulse=0.5,perturbation=0.01):
		self.lower=lower_bound 
		self.upper=upper_bound 
		self.freq_min=min_freq
		self.freq_max=max_freq
		self.D=len(matrix)*3
		self.matrix=matrix
		self.generations=int(generations) #how many iterations of the algorithm
		self.perturbation=perturbation # this dictates the size of a bats random walk after pulsing
		# self.loss=func
		self.loudness=volume
		self.num_bats=int(num_bats)
		self.best_fit=np.Inf

		self.pulse=pulse
		self.round=lambda x:0 if x>=self.pulse else 1
		self.freq=np.zeros(self.num_bats)
		self.sols=self.lower+(self.upper-self.lower)*np.random.rand(self.num_bats,self.D)
		self.velo=np.zeros((self.num_bats,self.D)) #bats velocity vector with shape (num_bats,D), where D is the number of dimensions 
		self.fitness=np.apply_along_axis(self.loss,1,self.sols)
		self.evalBats()

	#this is the loss function
	#It is a sum of the differnce between the target distance matrix and the once induced from
	# the proposed solution xyz coordinates
	def loss(self,sol):
		key=self.matrix
		sum_matrix=sol2distVec(sol)
		loss_matrix=key-sum_matrix
		loss_matrix=np.nan_to_num(loss_matrix,copy=True,posinf=0)
		return np.sum(np.absolute(np.tril(loss_matrix)))
	#this method assumes self.fitness corrospondes to self.sol
	#sets self.best_sol,self.best_fitness accordingly
	def evalBats(self):
		best_bat=np.argmin(self.fitness)
		if self.fitness[best_bat]<=self.best_fit:
			self.best_sol=self.sols[best_bat]
			self.best_fit=self.fitness[best_bat]
	#this is the actual bat algorithm
	#vectorized implementation has some redunant computation, but almost all computation time is spent in oracle calls
	def fly(self):
		for g in range(self.generations):
			#bat's mobility logic
			self.freq=self.freq_min+(self.freq_max-self.freq_min)*np.random.rand(self.num_bats)
			self.velo=self.velo+(self.sols-self.best_sol)*self.freq[:,np.newaxis]
			self.temp_sols=self.sols+self.velo
			self.temp_sols=np.clip(self.temp_sols,self.lower,self.upper)

			#bat's pulse logic
			self.pulse_array=(np.array([self.pulse]*self.num_bats)>=np.random.rand(self.num_bats)).astype(int)
			self.flip_pulse_array=(self.pulse_array+1)%2
			self.temp_sols=(self.temp_sols*self.flip_pulse_array[:,np.newaxis])+(np.tile(self.best_sol,(self.num_bats,1))*self.pulse_array[:,np.newaxis]+self.perturbation*np.random.normal(0,1,(self.num_bats,self.D)))
			self.temp_sols=np.clip(self.temp_sols,self.lower,self.upper)

			#oracle calls
			self.temp_fitness=np.apply_along_axis(self.loss,1,self.temp_sols)

			#accepting soln logic
			self.accepted_sols=(self.temp_fitness<=self.fitness).astype(int)*(np.array([self.loudness]*self.num_bats)>=np.random.rand(self.num_bats)).astype(int)
			self.declined_sols=(self.accepted_sols+1)%2
			self.sols=(self.sols*self.declined_sols[:,np.newaxis])+(self.temp_sols*self.accepted_sols[:,np.newaxis])
			self.fitness=(self.fitness*self.declined_sols)+(self.temp_fitness*self.accepted_sols)
			self.evalBats()
		return self.best_sol


#accepts symmetric square matrix a
#returns a with zero adjancent bins replaced with an average adjancency value
def adjancenyPreprocess(a):
	a=a.copy() #method is nondestructive
	adj_diagonal=np.diagonal(a,1).copy() #this is the diagonal of the matrix off by 1
	mean_adj=sum(adj_diagonal)/len(adj_diagonal[adj_diagonal != 0]) #mean of all values on off diagonal not including 0s
	adj_diagonal[adj_diagonal == 0] = mean_adj #this is corrected off diagonal
	a.flat[1::len(a)+1]=adj_diagonal #this fixes the diagonal to the right of the main diagonal
	a.flat[len(a)::len(a)+1]==adj_diagonal #this fixes the diagonal to the left of the main diagonal
	return a

#used for making the pdb files
#pads the input string so it is size long
#left=True pads on the left False on the right
def pad(value,size,left=False):
	value=str(value)
	space_to_add=size-len(value)
	if left:
		return space_to_add*" "+value
	else:
		return value+space_to_add*" "
#writes pdb based on input solution and an outfile name
def outputPdb(input_sol,outfile=None):
	if outfile is None:
		outfile="bat.pdb"
	else:
		outfile+=".pdb"
	sol=np.array(input_sol)
	if len(np.shape(sol))==1:
		sol=np.array(np.split(sol,len(sol)/3))
	if len(np.shape(sol))!=2:
		raise ValueError(f"Invalid solution shape {np.shape(input_sol)}")

	out_string="pdb carefully constructed by Brandon\n"
	for i in range(len(sol)):
		x,y,z=sol[i]
		out_string+="ATOM  "+pad(i+1,5,True)+"   " \
			+"CA MET "+pad("B"+str(i+1),6)+"   " \
			+pad(format(x,".3f"),8,True)+pad(format(y,".3f"),8,True) \
			+pad(format(z,".3f"),8,True)+"  " \
			+"0.20 10.00\n"
	connect_string="\n".join([f"CONECT"+pad(i+1,5,True)+pad(i+2,5,True) for i in range(len(sol)-1)])
	out_string+=connect_string+"\n"+"END\n"

	f=open(outfile,"w")
	f.write(out_string)
	f.close()
#this does all the command line processing
#return (contact matrix,alpha,bat algo parameters,outfile,structs)
def processInput(input_fname,parameter_fname=None):
	contact=np.genfromtxt(input_fname)
	alpha=None
	structs=1
	outfile=None
	if parameter_fname is None:
		return (contact,alpha,dict(),outfile,structs)
	pfile=open(parameter_fname,"r")
	plines=pfile.readlines()
	param_dict=dict()
	
	for line in plines:
		line=re.sub(comment_re,"",line)
		line="".join(line.split())
		line=line.replace("\t","")
		if line=="":
			continue
		param,arg=line.split("=")
		if param=="output_file":
			outfile=arg
			continue
		elif param not in Bat.params:
			raise ValueError(f"{param} is not a valid parameter")
		elif param=="structs":
			structs=int(arg)
		elif numeric_re.match(arg) is None:
			raise ValueError(f"{arg} must be numeric")
		elif param=="alpha":
			alpha=float(arg)
		else:
			param_dict[param]=float(arg)

	# distance_matrix=contactToDistance(contact,alpha)
	return (contact,alpha,param_dict,outfile,structs)

#accepts the proposed solution and the target distance matrix
#returns a string with PCC,SCC,RMSE
#string=False makes it return the tuple (rmse,scc,pcc)
def formatMetrics(sol,key,string=True):
	distance_matrix=sol2distVec(sol)
	# key=np.nan_to_num(key,copy=True,posinf=0)
	distance_list=[]
	key_list=[]
	for i in range(len(distance_matrix)):
		for j in range(i+1):
			if key[i,j]==np.inf:
				continue
			distance_list.append(distance_matrix[i,j])
			key_list.append(key[i,j])
	pearson=scipy.stats.pearsonr(distance_list,key_list)[0]
	spearman=scipy.stats.spearmanr(distance_list,key_list).correlation
	rmse=sklearn.metrics.mean_squared_error(key_list,distance_list)**0.5
	if string:
		metrics=f"AVG RMSE: {rmse}\nAVG Spearman correlation Dist vs. Reconstructed Dist: {spearman}\nAVG Pearson correlation Dist vs. Reconstructed Dist: {pearson}"
	else:
		metrics=(rmse,spearman,pearson)
	return metrics
#this writes the log file
#accepts the metric string from formatMetrics,input file name, and an output file name
#writes .log file accordingly
def outputLog(metrics,alpha,input_fname,out_fname=None,bat_params=None,runtime=None,searched_alphas=None,structs=None):
	if out_fname is None:
		out_fname="bat.log"
	else:
		out_fname+=".log"
	outstring=f"Input file: {input_fname}\nConvert factor: {alpha}\n"+metrics
	if runtime is not None:
		outstring+=f"\nRuntime: {runtime:.2f} seconds"
	if searched_alphas is not None:
		outstring+=f"\nPerformed Alpha Search: {searched_alphas}"
	if structs is not None:
		outstring+=f"\nGenerated Stuctures Count: {structs}"
	if bat_params is not None:
		outstring+=f"\n"+"\n".join([f"{key}={value}" for key,value in bat_params.items()])
	f=open(out_fname,"w")
	f.write(outstring)
	f.close()

#simple function for multiprocessing
def optimize(distance,params):
	bats=Bat(distance,**params)
	return bats.fly()


#driver code
if __name__=="__main__":


	#command line/input processing
	start_time=time.time()
	parser=argparse.ArgumentParser()
	parser.add_argument("contact_matrix",help="File name of a square contact matrix")
	parser.add_argument("params",nargs="?",default=None,help="File name of the parameters file")
	# parser.add_argument("parameters",default="parameters.txt",)	
	args=parser.parse_args()
	contact,alpha,param_dict,outfile,structs=processInput(args.contact_matrix,args.params)

	print("Running ChromeBat")

	#perform alpha search if no alpha found in the parameter file
	#sets alpha to best alpha found
	PROC_COUNT=cpu_count()
	searched_alphas=False
	if alpha is None:
		print("Performing Alpha Search...")

		#multiprocessing to search the alpha values
		searched_alphas=True
		alphas=[0.1,0.3,0.5,0.7,0.9,1]
		distance_m_list=[adjancenyPreprocess(contactToDistance(contact,a)) for a in alphas]
		pool = Pool(processes=PROC_COUNT)
		swarms=pool.starmap(optimize, zip(distance_m_list,repeat(param_dict)))
		pool.close()
		pool.join()

		#metrics computations and alpha selection
		spearmans=[formatMetrics(sol,distance_m_list[i],string=False)[1] for i,sol in enumerate(swarms)]
		best_index=np.argmax(spearmans)
		alpha=alphas[best_index]
		searched_sol=swarms[best_index]
		end_alphas_time=time.time()
		pad_len=len(str(max(alphas,key=lambda x:len(str(x)))))+7
		alpha_scc_string="\n".join([pad(alphas[i],pad_len,left=False)+str(spearmans[i]) for i in range(len(spearmans))])
		print("Search Results:\n"+pad("alpha",pad_len)+"SCC\n"+alpha_scc_string)
		print(f"Best Alpha={alpha}, found in {end_alphas_time-start_time:.2f} seconds")

	#geneterates structs structures 
	distance_m=contactToDistance(contact,alpha)
	if structs>0: 
		print(f"Generating {structs} more structures using Alpha={alpha}")
		distance_m_prepro=adjancenyPreprocess(distance_m)
		pool = Pool(processes=PROC_COUNT)
		swarms=pool.starmap(optimize, zip(repeat(distance_m_prepro,structs),repeat(param_dict)))
		pool.close()
		pool.join()
	else: # reset swarms variable if structs happened to be 0
		swarms=[]

	#metric computations on generated structures
	if searched_alphas: #reuse structure generated in the search portion
		swarms.append(searched_sol)

	spearmans=[formatMetrics(sol,distance_m,string=False)[1] for sol in swarms]
	best_index=np.argmax(spearmans)
	
	best_sol=swarms[best_index]*10 # output structures tended to be very small so now they are not
	metrics=formatMetrics(best_sol,distance_m)

	final_end_time=time.time()
	print(f"Done in {final_end_time-start_time:.2f} seconds total")

	#output results
	print("Writing best structure....")

	print(metrics)
	

	
	outputLog(metrics,alpha,args.contact_matrix,outfile,bat_params=param_dict,runtime=final_end_time-start_time,searched_alphas=searched_alphas,structs=structs)
	outputPdb(best_sol,outfile)
	print(f"{outfile}.log and {outfile}.pdb written!")

	





