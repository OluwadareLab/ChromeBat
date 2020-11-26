#GenomeBat
#Written By Brandon Collins for CS5435
#11/24/2020


import argparse
import numpy as np
import random
import re
import scipy.stats
import sklearn.metrics

#regular expressions
comment_re=re.compile("#.*")
numeric_re=re.compile('\d+(\.\d+)?')

#accepts a canidate solution(xyz coordinates of all loci)
#returns the distance matrix between all of said loci
def sol2dist(sol):
	xyz_col=np.array(np.split(sol,len(sol)/3))
	sum_matrix=np.zeros((len(xyz_col),len(xyz_col)))
	for i in [1,2,3]:
		b=np.tile(xyz_col[:,i-1:i].T,(len(xyz_col),1))
		c=np.tile(xyz_col[:,i-1:i],(1,len(xyz_col)))
		d=b-c
		sum_matrix=sum_matrix+(d*d)
	sum_matrix=sum_matrix **0.5
	return sum_matrix

#accepts contact matrix and scaling parameter alpha
#returns the distance matrix
def contactToDistance(contact,alpha=.5):
	with np.errstate(divide="ignore"):
		matrix=1/(contact**alpha)
	return matrix



#this is the optimization algorithm
#it is based on the paper and https://github.com/buma/BatAlgorithm
#call it by using fly()
class Bat:
	params=set(["matrix","alpha","func","num_bats","upper_bound","lower_bound","min_freq","max_freq","volume","generations","pulse"])
	def __init__(self,matrix,num_bats=20,lower_bound=0,upper_bound=3,min_freq=0,max_freq=2,volume=0.5,generations=100,pulse=0.5):
		self.lower=lower_bound
		self.upper=upper_bound
		self.freq_min=min_freq
		self.freq_max=max_freq
		self.D=len(matrix)*3
		self.matrix=matrix
		self.generations=int(generations)
		# self.loss=func
		self.loudness=volume
		self.num_bats=int(num_bats)
		self.best_fit=np.Inf

		self.pulse=pulse
		self.round=lambda x:0 if x>=self.pulse else 1
		self.freq=np.zeros(self.num_bats)
		self.sols=self.lower+(self.upper-self.lower)*np.random.rand(self.num_bats,self.D)
		self.velo=np.zeros((self.num_bats,self.D))
		self.fitness=np.apply_along_axis(self.loss,1,self.sols)
		self.evalBats()

	#this is the loss function
	#It is a sum of the differnce between the target distance matrix and the once induced from
	# the proposed solution xyz coordinates
	def loss(self,sol):
		key=self.matrix
		sol=np.array(sol)
		xyz_col=np.array(np.split(sol,len(sol)/3))
		sum_matrix=np.zeros(len(xyz_col))
		for i in [1,2,3]:
			b=np.tile(xyz_col[:,i-1:i].T,(len(xyz_col),1))
			c=np.tile(xyz_col[:,i-1:i],(1,len(xyz_col)))
			d=np.tril(b-c)
			sum_matrix=sum_matrix+(d*d)
		sum_matrix=sum_matrix **0.5
		loss_matrix=key-sum_matrix
		loss_matrix=np.nan_to_num(loss_matrix,copy=True,posinf=0)
		return np.sum(np.absolute(np.tril(loss_matrix)))
	def evalBats(self):
		best_bat=np.argmin(self.fitness)
		if self.fitness[best_bat]<self.best_fit:
			self.best_sol=self.sols[best_bat].copy()
			self.best_fit=self.fitness[best_bat]
	def fly(self):
		for g in range(self.generations):
			self.freq=self.freq_min+(self.freq_max-self.freq_min)*np.random.rand(self.num_bats)
			self.velo=self.velo+(self.sols-self.best_sol)*self.freq[:,np.newaxis]
			temp_sols=self.sols+self.velo
			temp_sols=np.clip(temp_sols,self.lower,self.upper)

			#this could become individuals bats pulse rate
			pulse_array=(np.array([self.pulse]*self.num_bats)>=np.random.rand(self.num_bats)).astype(int)
			flip_pulse_array=(pulse_array+1)%2
			temp_sols=(temp_sols*flip_pulse_array[:,np.newaxis])+(np.tile(self.best_sol,(self.num_bats,1))*pulse_array[:,np.newaxis]+0.001*np.random.normal(0,1,(self.num_bats,self.D)))
			temp_sols=np.clip(temp_sols,self.lower,self.upper)

			temp_fitness=np.apply_along_axis(self.loss,1,temp_sols,args=self.matrix)

			#this could become individual bats volume
			accepted_sols=(temp_fitness<=self.fitness).astype(int)*(np.array([self.loudness]*self.num_bats)>=np.random.rand(self.num_bats)).astype(int)
			declined_sols=(accepted_sols+1)%2
			self.sols=(self.sols*declined_sols[:,np.newaxis])+(temp_sols*accepted_sols[:,np.newaxis])
			self.fitness=(self.fitness*declined_sols)+(temp_fitness*accepted_sols)
			self.evalBats()
		return self.best_sol
	def fly2(self):
		S=np.zeros((self.num_bats,self.D))
		for t in range(self.generations):
			self.freq = np.array([self.freq_min]*self.num_bats) + (np.array([self.freq_max]*self.num_bats) - np.array([self.freq_min]*self.num_bats)) * np.random.rand(self.num_bats)
			self.velo = self.velo + (self.sols - self.best_sol) * self.freq[:,np.newaxis]
			S = self.sols + self.velo
			S=np.clip(S,self.lower,self.upper)
			for i in range(self.num_bats):
				if np.random.random_sample() >= self.pulse:
					for j in range(self.D):
						S[i][j] = self.best_sol[j] + 0.01 * np.random.normal(0,1)
						S[i][j] = max(self.lower,S[i][j])
						S[i][j] = min(self.upper,S[i][j])
				Fnew = self.loss(S[i])
				rnd = np.random.random_sample()
				if (Fnew <= self.fitness[i]) and (rnd < self.loudness):
					for j in range(self.D):
						self.sols[i][j] = S[i][j]
					self.fitness[i] = Fnew

				if Fnew <= self.best_fit:
					for j in range(self.D):
						self.best_sol[j] = S[i][j]
					self.best_fit = Fnew
		return self.best_sol
#used for making the pdb files
#pads the input string so it is size long
#left=True pads on the left False on the Right
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
		out_string+="ATOM  "+pad(i+1,5,True)+" " \
			+"CA MET "+pad("A"+str(i+1),5) \
			+pad(format(x,".3f"),8,True)+pad(format(y,".3f"),8,True) \
			+pad(format(z,".3f"),8,True) \
			+" 0.20 10.00\n"
	connect_string="\n".join([f"CONECT"+pad(i+1,5,True)+pad(i+2,5,True) for i in range(len(sol)-1)])
	out_string+=connect_string+"\n"+"END\n"

	f=open(outfile,"w")
	f.write(out_string)
	f.close()
#this does all the command line processing
#return (contact matrix,alpha,bat algo parameters,outfile)
def processInput(input_fname,parameter_fname=None):
	contact=np.genfromtxt(input_fname)
	alpha=0.5
	outfile=None
	if parameter_fname is None:
		return (contact,alpha,dict(),outfile)
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
		if param not in Bat.params:
			raise ValueError(f"{param} is not a valid parameter")
		if numeric_re.match(arg) is None:
			raise ValueError(f"{arg} must be numeric")
		if param=="alpha":
			alpha=float(arg)
		else:
			param_dict[param]=float(arg)

	distance_matrix=contactToDistance(contact,alpha)
	return (distance_matrix,alpha,param_dict,outfile)

#accepts the proposed solution and the target distance matrix
#returns a string with PCC,SCC,RMSE
def formatMetrics(sol,key):
	distance_matrix=sol2dist(sol)
	key=np.nan_to_num(key,copy=True,posinf=0)
	pearson=scipy.stats.pearsonr(distance_matrix.flatten(),key.flatten())[0]
	spearman=scipy.stats.spearmanr(distance_matrix.flatten(),key.flatten()).correlation
	rmse=sklearn.metrics.mean_squared_error(key.flatten(),distance_matrix.flatten())**0.5
	metrics=f"AVG RMSE: {rmse}\nAVG Spearman correlation Dist vs. Reconstructed Dist: {spearman}\nAVG Pearson correlation Dist vs. Reconstructed Dist: {pearson}"
	return metrics
#this writes the log file
#accepts the metric string from formatMetrics,input file name, and an output file name
def outputLog(metrics,alpha,input_fname,out_fname=None):
	if out_fname is None:
		out_fname="bat.log"
	else:
		out_fname+=".log"
	outstring=f"Input file: {input_fname}\nConvert factor: {alpha}\n"+metrics
	f=open(out_fname,"w")
	f.write(outstring)
	f.close()



#driver code
if __name__=="__main__":
	parser=argparse.ArgumentParser()
	parser.add_argument("contact_matrix",help="File name of a square contact matrix")
	parser.add_argument("params",nargs="?",default=None,help="File name of the parameters file")
	# parser.add_argument("parameters",default="parameters.txt",)	
	args=parser.parse_args()
	distance,alpha,param_dict,outfile=processInput(args.contact_matrix,args.params)
	bats=Bat(distance,**param_dict)
	best_sol=bats.fly2()
	metrics=formatMetrics(best_sol,distance)
	print(metrics)
	outputLog(metrics,alpha,args.contact_matrix,outfile)
	outputPdb(best_sol,outfile)






