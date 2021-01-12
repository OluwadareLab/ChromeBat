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

#regular expressions
comment_re=re.compile("#.*")
numeric_re=re.compile('\d+(\.\d+)?')

#accepts a canidate solution(xyz coordinates of all loci)
#returns the distance matrix between all of said loci
def sol2dist(sol):
	if len(np.shape(sol))==1:
		xyz_col=np.array(np.split(sol,len(sol)/3))
	else:
		xyz_col=sol
	
	matrix=np.zeros((len(xyz_col),len(xyz_col)))
	for i in range(len(xyz_col)):
		x_i,y_i,z_i=xyz_col[i]
		for j in range(len(xyz_col)):
			x_j,y_j,z_j=xyz_col[j]
			matrix[i,j]=((x_i-x_j)**2+(y_i-y_j)**2+(z_i-z_j)**2)**.5
	return matrix

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



#this is the optimization algorithm
#it is based on the paper and https://github.com/buma/BatAlgorithm
#call it by using fly()
class Bat:
	params=set(["matrix","fly_ver","alpha","func","num_bats","upper_bound","lower_bound","min_freq","max_freq","volume","generations","pulse","perturbation"])
	def __init__(self,matrix,num_bats=20,lower_bound=0,upper_bound=3,min_freq=0,max_freq=2,volume=0.5,generations=100,pulse=0.5,perturbation=0.01):
		self.lower=lower_bound
		self.upper=upper_bound
		self.freq_min=min_freq
		self.freq_max=max_freq
		self.D=len(matrix)*3
		self.matrix=matrix
		self.generations=int(generations)
		self.perturbation=perturbation
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
		sum_matrix=sol2distVec(sol)
		loss_matrix=key-sum_matrix
		loss_matrix=np.nan_to_num(loss_matrix,copy=True,posinf=0)
		return np.sum(np.absolute(np.tril(loss_matrix)))
	# def loss(self,sol):
	# 	key=self.matrix
	# 	sol=np.array(sol)
	# 	xyz_col=np.array(np.split(sol,len(sol)/3))
	# 	sum_matrix=np.zeros(len(xyz_col))
	# 	for i in [1,2,3]:
	# 		b=np.tile(xyz_col[:,i-1:i].T,(len(xyz_col),1))
	# 		c=np.tile(xyz_col[:,i-1:i],(1,len(xyz_col)))
	# 		d=np.tril(b-c)
	# 		sum_matrix=sum_matrix+(d*d)
	# 	sum_matrix=sum_matrix **0.5
	# 	loss_matrix=key-sum_matrix
	# 	loss_matrix=np.nan_to_num(loss_matrix,copy=True,posinf=0)
	# 	return np.sum(np.absolute(np.tril(loss_matrix)))
	def evalBats(self):
		best_bat=np.argmin(self.fitness)
		if self.fitness[best_bat]<=self.best_fit:
			self.best_sol=self.sols[best_bat]
			self.best_fit=self.fitness[best_bat]
	def fly(self):
		for g in range(self.generations):
			self.freq=self.freq_min+(self.freq_max-self.freq_min)*np.random.rand(self.num_bats)
			self.velo=self.velo+(self.sols-self.best_sol)*self.freq[:,np.newaxis]
			self.temp_sols=self.sols+self.velo
			self.temp_sols=np.clip(self.temp_sols,self.lower,self.upper)

			#this could become individuals bats pulse rate
			self.pulse_array=(np.array([self.pulse]*self.num_bats)>=np.random.rand(self.num_bats)).astype(int)
			self.flip_pulse_array=(self.pulse_array+1)%2
			self.temp_sols=(self.temp_sols*self.flip_pulse_array[:,np.newaxis])+(np.tile(self.best_sol,(self.num_bats,1))*self.pulse_array[:,np.newaxis]+self.perturbation*np.random.normal(0,1,(self.num_bats,self.D)))
			self.temp_sols=np.clip(self.temp_sols,self.lower,self.upper)

			self.temp_fitness=np.apply_along_axis(self.loss,1,self.temp_sols)

			#this could become individual bats volume
			self.accepted_sols=(self.temp_fitness<=self.fitness).astype(int)*(np.array([self.loudness]*self.num_bats)>=np.random.rand(self.num_bats)).astype(int)
			self.declined_sols=(self.accepted_sols+1)%2
			self.sols=(self.sols*self.declined_sols[:,np.newaxis])+(self.temp_sols*self.accepted_sols[:,np.newaxis])
			self.fitness=(self.fitness*self.declined_sols)+(self.temp_fitness*self.accepted_sols)
			self.evalBats()
		return self.best_sol
	def fly2(self):
		self.S=np.zeros((self.num_bats,self.D))
		for t in range(self.generations):
			self.freq = np.array([self.freq_min]*self.num_bats) + (np.array([self.freq_max]*self.num_bats) - np.array([self.freq_min]*self.num_bats)) * np.random.rand(self.num_bats)
			self.velo = self.velo + (self.sols - self.best_sol) * self.freq[:,np.newaxis]
			self.S = self.sols + self.velo
			self.S=np.clip(self.S,self.lower,self.upper)
			for i in range(self.num_bats):
				if np.random.random_sample() >= self.pulse:
					for j in range(self.D):
						self.S[i][j] = self.best_sol[j] + 0.01 * np.random.normal(0,1)
						self.S[i][j] = max(self.lower,self.S[i][j])
						self.S[i][j] = min(self.upper,self.S[i][j])
				self.Fnew = self.loss(self.S[i])
				self.rnd = np.random.random_sample()
				if (self.Fnew <= self.fitness[i]) and (self.rnd < self.loudness):
					for j in range(self.D):
						self.sols[i][j] = self.S[i][j]
					self.fitness[i] = self.Fnew

				if self.Fnew <= self.best_fit:
					for j in range(self.D):
						self.best_sol[j] = self.S[i][j]
					self.best_fit = self.Fnew
		return self.best_sol
	def PSO(self):
		lr=0.2 #learning rate
		#temp_sols is bats current location
		#self.sols is a bats best known location
		temp_sols=self.sols.copy()
		self.velo=np.random.uniform(self.lower-self.upper,self.upper-self.lower,(self.num_bats,self.D))
		for g in range(self.generations):
			global_mult_m=np.random.uniform(0,1,(self.num_bats,self.D))#r_g from wiki
			local_mult_m=np.random.uniform(0,1,(self.num_bats,self.D))#r_p from wiki
			#we omit omega,phi_p,phi_g from wiki for now
			self.velo=self.velo+local_mult_m*(self.sols-temp_sols)+global_mult_m*(np.tile(self.best_sol,(self.num_bats,1))-temp_sols)
			temp_sols=temp_sols+self.velo*lr
			temp_fitness=np.apply_along_axis(self.loss,1,temp_sols)
			accepted_sols=(temp_fitness<self.fitness).astype(int)
			declined_sols=(accepted_sols+1)%2
			self.sols=(self.sols*declined_sols[:,np.newaxis])+(temp_sols*accepted_sols[:,np.newaxis])
			self.fitness=(self.fitness*declined_sols)+(temp_fitness*accepted_sols)
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
	a.flat[len(a)::len(a)+1]==adj_diagonal #this fixes the diagonal to theleft of the main diagonal
	return a

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
#return (contact matrix,alpha,bat algo parameters,outfile,fly_ver)
def processInput(input_fname,parameter_fname=None):
	contact=np.genfromtxt(input_fname)
	alpha=None
	fly_ver="1"
	outfile=None
	if parameter_fname is None:
		return (contact,alpha,dict(),outfile,fly_ver)
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
		elif param=="fly_ver":
			fly_ver=arg
		elif numeric_re.match(arg) is None:
			raise ValueError(f"{arg} must be numeric")
		elif param=="alpha":
			alpha=float(arg)
		else:
			param_dict[param]=float(arg)

	# distance_matrix=contactToDistance(contact,alpha)
	return (contact,alpha,param_dict,outfile,fly_ver)

#accepts the proposed solution and the target distance matrix
#returns a string with PCC,SCC,RMSE
def formatMetrics(sol,key):
	distance_matrix=sol2dist(sol)
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
	metrics=f"AVG RMSE: {rmse}\nAVG Spearman correlation Dist vs. Reconstructed Dist: {spearman}\nAVG Pearson correlation Dist vs. Reconstructed Dist: {pearson}"
	return metrics
#this writes the log file
#accepts the metric string from formatMetrics,input file name, and an output file name
def outputLog(metrics,alpha,input_fname,out_fname=None,bat_params=None,runtime=None):
	if out_fname is None:
		out_fname="bat.log"
	else:
		out_fname+=".log"
	outstring=f"Input file: {input_fname}\nConvert factor: {alpha}\n"+metrics
	if runtime is not None:
		outstring+=f"\nRuntime: {runtime:.2f} seconds"
	if bat_params is not None:
		outstring+=f"\n"+"\n".join([f"{key}={value}" for key,value in bat_params.items()])
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
	contact,alpha,param_dict,outfile,fly_ver=processInput(args.contact_matrix,args.params)
	if alpha is None:
		for alpha in [0.1,0.3,0.5,0.7,0.9,1.0]:
			start_time=time.time()
			distance=contactToDistance(contact,alpha)
			distance_prepro=adjancenyPreprocess(distance)
			bats=Bat(distance_prepro,**param_dict)
			if fly_ver=="2":
				best_sol=bats.fly2()
			elif fly_ver=="PSO":
				best_sol=bats.PSO()
			else:
				best_sol=bats.fly()
			best_sol*=10 #resizing solotion as it is quite small
			metrics=formatMetrics(best_sol,distance)
			end_time=time.time()
			print(f"alpha={alpha}\n"+metrics)
			print(f"Execution Time: {end_time-start_time:.2f} seconds using fly_ver={fly_ver}")

			outfile_suffix="a"+(str(alpha)).strip(".")
			outputLog(metrics,alpha,args.contact_matrix,outfile+outfile_suffix,bat_params=param_dict,runtime=end_time-start_time)
			outputPdb(best_sol,outfile+outfile_suffix)
	else:
		start_time=time.time()
		distance=contactToDistance(contact,alpha)
		distance_prepro=adjancenyPreprocess(distance)
		bats=Bat(distance_prepro,**param_dict)
		if fly_ver=="2":
			best_sol=bats.fly2()
		elif fly_ver=="PSO":
			best_sol=bats.PSO()
		else:
			best_sol=bats.fly()
		best_sol*=10
		metrics=formatMetrics(best_sol,distance)
		end_time=time.time()
		print(metrics)
		print(f"Execution Time: {end_time-start_time:.2f} seconds using fly_ver={fly_ver}")
		outputLog(metrics,alpha,args.contact_matrix,outfile,bat_params=param_dict,runtime=end_time-start_time)
		outputPdb(best_sol,outfile)






