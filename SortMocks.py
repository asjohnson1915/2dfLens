#! /usr/bin/env python

# Purpose of code:
# Take Kids source and lens simulations and put into format that can be used with athena..
# Take redshift slices of the simulations
# write out new files, and write list of new file names.. 

## Check trying to updo file modifications

# Notes:
# Format wanted for sources: standard  x y e1 e2 w. 
# Format wanted for lens : standard    x y e1 e2 w -- set e1 e2 and w to zero.

# Current format for Kids sources: x, y, eobs1, eobs2, weight, z
# Current format for lenses DA,REC,z

## Run using: mpiexec -n 8 python SortMocks.py

from __future__ import print_function

import os
import pyfits
import re
import pandas as pd 
from optparse import OptionParser 
import numpy as np 
import subprocess
import math
import errno
import shlex
import sys
from mpi4py import MPI

class Kids_mock_S:
	""" Used to store all mock data for sources.
	"""
	def __init__(self,file_name):
		""" Read in kids simulation data.. 
		"""
		self.error_sd = 0.29  					  ## Size of error to add to simulation measurement of shear.
		self.infile = file_name.replace('\n','') # remove end of line indicator
		self.line_num  = 0 

		for line in open(self.infile):
	    		if not line.strip().startswith("#"):
	    			self.line_num = int(line)
	    			print('line number is',self.line_num)
	    			break
     		  	else:
     		   		print('Initial file content:',line)
     		   	
		self.allocate_arrays()
		
		data = np.loadtxt(self.infile, usecols=range(0,6),skiprows = 3,comments='#')
		
		for inum in range(self.line_num):
		
			self.RA[inum] = data[inum][0]
			self.DEC[inum] = data[inum][1] 
			self.e1[inum] = data[inum][2]
			self.e2[inum] = data[inum][3]
			self.weight[inum] = data[inum][4]
			self.z[inum] = data[inum][5]
		
		print('Loaded Kids simulation',file_name)
	
	def Gen_Mock_z(self,z_index,outfile_z_bin,num_bin):
		""" Write out tomographic mocks (in an new directory) using z_bin index arrays..
		"""

		full_outfile = outfile_z_bin.replace('.dat','_z_bin_Num_' + str(num_bin) + '.dat')
		
		size = len(z_index)
		
		outdata = np.empty([size,5],dtype=float)
		
		for jnum,z_pos in enumerate(z_index):
	
			outdata[jnum][0] = self.RA[z_pos]
			outdata[jnum][1] = self.DEC[z_pos]
			outdata[jnum][2] = self.e1[z_pos]
			outdata[jnum][3] = self.e2[z_pos]
			outdata[jnum][4] = self.weight[z_pos]
		
		check = np.savetxt(full_outfile,outdata)
						
		if check == 1:
			print('write z bin num',num_bin,'complete')
		else: 
			print('problem with save',num_bin)
	
		return 1 
	
	def AddShearError(self):	

		mu,sd = 0.00,self.error_sd
		
		error_array1 = np.random.normal(mu,sd,self.line_num)
		error_array2 = np.random.normal(mu,sd,self.line_num)

		for inum in range(self.line_num):
		
			e1 = self.e1[inum] 
			e2 = self.e2[inum] 

			n1 = error_array1[inum]
			n2 = error_array2[inum]

			scatter_a1 = e1 + n1
			scatter_a2 = e2 + n2
			scatter_a3 = 1.0 + e1*n1 + e2*n2
			scatter_a4 = e1*n2 - e2*n1

			self.e1[inum]=(scatter_a1*scatter_a3 + scatter_a2*scatter_a4)/(scatter_a3*scatter_a3 + scatter_a4*scatter_a4)
			self.e2[inum]=(scatter_a2*scatter_a3 - scatter_a1*scatter_a4)/(scatter_a3*scatter_a3 + scatter_a4*scatter_a4)

	def allocate_arrays(self):
		
		self.RA = np.empty(self.line_num,dtype=float)
		self.DEC = np.empty(self.line_num, dtype=float)
		self.z = np.empty(self.line_num, dtype=float)
		self.e1 = np.empty(self.line_num, dtype=float)
		self.e2 = np.empty(self.line_num, dtype=float)
		self.weight = np.empty(self.line_num, dtype=float)
	
	
class Kids_mock_L:
	""" Used to store all mock data for sources.
	"""
	def __init__(self,file_name):
		""" Read in kids simulation data.. 
		"""
		self.infile = file_name.replace('\n','') # remove end of line indicator
		self.line_num  = 0 

		for line in open(self.infile):
	    		if not line.strip().startswith("#"):
	    			self.line_num = int(line)
	    			print('line number is',self.line_num)
	    			break
     		  	else:
     		   		print('Initial file content:',line)
     		   	
		self.allocate_arrays()
		
		data = np.loadtxt(self.infile, usecols=range(0,3),skiprows = 3,comments='#')
		
		for inum in range(self.line_num):
		
			self.RA[inum] = data[inum][0]
			self.DEC[inum] = data[inum][1] 
			self.z[inum] = data[inum][2]
		
		print('Loaded Kids simulation',file_name)
	
	def Gen_Mock_z(self,z_index,outfile_z_bin,num_bin):
		""" Write out tomographic mocks (in an new directory) using z_bin index arrays..
		"""

		outfile = outfile_z_bin.replace('.dat','_z_bin_Num_')
		z_index =  str(num_bin)
		filename_suffix = '.dat' + '\n'

		full_outfile = os.path.join(outfile, z_index + filename_suffix)

		size = len(z_index)
		
		outdata = np.empty([size,5],dtype=float)
		
		for jnum,z_pos in enumerate(z_index):
	
			outdata[jnum][0] = self.RA[z_pos]
			outdata[jnum][1] = self.DEC[z_pos]
			outdata[jnum][2] = 0
			outdata[jnum][3] = 0
			outdata[jnum][4] = 1.0
		
		check = np.savetxt(full_outfile,outdata)
						
		if check == 1:
			print('write z bin num',num_bin,'complete')
		else: 
			print('problem with save',num_bin)
	
		return 1 	


	def allocate_arrays(self):
		
		self.RA = np.empty(self.line_num,dtype=float)
		self.DEC = np.empty(self.line_num, dtype=float)
		self.z = np.empty(self.line_num, dtype=float)



def GetSourceFileNamesFull(dir_main,file_main,path_obs):
	
   FileNames = []	
   txtfile = dir_main + '/' + file_main
   count = 0
   for line in open(txtfile,'r').readlines():
       FileNames.append(path_obs + '/' + line)
       count += 1
  	
   return FileNames
  
def GetLensFileNamesFull(dir_main,file_main,path_obs):
	
   FileNames = []	
   txtfile = dir_main + '/' + file_main
   count = 0
   for line in open(txtfile,'r').readlines():
       FileNames.append(path_obs + '/' + line)
       count += 1
  	
   return FileNames  
   
def index_z_bins(z_array):

	""" Sort the z positions into z-bins..
	"""
	
	num_bins = 5	## note needs to be changed!
	
	z_index1 = []
	z_index2 = []
	z_index3 = []
	z_index4 = []
	z_index5 = []
			
	for inum,z_value in enumerate(z_array):

		if  0.0 < z_value < 0.3:
			z_index1.append(inum)
			
		elif 0.3 <= z_value < 0.5:
			z_index2.append(inum)
			
		elif 0.5 <= z_value < 0.7:
			z_index3.append(inum)				
		
		elif 0.7 <= z_value < 0.9:
			z_index4.append(inum)
		
		elif 0.9 <= z_value < 3:
			z_index5.append(inum)
		else:
			print('Problem! Check z value', z_value)
				
	z_location = [z_index1,z_index2,z_index3,z_index4,z_index5]			
				
	return z_location,num_bins

def main():
	
	### Calculation options..

	Tomography_source = False
	Tomography_lens = False
	Reformat_sources_files = True
	Reformat_lens_files = True

	### Default mock dir
	mock_dir = '/lustre/projects/p018_swin/ajohnson/2dfLens/Mocks/mockcats'

	### Default Mock files
	mock_source_files = 'Mock_source_files.txt'
	mock_lens_files = 'Mock_lenses_files.txt'

	### New Mock files (write list of new file names here)
	new_mock_source_files = 'Mock_source_files_Group1.txt'
	new_mock_lens_files = 'Mock_lenses_files_Group1.txt'
	
	### Directory of list of file names, assume current dir..
	dir_list = os.getcwd()
	
	### array of all mock file names
	source_list= []
	lens_list = []
	
	# New array of all mock file names
	new_source_list= []
	new_lens_list = []
	
	## Write directoy..
	write_dir = '/lustre/projects/p018_swin/ajohnson/2dfLens/Mocks/mock_group1'
	
	parser = OptionParser(usage='Compute shear-shear and shear-galaxy correlation functions using athena')
   	
	parser.add_option('-L', '--mocks_lens', dest='lens_files', type='string', default=mock_lens_files,help='lens mocks names here')
	parser.add_option('-S', '--mocks_sources', dest='source_files', type='string', default=mock_source_files, help='source mocks names here')    
                      
	options, args = parser.parse_args()
   		
	### Find mock dirs
	
	source_list = GetSourceFileNamesFull(dir_list,options.source_files,mock_dir)   	
   	lens_list = GetLensFileNamesFull(dir_list,options.lens_files,mock_dir)

	if Tomography_source:

		### Create redshift bins and write out new simulations..
		
		rank = MPI.COMM_WORLD.Get_rank()
		size = MPI.COMM_WORLD.Get_size()
   			
		num_jobs = len(source_list)
		
		jobs_per_core = num_jobs/size
		
		if rank == 0:
			print('number of simulations to process is',num_jobs)
			print('Number of MPI cores in use is',size)
			print('jobs per core is',jobs_per_core)
		
		MPI_core_source_list = []
		
		for MPI_RANK in range(size):
		
			if rank == MPI_RANK:
				
				index_1 = MPI_RANK*jobs_per_core
				index_2 = MPI_RANK*jobs_per_core + 1
				
				MPI_core_source_list.append(source_list[index_1])
				MPI_core_source_list.append(source_list[index_2])
				
				print('MPI Rank',rank,'processing files',MPI_core_source_list[0],'and',MPI_core_source_list[1])
		
		
		for infile in MPI_core_source_list:
		
			mock_S = Kids_mock_S(infile)
		
			z_location, num_bins = index_z_bins(mock_S.z)

			for zbin in range(num_bins):
					
				mock_S.Gen_Mock_z(z_location[zbin],infile,zbin)			

		return 1

	if Tomography_lens: 

		### Create redshift bins and write out new simulations..
		
		rank = MPI.COMM_WORLD.Get_rank()
		size = MPI.COMM_WORLD.Get_size()
   			
		num_jobs = len(source_list)
		
		jobs_per_core = num_jobs/size
		
		if rank == 0:
			print('number of simulations to process is',num_jobs)
			print('Number of MPI cores in use is',size)
			print('jobs per core is',jobs_per_core)
		
		MPI_core_lens_list = []
		
		for MPI_RANK in range(size):
		
			if rank == MPI_RANK:
				
				index_1 = MPI_RANK*jobs_per_core
				index_2 = MPI_RANK*jobs_per_core + 1
				
				MPI_core_lens_list.append(lens_list[index_1])
				MPI_core_lens_list.append(lens_list[index_2])
				
				print('MPI Rank',rank,'Processing Files',MPI_core_lens_list[0],'and',MPI_core_lens_list[1])
		
		for infile in MPI_core_lens_list:
		
			mock_L = Kids_mock_L(infile)
		
			z_location, num_bins = index_z_bins(mock_L.z)

			for zbin in range(num_bins):
					
				mock_L.Gen_Mock_z(z_location[zbin],infile,zbin)			

		return 1
	
	if Reformat_sources_files:

		### Load mocks and save reformatted files..
		   
		for inum in range(len(source_list)):
			new_source_list.append(source_list[inum].replace('.dat','_group1.dat'))
			
		for inum,infile in enumerate(source_list):
			
			Mocks_S = Kids_mock_S(infile)

			Mocks_S.AddShearError()

			outdata_S = np.empty([Mocks_S.line_num,5],dtype=float)

			outfile = new_source_list[inum]
			
			for jnum in range(Mocks_S.line_num):
		
				outdata_S[jnum][0] = Mocks_S.RA[jnum]
				outdata_S[jnum][1] = Mocks_S.DEC[jnum]
				outdata_S[jnum][2] = Mocks_S.e1[jnum]
				outdata_S[jnum][3] = Mocks_S.e2[jnum]
				outdata_S[jnum][4] = Mocks_S.weight[jnum]
			
			print('writing out',outfile)
			
			np.savetxt(outfile,outdata_S)

	if Reformat_lens_files:
		### Load mocks and save reformatted files..
		   	
		rank = MPI.COMM_WORLD.Get_rank()
		size = MPI.COMM_WORLD.Get_size()   			
		num_jobs = len(lens_list)
	
		jobs_per_core = num_jobs/size
	
		if rank == 0:
			print('number of simulations to process is',num_jobs)
			print('Number of MPI cores in use is',size)
			print('jobs per core is',jobs_per_core)
	
		MPI_core_lens_list = []
	
		for MPI_RANK in range(size):
	
			if rank == MPI_RANK:
			
				index_1 = MPI_RANK*jobs_per_core
				index_2 = MPI_RANK*jobs_per_core + 1
			
				MPI_core_lens_list.append(lens_list[index_1])
				MPI_core_lens_list.append(lens_list[index_2])
			
				print('MPI Rank',rank,'processing files',MPI_core_lens_list[0],'and',MPI_core_lens_list[1])

		for inum in range(len(lens_list)):
			new_lens_list.append(lens_list[inum].replace('.dat','_group1.dat'))
		
		for inum,infile in enumerate(lens_list):
			
			Mocks_L = Kids_mock_L(infile)
			outdata_L = np.empty([Mocks_L.line_num,5],dtype=float)
			outfile = new_lens_list[inum]
			
			for jnum in range(Mocks_L.line_num):
		
				outdata_L[jnum][0] = Mocks_L.RA[jnum]
				outdata_L[jnum][1] = Mocks_L.DEC[jnum]
				outdata_L[jnum][2] = 0
				outdata_L[jnum][3] = 0
				outdata_L[jnum][4] = 1.0
			
			print('writing out',outfile)
			
			np.savetxt(outfile,outdata_L)

	return 1
	
if __name__ == "__main__":

   main()