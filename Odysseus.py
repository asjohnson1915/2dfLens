#! /usr/bin/env python
## Check branch structure
## Checkcheck
# Author: Andrew Johnson -- with some code taken from Athena 

from __future__ import print_function

import os
import pyfits
import re
from optparse import OptionParser 
import numpy as np 
import subprocess
import math
import errno
import shlex
import sys
from mpi4py import MPI

def run_cmd(cmd, run=True, verbose=True, stop=False):
    """ Borrowed from mkstuff.py, runs shell command  
    """
    cmds = shlex.split(cmd)

    if run is True:
        if verbose is True:
            print('Running command \'{0}\''.format(cmd))
        try:
            ex = subprocess.call(cmds)
        except OSError as e:
            print('Error: {0}'.format(e.strerror))
            ex = e.errno
    else:
        if verbose is True:
            print('Not running command \'{0}\''.format(cmd))
        return 0

    if ex != 0:
        if verbose is True:
            print('Last command returned {0}'.format(ex), end='')
            if stop is True:
                print(', stopping')
            else:
                print(', continuing')

        if stop is True:
            sys.exit(ex)

    return ex


class athena:
  """Run athena on mocks via shell commands
  """
  def __init__(self, dir, func, args, bindir):
    self.dir      = dir                 # Directory of test file (mocks)
    self.func     = func                # Program to be called (athena)
    self.args     = args                # Arguments with function call -- give configure file + output file 
    self.bindir   = bindir              # Location of test code

  def chdir(self):
    pwd = os.getcwd()
    try:
      os.chdir(self.dir)
    except:
      print('chdir: Could not change into directory \'{0}\''.format(self.dir))
      return None

    return pwd

  def run(self):
    """Run test via shell command
    """
    pwd = self.chdir()
    if pwd is None: return -1
    res = run_cmd(self.bindir + '/' + self.func + ' ' + self.args)
    os.chdir(pwd)
    return res
        

class Config_shear_shear:
    """Use to generate configure files for athena... Set based on simulations being used
    Input file names, and output name for configure file. 
    Kids simulations are in Deg.
    """ 

    def __init__(self, galcat1, galcat2,output):

        self.galcat1 = galcat1       # Catalogue 1
        self.galcat2 = galcat2       # Catalogue 2
        self.output   = output       # Output file written to current dir then read by anethna.. 
        self.WCORR = '1'            # 1: shear-shear, 2: shear-position, 4: position-position
        self.SFORMAT = 'standard'        # One of standard, hamana, position
        self.SCOORD_INPUT = 'deg'      # Input catalogue coordinates, {arcsec|arcmin|rad|deg}
        self.SCOORD_OUTPUT = 'arcmin'      # Output coordinates
        self.THMIN = '0.02'            # Smallest scale  in units of 'SCOORD_OUTPUT' 
        self.THMAX = '45'              # Largest scale in units of 'SCOORD_OUTPUT'
        self.NTH = '10'              # Number of bins
        self.BINTYPE = 'LOG'          # LIN or LOG
        self.RADEC = '0'            # 0: Cartesian, 1: spherical coordinates
        self.OATH = '0.02'            # Open angle threshold [rad]
        self.SERROR = 'jackknife'        # Error type ('none', 'bootstrap', 'jackknife')
        self.NRESAMPLE = '20'    

    def write_config(self):
        """Write Configure File Out"""
      
        fout = open(self.output, 'w')
        format = []
        format = self.gen_config()        
        for inum in range(len(format)):
            print(format[inum],file = fout)
        print(file=fout, end='\n')

        return 1  

    def gen_config(self):
        """build Configure File Out for athena"""
      
        format = []
        format.append("GALCAT1" + ' '+ self.galcat1)
        format.append("GALCAT2" + ' '+ self.galcat2)
        format.append("WCORR" + ' '+ self.WCORR)
        format.append("SFORMAT" + ' '+ self.SFORMAT)
        format.append("SCOORD_INPUT" + ' '+ self.SCOORD_INPUT)
        format.append("SCOORD_OUTPUT" + ' '+ self.SCOORD_OUTPUT)
        format.append("THMIN" + ' '+ self.THMIN)
        format.append("THMAX" + ' '+ self.THMAX)
        format.append("NTH" + ' '+ self.NTH)
        format.append("BINTYPE" + ' '+ self.BINTYPE)
        format.append("RADEC" + ' '+ self.RADEC)
        format.append("OATH" + ' '+ self.OATH)
        format.append("SERROR" + ' '+ self.SERROR)
        format.append("NRESAMPLE" + ' '+ self.NRESAMPLE)

        return format  

class Config_shear_position:
    """Use to generate configure files for athena... Set based on simulations being used
    Input file names, and output name for configure file. 
    Kids simulations are in Deg.
    """ 

    def __init__(self, galcat1, galcat2,output):

        self.galcat1 = galcat1       # Catalogue 1
        self.galcat2 = galcat2       # Catalogue 2
        self.output   = output       # Output file written to current dir then read by anethna.. 
        self.WCORR = '2'            # 1: shear-shear, 2: shear-position, 4: position-position
        self.SFORMAT = 'standard'        # One of standard, hamana, position
        self.SCOORD_INPUT = 'deg'      # Input catalogue coordinates, {arcsec|arcmin|rad|deg}
        self.SCOORD_OUTPUT = 'arcmin'      # Output coordinates
        self.THMIN = '0.5' #'0.02'            #Smallest scale  in units of 'SCOORD_OUTPUT' 
        self.THMAX = '60' #'45'              # Largest scale in units of 'SCOORD_OUTPUT'
        self.NTH = '20' #'15'              # Number of bins
        self.BINTYPE = 'LOG'          # LIN or LOG
        self.RADEC = '0'            # 0: Cartesian, 1: spherical coordinates
        self.OATH = '0.05'            # Open angle threshold [rad]
        self.SERROR = 'jackknife'        # Error type ('none', 'bootstrap', 'jackknife')
        self.NRESAMPLE = '50 50'    

    def write_config(self):
        """Write Configure File Out"""
      
        fout = open(self.output, 'w')
        format = []
        format = self.gen_config()        
        for inum in range(len(format)):
            print(format[inum],file = fout)
        print(file=fout, end='\n')

        return 1  

    def gen_config(self):
        """build Configure File Out for athena"""
      
        format = []
        format.append("GALCAT1" + ' '+ self.galcat1)
        format.append("GALCAT2" + ' '+ self.galcat2)
        format.append("WCORR" + ' '+ self.WCORR)
        format.append("SFORMAT" + ' '+ self.SFORMAT)
        format.append("SCOORD_INPUT" + ' '+ self.SCOORD_INPUT)
        format.append("SCOORD_OUTPUT" + ' '+ self.SCOORD_OUTPUT)
        format.append("THMIN" + ' '+ self.THMIN)
        format.append("THMAX" + ' '+ self.THMAX)
        format.append("NTH" + ' '+ self.NTH)
        format.append("BINTYPE" + ' '+ self.BINTYPE)
        format.append("RADEC" + ' '+ self.RADEC)
        format.append("OATH" + ' '+ self.OATH)
        format.append("SERROR" + ' '+ self.SERROR)
        format.append("NRESAMPLE" + ' '+ self.NRESAMPLE)

        return format  

def GetSourceFileNames(dir_main,file_main):
    
    FileNames = []  
    txtfile = dir_main + '/' + file_main
    count = 0
    for line in open(txtfile,'r').readlines():
        FileNames.append(line)
        count += 1
    
    return FileNames
  
def GetLensFileNames(dir_main,file_main):
    
    FileNames = []  
    txtfile = dir_main + '/' + file_main
    count = 0
    for line in open(txtfile,'r').readlines():
        FileNames.append(line)
        count += 1
    
    return FileNames    


# def Get_All_shear_shear_files(source_files,NUM_zBINS): ## For sim 1 now..
    
#     NumSims = 16

#     list_1 = []
#     list2 = []

#     Simulations = []

#     suffix = '_z_bin_Num'

#     for inum in range(NumSims):
#         Simulations.append([])


#     for num,name in enumerate(source_files):
#         for inum in length(NUM_zBINS):

#             extrapath = '_z_bin_Num' + str(inum) + '.dat'
#             file_name = name.replace('.dat',extrapath)
#             Simulations[num].append(file_name)


#    return list_1 list_2 


def main():

    Calc_Shear_Shear = False

    Calc_Pos_Shear = True

    Tomography = False

    # Default mock dir
    mock_dir = '/lustre/projects/p018_swin/ajohnson/2dfLens/Mocks/mockcats'

    # Default Mock files
    mock_source_files = 'Mock_source_files_updated.txt'
    mock_lens_files = 'Mock_lenses_files_updated.txt'

    # Directory of list of file names, assume current dir..
    dir_list = os.getcwd()

    # Athena dir 
    athena_dir = '/lustre/projects/p018_swin/ajohnson/2dfLens/athena_1.7/bin'

    #results Dir 
    results_dir = '/lustre/projects/p018_swin/ajohnson/2dfLens/mock_results'

    # array of all mock files
    source_list= []

    lens_list = []

    parser = OptionParser(usage='Compute shear-shear and shear-galaxy correlation functions using athena')

    parser.add_option('-T', '--path_test', dest='path_to_test', type='string', default=mock_dir, help='Absolute path to mock directory')
    parser.add_option('-L', '--mocks_lens', dest='lens_files', type='string', default=mock_lens_files,help='lens mocks names here')
    parser.add_option('-S', '--mocks_sources', dest='source_files', type='string', default=mock_source_files, help='source mocks names here')
    parser.add_option('-A', '--path_athena', dest='path_to_athena', type='string', default=athena_dir, help='Absolute path to programs athena')        

    options, args = parser.parse_args()

    source_list = GetSourceFileNames(dir_list,options.source_files)

    lens_list = GetLensFileNames(dir_list,options.lens_files)
    
    out_config = options.path_to_test + '/' + 'configure_tree_temp' 

    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
  
    num_jobs = len(source_list)
        
    jobs_per_core = num_jobs/size
        
    if rank == 0:

        print('number of simulations to process is',num_jobs)
        print('Number of MPI cores in use is',size)
        print('jobs per core is',jobs_per_core)
        
    MPI_core_source_list = []

    MPI_core_lens_list = []
    
    if not (Tomography):

        if(Calc_Shear_Shear):

            for MPI_RANK in range(size):
                
                if rank == MPI_RANK:
                        
                    index_1 = MPI_RANK*jobs_per_core
                    index_2 = MPI_RANK*jobs_per_core + 1
                        
                    MPI_core_source_list.append(source_list[index_1])
                    MPI_core_source_list.append(source_list[index_2])
                        
                    print('MPI Rank',rank,'processing files',MPI_core_source_list[0],'and',MPI_core_source_list[1])
              
            for inum, mockname in enumerate(MPI_core_source_list):
                
                outpt = str(inum + jobs_per_core*rank) ## To distingush between output files 

                out_config += str(rank)

                config = Config_shear_shear(mockname,mockname,out_config)

                config.write_config()

                arg1 = '-c configure_tree_temp' + str(rank)

                arg2 = '--out_xi' + ' ' + results_dir + '/' + 'xi_xi_kidss_final_' + outpt

                athena_mock = athena(options.path_to_test, 'athena',arg1 + ' ' + arg2 , options.path_to_athena)

                athena_mock.run()

        if(Calc_Pos_Shear):

            for MPI_RANK in range(size):
                
                if rank == MPI_RANK:
                        
                    index_1 = MPI_RANK*jobs_per_core
                    index_2 = MPI_RANK*jobs_per_core + 1
                        
                    MPI_core_lens_list.append(lens_list[index_1])
                    MPI_core_lens_list.append(lens_list[index_2])

                    MPI_core_source_list.append(source_list[index_1])
                    MPI_core_source_list.append(source_list[index_2])
                        
                    print('MPI Rank',rank,'processing files',MPI_core_lens_list,'and',MPI_core_lens_list)

            if(len(MPI_core_lens_list) != len(MPI_core_source_list)):

                print('Problem source lens list differnt sizes')

                exit()

            for inum,lens_name in enumerate(MPI_core_lens_list):
                
                out_config += str(rank)

                outpt = str(inum + jobs_per_core*rank) ## To distingush between output files 

                source_name = MPI_core_source_list[inum]

                config_g = Config_shear_position(source_name,lens_name,out_config)

                config_g.write_config()

                arg1 = '-c configure_tree_temp' + str(rank)

                arg2 = '--out_gl' + ' ' + results_dir + '/' + 'galaxy_galaxy_lens_kidss_final_UPDATED2' + outpt

                athena_mock = athena(options.path_to_test, 'athena',arg1 + ' ' + arg2 , options.path_to_athena)

                athena_mock.run()


    exit()

    # if(Tomography):

    #     NUM_zBINS = 5

    #     source_list_z_1,source_list_z_2 = 1, 2 #Get_All_shear_shear_files(source_list,NUM_zBINS)

    #     if(Calc_Shear_Shear):

    #         for MPI_RANK in range(size):
                
    #             if rank == MPI_RANK:
                        
    #                 index_1 = MPI_RANK*jobs_per_core
    #                 index_2 = MPI_RANK*jobs_per_core + 1
                        
    #                 MPI_core_source_list_z_bin.append(   [index_1])
    #                 MPI_core_source_list_z_bin.append(   [index_2])
                        
    #                 print('MPI Rank',rank,'processing files',MPI_core_source_list_z_bin)

    #         for inum, mockname in enumerate(MPI_core_source_list):

    #             for zbin1 in range(NUM_zBINS)

    #                 for zbin2 in range(NUM_zBINS)

    #                     outpt = str(inum + jobs_per_core*rank)

    #                     out_config += str(rank)

    #                     config = Config_shear_shear(mockname,mockname,out_config)

    #                     config.write_config()

    #                     arg1 = '-c configure_tree_temp' + str(rank)

    #                     arg2 = '--out_xi' + ' ' + results_dir + '/' + 'xi_xi_kidss_final_' + outpt

    #                     athena_mock = athena(options.path_to_test, 'athena',arg1 + ' ' + arg2 , options.path_to_athena)

    #                     athena_mock.run()

    #     if(Calc_Pos_Shear):

    #         for MPI_RANK in range(size):
                
    #             if rank == MPI_RANK:
                        
    #                 index_1 = MPI_RANK*jobs_per_core
    #                 index_2 = MPI_RANK*jobs_per_core + 1
                        
    #                 MPI_core_lens_list.append(lens_list[index_1])
    #                 MPI_core_lens_list.append(lens_list[index_2])

    #                 MPI_core_source_list.append(source_list[index_1])
    #                 MPI_core_source_list.append(source_list[index_2])
                        
    #                 print('MPI Rank',rank,'processing files',MPI_core_lens_list,'and',MPI_core_lens_list)

    #         if(len(MPI_core_lens_list) != len(MPI_core_source_list)):
    #             print('Problem source lens list differnt sizes')
    #             exit()

    #         for inum in len(MPI_core_lens_list):
                
    #             out_config += str(rank)
    #             outpt = str(inum + jobs_per_core*rank) ## To distingush between output files 

    #             lens_name = MPI_core_lens_list[inum]
    #             source_name = MPI_core_source_list[inum]

    #             config_g = Config_shear_position(source_name,lens_name,out_config)

    #             config_g.write_config()

    #             arg1 = '-c configure_tree_temp' + str(rank)

    #             arg2 = '--out_gl' + ' ' + results_dir + '/' + 'galaxy_galaxy_lens_kidss_final_' + outpt

    #             athena_mock = athena(options.path_to_test, 'athena',arg1 + ' ' + arg2 , options.path_to_athena)

    #             athena_mock.run()


if __name__ == "__main__":

  main()

