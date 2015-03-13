#! /usr/bin/env python

import os
import pyfits
import re
from optparse import OptionParser 
import numpy as np 
import subprocess

def GetFitsFilesNames(txtfile):
    
    fileNames=[]
    for line in open(txtfile,'r').readlines():
        fileNames.append(line)

    StackFiles(fileNames,len(fileNames))    

def SortFiles(imlist,numfits):

    if(numfits == 1): 
        GENFITS(imlist[0])
    else:          
        StackFiles(imlist,numfits)
                          
def StackFiles(fitsfiles,numfits):

    # lets keep the headers from the first file and just update the data.. might work?
    # size given in header files appears to be updated automatically given .data.shape..

    HDUarray,rownum,colnum = [],[],[]

    hdulist0 = pyfits.open(fitsfiles[0])

    for inum in range(len(fitsfiles)):

        hdulist = pyfits.open(fitsfiles[inum])
        HDUarray.append(hdulist)
        row_temp, col_temp = HDUarray[inum][0].data.shape
        rownum.append(row_temp)
        colnum.append(colnum)

    rowtotal = sum(rownum)
    dataUnitTemp = np.zeros((rowtotal,col_temp)) 
    row_index = np.arange(0,rowtotal + rownum[0],rownum[0])

    for inum in range(len(fitsfiles)):
        row_i = row_index[inum]
        row_f = row_index[inum + 1]
        dataUnitTemp[row_i:row_f][:] = HDUarray[inum][0].data[:][:]
    
    newpri=pyfits.PrimaryHDU(data=dataUnitTemp, header=hdulist0[0].header)  
    newhdulist=pyfits.HDUList([newpri])

    next=len(hdulist0) 
   
    for iext in range(1,next):

        if hdulist0[iext].header['XTENSION']=='IMAGE':
            cmdstr='pyfits.ImageHDU'

        elif hdulist0[iext].header['XTENSION']=='TABLE':
            cmdstr='pyfits.TableHDU'

        elif hdulist0[iext].header['XTENSION']=='BINTABLE':
            cmdstr='pyfits.BinTableHDU'

        else:
            print "Does not recognize extension type %s."%hdulist0[iext].header['XTENSION']

        if hdulist0[iext].header['EXTNAME'] in ['VARIANCE','RWSS']:

            row1,col1,row_index1 =[],[],[]

            for inum in range(len(fitsfiles)):    
                row_temp, col_temp = HDUarray[inum][iext].data.shape
                row1.append(row_temp)
                col1.append(col_temp)

            totalnumrow = sum(row1)
            dataUnitTemp = np.zeros((totalnumrow,col1[0]))
            row_index1 = np.arange(0,totalnumrow + row1[0],row1[0])

            for inum in range(len(fitsfiles)):
                row_i = row_index1[inum]
                row_f = row_index1[inum + 1]
                dataUnitTemp[row_i:row_f][:] = HDUarray[inum][iext].data[:][:]

            newhdu=eval(cmdstr+'(data=dataUnitTemp, header=hdulist0[iext].header)')

        elif  hdulist0[iext].header['EXTNAME'] in ['FIBRES']:
            ## Fibers are annyoing, need custom datatype for arrays. Hence seperated from above, much to my annoyance.
            
            row1,col1,row_index1 =[],[],[]

            for inum in range(len(fitsfiles)):    
                col_temp = HDUarray[inum][iext].header['TFIELDS']
                row_temp = HDUarray[inum][iext].header['NAXIS2']
                row1.append(row_temp)
                col1.append(col_temp)

            totalrow = sum(row1)

            dt = np.dtype([('NAME',np.str_, 80),('RA',np.float64),('DEC',np.float64),('X',np.int32),('Y',np.int32),('XERR',np.int32),\
            ('YERR',np.int32),('THETA',np.float64),('TYPE',np.str_, 1),('PIVOT',np.int32),('MAGNITUDE',np.float64),('PID',np.int32),\
            ('COMMENT',np.str_, 80),('RETRACTOR',np.str_, 10),('WLEN',np.float64),('PMRA',np.float64),('PMDEC',np.float64)])

            dataUnitTemp = np.zeros(totalrow,dtype=dt)
            row_index1 = np.arange(0,totalrow + row1[0],row1[0])

            for inum in range(len(fitsfiles)):
                for pnum in range(row1[0]):
                    for qnum in range(col1[0]):
                        row_i = row_index1[inum]
                        dataUnitTemp[pnum + row_i][qnum] = HDUarray[inum][iext].data[pnum][qnum]
                        
            newhdu=eval(cmdstr+'(data=dataUnitTemp, header=hdulist0[iext].header)')

        elif  hdulist0[iext].header['EXTNAME'] in ['NDF_CLASS','SKY','TELCOR','REDUCTION_ARGS','REDUCED']:
            # Dont stack these HDUs (use values for first file), should not effect runz (hopefully)

            newhdu=eval(cmdstr+'(data=hdulist0[iext].data, header=hdulist0[iext].header)')

        elif  hdulist0[iext].header['EXTNAME'] in ['SHIFTS']:     

            row1,col1,row_index1 =[],[],[]

            for inum in range(len(fitsfiles)):    
                col_temp,row_temp = HDUarray[inum][iext].data.shape
                row1.append(row_temp)
                col1.append(col_temp)

            totalrow = sum(row1)
            dataUnitTemp = np.zeros((col1[0],totalrow))
            row_index1 = np.arange(0,totalrow + row1[0],row1[0])

            for inum in range(len(fitsfiles)):
                for qnum in range(row1[0]):
                    for pnum in range(col1[0]):
                        row_i = row_index1[inum]
                        dataUnitTemp[pnum][qnum + row_i] = HDUarray[inum][iext].data[pnum][qnum]
                   

            newhdu=eval(cmdstr+'(data=dataUnitTemp, header=hdulist0[iext].header)')

        elif  hdulist0[iext].header['EXTNAME'] in ['THPUT',]:    

            length =[]

            for inum in range(len(fitsfiles)):    
                length_temp = HDUarray[inum][iext].data.shape[0] 
                length.append(length_temp)

            totalsize = sum(length)
            dataUnitTemp = np.zeros((totalsize))
            len_index = np.arange(0,totalsize + length[0],length[0])

            for inum in range(len(fitsfiles)):
                for qnum in range(length[0]):
                    size_i = len_index[inum]
                    dataUnitTemp[pnum + size_i] = HDUarray[inum][iext].data[pnum]

            newhdu=eval(cmdstr+'(data=dataUnitTemp, header=hdulist0[iext].header)')

        else:
            print 'Unknown HDU -- need to update code ?'
            return 1
           
        newhdulist.append(newhdu)

    outfile = 'TempStacked.fits'

    newhdulist.writeto(outfile,clobber='True')

    GENFITS(outfile)

def GENFITS(newfile):

    outfile1 = '2df_spectra.fits'

    newoutputG1=os.path.basename(outfile1).replace('.fits','G_1.fits')
    newoutputG2=os.path.basename(outfile1).replace('.fits','G_2.fits')
    newoutputG3=os.path.basename(outfile1).replace('.fits','G_3.fits')
    newoutputG4=os.path.basename(outfile1).replace('.fits','G_4.fits') 
    newoutputG5=os.path.basename(outfile1).replace('.fits','G_5.fits') 

    file_names=[newoutputG1,newoutputG2,newoutputG3,newoutputG4,newoutputG5]

    data,header = pyfits.getdata(newfile,extname='FIBRES',header=True)

    numG1=[] # Direct photo-z plus bright (m < 19) galaxies 
    numG2=[] # CMASS galaxies
    numG3=[] # eBOSS galaxies
    numG4=[] # red nugget
    numG5=[] # (m > 19) low-z galaxies plus LENS-SDSS galaxies

    magcut=19 

    groupindex = [numG1,numG2,numG3,numG4,numG5]

    for inew,newrow in enumerate(data):
        if newrow['COMMENT']=='DIRECTPHOT' or newrow['MAGNITUDE'] < magcut and newrow['COMMENT']!='Parked':
            numG1.append(inew)
        elif newrow['COMMENT']=='LENSHIZ':
            numG2.append(inew)
        elif newrow['COMMENT']=='CROSSLRG':
            numG3.append(inew)
        elif newrow['COMMENT']=='REDNUGGET':
            numG4.append(inew)
        elif newrow['COMMENT']=='LENSLOZ' or newrow['COMMENT']=='LENSSDSS':   
            numG5.append(inew)
        else: print 'Ingnoring type',newrow['COMMENT']
            
    hdulist = pyfits.open(newfile)

    for inum in range(len(file_names)):

        newpri=pyfits.PrimaryHDU(data=hdulist[0].data[groupindex[inum]], header=hdulist[0].header)
        newhdulist=pyfits.HDUList([newpri])
        next=len(hdulist) 

        for iext in range(1,next):
            if hdulist[iext].header['XTENSION']=='IMAGE':
                cmdstr='pyfits.ImageHDU'
            elif hdulist[iext].header['XTENSION']=='TABLE':
                cmdstr='pyfits.TableHDU'
            elif hdulist[iext].header['XTENSION']=='BINTABLE':
                cmdstr='pyfits.BinTableHDU'
            else:
                print "Does not recognize extension type %s."%hdulist[iext].header['XTENSION']

            if hdulist[iext].header['EXTNAME'] in ['VARIANCE','FIBRES','RWSS']:
                newhdu=eval(cmdstr+'(data=hdulist[iext].data[groupindex[inum]], header=hdulist[iext].header)')
            else:
                newhdu=eval(cmdstr+'(data=hdulist[iext].data, header=hdulist[iext].header)')
        
            newhdulist.append(newhdu)

        outfile = file_names[inum]

        newhdulist.writeto(outfile,clobber='True')
    
def main():

    "Options: Pass file names from command line or supply a text file with a list with flag  --file"

    parser = OptionParser(usage='Sort 2dF Spectra into groups')
    
    parser.add_option("-f", "--file", dest="filename", help="txt file with list of .fits files to use", metavar="FILE")

    (opts, args) = parser.parse_args()

    numinputfits = len(args)

    if len(args)<1:
        GetFitsFilesNames(opts.filename)
    elif len(args)>1:

        SortFiles(args,numinputfits)
    else:
        print 'Need to supply .fits files'

if __name__ == "__main__":

    main()
    

# Comments on files produced..  
