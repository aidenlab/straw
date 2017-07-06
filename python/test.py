#!/usr/bin/python
from numpy import *
from numpy import shape
from sklearn import svm
from scipy.cluster.vq import kmeans,vq
from sklearn import decomposition
from scipy import stats

import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist
import scipy.ndimage
import time
import sys
import os
import struct
from scipy import sparse
import matplotlib
#import matplotlib.pyplot as plt
import glob
import argparse
import random
import straw
import timeit

jb_jar_path = '/home/jaeweon/Downloads/juicer_tools.1.6.2_jcuda.0.7.0.jar' 
norm = "KR"
# hic_file = '/home/jaeweon/research/imr90_intra_nofrag_30.hic'
# hic_file = "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hicfiles"
# hic_file = "/home/jaeweon/research/gm12878_intra_nofrag_30.hic"
hic_file = "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/primary.hic"

chr1 = "22"	
chr2 = "22"
units = "BP"
res = 100000

result = straw.straw(norm, hic_file, chr1, chr2, units, res)
# print result, shape(result)

# find the matrix size of result 

# call_juicer_jar = 'java -Xms512m -Xmx2048m -jar ' + jb_jar_path + ' dump observed '+ norm + ' ' + hic_file + ' ' + chr1 + ' ' + chr2 + ' ' + units + ' ' + str(res)
# output_file = "local_results.txt"

# print( call_juicer_jar  + ' ' +  output_file)
# os.system(call_juicer_jar + ' ' + output_file)

# data = loadtxt(output_file)
# os.system('rm '+ output_file)

# c_temp = zeros((1000, 1000), dtype=float32)

# I = array(data[:, 0]/res, dtype=int)
# J = array(data[:, 1]/res, dtype=int)
# V = data[:, 2]
# V[isnan(V)] = 0
# print( "csize1 ", shape(c_temp))
# c_temp[I, J] = V
# if int(chr1) == int(chr2):
# 	c_temp[J, I] = V

# print c_temp, shape(c_temp)
# compare the matrix size with straw
