import pylab as pl
from colorSpace import *
import numpy as np
import math
from skimage import img_as_float
from getColorExact import getColorExact

inputFile = 'example.bmp'
markedFile = 'example_m.bmp'
outputFile = 'example_rp.bmp'

gI = img_as_float(pl.imread(inputFile))
cI = img_as_float(pl.imread(markedFile))
print gI[113:120,112:120,1]
colorIm = np.abs(gI-cI)
colorIm = colorIm.sum(axis=2)
colorIm = colorIm>0.01
colorIm = colorIm.astype(float)

print "colorIm.shape=\n"
print colorIm.shape
print colorIm[113:120,112:120]
sgI = rgb2ntsc(gI)
scI = rgb2ntsc(cI)
pl.imshow(scI)
pl.show()
'''
ntscIm = np.zeros(np.shape(sgI),float)

ntscIm[:,:,0] = sgI[:,:,0]
ntscIm[:,:,1] = scI[:,:,1]
ntscIm[:,:,2] = scI[:,:,2] 

max_d = math.floor(math.log(min(np.size(ntscIm,0),np.size(ntscIm,1)))/math.log(2)-2)
iu = math.floor(np.size(ntscIm,0)/(2**(max_d-1)))*(2**(max_d-1))
ju = math.floor(np.size(ntscIm,1)/(2**(max_d-1)))*(2**(max_d-1))
id = 0
jd = 0
print "max_d=\n"
print max_d
print "iu=\n"
print iu
print "ju=\n"
print ju
colorIm = colorIm[id:iu+1,jd:ju+1]
#pl.imshow(colorIm,origin='lower')
#pl.show()
print colorIm.shape
ntscIm = ntscIm[id:iu,jd:ju,:]
#pl.imshow(colorIm,origin='lower')
#pl.show()

nI = getColorExact(colorIm,ntscIm)

pl.imshow(nI)

pl.show()
'''
