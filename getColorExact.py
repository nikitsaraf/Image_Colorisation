import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import math
from colorSpace import ntsc2rgb
import pdb
import pickle

def getColorExact(colorIm,ntscIm):
	n = np.size(ntscIm,0)
	m = np.size(ntscIm,1)
	imgSize = n*m
	
	nI = np.zeros(ntscIm.shape)
	
	nI[:,:,0] = ntscIm[:,:,0]

	indsM = np.arange(0,imgSize).reshape((n,m),order='F')
	print "indsM=\n"
	print indsM
	test = np.ravel(colorIm,order='F')
	test = test.nonzero()
	lblInds = np.array(test)
	lblInds = lblInds.transpose()
	print "lblInds.shape=\n"
	print lblInds.shape
	
	wd = 1;
	
	len = 0;
	consts_len = 0;
	col_inds = np.zeros((imgSize*(2*wd+1)**2,1),float)
	row_inds = np.zeros((imgSize*(2*wd+1)**2,1),float)
	vals = np.zeros((imgSize*(2*wd+1)**2,1),float)
	gvals = np.zeros((1,(2*wd+1)**2))
	colorIm = colorIm.astype(int)
	#pdb.set_trace()
	for j in range(m):
		for i in range(n):
			if (~colorIm[i,j]):
				tlen = 0;
				ii = max(0,i-wd)
				while (ii<=min(i+wd,n-1)):
					jj = max(0,j-wd)
					while (jj<=min(j+wd,m-1)):
						if (ii!=i) or (jj!=j):
							row_inds[len] = consts_len
							col_inds[len] = indsM[ii,jj]
							gvals[0,tlen] = ntscIm[ii,jj,0]
							len=len+1
							tlen=tlen+1
						jj = jj + 1
					ii = ii + 1
				t_val = ntscIm[i,j,0]
				gvals[0,tlen] = t_val
				c_var = np.mean((gvals[0,0:tlen+1] - np.mean(gvals[0,0:tlen+1]))**2)
				csig = c_var*0.6
				mgv = np.min((gvals[0,0:tlen]-t_val)**2)
				if (csig<(-mgv/math.log(0.01))):
					csig = -mgv/math.log(0.01)
				if (csig<0.000002):
					csig = 0.000002
			
				gvals[0,0:tlen] = np.exp(-(((gvals[0,0:tlen]-t_val)**2)/csig))
				gvals[0,0:tlen] = gvals[0,0:tlen]/np.sum(gvals[0,0:tlen])
				vals[(len-1)-(tlen-1):len,0] = -(gvals[0,0:tlen].T)

			row_inds[len] = consts_len
			col_inds[len] = indsM[i,j]
			vals[len] = 1
			consts_len = consts_len + 1

		
	vals = vals[0:len]
	col_inds = col_inds[0:len]
	row_inds = row_inds[0:len]
	#print "vals=\n" + vals
	#print "col_inds=\n" + col_inds
	#print "row_inds=\n" + row_inds
	vals = vals.transpose()
	vals.shape = (np.size(vals,1),)
	col_inds = col_inds.transpose()
	col_inds.shape = (np.size(col_inds,1),)
	row_inds = row_inds.transpose()
	row_inds.shape = (np.size(row_inds,1),)
	print "row_inds.max="
	print row_inds.max()
	print "consts_len="
	print consts_len
	
	A = coo_matrix((vals,(row_inds,col_inds)),shape=(consts_len,imgSize)).tocsc()
	#A = np.squeeze(np.asarray(A))
	#A.shape = (consts_len,imgSize)
	#print type(A)
	#print A
	#print A.shape
	b = np.zeros((consts_len,1))
	#print "A=\n" + A
	#print "b=\n" + b
	
	with open('imgcolor.pkl','wb') as f:
		pickle.dump(vals,f)
		pickle.dump(row_inds,f)
		pickle.dump(col_inds,f)
		pickle.dump(consts_len,f)
		pickle.dump(imgSize,f)
		pickle.dump(b,f)
		
	t = 1
	while (t<=2):
		curIm = ntscIm[:,:,t]
		for i in range(np.size(lblInds,0)):
			b = np.ravel(b,order='F')
			curIm = np.ravel(curIm,order='F')
			b[lblInds[i]]=curIm[lblInds[i]]
		b = b.reshape((consts_len,1),order='F')
		new_vals = spsolve(A,b)
		print "Pass"
		print t
		print "Compeleted"
		nI[:,:,t] = new_vals.reshape((n,m),order='F')
		t = t+1
	
	snI = nI
	#nI = ntsc2rgb(nI)

	return snI

			
	
