import numpy as np

def rgb2ntsc(rgb):
	trans = np.array([[0.299,0.596,0.211],[0.587,-0.274,-0.523],[0.114,-0.322,0.312]])
	#trans = np.array([[0.299,0.587,0.114],[-0.14713,-0.28886,0.436],[0.615,-0.51498,-0.10001]])
	#trans = np.array([[0.299,-0.14713,0.614],[0.587,-0.28886,-0.51498],[0.114,0.436,-0.10001]])	
	shape = np.shape(rgb)

	r = rgb[:,:,0]
	g = rgb[:,:,1]
	b = rgb[:,:,2]

	
	r = np.reshape(r,(r.size,1))
	g = np.reshape(g,(g.size,1))
	b = np.reshape(b,(b.size,1))
	'''	
	r = np.reshape(r,(1,r.size))
	g = np.reshape(g,(1,g.size))
	b = np.reshape(b,(1,b.size))
	'''
	#rgb = np.concatenate((r,g,b))
	rgb = np.concatenate((r,g,b),1)
	yiq = np.dot(rgb,trans)
	yiq = np.reshape(yiq,shape)
	return yiq

def ntsc2rgb(yiq):
        trans = np.array([[1.0,1.0,1.0],[0.95617,-0.27269,-1.10374],[0.62143,-0.64681,1.70062]])

        shape = np.shape(yiq)

        y = yiq[:,:,0]
        i = yiq[:,:,1]
        q = yiq[:,:,2]

        y = np.reshape(y,(y.size,1))
        i = np.reshape(i,(i.size,1))
        q = np.reshape(q,(q.size,1))

        yiq = np.concatenate((y,i,q),1)
        rgb = np.dot(yiq,trans)
        rgb = np.reshape(rgb,shape)
        return rgb

#a = np.array([[[1,2,3],[4,5,6]],[[7,8,9],[10,11,12]]])
#print rgb2ntsc(a)
	

