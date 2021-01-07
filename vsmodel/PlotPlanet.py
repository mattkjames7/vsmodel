import numpy as np


def PlotPlanetXY(fig,R=1.0,Center=[0.0,0.0,0.0],zorder=10,NoBlack=False,NoonTop=True):
	a = 2*np.pi*np.arange(361,dtype='float32')/360
	x = R*np.sin(a) + Center[0]
	y = R*np.cos(a) + Center[1]
	
	if NoonTop:	
		fig.fill(y,x,color=[1.0,1.0,1.0],zorder=zorder)
		fig.plot(y,x,color=[0,0,0],zorder=zorder+1)
		if NoBlack == False:
			fig.fill(y[180:360],x[180:360],color=[0.0,0.0,0.0],zorder=zorder+1)
	else:
		fig.fill(x,y,color=[1.0,1.0,1.0],zorder=zorder)
		fig.plot(x,y,color=[0,0,0],zorder=zorder+1)
		if NoBlack == False:
			fig.fill(x[180:360],y[180:360],color=[0.0,0.0,0.0],zorder=zorder+1)


