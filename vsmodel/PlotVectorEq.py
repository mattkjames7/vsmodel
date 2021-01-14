import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from .PlotPlanet import PlotPlanetXY
from .ModelE import ModelECart
from .VExBModel import VExBModel


zlabs = {	'E':	r'Electric Field, $\mathbf{E}$, (mV m$^{-1}$)',
			'V':	r'Velocity, $\mathbf{V}$, (10,000 $\times$ m s$^{-1}$)',
			}


def _ModelGrid(Model,Rmax,dR,Kp,Esw,Vsw,Bz):
	#create a grid of points to calculate the model at
	nR = 2*np.int32(Rmax/dR)
	xe = np.linspace(-Rmax,Rmax,nR+1)
	ye = xe
	xc = xe[:-1] + 0.5*dR
	yc = ye[:-1] + 0.5*dR
	xe,ye = np.meshgrid(xe,ye)
	xc,yc = np.meshgrid(xc,yc)
	zc = np.zeros(xc.shape)
	R = np.sqrt(xc**2 + yc**2 + zc**2)
	
	#calculate the model
	if  'E' in Model:
		mx,my,mz = ModelECart(xc,yc,Kp,Esw,Vsw,Bz)
	elif 'V' in Model:
		mx,my,mz = VExBModel(xc,yc,zc,Kp,Esw,Vsw,Bz)
	m = np.sqrt(mx**2 + my**2 + mz**2)

	
	#remove stuff from within the planet
	bad = np.where(R < 1.0)
	m[bad] = np.nan
	mx[bad] = np.nan
	my[bad] = np.nan
	mz[bad] = np.nan
	
	return xc,yc,mx,my,m	


def PlotVectorEq(Model='V',Rmax=10.0,dR=1.0,Kp=1.0,Esw=None,Vsw=None,Bz=None,fig=None,maps=[1,1,0,0],
		ShowContour=True,zlog=True,cmap='gnuplot',scale=None,Verbose=False,fmt='%4.2f'):
	'''
	Plots various models in the equatorial plane.
	
	Inputs
	======
	Model : str
		'E' - Electric field
		'V' - ExB drift
	Rmax : float
		Limits of the plot
	dR : float
		Size of bins on plot
	Kp : float
		Kp index
	Esw : float
		Electric field due to the solar wind propagating past the 
		magnetosphere in mV/m.
	Vsw : float
		x-component (positive towards the Sun) of the solar wind 
		velocity in km/s.
	Bz : float
		Z-component of the IMF in nT.
	fig : object
		This should be either None (to create a new figure), a matplotlib.pyplot
		object (to add a new set of axes to an existing plot) or an 
		instance of matplotlib.pyplot.Axes in order to add to an existing
		set of axes.
	maps : list/tuple
		This will determine the position of the subplot on the figure.
	zlog : bool
		If True - the color scale will be logarithmic
	cmap : str
		String to say which colormap to use.
	scale : None or list/tuple
		If None - scale will be generated automatically, otherwise set 
		scale = [minimum,maximum]
	Verbose : bool
		If True, model calculation progress will be displayed	
	
	
	
	'''
	#get the model grid
	xc,yc,mx,my,m = _ModelGrid(Model,Rmax,dR,Kp,Esw,Vsw,Bz)
	
	#get the scale limits
	if scale is None:
		if zlog:
			scale = [np.nanmin(m[m > 0]),np.nanmax(m)]
		else:
			scale = [np.nanmin(m),np.nanmax(m)]
			
	zlabel = zlabs[Model]
	if zlog:
		norm = colors.LogNorm(vmin=scale[0],vmax=scale[1])
		lvl = 10**np.linspace(np.log10(scale[0]),np.log10(scale[1]),10)
	else:
		norm = colors.Normalize(vmin=scale[0],vmax=scale[1])	
		lvl = 10**np.linspace(scale[0],scale[1],10)
	

	

	#create the plot window
	if fig is None:
		fig = plt
		fig.figure()
	if hasattr(fig,'Axes'):
		ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
		ax.set_xlim([Rmax,-Rmax])
		ax.set_ylim([-Rmax,Rmax])
		ax.set_aspect(1.0)
	else:
		ax = fig		
		
	#plot the arrows
	if Model == 'V':
		qkl = 10000.0
	else:
		qkl = 1
	Q = ax.quiver(yc,xc,my,mx,m,cmap='jet',norm=norm)
	qk = ax.quiverkey(Q, 0.0, 1.02, qkl, zlabel + ', $K_p$ = {:3.1f}'.format(Kp), labelpos='E',
					   coordinates='axes')
	if ShowContour:
		#get a dense version of the model grid
		xcc,ycc,mxc,myc,mc = _ModelGrid(Model,Rmax,Rmax/100.0,Kp,Esw,Vsw,Bz)
		
		cs = ax.contour(ycc,xcc,mc,lvl,cmap='Greys',norm=norm)
		ax.clabel(cs, inline=1, fontsize=10,fmt=fmt)
		
	#plot the planet
	PlotPlanetXY(ax)


	#sort the axis labels
	ax.set_ylabel('$x_{SM}$ ($R_E$)')
	ax.set_xlabel('$y_{SM}$ ($R_E$)')	
	
	return ax
