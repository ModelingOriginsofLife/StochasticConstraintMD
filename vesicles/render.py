import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc
import matplotlib.cm as cm
import sys

def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence 
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circle in data unit.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
    plt.colorbar()

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None
    if 'fc' in kwargs: kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs: kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs: kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs: kwargs.setdefault('linewidth', kwargs.pop('lw'))

    patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(x, y, s)]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        collection.set_array(np.asarray(c))
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    if c is not None:
        plt.sci(collection)

    return collection

pts = np.loadtxt("parts.out")
bonds = np.loadtxt("bonds.out")

lbonds = []
for i in range(bonds.shape[0]):
	lbonds.append([(bonds[i,0],bonds[i,1]),(bonds[i,2],bonds[i,3])])

fig, ax = plt.subplots()

#~ out2 = circles(pts[:,0],pts[:,1],1.75,c=pts[:,2],cmap=cm.Dark2,alpha=0.1,ec=[0,0,0,0])

#~ out = circles(pts[:,0],pts[:,1],1.0,c=pts[:,2],cmap=cm.Dark2,ec=[0,0,0,0])
#~ lc = mc.LineCollection(lbonds,colors='r',alpha=0.5)
#~ ax.add_collection(lc)
#~ plt.xlim(0,5*20)
#~ plt.ylim(0,5*20)

data = np.loadtxt("avg.out")
data = data.reshape( (200,200,3) )/50.0
data[data>1]=1
data[data<0]=0

plt.xlim(0,200)
plt.ylim(0,200)
plt.imshow(data,interpolation='nearest')
frame = plt.gca()
frame.axes.get_xaxis().set_ticks([])
frame.axes.get_yaxis().set_ticks([])

plt.savefig(sys.argv[1],bbox_inches='tight')
