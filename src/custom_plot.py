"""
A module for customizing matplotlib plots.
"""

__author__ = ("KT")
__version__ = "3.3, September 2014"

# python modules
import sys

# other modules
import numpy
import scipy.stats

# graphics modules
import matplotlib.pyplot as pl
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.collections import LineCollection
import matplotlib as mpl
import matplotlib.transforms as mtransforms
from matplotlib.artist import Artist
from matplotlib.offsetbox import AnchoredOffsetbox

from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes


colors = ['#0000FF','#FF0000','#008000','#BF00BF','#000000','#00FFFF','#EE82EE',
          '#808000','#800080','#FF6347','#FFFF00','#9ACD32','#4B0082',
          '#FFFACD','#C0C0C0','#A0522D','#FA8072','#FFEFD5','#E6E6FA',
          '#F1FAC1','#C5C5C5','#FF00FF','#A152ED','#FADD72','#F0EFD0','#EEE6FF',
          '#01FAC1','#F5F5F5','#A152FF','#FAFD72','#F0EFDF','#EEEFFF',
          '#F1FA99','#C9C9C9','#A152DD','#FA5572','#FFFFD0','#EDD6FF']

pretty_colors_set1 = [(0.8941176470588236, 0.10196078431372549, 0.10980392156862745),
 (0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
 (0.30196078431372547, 0.6862745098039216, 0.2901960784313726),
 (0.596078431372549, 0.3058823529411765, 0.6392156862745098),
 (1.0, 0.4980392156862745, 0.0),
 (1.0, 1.0, 0.2),
 (0.6509803921568628, 0.33725490196078434, 0.1568627450980392),
 (0.9686274509803922, 0.5058823529411764, 0.7490196078431373),
 (0.6, 0.6, 0.6)]

pretty_colors_set2 = [(0.4, 0.7607843137254902, 0.6470588235294118),
 (0.9882352941176471, 0.5529411764705883, 0.3843137254901961),
 (0.5529411764705883, 0.6274509803921569, 0.796078431372549),
 (0.9058823529411765, 0.5411764705882353, 0.7647058823529411),
 (0.6509803921568628, 0.8470588235294118, 0.32941176470588235),
 (1.0, 0.8509803921568627, 0.1843137254901961),
 (0.8980392156862745, 0.7686274509803922, 0.5803921568627451),
 (0.7019607843137254, 0.7019607843137254, 0.7019607843137254)]
 # (0.3, 0.8, 0.8),
 # (0.3, 1.0, 1.0),
 # (0.0, 0.38, 1.0),
 # (0.0, 1.0, 0.0),
 # (1.0, 0.0, 1.0)]

monochromes = ['#111111', '#444444', '#777777', '#AAAAAA']


hellgrau = numpy.array([1,1,1])*.5
dunkelgrau = numpy.array([1,1,1])*.35
grau = numpy.array([1,1,1])*.6
grau1 = numpy.array([1,1,1])*.85
grau2 = numpy.array([1,1,1])*.95
grau3 = numpy.array([1,1,1])*.7
grau4 = numpy.array([1,1,1])*.65

Bildformat='pdf'

fontsize = 16.0
pValue_fontsize = fontsize - 4
linewidth = 2
mpl.rcParams['lines.linewidth'] = linewidth
mpl.rcParams['font.size'] = fontsize
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.interactive(1)



###################################################### FUNCTIONS


def add_scalebar(ax, matchx=True, matchy=True, hidex=True, hidey=True, labelExtra='', **kwargs):
    # Adapted from mpl_toolkits.axes_grid2
    # LICENSE: Python Software Foundation (http://docs.python.org/license.html)
    """ Add scalebars to axes

    Adds a set of scale bars to *ax*, matching the size to the ticks of the plot
    and optionally hiding the x and y axes

    - ax : the axis to attach ticks to
    - matchx,matchy : if True, set size of scale bars to spacing between ticks
                    if False, size should be set using sizex and sizey params
    - hidex,hidey : if True, hide x-axis and y-axis of parent
    - **kwargs : additional arguments passed to AnchoredScaleBars

    Returns created scalebar object
    """
    def f(axis):
        l = axis.get_majorticklocs()
        return len(l)>1 and (l[1] - l[0])

    if matchx:
        kwargs['sizex'] = f(ax.xaxis)
        kwargs['labelx'] = str(kwargs['sizex'])+labelExtra
    if matchy:
        kwargs['sizey'] = f(ax.yaxis)
        kwargs['labely'] = str(kwargs['sizey'])+labelExtra

    sb = AnchoredScaleBar(ax.transData, **kwargs)
    ax.add_artist(sb)

    if hidex : ax.xaxis.set_visible(False)
    if hidey : ax.yaxis.set_visible(False)

    return sb


def allOff(ax):
    turnOffAxes(ax, ['left', 'right', 'bottom', 'top'])


def breakAxes(ax1, ax2, axes='x', d=.01):
        # d ... how big to make the diagonal lines in axes coordinates

        if axes == 'x':
            d1 = d / ax1.get_position().width
            d2 = d / ax2.get_position().width
            print d1, d2
            points = numpy.array([[(1-d1,1+d1), (-d1,+d1), (1-d1,1+d1),(1-d1,1+d1)],\
                [(-d2,+d2), (-d,+d2), (-d2,+d2), (1-d2,1+d2)]])
        elif 'y':
            points = numpy.array([[(-d,+d), (-d,+d), (1-d,1+d),(-d,+d)],\
                [(-d,+d), (1-d,1+d), (1-d,1+d), (1-d,1+d)]])

        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False,\
            linewidth=ax1.spines['bottom'].get_linewidth())
        ax1.plot(points[0][0], points[0][1], **kwargs)      # top-left diagonal
        ax1.plot(points[0][2], points[0][3], **kwargs)    # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot(points[1][0], points[1][1], **kwargs)      # top-left diagonal
        ax2.plot(points[1][2], points[1][3], **kwargs)    # top-right diagonal


def add_rotated_axes(fig, rect=111, angle=-45):
    """ Add a rotated axes to figure fig.
    """
    tr = Affine2D().scale(1, 1).rotate_deg(angle)

    grid_helper = floating_axes.GridHelperCurveLinear(tr, extremes=(0, 1, 0, 1))

    ax = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_axes(ax)

    aux_ax = ax.get_aux_axes(tr)

    #grid_helper.grid_finder.grid_locator1._nbins = 4
    #grid_helper.grid_finder.grid_locator2._nbins = 4

    return ax, aux_ax


def turnOffAxes(ax, spines=[]):

    # remove all axis spines
    for spine in spines:
        ax.spines[spine].set_visible(False)

        # run through all lines drawn for xticks and yticks and
        # remove ticks
        if spine == 'left':
            for i, line in enumerate(ax.get_yticklines()):
                if i%2 == 0:   # even indices
                    line.set_visible(False)
        elif spine == 'right':
            for i, line in enumerate(ax.get_yticklines()):
                if i%2 == 1:   # odd indices
                    line.set_visible(False)
        elif spine == 'bottom':
            for i, line in enumerate(ax.get_xticklines()):
                if i%2 == 0:   # even indices
                    line.set_visible(False)
        elif spine == 'top':
            for i, line in enumerate(ax.get_xticklines()):
                if i%2 == 1:   # odd indices
                    line.set_visible(False)

        # remove remaining tick labels
        if spine == 'bottom':
            ax.set_xticklabels([])
            ax.set_xlabel('')
        elif spine == 'left':
            ax.set_yticklabels([])
            ax.set_ylabel('')


def turnOnAxes(ax, spines=[]):

    # remove all axis spines
    for spine in spines:
        ax.spines[spine].set_visible(True)

        # run through all lines drawn for xticks and yticks and
        # remove ticks
        if spine == 'left':
            for i, line in enumerate(ax.get_yticklines()):
                if i%2 == 0:   # even indices
                    line.set_visible(True)
        elif spine == 'right':
            for i, line in enumerate(ax.get_yticklines()):
                if i%2 == 1:   # odd indices
                    line.set_visible(True)
        elif spine == 'bottom':
            for i, line in enumerate(ax.get_xticklines()):
                if i%2 == 0:   # even indices
                    line.set_visible(True)
        elif spine == 'top':
            for i, line in enumerate(ax.get_xticklines()):
                if i%2 == 1:   # odd indices
                    line.set_visible(True)

##        # remove remainig tick labels
##        if spine == 'bottom':
##            ax.set_xticklabels([])
##        elif spine == 'left':
##            ax.set_yticklabels([])




def avgData(*args):
    """ Average data from array y (might be of different length).

    Parameters
    ----------
    x : array_like
        Array of x values corresponding to y.
    y : array_like
        Array that contains data ordered like the values in x. It is expected to be a list containing several datasets.

    Returns
    -------
    avg, std, num [, x[index_maxLength]]
    """

    if args.__len__() == 1:
        y = args[0]
    elif args.__len__() == 2:
        x = args[0]
        y = args[1]
    else:
        print "ERROR: Too many arguments for custom_plot.avgData()!"
        return 0

    maxLength = 0
    index_maxLength = 0
    for i, v in enumerate(y):
        dummy = maxLength
        maxLength = max(v.size, maxLength)
        if dummy != maxLength:
            index_maxLength = i

    # average
    avg = numpy.zeros(maxLength)
    std = numpy.zeros(maxLength)
    num = numpy.zeros(maxLength)
    for v in y:
        for i, value in enumerate(v):
            if not numpy.isnan(value):
                avg[i] += value
                std[i] += value**2
                num[i] += 1
    avg /= num
    std /= num
    std -= avg**2
    std = numpy.sqrt(std)

    if args.__len__() == 2:
        return avg, std, num, x[index_maxLength]
    else:
        return avg, std, num


def avgPlot(ax, xvalues, avg, err=[], avg_col=[0., 0., 0.], err_col=[.75, .75, .75], shadow=False, linewidth=2.5, zorder=0):

    line, = ax.plot(xvalues, avg, '-', color=avg_col, linewidth=linewidth, zorder=zorder)
    if shadow:
        drop_shadow_line(ax, line)

    if err.__len__():
        ax.fill_between(xvalues, avg-err, avg+err, color=err_col, alpha=.25, interpolate=True, zorder=zorder)
        # line, = ax.plot(xvalues, avg+err, '-', color=err_col, linewidth=1, zorder=zorder)
        if shadow:
            drop_shadow_line(ax, line)
        # line, = ax.plot(xvalues, avg-err, '-', color=err_col, linewidth=1, zorder=zorder)
        if shadow:
            drop_shadow_line(ax, line)

    huebschMachen(ax)


def plotMarkerLine(ax, data, x=[], color='b', markersize=10, zorder=0):
    if not len(x):
        x = numpy.arange(len(data))
    line, = ax.plot(x, data, 'o-',
                markerfacecolor='w',
                color=color,
                markeredgewidth=mpl.rcParams['lines.linewidth']-1,
                markersize=markersize, zorder=zorder)
    line.set_markeredgecolor(line.get_color())
    return line


def plotMarkerLinewErrors(ax, data, x=[], yerr=[], color='b', markersize=10, zorder=2):
    if not len(x):
        x = numpy.arange(len(data))
    if color:
        line = ax.errorbar(x, data, yerr=yerr,
                marker='o', linestyle='-', elinewidth=mpl.rcParams['lines.linewidth']-1,
                color=color, capsize=5, capthick=mpl.rcParams['lines.linewidth']-1, markerfacecolor='w',
                markeredgewidth=mpl.rcParams['lines.linewidth']-1,
                markersize=markersize, zorder=zorder)
    else:
        line = ax.errorbar(x, data, yerr=yerr,
                marker='o', linestyle='-', elinewidth=mpl.rcParams['lines.linewidth']-1,
                capsize=5, capthick=mpl.rcParams['lines.linewidth']-1, markerfacecolor='w',
                markeredgewidth=mpl.rcParams['lines.linewidth']-1,
                markersize=markersize, zorder=zorder)
    line = line[0]
    line.set_markeredgecolor(line.get_color())
    return line


def plotMarkerLinewpValues(ax, data, pValues, x=[], color='b', stars=False, pValue_fontsize=14):

    if not len(x):
        x = numpy.arange(len(data))
    line = plotMarkerLine(ax, data, x)
    line.set_marker('s')
    markersize = line.get_markersize()
    line.set_markersize(markersize/2.)
    line.set_alpha(.75)
    line.set_color(color)
    line.set_markeredgecolor(color)
    xticks = x
    for j, p in enumerate(pValues):
        if stars:
            ax.text(xticks[j], data[j]+1, pValues2stars(p), fontsize=pValue_fontsize, horizontalalignment='center')
        else:
            if p <= 0.05:
##                col = numpy.array(mpl.colors.colorConverter.to_rgba(line.get_color()))
##                col += pValues2rgba(p)
##                col = col.clip(0, 1)
##                ax.plot(xticks[j], data[j], 'o-', color = col,\
##                    markeredgecolor=col,\
##                    markerfacecolor='w',\
##                    markeredgewidth=mpl.rcParams['lines.linewidth']-1,
##                    markersize=10)
                l1, = ax.plot(xticks[j], data[j], 'o-', color = line.get_color(),
                    markeredgecolor = line.get_color(),
                    markerfacecolor = 'w',
                    markeredgewidth=mpl.rcParams['lines.linewidth']-1,
                    markersize=markersize*pValues2markersize(p))
                #drop_shadow_line(ax, l1)

    return line


def plotPooledData(x, y, fig=None, ax=None, showError=None, highlight_index=None, color=None, zorder=0):
    """ Plot a learning curve for a certain learning parameter.

    Parameters
    ----------
    x : ndarray
        Array of x values corresponding to y.
    y : list
        Array that contains data ordered like the values in x. It is expected to be a list containing several datasets.
    fig, ax : matplotlib figure and axes, optional
        For plotting the data into.
    showError : bool, optional
        Show not only avg but also some error measure (i.e. SEM).
    highlight_index : int, optional
        Index into arr for the data set to be highlighted (plotted in a special way).

    Returns
    -------
    fig, ax : figure and axes handle for the matplotlib axes into which the data was plotted.
    """

    # initialization
    if not fig:
        fig = pl.figure()
    if not ax:
        ax = fig.add_subplot(111)

    if not x.__class__ == list:
        x = [x]
    if not y.__class__ == list:
        y = [y]

    if not color:
        color = grau


    for i, v in enumerate(y):
        ax.plot(x[i], v, '-', color=color, markeredgecolor=color, linewidth=linewidth/2, zorder=zorder)
        if i == highlight_index:
            ax.plot(x[i], v, 'k--', linewidth=linewidth/2, zorder=zorder)


    avg, std, num = avgData(x, y)

    maxX = numpy.argmax([i.__len__() for i in x])
    xvalues = x[maxX]#numpy.arange(avg.size)
    ax.plot(xvalues, avg, 'k-', linewidth=linewidth/2, zorder=zorder)
    if showError:
        if showError == 'std':
            error = std
        elif showError == 'sem':
            error = std/numpy.sqrt(num)
        ax.fill_between(xvalues, avg-error, avg+error, color=color, alpha=.25, interpolate=True, zorder=zorder)
        ax.plot(xvalues, avg+error, '-', color=scipy.stats.threshold(numpy.array(color) - numpy.ones(3)*.3, threshmin=0),
                linewidth=1, zorder=zorder)
        ax.plot(xvalues, avg-error, '-', color=scipy.stats.threshold(numpy.array(color) - numpy.ones(3)*.3, threshmin=0),
                linewidth=1, zorder=zorder)

    huebschMachen(ax)

    return avg, error


def huebschMachen(ax):
    turnOffAxes(ax, ['right', 'top'])


def drop_shadow_line(ax, line, sigma=3, xOffset=1., yOffset=-1.):
    # adapted from examples/misc/svg_filter_line.py

    gauss = DropShadowFilter(sigma)

    # draw shadows with same lines with slight offset.
    xx = line.get_xdata()
    yy = line.get_ydata()
    shadow, = ax.plot(xx, yy)
    shadow.update_from(line)

    # offset transform
    ot = mtransforms.offset_copy(line.get_transform(), ax.figure,\
                                    x=xOffset, y=yOffset, units='points')
    shadow.set_transform(ot)


    # adjust zorder of the shadow lines so that it is drawn below the
    # original lines
    shadow.set_zorder(line.get_zorder()-0.75)
    shadow.set_agg_filter(gauss)
    shadow.set_rasterized(True) # to support mixed-mode renderers


def drop_shadow_patches(ax, patches, sigma=5, offsets=(1,1)):
    # adapted from barchart_demo.py

    if not type(patches) == list:
        patches = [patches]
    gauss = DropShadowFilter(sigma, offsets=offsets)

    for patch in patches:
        shadow = FilteredArtistList([patch], gauss)
        ax.add_artist(shadow)
        shadow.set_zorder(patch.get_zorder()-.1)



###################################################### TOOLS


def brighten(color, value=.1):
    if not type(color) == tuple:
        sys.exit('Color style not implemented!')
    color = numpy.array(color) + numpy.ones(3)*value

    return color.clip(0, 1)



def pValues2stars(pValues):
    pValues_stars = []
    if not type(pValues) == list and not type(pValues) == numpy.ndarray:
        pValues = [pValues]
    for p in pValues:
        star = ''
        if p <= .05:
            star = '*'
        if p <= .01:
            star = '**'
        if p <= .001:
            star = '***'
        pValues_stars.append(star)

    # if pValues_stars.__len__() == 1:
    #     return pValues_stars[0]
    # else:
    return numpy.array(pValues_stars)


def pValues2rgba(pValues):
    pValues_col = []
    col = [0, 0, 0, 0]
    if not type(pValues) == list or not type(pValues) == numpy.ndarray:
        pValues = [pValues]
    for p in pValues:
        if p <= .05:
            col = [0, 0, 0, 1]
        if p <= .01:
            col = [-.33, -.33, -.33, 1]
        if p <= .001:
            col = [-.75, -.75, -.75, 1]
        pValues_col.append(col)

    if pValues_col.__len__() == 1:
        return pValues_col[0]
    else:
        return numpy.array(pValues_col)


def pValues2markersize(pValues):
    pValues_m = []
    if not type(pValues) == list or not type(pValues) == numpy.ndarray:
        pValues = [pValues]
    for p in pValues:
        if p <= .05:
            m = 1.2
        if p <= .01:
            m = 1.5
        if p <= .001:
            m = 1.9
        pValues_m.append(m)

    if pValues_m.__len__() == 1:
        return pValues_m[0]
    else:
        return numpy.array(pValues_m)


def smooth1d(x, window_len):
    # copied from http://www.scipy.org/Cookbook/SignalSmooth

    s = numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    w = numpy.hanning(window_len)
    y = numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def smooth2d(A, sigma=3):

    window_len = max(int(sigma), 3)*2+1
    A1 = numpy.array([smooth1d(x, window_len) for x in numpy.asarray(A)])
    A2 = numpy.transpose(A1)
    A3 = numpy.array([smooth1d(x, window_len) for x in A2])
    A4 = numpy.transpose(A3)

    return A4


###################################################### CLASSES

class BaseFilter(object):
    def prepare_image(self, src_image, dpi, pad):
        ny, nx, depth = src_image.shape
        #tgt_image = numpy.zeros([pad*2+ny, pad*2+nx, depth], dtype="d")
        padded_src = numpy.zeros([pad*2+ny, pad*2+nx, depth], dtype="d")
        padded_src[pad:-pad, pad:-pad,:] = src_image[:,:,:]

        return padded_src#, tgt_image

    def get_pad(self, dpi):
        return 0

    def __call__(self, im, dpi):
        pad = self.get_pad(dpi)
        padded_src = self.prepare_image(im, dpi, pad)
        tgt_image = self.process_image(padded_src, dpi)
        return tgt_image, -pad, -pad


class OffsetFilter(BaseFilter):

    def __init__(self, offsets=None):
        if offsets is None:
            self.offsets = (0, 0)
        else:
            self.offsets = offsets

    def get_pad(self, dpi):
        return int(max(*self.offsets)/72.*dpi)

    def process_image(self, padded_src, dpi):
        ox, oy = self.offsets
        a1 = numpy.roll(padded_src, int(ox/72.*dpi), axis=1)
        a2 = numpy.roll(a1, -int(oy/72.*dpi), axis=0)
        return a2


class GaussianFilter(BaseFilter):
    """ Simple gaussian filter.
    """

    def __init__(self, sigma, alpha=0.5, color=None):
        self.sigma = sigma
        self.alpha = alpha
        if color is None:
            self.color=(0, 0, 0)
        else:
            self.color=color

    def get_pad(self, dpi):
        return int(self.sigma*3/72.*dpi)


    def process_image(self, padded_src, dpi):
        #offsetx, offsety = int(self.offsets[0]), int(self.offsets[1])
        tgt_image = numpy.zeros_like(padded_src)
        aa = smooth2d(padded_src[:,:,-1]*self.alpha,
                      self.sigma/72.*dpi)
        tgt_image[:,:,-1] = aa
        tgt_image[:,:,:-1] = self.color
        return tgt_image


class DropShadowFilter(BaseFilter):

    def __init__(self, sigma, alpha=0.3, color=None, offsets=None):
        self.gauss_filter = GaussianFilter(sigma, alpha, color)
        self.offset_filter = OffsetFilter(offsets)

    def get_pad(self, dpi):
        return max(self.gauss_filter.get_pad(dpi),
                   self.offset_filter.get_pad(dpi))

    def process_image(self, padded_src, dpi):
        t1 = self.gauss_filter.process_image(padded_src, dpi)
        t2 = self.offset_filter.process_image(t1, dpi)
        return t2



class FilteredArtistList(Artist):
    """ A simple container to draw filtered artist.
    """
    def __init__(self, artist_list, filter):
        self._artist_list = artist_list
        self._filter = filter
        Artist.__init__(self)

    def draw(self, renderer):
        renderer.start_rasterizing()
        renderer.start_filter()
        for a in self._artist_list:
            a.draw(renderer)
        renderer.stop_filter(self._filter)
        renderer.stop_rasterizing()



class AnchoredScaleBar(AnchoredOffsetbox):
    """ Scale bar class.

    Adapted from mpl_toolkits.axes_grid2
    LICENSE: Python Software Foundation (http://docs.python.org/license.html)
    """

    def __init__(self, transform, sizex=0, sizey=0, labelx=None, labely=None, loc=4,
                 pad=0.1, borderpad=0.1, sep=2, thickness=0, prop=None, **kwargs):
        """
        Draw a horizontal and/or vertical  bar with the size in data coordinate
        of the given axes. A label will be drawn underneath (center-aligned).

        - transform : the coordinate frame (typically axes.transData)
        - sizex,sizey : width of x,y bar, in data units. 0 to omit
        - labelx,labely : labels for x,y bars; None to omit
        - loc : position in containing axes
        - pad, borderpad : padding, in fraction of the legend font size (or prop)
        - sep : separation between labels and bars in points.
        - **kwargs : additional arguments passed to base class constructor
        """
        from matplotlib.patches import Rectangle
        from matplotlib.offsetbox import AuxTransformBox, VPacker, HPacker, TextArea
        bars = AuxTransformBox(transform)
        if sizex:
            bars.add_artist(Rectangle((0,0), sizex, thickness, fc="k"))
        if sizey:
            bars.add_artist(Rectangle((0,0), thickness, sizey, fc="k"))

        if sizex and labelx:
            bars = VPacker(children=[bars, TextArea(labelx, minimumdescent=False)],
                           align="center", pad=0, sep=sep)
        if sizey and labely:
            bars = HPacker(children=[TextArea(labely), bars],
                            align="center", pad=0, sep=sep)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=bars, prop=prop, frameon=False, **kwargs)


