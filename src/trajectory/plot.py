"""
trajectory.plot
===============

A module for plotting trajectories. 
"""

__author__ = ("KT")
__version__ = "6.11, May 2013"



# python modules

# other modules
import numpy

import matplotlib.pyplot as pl

# custom made modules
import custom_plot

# package modules




grau = numpy.ones(3) * .4
grau2 = numpy.array([1, 1, 1]) * .95
fontsize = 18.0

def plotTrajectories_monochrome(traj, collTraj, rewardTraj, fig=None, ax=[], title='', scale_bars=True, offset=0):
    """ Plot trajectory data for paper figures.
    
    NOTE: Not flexibly programmed!
    """
    
    if not fig:
        fig = pl.figure(figsize=(8, 6))
    if not ax:
        ax = fig.add_axes([.2, .15, .81, .75])
    
    if not isinstance(ax, list):
        ax = [ax]

    # orient trajectory
    if traj.mazeType.find('linearMaze') + 1 and traj.euler[0] != 0:
        if traj:
            traj.orient()
        if collTraj:
            collTraj.orient()
        if rewardTraj:
            rewardTraj.orient()

    #--- plot trajectory
    if traj.mazeType.find('yDecision') + 1:
        rewardTraj.purgeRewards(traj.dt)  # sometimes reward postitions appear several times
                                            # due to higher sampling rate for rewards than for position
        if not traj.has_timeAxis:  # if traj does not include a recorded time axis, recalculate it according to reward times
            traj.recalcTimesAxis(rewardTraj.places, rewardTraj.times, 5)
        traj.plot(rewardTraj.times, fig, ax[0], offset, language='e', chic=True, monochrome=True)
    else:     
        ax[0].plot(traj.places[offset:, 0], traj.places[offset:, 1], '-', linewidth=1, color=grau)
    

    # add collisions
    collTraj_h = ax[0].plot(collTraj.places[:, 0], collTraj.places[:, 1], 'o', \
                         markerfacecolor=numpy.ones(3) * .65, markeredgecolor=numpy.ones(3) * .75)
    for c in collTraj_h:
        c.set_zorder(c.get_zorder() - .5)

    # add rewards
    if rewardTraj.places.size:
        rewardTraj_h = ax[0].plot(rewardTraj.places[:, 0], rewardTraj.places[:, 1], 'o', \
                markerfacecolor=numpy.ones(3) * 0, markeredgecolor=numpy.ones(3) * .1)
    
    # huebsch machen
    custom_plot.allOff(ax[0])
    ax[0].set_xlabel('')
    ax[0].set_ylabel('')
    ax[0].set_title(title, fontsize=fontsize)
    
    xoffset = traj.xWidth * .08
    yoffset = traj.yWidth * .2

    xmin = -xoffset
    xmax = traj.xWidth + xoffset
    ymin = traj.ylim[0] - yoffset
    ymax = traj.ylim[1] + yoffset

    ax[0].set_xlim(xmin, xmax)
    ax[0].set_ylim(ymin, ymax)

    # add scale bars
    if scale_bars:
        custom_plot.add_scalebar(ax[0], matchx=True, matchy=True, labelExtra=' m', \
            loc=3, borderpad=-1.5, sep=5)
    
    # draw rewardsPos
    for rPos in traj.rewardsPos:
        if traj.rewardsAreaType == 'circular':
            circ = pl.Circle(rPos, radius=traj.rewardsArea[0], color=custom_plot.grau2)
            ax[0].add_artist(circ)
            custom_plot.drop_shadow_patches(ax[0], circ, sigma=5, offsets=(2, 2))
    
    if traj.mazeType.find('linearMaze') + 1:
        #--- plot x-components
        xComponents = traj.getXComponents()
        ax[1].plot(traj.times, xComponents, '-', color=grau)
        # add rewards
        if rewardTraj.places.size:
            rewardTraj_h = ax[1].plot(rewardTraj.times, rewardTraj.getXComponents(), 'o', \
                markerfacecolor=numpy.ones(3) * 0, markeredgecolor=numpy.ones(3) * .1)

        # huebsch machen
        custom_plot.huebschMachen(ax[1])
        ax[1].set_xlabel('Time (' + traj.timeUnit + ')')
        ax[1].set_ylabel('x position (' + traj.spaceUnit + ')')
        
        minY = xComponents.min()
        maxY = xComponents.max()
        dY = 2.
        ax[1].set_yticks(numpy.arange(round(minY), round(maxY) + dY, dY))
        ax[1].set_ylim(minY - dY / 10, maxY + dY / 10)
        dx = numpy.round((int(traj.t_stop) / 60 * 60 - int(traj.t_start) / 60 * 60) / 4.)
        if not dx:
            dx = 100
        xticks = numpy.arange(round(traj.t_start), round(traj.t_stop) + 1, dx)
        ax[1].set_xticks(xticks)
        ax[1].set_xlim(traj.t_start, traj.t_stop)




def plotDatawpValues(arr, p=[], fig=None, ax=None):
    
    if not fig:
        fig = pl.figure()
    if not ax:
        ax = fig.add_subplot(111)
    
    if not len(p):
        p = numpy.ones([arr.__len__(), arr[0].__len__()])
    
    maxi = -10 ** 10
    maxLength = 0
    for i, v in enumerate(arr):        
        line = custom_plot.plotMarkerLinewpValues(ax, v, p[i], \
            color=custom_plot.colors[i], pValue_fontsize=custom_plot.pValue_fontsize)
        custom_plot.drop_shadow_line(ax, line)
    
        maxi = max(v.max(), maxi)
        maxLength = max(v.size, maxLength)

    return fig, ax


def plotLearningCurve(arr, fig=None, ax=None, showError=False, highlight_index=None):
    """Plot a learning curve for a certain learning parameter.
    
        arr ... array that contains data of the learning parameter ordered in time
                expected to be a list containing several datasets
        fig, ax ... matplotlib figure and axes to plot into
        showError ... show not only avg but also some error measure (i.e. SEM)
        highlight_index ... index into arr for the data set to be highlighted (plotted 
                            in a special way)

        Returns:
        fig, ax ... figure and axes handle for the matplotlib axes which was plotted into
    """
    
    if not fig:
        fig = pl.figure()
    if not ax:
        ax = fig.add_subplot(111)
    
    if not arr.__class__ == list:
        # in case
        arr = [arr]
    
    maxi = -10 ** 10
    maxLength = 0
    for i, v in enumerate(arr):
        ax.plot(v, '-', color=custom_plot.grau, markeredgecolor=custom_plot.grau, linewidth=custom_plot.linewidth / 2)
        if i == highlight_index:
            ax.plot(v, 'k--', linewidth=custom_plot.linewidth / 2)
        maxi = max(v.max(), maxi)
        maxLength = max(v.size, maxLength)

    # average
    avg = numpy.zeros(maxLength)
    std = numpy.zeros(maxLength)
    num = numpy.zeros(maxLength)
    for v in arr:
        for i, value in enumerate(v):
            if not numpy.isnan(value):
                avg[i] += value
                std[i] += value ** 2
                num[i] += 1
    avg /= num
    std /= num
    std -= avg ** 2
    std = numpy.sqrt(std)

    xvalues = numpy.arange(avg.size)
    ax.plot(xvalues, avg, 'k-')
    if showError:
        error = std / numpy.sqrt(num)
        ax.fill_between(xvalues, avg - error, avg + error, color=custom_plot.grau1, interpolate=True)
        ax.plot(xvalues, avg + error, '-', color=numpy.ones(3) * .3, linewidth=1)
        ax.plot(xvalues, avg - error, '-', color=numpy.ones(3) * .3, linewidth=1)
    
    custom_plot.huebschMachen(ax)
    ax.set_xlabel('Session #')
    
    return fig, ax

