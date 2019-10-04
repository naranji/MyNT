"""
signale
=======

A collection of classes and functions to create, manipulate and play with 
analog signals and spike trains. 

Classes
-------


Functions
---------

"""

__author__ = ("KT", "Moritz Dittmeyer", "Franziska Hellmundt", "Christian Leibold")
__version__ = "4.0, May 2013"

# append parent folder to package path,
# such that custom modules located in parent folder can be found
# NOTE: not really nice, should maybe changed into global package
__path__.append(__path__[0].rpartition('/')[0])

# python modules

# other modules

# custom made modules

# package modules
from cscs import *
from io import *
from place_cells import *
from signals import *
from spikes import *
from tools import *




