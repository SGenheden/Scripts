# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Helper routines for plotting programs
"""

import numpy as np

def color(idx) :
  """
  Returns a color of index

  For instances when the index is larger than the number of defined colors,
  this routine takes care of this by periodicity, i.e.
  color at idx=0 is the same color as idx=n

  Parameters
  ----------
  idx : int
    the index of the color

  Returns
  -------
  list of floats
    the color
  """
  colors = []
  colors.append((0.0/255.0,69.0/255.0,134.0/255.0))
  colors.append((255.0/255.0,66.0/255.0,14.0/255.0))
  colors.append((255.0/255.0,211.0/255.0,32.0/255.0))
  colors.append((87.0/255.0,157.0/255.0,28.0/255.0))
  colors.append((126.0/255.0,0.0/255.0,33.0/255.0))
  colors.append((131.0/255.0,202.0/255.0,255.0/255.0))
  colors.append((49.0/255.0,64.0/255.0,4.0/255.0))
  colors.append((174.0/255.0,207.0/255.0,0.0/255.0))
  colors.append((74.0/255.0,35.0/255.0,109.0/255.0))
  colors.append((253.0/255.0,148.0/255.0,43.0/255.0))
  colors.append((197.0/255.0,0.0/255.0,12.0/255.0))
  colors.append((1.0/255.0,132.0/255.0,209.0/255.0))  
  d = int(len(colors)*np.floor(idx / float(len(colors))))
  return colors[idx-d]

def darker(color, size=0.90):

    return (size*color[0],
            size*color[1],
            size*color[2])

def style(idx) :

  styles = "-. : - --".split()
  d = int(len(styles)*np.floor(idx / float(len(styles))))
  return styles[idx-d]
