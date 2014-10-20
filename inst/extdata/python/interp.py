#!/usr/bin/env python

import numpy as np
import scipy as sc
import scipy.interpolate

def interpolaten(points, values, newpoints, dim):
  npoints = int(len(points) / dim)
  nnewpoints = int(len(newpoints) / dim)

  points = np.array(points, dtype = float)
  newpoints = np.array(newpoints, dtype = float)
  values = np.array(values, dtype = float)

  points = np.ndarray(buffer = points, shape = (npoints, dim), order = 'F')
  values = np.ndarray(buffer = values, shape = (npoints, 1), order = 'F')
  newpoints = np.ndarray(buffer = newpoints, shape = (nnewpoints, dim), order = 'F')

  interpol = scipy.interpolate.LinearNDInterpolator(points, values)

  newvalues = interpol(newpoints)

  # Convert to list of doubles
  # This is necessary for correct transfer to R
  newvalues = map(float, newvalues)

  # Eval map (needed for Python 3)
  newvalues = list(newvalues)

  return newvalues
#
#
# dim = 3
# npoints = 1000
# mu = 0
# sigma = 1
# points = np.random.normal(mu, sigma, size = (npoints, dim))
#
# coef = np.random.normal(mu, sigma, (dim, 1))
#
# values = np.matrix(points, copy = False) * coef
# interpolaten(points, values, points, dim)

# dim = 3
# npoints = 1000000
# mu = 0
# sigma = 1
# points = np.random.normal(mu, sigma, size = (npoints, dim))
#
# coef = np.random.normal(mu, sigma, (dim, 1))
#
# values = np.matrix(points, copy = False) * coef
#
# interpol = scipy.interpolate.LinearNDInterpolator(points, values)
#
#
# new_points = np.random.normal(mu, sigma, size = (npoints, dim))
#
# new_values = np.matrix(new_points, copy = False) * coef
#
#
#
# np.allclose(values, interpol(points))
#
# int_new_values = interpol(new_points)
# np.allclose(new_values[~np.isnan(int_new_values)], int_new_values[~np.isnan(int_new_values)])
# new_values = np.matrix(new_points, copy = False) * coef

