#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 18:17:02 2022

@author: 
"""

import numpy as np
import matplotlib.pyplot as plt



def find_connectivity(points, low_th = 0.95, up_th = 1.05):
    """
    

    Parameters
    ----------
    points : 2-d numpy.array
        X and Y coordinates of points in 2D space.
    low_th : float, optional
        Lower accepted threshold values for the Euclidean distance squared value.
        The default is 0.95.
    up_th : float, optional
        Upper accepted threshold values for the Euclidean distance squared value. 
        The default is 1.05.

    Returns
    -------
    DistMat : 2-d numpy.array
        Euclidean distance squared of two points to decide the connectivity. If the distance is outside the upper and
        lower thresholds, this value is set to zero.

    """
    DistMat = np.zeros((points.shape[0], points.shape[0]))                      # initializes a square array based on the dimensions of the point arrays
    
    for i in range(len(points)):
        for j in range(len(points)):
            DistMat[i, j] = np.sum((points[i,:] - points[j,:])**2)              # finds the sum of the distances
            if DistMat[i, j] <= low_th or DistMat[i,j] >= up_th:                # checks if the sum of the distances are within the thresholds
                DistMat[i,j] = 0                                                # sets to zero if distance under or above threshold
                    
    return DistMat


def get_connectivity_ind(points, low_th = 0.95, up_th = 1.05):
    """
    

    Parameters
    ----------
    points : 2-d numpy.array
        X and Y coordinates of points in 2D space.
    low_th : float, optional
        Lower accepted threshold values for the Euclidean distance squared value.
        The default is 0.95.
    up_th : float, optional
        Upper accepted threshold values for the Euclidean distance squared value. 
        The default is 1.05.

    Returns
    -------
    index_points_limits : numpy.array
        Two numpy arrays which are indexes of points where the distance squared is within the defined limit. 
        First numpy array contains the indexes for the first points and the second array contains the indexes
        for the second set of points.

    """
    d = find_connectivity(points, low_th, up_th)                                # calls find_connectivity and assigns DistMat as variable d
    index_points_limits = np.where(d > 0)                                       # returns the indexes of where the distance is greater than 0  
    
    return index_points_limits
    
    


def plot_mesh(xx,yy,points):
    """
    

    Parameters
    ----------
    xx : numpy.array
        Indexes for the first points.
    yy : numpy.array
        Indexes for the second points.
    points : numpy.array
        X and Y coordinates of carbons in 2D space.

    Returns
    -------
    None.

    """
    
    for i in range(len(xx)): 
        p1 = points[xx[i],:]                                                    # assigns xx values at index i, as the range at which the points are recorded as p1 
        p2 = points[yy[i], :]                                                   # assigns yy values at index i, as the range at which the points are recorded as p1 
        plot_carbon_carbon(p1, p2)                                              # calls the plot_carbon_carbon
    plt.show()                                                                  # reveals the plot on a window

            
            

def plot_carbon_carbon(p1, p2):
    """
    

    Parameters
    ----------
    p1 : numpy.array
        Coordinates of point 1.
    p2 : numpy.array
        Coordinates of point 2.

    Returns
    -------
    None.

    """
    x = [p1[0], p2[0]]                                                          # extracts x values from both numpy coordinates of point 1 and point 2
    y = [p1[1], p2[1]]                                                          # extracts y values from both numpy coordinates of point 1 and point 2 
    plt.scatter(x, y, s = 75, c = 'red', marker = 'o')                          # plots two red points 
    plt.plot(x, y, linestyle = "-", color = "b")                                # plots the line between the two points 

if __name__ == '__main__':
    pass
    
    
  
  