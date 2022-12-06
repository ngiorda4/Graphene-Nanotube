#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 16:14:27 2022

@author: dimitris
"""
import math
from plot_nano import *
from misc import include_atoms_gr,  grid_pq, include_atoms_nt
import numpy as np


class CarbonFiller:
    """
    This is a class to construct carbon nanofillers, e.g. Graphene sheets and/
    or carbon nanotubes.
    """
    
    def __init__(self, n, m, l):
        """
        

        Parameters
        ----------
        n : int
            Number of hops along the a1 direction.
        m : int
            Number of hops along the a2 direction.
        l : int
            Length of the sheet or tube as measured in carbon-carbon bond lengths.

        Returns
        -------
        None.

        """
        self.name = 'Carbon Filler'
        self.n = n
        self.m = m
        self.l = l
        
    def vector(self):
        """
        

        Returns
        -------
        vec : numpy.array
            Coordinates of the vector Ch, which extends from the carbon
            at the origin to the final carbon after n hops in the a1 direction
            and m hops in the a2 direction.

        """
        x = round(self.n*np.sqrt(3) + (self.m*np.sqrt(3))/2, 3)                 # Must accomodate for half jumps; half jumps are determined by self.m 
        y = round(self.m*-1.5, 3) 
        vec = np.array([x, y])                                                  # turns x and y values into numpy array
        print(type(vec))
        return vec
    
    
    def TVector(self, Ch):   
        """
        

        Parameters
        ----------
        Ch : numpy.array
            Coordinates of vector Ch.

        Returns
        -------
        T : numpy.array
            Coordinates of vector T, which is perpendicular to vector Ch.

        """
        vector = self.normTvector(self.normVector(Ch)[0])                       # Only obtains the returned norm_vec from normVector
        T =  vector*self.l                                                      # Mutliplies each coordinate of the vector by the length
        return T
    
    @staticmethod # Does not use self
    def normVector(vec):
        """
        

        Parameters
        ----------
        vec : numpy.array
            Coordinates of vector Ch.

        Returns
        -------
        norm_vec : numpy.array
            Coordinates of the normalized vector c^.
        norm : float
            Length of vector Ch.

        """
        norm = round(np.sqrt(vec[0]**2 + vec[1]**2), 3)                         # Calculation of length from formula
        norm_vec = np.round(vec/norm, 3)                                        # Calculated by dividing the vector by the length of the vector
            
        return norm_vec, norm 
    
    @staticmethod # Does not use self
    def normTvector(c_hat):
        """
        

        Parameters
        ----------
        c_hat : numpy.array
            Coordinates of normalized vector c^.

        Returns
        -------
        t_hat : numpy.array
            Coordinates of vector t^ (-Ny, Nx).

        """
        t_hat = np.array([-c_hat[1], c_hat[0]])                                 # inverses x and y and makes y negative
        return t_hat
        
    @staticmethod # Does not use self  
    def pq(Ch, T):
        """
        

        Parameters
        ----------
        Ch : numpy.array
            Coordinates of vector Ch.
        T : numpy.array
            Coordinates of vector T.

        Returns
        -------
        array : numpy.array
            Two numpy arrays holding the values of p and q (of type int), respectively.

        """
        p_min = 0                                                               # Return to the numbering system 
        p_max = math.ceil((Ch[0] + T[0]) / (math.sqrt(3)*0.5))                  # Removing the scaling factor
        q_min = math.floor((Ch[1] / 1.5))                                       # qs are the Y coordinates of the first vector which is why they use the [1]; divided by 1.5 to revert back to numbering system
        q_max = math.ceil((T[1] / 1.5))                                         # qs are the Y coordinates of the second vector which is why they use the [1]; divided by 1.5 to revert back to numbering system
        array = np.array([[ p_min , p_max ], [ q_min , q_max ]])
        return array
    
    def grid_pq(p, q):
        """
        

        Parameters
        ----------
        p : numpy.array
            Minimum and maximum values of p.
        q : numpy.array
            Minimum and maximum values of q.

        Returns
        -------
        pg : 2-d numpy.array
            Integers used to identify each atom along the p-direction.
        qg : 2-d numpy.array
            Integers used to identify each atom along the q-direction.

        """
        grid_pq()                                                               # Simply calls the function from the misc.py file
        pg = pgrid
        qg = qgrid
        return pg, qg                                                           # returns two 2-d nump arrays holding the integers used identify each atom along the p- and q−direction
        
    
    @staticmethod  # Does not use self  
    def coordinates(pg, qg):
        """
        

        Parameters
        ----------
        pg : 2-d numpy.array
            Integers used to identify each atom along the p-direction.
        qg : 2-d numpy.array
            Integers used to identify each atom along the q-direction.

        Returns
        -------
        x : 2-d numpy.array
            X-coordinates of each atom.
        y : 2-d numpy.array
            Y-coordinates of each atom.

        """
        pg = np.array(pg)
        qg = np.array(qg)
        x =  np.zeros(pg.shape)                                                 # initializes an array of the same size with zeros
        y = np.zeros(qg.shape)                                                  # initializes an array of the same size with zeros
        
        totalRows, totalColumns = pg.shape
        
        for i in range(totalRows):
            for j in range(totalColumns):
                xCoord = pg[i,j] * (math.sqrt(3)/2)
                x[i, j] = np.round(xCoord, 3)                                   # rounds the array to 3 decimal places
                if (pg[i,j] + qg[i,j]) % 2 == 0:                                # checks if the sum is even, multiplies times 1.5
                    y[i,j] = qg[i,j] * 1.5      
                else:                                                           # if sum is odd, substracts 0.5 
                    y[i,j] = (qg[i,j] * 1.5) - 0.5
            
        return x, y
        
        

        
    @staticmethod # Does not use self    
    def distance(x, y, c_hat):
        """
        

        Parameters
        ----------
        x : 2-d numpy.array
            X-coordinates of each atom.
        y : 2-d numpy.array
            Y-coordinates of each atom.
        c_hat : numpy.array
            Coordinates of vector Ch.

        Returns
        -------
        s : 2-d numpy.array
            Distance of atoms from the mouth of the tube (x−axis).
        t : 2-d numpy.array
            Distance of atoms along the tube axis (y−axis).

        """
        array1 = (np.multiply(-c_hat[1], x)) + (np.multiply(c_hat[0], y))       # formula obtained from the project instructions # c_hat[0] = Nx, c_hat[1]
        array2 = (np.multiply(c_hat[0], x)) + (np.multiply(c_hat[1],y))         # formula obtained from the project instructions
        t = np.round(array1, 3)                                                 # rounds the array to 3 decimal places
        s = np.round(array2, 3)                                                 # rounds the array to 3 decimal places
        return s, t
    
    def include_atoms_gr(x, y, s, t, arclen, l):
        """
        

        Parameters
        ----------
        x : 2-d numpy.array
            X−coordinates of the atoms.
        y : 2-d numpy.array
            Y−coordinates of the atoms.
        s : 2-d numpy.array
            Distance of atoms from the mouth of the tube (x−axis).
        t : 2-d numpy.array
            Distance of atoms along the length direction (y−axis).
        arclen : float
            Length of vector Ch.
        l : float
            Length of the sheet, as measured in carbon-carbon bond lengths.

        Returns
        -------
        pos_gr : numpy.array
            Position of the atoms in the graphene sheet .

        """
        include_atoms_gr()                                                      # calls the function from the misc.py file
        pos_gr = np.hstack([x[include].reshape(-1, 1), y[include].reshape(-1, 1)]) # renames the output from the function in misc.py
        return pos_gr
    
    def include_atoms_nt(pos_gr, c_hat, arclen, tubrad):
        """
        

        Parameters
        ----------
        pos_gr : numpy.array
            X, Y coordinates of the atoms of the graphene sheet.
        c_hat : numpy.array
            Coordinates of the normal vector c^.
        arclen : float
            Length of vector Ch.
        tubrad : float
            Radius of the tube.

        Returns
        -------
        pos_nt : numpy.array
            Position of the atoms in the nanotube.

        """
        include_atoms_nt()                                                      # calls the function from the misc.py file
        pos_nt = np.vstack((pos_[0],pos_[1],pos_[2])).T                         # renames the output from the function in misc.py
        return pos_nt
    
def Graphene(n, m ,l):
    """
    

    Parameters
    ----------
    n : int
        Number of hops along the a1 direction.
    m : int
        Number of hops along the a2 direction.
    l : int
        Length of the sheet or tube as measured in carbon-carbon bond lengths.

    Returns
    -------
    pos_gr : numpy.array
        Position of the atoms in the graphene sheet.

    """
    Cf = CarbonFiller (n, m, l)                                                 # calls function and assigns its output to a variable
    Ch = Cf.vector ()                                                           # calls function and assigns its output to a variable
    T = Cf.TVector(Ch)                                                          # calls function and assigns its output to a variable
    p, q = Cf.pq(Ch , T)                                                        # calls function and assigns its output to variables
    Pgrid , Qgrid = grid_pq(p, q)                                               # calls function and assigns its output to variables
    x, y = Cf.coordinates (Pgrid , Qgrid)                                       # calls function and assigns its output to variables
    c_hat , arclen = Cf.normVector(Ch)                                          # calls function and assigns its output to variables
    s, t = Cf.distance(x, y, c_hat)                                             # calls function and assigns its output to variables
    pos_gr = include_atoms_gr(x, y, s, t, arclen , l)                           # calls function and assigns its output to a variable
    
    return pos_gr
    
    
    
if __name__ == '__main__':
     pass
      # pos_gr = Graphene(7,7,10)
      # low_th = 0.95
      # up_th = 1.05
      # xx,yy = get_connectivity_ind(pos_gr,low_th,up_th)
      # plot_mesh(xx,yy,pos_gr) 
    


