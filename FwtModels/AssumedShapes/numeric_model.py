import sympy as sym
from .fwt_parameters import FwtParameters
import numpy as np
import pandas as pd
from scipy.linalg import eig
from dataclasses import dataclass, InitVar, field

@dataclass
class NumericModel:
    M : np.array
    K : np.array
    B : np.array
    C : np.array
    G : np.array
    Z : any
    Alpha : any
    Params : FwtParameters

    def FreeVibrationVals(self):
        """Calculate the eigen vectors and values of the free vibration 
        problem
        """
        return eig(self.K,self.M)

    def FreeVibrationValsPd(self):
        """ takes the results of the free vibration problem and places them in
        a data frame"""
        eVals,eVecs = self.FreeVibrationVals()

        # convert vectors into a 2D array and append the Freq's on the front
        vecs = np.array(eVecs)
        data = np.insert(vecs,0,0,axis=1)
        data[:,0] = np.sqrt(np.real(eVals))

        # Make column names for the data frame 
        coords = self.M.shape[1]
        cols = ['Freq']
        for i in range(0,coords):
            cols.append(f'Comp {i}')
        # Create the Data Frame
        df = pd.DataFrame(data,columns=cols)

        # return a sorted list of modes with Nans removed
        return df.dropna().sort_values(by='Freq').reset_index(drop=True)

    def AeroVibrations(self):
        Mp = self.Mprime()
        Kp = self.Kprime()
        ev,eve = eig(Kp,Mp)
        return (ev,eve)

    def Mprime(self):
        coords = self.M.shape[1]
        Mprime = np.eye(coords*2)
        Mprime[coords:,coords:]=self.M
        return Mprime

    def Kprime(self):
        coords = self.M.shape[1]
        Kprime = np.zeros((coords*2,coords*2))
        Kprime[coords:,coords:]=self.B
        Kprime[0:coords,coords:]=np.eye(coords)
        Kprime[coords:,0:coords]=self.C-self.K
        return Kprime

    def MprimeSym(self):
        coords = self.M.shape[1]
        Mprime = sym.eye(coords*2)
        Mprime[coords:,coords:]=self.M
        return Mprime

    def KprimeSym(self):
        coords = self.M.shape[1]
        Kprime = sym.zeros(coords*2)
        Kprime[coords:,coords:]=self.B
        Kprime[0:coords,coords:]=np.eye(coords)
        Kprime[coords:,0:coords]=self.C-self.K
        return Kprime
