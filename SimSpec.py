# -*- coding: utf-8 -*-
"""Basic spectrum simulation.

This module simulates a spectrum with a slowly varying background and random
peaks of various height, width, and location. Once the signal expectation is
generated, one can add noise using a Gaussian noise model.
"""

__version__ = 1.0

import numpy as np
from scipy.stats import norm

def genBkg(x, nperiod, ampl=1., offset=0.):
    """Generate a sinusoidal background as a function of x.

    Args:
        x (numpy.array): independent ordinate.
        nperiod (float): number of periods of sinusoid in x range.
        ampl (float): amplitude of sinusoid.
        offset (float): DC offset of background.

    Returns:
        numpy.array: sinusoidal background.
    """
    freq = 2*np.pi*nperiod/(x[-1] - x[0])
    return offset + ampl * np.sin(freq*x)

def genPeaks(x, npeaks, ampl, width, loc=[0.,1.]):
    """Generate a spectrum of Gaussian peaks.

    Args:
        x (numpy.array): indepdent ordinate.
        npeaks (int): number of peaks to simulate.
        ampl (list): range of peak amplitudes, sampled uniformly.
        width (list): range of peak widths, sampled uniformly.
        loc (list): range of location, in fraction of x range [0..1].

    Returns:
        numpy.array: spectrum of peaks.
    """
    dx = x[-1] - x[0]
    y = np.zeros(x.shape, dtype=float)
    for peak in range(npeaks):
        A = np.random.uniform(ampl[0], ampl[1])
        mu = np.random.uniform(x[0]+loc[0]*dx, x[0]+loc[1]*dx)
        sigma = np.random.uniform(width[0], width[1])
        y += A * norm.pdf(x, loc=mu, scale=sigma)
    return y

def genNoise(x, ampl=1.):
    """Generate white noise (flat amplitude as a function of x).

    Args:
        x (numpy.array): independent ordinate.
        ampl (float): noise amplitude.

    Returns:
        numpy.array: noise spectrum as a function of x.
    """
    return ampl*np.random.randn(len(x))

