import numpy as np
import pandas

class Region:
    def __init__(self, eta_min=-2.5, eta_max=2.5, phi_min = -np.pi, phi_max = np.pi):
        self.eta = (eta_min, eta_max)
        self.phi = (phi_min, phi_max)

    def in_region(self, particles):
        # Returns a boolean array for each particle in list of particles
        is_in = np.logical_and(particles.eta >= self.eta[0], particles.eta < self.eta[1])
        is_in = np.logical_and(is_in, np.logical_and(particles.phi >= self.phi[0], particles.phi < self.phi[1]))
        return is_in

    def __str__(self):
        return "Region, (eta_min, eta_max), (phi_min, phi_max): ({}, {}), ({}, {})".format(self.eta[0], self.eta[1], self.phi[0], self.phi[1])

    def __repr__(self):
        return self.__str__()


def regionize(particles, regions):
    '''Return a list of DataFrames of the supplied particles partioned into the specified regions'''
    if len(regions.shape) == 1:
        return [particles[reg.in_region(particles)] for reg in regions]
    elif len(regions.shape) == 2:
        return [regionize(particles, reg1D) for reg1D in regions]

def merge1D(regionized, truncate=False):
    ''' Merge a list of DataFrame PF regions
        If an integer is supplied for truncate, the first 'truncate' elements only will be returned.'''
    merged = pandas.concat(regionized)
    if truncate:
        merged = merged.iloc[:truncate]
    return merged

def merge2D(regionized, truncate=False, axis=0):
    ''' Merge a 2D list of DataFrame PF regions.
        Two steps of merge1D merges are used, in the order specified by the axis argument.
        The truncate parameter can be either 'False' for no truncation, an int, or a tuple of separate
        truncation to apply at each truncation step.
    '''
    t0 = False
    t1 = False
    if truncate:
        if isinstance(truncate, int):
            t0, t1 = truncate, truncate
        elif isinstance(truncate, tuple) and len(truncate) == 2:
            t0, t1 = truncate[0], truncate[1]
        else:
            print("Invalid truncation provided, options are: False, integer, 2-tuple of integer")
    if axis==0:
        merged = [merge1D(reg, truncate=t0) for reg in regionized]
    else:
        merged = [merge1D(regionized[:,i], truncate=t0) for i in range(regionized.shape[0])]
    merged = merge1D(merged, truncate=t1)
    return merged



# Spec for 36 board Layer 1, ignoring forward
eta_edges = [-3, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3]
phi_edges = np.linspace(-np.pi, np.pi, 10)

PFRegions = np.array([[Region(eta_edges[i], eta_edges[i+1], phi_edges[j], phi_edges[j+1]) for j in range(len(phi_edges) - 1)] for i in range(len(eta_edges) - 1)])