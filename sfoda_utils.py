# -*- coding: utf-8 -*-
"""

Tools for handling SUNTANS output data

Created on Mon Sep 24 16:55:45 2012

@author: mrayson
"""

from netCDF4 import MFDataset, Dataset

import numpy as np
import os

# Constants
GRAV = 9.81
FILLVALUE = -999999

###############################################################
# Dictionary with lookup table between object variable name and netcdf file
# variable name
suntans_gridvars = {'xp': 'xp',
                    'yp': 'yp',
                    'xv': 'xv',
                    'yv': 'yv',
                    'xe': 'xe',
                    'ye': 'ye',
                    'lonv': 'lonv',
                    'latv': 'latv',
                    'lonp': 'lonp',
                    'latp': 'latp',
                    'cells': 'cells',
                    'face': 'face',
                    'nfaces': 'nfaces',
                    'edges': 'edges',
                    'neigh': 'neigh',
                    'grad': 'grad',
                    # 'gradf':'gradf',
                    'mark': 'mark',
                    'normal': 'normal',
                    'mnptr': 'mnptr',
                    'eptr': 'eptr',
                    'n1': 'n1',
                    'n2': 'n2',
                    'df': 'df',
                    'dg': 'dg',
                    'def': 'def',
                    'Ac': 'Ac',
                    'dv': 'dv',
                    'dz': 'dz',
                    'z_r': 'z_r',
                    'z_w': 'z_w',
                    'Nk': 'Nk',
                    'Nke': 'Nke',
                    'time': 'time'
                    }
suntans_dimvars = {'Np': 'Np',
                   'Ne': 'Ne',
                   'Nc': 'Nc',
                   'Nkmax': 'Nk',
                   'Nk': 'Nk',
                   'maxfaces': 'numsides'
                   }


########
# Utility functions
########
def calc_z(dz):
    z_bot = np.cumsum(dz)
    z_top = np.hstack((0.0, z_bot[:-1]))
    return 0.5 * (z_bot + z_top)


##########
# Classes
##########
class Grid(object):
    """ Class for handling SUNTANS grid data"""

    MAXFACES = 3  # Default number of faces
    _FillValue = FILLVALUE
    gridvars = suntans_gridvars
    griddims = suntans_dimvars

    # Some grid properties
    DEF = None

    VERBOSE = True

    ###
    # Grid projection details
    projstr = None
    utmzone = None
    isnorth = None

    def __init__(self, infile, **kwargs):

        self.__dict__.update(kwargs)

        if isinstance(infile, list):
            infile2 = infile[0]
        else:
            infile2 = infile

        if os.path.isdir(infile2):
            # Load ascii grid file
            self.infile = infile
            self.__loadascii()

        else:
            # Load the grid fromm a netcdf file
            self.infile = infile
            self.ncfile = infile
            self.__loadGrdNc()

        # Find the grid limits
        self.xlims = [self.xp.min(), self.xp.max()]
        self.ylims = [self.yp.min(), self.yp.max()]

        # Cell polygon attribute for plotting
        # self.cells[self.cells.mask]=0
        # self.grad[self.grad.mask]=0
        # self.face[self.face.mask]=0
        self.xy = self.cellxy(self.xp, self.yp)

    def __loadascii(self):
        """
        Load the grid variables from the ascii files: points.dat, edges.dat, cells.dat
        """
        pointdata = readTXT(self.infile + '/points.dat')
        celldata = readTXT(self.infile + '/cells.dat')
        edgedata = readTXT(self.infile + '/edges.dat')

        self.xp = pointdata[:, 0]
        self.yp = pointdata[:, 1]
        # self.dv = pointdata[:,2] # zero to start
        self.Np = len(self.xp)

        # Work out if cells.dat is in the quad grid format based on number of
        # columns
        self.Nc = celldata.shape[0]
        if celldata.ndim == 2:
            ncols = celldata.shape[1]
            # print '!!!cells.dat has %d columns!!!'%ncols
            if ncols == 8:  # Old format
                self.xv = celldata[:, 0]
                self.yv = celldata[:, 1]
                self.cells = np.asarray(celldata[:, 2:5], int)
                self.neigh = np.asarray(celldata[:, 5:8])
                self.nfaces = 3 * np.ones((self.Nc,), np.int)
                self.maxfaces = self.MAXFACES
            elif ncols == 9:  # New format
                nfaces = celldata[:, 0]
                self.nfaces = np.zeros((self.Nc,), np.int)
                self.nfaces[:] = nfaces
                self.xv = celldata[:, 1]
                self.yv = celldata[:, 2]
                self.cells = np.asarray(celldata[:, 3:6], int)
                self.neigh = np.asarray(celldata[:, 6:9])
                self.nfaces = 3 * np.ones((self.Nc,), np.int)
                self.maxfaces = self.MAXFACES

            elif ncols == 11:  # Quad grid format
                nfaces = celldata[:, 0]
                self.nfaces = np.zeros((self.Nc,), np.int)
                self.nfaces[:] = nfaces
                self.xv = celldata[:, 1]
                self.yv = celldata[:, 2]
                self.cells = np.asarray(celldata[:, 3:7], int)
                self.neigh = np.asarray(celldata[:, 7:11])
                self.maxfaces = 4
        else:  # Uneven number of cells
            celldata = celldata.tolist()
            nfaces = [ll[0] for ll in celldata]
            self.nfaces = np.array(nfaces, np.int)
            self.maxfaces = self.nfaces.max()
            self.cells = self._FillValue * np.ones((self.Nc, self.maxfaces), int)
            self.neigh = self._FillValue * np.ones((self.Nc, self.maxfaces), int)
            self.xv = np.zeros((self.Nc,))
            self.yv = np.zeros((self.Nc,))
            for ii in range(self.Nc):
                nf = self.nfaces[ii]
                self.xv[ii] = celldata[ii][1]
                self.yv[ii] = celldata[ii][2]
                self.cells[ii, 0:nf] = celldata[ii][3:3 + nf]
                self.neigh[ii, 0:nf] = celldata[ii][3 + nf:3 + 2 * nf]

        self.edges = np.asarray(edgedata[:, 0:2], int)
        self.Ne = self.edges.shape[0]
        self.mark = np.asarray(edgedata[:, 2], int)
        self.grad = np.asarray(edgedata[:, 3:5], int)
        if np.size(edgedata, 1) == 6:
            self.edge_id = np.asarray(edgedata[:, 5], int)
        else:
            self.edge_id = np.zeros((self.Ne,), int)

        # Load the vertical grid info from vertspace.dat if it exists
        try:
            vertspace = readTXT(self.infile + '/vertspace.dat')
        except:
            if self.VERBOSE:
                print('Warning could not find vertspace.dat in folder, setting Nkmax=1')
            vertspace = 0.0

        self.setDepth(vertspace)

        self.maskgrid()

    def __loadGrdNc(self):

        """
        Load the grid variables into the object from a netcdf file

        Try statements are for backward compatibility

        Variables loaded are presently:
        'xp','yp','xv','yv','xe','ye','cells','face','nfaces','edges','neigh','grad',
        'gradf','mark','normal','n1','n2','df','dg','def','Ac','dv','dz','z_r','z_w','Nk','Nke'
        """

        self.__openNc()
        nc = self.nc

        # Get the dimension sizes
        for vv in list(self.griddims.keys()):
            try:
                setattr(self, vv, nc.dimensions[self.griddims[vv]].__len__())
            except:
                if self.VERBOSE:
                    print('Cannot find dimension: %s' % self.griddims[vv])

        for vv in list(self.gridvars.keys()):
            try:
                if vv == 'def':  # Cannot have this attribute name in python!
                    setattr(self, 'DEF', nc.variables[self.gridvars[vv]][:])
                else:
                    setattr(self, vv, nc.variables[self.gridvars[vv]][:])
            except:
                if self.VERBOSE:
                    print('Cannot find variable: %s' % self.gridvars[vv])

        if 'Nk' in self.__dict__:
            self.Nk -= 1  # These need to be zero based

        if 'nfaces' not in self.__dict__:
            self.MAXFACES = self.cells.shape[1]
            self.nfaces = self.MAXFACES * np.ones((self.Nc,), np.int)
            self.maxfaces = self.MAXFACES

        # If edges, grad or neigh have not been stored then calculate them
        if 'edges' not in self.__dict__:
            self.reCalcGrid()
        elif 'grad' not in self.__dict__:
            self.reCalcGrid()
        # elif not self.__dict__.has_key('neigh'):
        #    self.reCalcGrid()

        # Set the mark equal zero if doesn't exist
        if 'mark' not in self.__dict__:
            self.mark = np.zeros((self.Ne))

        # Check the _FillValue attribute is consistent with the grid
        # If the maximum cells value exceeds the number of points the fillvalue
        # is likely to be  999999
        if self.cells.max() > self.Np:
            if self.VERBOSE:
                print('Changing the _FillValue from {} to {}'.format(
                    self._FillValue, self.cells.max()))
            self._FillValue = self.cells.max()

        if type(self.cells) != type(np.ma.MaskedArray()):
            self.maskgrid()
        else:
            self.cellmask = self.cells.mask

        # if type(self.DEF) == type(np.ma.MaskedArray()):
        #    if np.all(self.DEF.mask):
        #        self.calc_def()
        try:
            self.calc_def()
        except:
            print('No def array...')

    def maskgrid(self):
        """
        Mask the cells, face and neigh arrays
        """
        self.cellmask = self.cells == int(self._FillValue)
        # for ii in range(self.Nc):
        #    self.cellmask[ii,self.nfaces[ii]::]=True

        self.cells[self.cellmask] = 0
        self.cells = \
            np.ma.masked_array(self.cells, mask=self.cellmask, fill_value=0)

        if 'face' in self.__dict__:
            self.face[self.cellmask] = 0
            self.face = \
                np.ma.masked_array(self.face, mask=self.cellmask, fill_value=0)

        if 'neigh' in self.__dict__:
            self.neigh = \
                np.ma.masked_array(self.neigh, mask=self.cellmask, fill_value=0)

    def calc_def(self):
        """
        Recalculate the edge to face distance
        """
        ne = np.array(self.face)

        # try:
        #    mask = ne.mask.copy()
        # except:
        #    mask = np.zeros(self.DEF.shape,np.bool)
        ne[self.cellmask] = 0

        def dist(x0, x1, y0, y1):
            return np.sqrt((x0 - x1) ** 2. + (y0 - y1) ** 2.)

        self.DEF = dist(self.xv, self.xe[ne].T, self.yv, self.ye[ne].T).T

        self.DEF = np.ma.masked_array(self.DEF, mask=self.cellmask)

    def cellxy(self, xpin, ypin):
        """
        Returns a list of Nx2 vectors containing the grid cell node coordinates

        Used by spatial ploting routines
        """
        xp = np.zeros((self.Nc, self.maxfaces + 1))
        yp = np.zeros((self.Nc, self.maxfaces + 1))

        cells = self.cells.copy()
        # cells[self.cells.mask]=0
        cellmask = cells == int(self._FillValue)
        cells[cellmask] = 0
        # print(self._FillValue)

        xp[:, :self.maxfaces] = xpin[cells]
        xp[list(range(self.Nc)), self.nfaces] = xpin[cells[:, 0]]
        yp[:, :self.maxfaces] = ypin[cells]
        yp[list(range(self.Nc)), self.nfaces] = ypin[cells[:, 0]]

        # xp[self.cells.mask]==0
        # yp[self.cells.mask]==0

        xy = np.zeros((self.maxfaces + 1, 2))

        def _closepoly(ii):
            nf = self.nfaces[ii] + 1
            xy[:nf, 0] = xp[ii, :nf]
            xy[:nf, 1] = yp[ii, :nf]
            return xy[:nf, :].copy()

        return [_closepoly(ii) for ii in range(self.Nc)]

        # Old Method
        # return [closePoly(self.xp[self.cells[ii,0:self.nfaces[ii]]],\
        #    self.yp[self.cells[ii,0:self.nfaces[ii]]]) for ii in range(self.Nc)]

    def saveBathy(self, filename):
        """
            Saves the grid bathymetry to an xyz ascii file
        """
        f = open(filename, 'w')

        for x, y, z in zip(self.xv, self.yv, self.dv):
            f.write('%10.6f %10.6f %10.6f\n' % (x, y, z))

        f.close()

    def loadBathy(self, filename):
        """
        Loads depths from a text file into the attribute 'dv'
        """
        depths = readTXT(filename)
        if len(depths) != self.Nc:
            print('Error - number of points in depth file (%d) does not match Nc (%d)' % (len(depths), self.Nc))
        else:
            dv = depths[:, 2]

        self.dv = dv
        return dv

    def calc_dg(self):
        """
        Manually calculate the distance between voronoi points, 'dg'
        """
        if self.VERBOSE:
            print('Calculating dg...')
            print(np.shape(self.grad))

        grad = self.grad
        Ne = len(grad)
        for ii in range(Ne):
            if grad[ii, 0] == -1:
                grad[ii, 0] = grad[ii, 1]
            elif grad[ii, 1] == -1:
                grad[ii, 1] = grad[ii, 0]

        x1 = self.xv[grad[:, 0]]
        x2 = self.xv[grad[:, 1]]
        y1 = self.yv[grad[:, 0]]
        y2 = self.yv[grad[:, 1]]

        dx = x1 - x2
        dy = y1 - y2

        self.dg = np.sqrt(dx * dx + dy * dy)

    def count_cells(self):
        """
        Count the total number of 3-D cells
        """
        return np.sum(self.Nk + 1)

    def calc_tangent(self):
        """
        Calculate the tangential vector for the edges of each cell
        """
        if '_tx' not in self.__dict__:
            dx = np.zeros(self.cells.shape)
            dy = np.zeros(self.cells.shape)

            dx[:, 0:-1] = self.xp[self.cells[:, 1::]] - self.xp[self.cells[:, 0:-1]]
            dy[:, 0:-1] = self.yp[self.cells[:, 1::]] - self.yp[self.cells[:, 0:-1]]

            for ii in range(self.Nc):
                dx[ii, self.nfaces[ii] - 1] = self.xp[self.cells[ii, 0]] - self.xp[self.cells[ii, self.nfaces[ii] - 1]]
                dy[ii, self.nfaces[ii] - 1] = self.yp[self.cells[ii, 0]] - self.yp[self.cells[ii, self.nfaces[ii] - 1]]

            mag = np.sqrt(dx * dx + dy * dy)

            self._tx = dx / mag
            self._ty = dy / mag
            self._mag = mag

        return self._tx, self._ty, self._mag

    def calc_edgecoord(self):
        """
        Manually calculate the coordinates of the edge points
        """
        self.xe = np.mean(self.xp[self.edges], axis=1)
        self.ye = np.mean(self.yp[self.edges], axis=1)

    def get_facemark(self):
        """
        Finds the cells next to type-2 or type-3 boundaries
        """
        mask = self.face.mask
        face = self.face.copy()
        face[mask] = 0
        facemark = self.mark[face]
        facemark[mask] = 0
        return np.min(np.max(facemark, axis=-1), 3)

    def setDepth(self, vertspace):
        """
        Calculates and sets the depth variables based on the vertspace vector
        """
        self.dz = vertspace
        self.Nkmax = np.size(self.dz)

        # Calculate the mid-point depth
        if not self.Nkmax == 1:
            # z_bot = np.cumsum(self.dz)
            # z_top = np.hstack((0.0,z_bot[:-1]))
            # self.z_r = 0.5*(z_bot+z_top)
            self.z_r = calc_z(self.dz)
        else:
            self.z_r = np.array([0.0])

    def calcVertSpace(self, Nkmax, r, depthmax):
        """
        Calculates the vertical spacing based on an exponential stretching function
        """

        vertspace = np.zeros((Nkmax,))

        if r < 1.0 or r > 1.1:
            print('r must be between 1.0 and 1.1')

        if Nkmax == 0:
            vertspace[0] = depthmax
        else:
            if r == 1.0:
                vertspace[0] = depthmax / Nkmax
            else:
                vertspace[0] = depthmax * (r - 1.0) / (r ** float(Nkmax) - 1.0)
            for k in range(1, Nkmax):
                vertspace[k] = r * vertspace[k - 1]

        return vertspace

    def pnt2cells(self, pnt_i):
        """
        Returns the cell indices for a point, pnt_i

        (Stolen from Rusty's TriGrid class)
        """
        if '_pnt2cells' not in self.__dict__:
            # build hash table for point->cell lookup
            self._pnt2cells = {}
            for i in range(self.Nc):
                for j in range(3):
                    if self.cells[i, j] not in self._pnt2cells:
                        # self._pnt2cells[self.cells[i,j]] = set()
                        self._pnt2cells[self.cells[i, j]] = []
                    # self._pnt2cells[self.cells[i,j]].add(i)
                    self._pnt2cells[self.cells[i, j]].append(i)
        return self._pnt2cells[pnt_i]

    def cell2node(self, cell_scalar):
        """
        Map a cell-based scalar onto a node

        This is calculated via a mean of the cells connected to a node(point)
        """
        # Simple mean
        # node_scalar = [np.mean(cell_scalar[self.pnt2cells(ii)]) for ii in range(self.Np)]

        # Area weighted interpolation
        node_scalar = [np.sum(cell_scalar[self.pnt2cells(ii)] * self.Ac[self.pnt2cells(ii)]) \
                       / np.sum(self.Ac[self.pnt2cells(ii)]) for ii in range(self.Np)]
        return np.array(node_scalar)

    def __del__(self):
        if 'nc' in self.__dict__:
            self.nc.close()

    def __openNc(self):
        # nc = Dataset(self.ncfile, 'r', format='NETCDF4')
        if self.VERBOSE:
            print('Loading: %s' % self.ncfile)
        try:
            self.nc = MFDataset(self.ncfile, aggdim='time')
        except:
            if type(self.ncfile) == list:
                self.ncfile = self.ncfile[0]
            self.nc = Dataset(self.ncfile, 'r')

    def __getitem__(self, y):
        x = self.__dict__.__getitem__(y)
        return x


####################################################################
#
# General functions to be used by all classes
#
####################################################################
def closePoly(x, y):
    """
    Returns an Nx2 closed polygon for vectors x and y

    This output is required for plotting by unsurf.
    """

    nx = len(x)
    ny = len(y)
    if nx != ny:
        print("Error: The lengths of vector x and y must be equal")
        return

    x = np.reshape(x, (nx, 1))
    y = np.reshape(y, (ny, 1))

    x = np.vstack((x, x[0]))
    y = np.vstack((y, y[0]))

    return np.hstack((x, y))


def readTXT(fname, sep=None):
    """
    Reads a txt file into an array of floats
    """

    fp = open(fname, 'rt')
    data = np.array([list(map(float, line.split(sep))) for line in fp])
    fp.close()

    return data
