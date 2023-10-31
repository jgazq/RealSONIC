# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2020-06-07 14:42:18
# @Last Modified by:   Theo Lemaire
# @Last Modified time: 2023-03-30 14:13:26

import numpy as np
import pandas as pd
from neuron import h

from .pyhoc import SquareMatrix, DiagonalMatrix, pointerMatrix, PointerVector
from ..utils import seriesGeq
from ..constants import *


class ConductanceMatrix(SquareMatrix):
    ''' Interface to an axial conductance matrix. '''

    def __new__(cls, Gvec, **_):
        ''' Instanciation. '''
        return super(ConductanceMatrix, cls).__new__(cls, Gvec.size)

    def __init__(self, Gvec, links=None):
        ''' Initialization.

            :param Gvec: vector of reference conductances for each element (S)
            :param links: list of paired indexes inicating links across nodes.
        '''
        self.Gvec = Gvec
        if links is not None:
            self.setLinks(links)

    def emptyClone(self):
        ''' Return empty matrix of identical shape. '''
        return ConductanceMatrix(self.Gvec)

    @property
    def Gvec(self):
        return self._Gvec

    @Gvec.setter
    def Gvec(self, value):
        assert value.size == self.nRow, 'conductance vector does not match number of rows'
        self._Gvec = value

    def Gij(self, i, j):
        ''' Half conductance in series. '''
        return 2 * seriesGeq(self.Gvec[i], self.Gvec[j])

    def addLink(self, i, j, w):
        ''' Add a bi-directional link between two nodes with a specific weight.

            :param i: first node index
            :param j: second node index
            :param w: link weight
        '''
        self.addVal(i, i, w)
        self.addVal(i, j, -w)
        self.addVal(j, j, w)
        self.addVal(j, i, -w)

    def link(self, i, j):
        ''' Add a link between two nodes. '''
        self.addLink(i, j, self.Gij(i, j))

    def unlink(self, i, j):
        ''' Remove a link between two nodes.

            :param i: first node index
            :param j: second node index
        '''
        self.addLink(i, j, -self.Gij(i, j))

    def setLinks(self, links):
        ''' Add cross-nodes links to the matrix.

            :param links: list of paired indexes inicating links across nodes.
        '''
        for i, j in links:
            self.link(i, j)
        self.checkNullRows()

    def checkNullRows(self):
        ''' Check that all rows sum up to zero (or close). '''
        for i in range(self.nRow):
            rsum = self.getrow(i).sum()
            assert np.isclose(rsum, .0, atol=1e-12), f'non-zero sum on line {i}: {rsum}'


class NormalizedConductanceMatrix(ConductanceMatrix):
    ''' Interface to an normalized axial conductance matrix. '''

    def __new__(cls, Gvec, *args, **kwargs):
        ''' Instanciation. '''
        return super(NormalizedConductanceMatrix, cls).__new__(cls, Gvec)

    def __init__(self, Gvec, rnorm=None, cnorm=None, **kwargs):
        ''' Initialization.

            :param rnorm: vector specifying the normalization factor for each row
            :param cnorm: vector specifying the normalization factor for each column
        '''
        self.rnorm = np.ones(Gvec.size) if rnorm is None else rnorm
        self.cnorm = np.ones(Gvec.size) if cnorm is None else cnorm
        super().__init__(Gvec, **kwargs)

    def emptyClone(self):
        ''' Return empty matrix of identical shape. '''
        return NormalizedConductanceMatrix(self.Gvec, rnorm=self.rnorm, cnorm=self.cnorm)

    @property
    def rnorm(self):
        return self._rnorm

    @rnorm.setter
    def rnorm(self, value):
        self._rnorm = value.copy()

    @property
    def cnorm(self):
        return self._cnorm

    @cnorm.setter
    def cnorm(self, value):
        self._cnorm = value.copy()

    def xnorm(self, i, j):
        ''' Normalizing factor for i-th row, j-th column. '''
        return self.rnorm[i] * self.cnorm[j]

    def setVal(self, i, j, x):
        ''' Set matrix element divided by the corresponding normalizing factor. '''
        super().setVal(i, j, x / self.xnorm(i, j))

    def getVal(self, i, j):
        ''' Get matrix element multiplied by the corresponding normalizing factor. '''
        return super().getVal(i, j) * self.xnorm(i, j)

    def setRNorm(self, value):
        ''' Set a new row-normalization vector. '''
        self.mulRows(self.rnorm / value)
        self.rnorm = value
        self.onFullUpdate()

    def setCNorm(self, value):
        ''' Set a new column-normalization vector. '''
        self.mulCols(self.cnorm / value)
        self.cnorm = value
        self.onFullUpdate()

    def scaled(self, rnorm=None, cnorm=None):
        ''' Return a "scaled" version of the matrix where each element is multiplied
            by the appropriate row and column normalizer.
        '''
        if rnorm is None:
            rnorm = self.rnorm
        if cnorm is None:
            cnorm = self.cnorm
        mout = self.fullClone()
        mout.mulRows(rnorm)
        mout.mulCols(cnorm)
        return mout


class HybridNetwork:
    ''' Interface used to build a hybrid voltage network amongst a list of sections.

        Consider a neuron model consisting of a list of sections connected in series.

        We define the following terms:

        - vi: intracelular voltage
        - vm: transmembrane voltage
        - vx: extracellular voltage
        - ex: imposed voltage outside of the surrounding extracellular membrane
        - is: stimulating current
        - cm: membrane capacitance
        - i(vm): transmembrane ionic current
        - ga: intracellular axial conductance between nodes
        - cx: capacitance of surrounding extracellular membrane (e.g. myelin)
        - gx: transverse conductance of surrounding extracellular membrane (e.g. myelin)
        - gp: extracellular axial conductance between nodes (e.g. periaxonal space)
        - j: index indicating the connection to a neighboring node

        Governing equations for internal and external nodes are, respectively:

        (1) cm * dvm/dt + i(vm) = is + ga_j * (vi_j - vi)
        (2) cx * dvx/dt + gx * (vx - ex) = cm * dvm/dt + i(vm) + gp_j * (vx_j - vx)

        Putting all voltage dependencies on the left-hand sides, and developing, we find:

        (1) cm * dvm/dt + ga_j * (vi - vi_j) = is - i(vm)
        (2) cx * dvx/dt - cm * dvm/dt + gx * vx + gp_j * (vx - vx_j) = i(vm) + gx * ex

        Re-expressing vi as (vm + vx), we find the matrix equation rows for the two nodes:

        (1) cm * dvm/dt + ga_j * (vm - vm_j) + ga_j * (vx - vx_j) = is - i(vm)
        (2) cx * dvx/dt - cm * dvm/dt + gx * vx + gp_j * (vx - vx_j) = i(vm) + gx * ex

        Special attention must be brought on the fact that we use NEURONS's v variable as an
        alias for the section's membrane charge density Qm. Therefore, the system must be
        adapted in consequence. To this effect, we introduce a change of variable:

        Qm = cm * vm; dQm/dt = cm * dvm/dt  (cm considered time-invariant)

        Recasting the system according to this change of variable thus yields:

        (1) dQm/dt + ga_j * (Qm/cm - Qm_j/cm_j) + ga_j * (vx - vx_j) = is - i(vm)
        (2) cx * dvx/dt - dQm/dt + gx * vx + gp_j * (vx - vx_j) = i(vm) + gx * ex

        Finally, we replace dQm/dt + i(vm) in (2) by equivalents intracellular axial and
        stimulation currents, in order to remove the need to access net membrane current.
        After re-arranging all linear voltage terms on the left-hand side, we have:

        (1) dQm/dt + ga_j * (Qm/cm - Qm_j/cm_j) + ga_j * (vx - vx_j) = is - i(vm)
        (2) cx * dvx/dt + ga_j * (Qm/cm - Qm_j/cm_j) + (ga_j + gp_j) * (vx - vx_j) + gx * vx
            = is + gx * ex

        Developing to isolate Qm elements, we have:

        (1) dQm/dt + ga_j/cm * Qm - ga_j/cm_j * Qm_j + ga_j * (vx - vx_j) = is - i(vm)
        (2) cx * dvx/dt + ga_j/cm * Qm - ga_j/cm_j * Qm_j + (ga_j + gp_j) * (vx - vx_j) + gx * vx
            = is + gx * ex

        Among these equation rows, 2 terms are automatically handled by NEURON, namely:
        - LHS-1: dQm/dt
        - RHS-1: is - i(vm)

        Hence, 8 terms remain to be handled as part of the additional network:
        - LHS-1: ga_j/cm * Qm - ga_j/cm_j * Qm_j
        - LHS-1: ga_j * vx - ga_j * vx_j
        - LHS-2: cx * dvx/dt
        - LHS-2: ga_j/cm * Qm - ga_j/cm_j * Qm_j
        - LHS-2: ga_j * vx - ga_j * vx_j
        - LHS-2: gp_j * (vx - vx_j)
        - LHS-2: gx * vx
        - RHS-2: is + gx * ex

        Note that axial conductance terms (ga and gp) must be normalized by the
        appropriate section membrane area in order to obtain a consistent system.

        How do we do this?

        The above equations describe a partial differential system of the form:

            C * dy/dt + G * y = I

        with C a capacitance matrix, G a conductance matrix, y a hybrid vector
        of transmembrane charge density and external voltage, and I a current vector.

        Thankfully, NEURON provides a so-called "LinearMechanism" class that allows
        to define such a system, where the first n elements of the y vector
        can be mapped to the "v" variable of a specific list of sections. This allows
        us to add a linear network atop of NEURON's native network, in which we define
        the additional terms.

        Concretely, a model of n sections connected in series can be represented by a "CGI"
        system of size 2*n, where the first n items correspond to membrane charge density
        nodes, and the following n items correspond to external voltage nodes.

                -------------------------------
        y =     |       Qm      |     vx      |
                -------------------------------

        The corresponding linear mechanism added atop of NEURON's native linear network
        would then consist of the following capacitance and conductance matrices C and G,
        and current vector I:

                -------------------------------
                |              |              |
                |       0      |       0      |
                |              |              |
        C  =    |-----------------------------|
                |              |              |
                |        0     |      Cx      |
                |              |              |
                -------------------------------

                -------------------------------
                |              |              |
                |     Ga/cm    |      Ga      |
                |              |              |
        G =     |-----------------------------|
                |              |              |
                |     Ga/cm    | Ga + Gp + Gx |
                |              |              |
                -------------------------------

                -------------------------------
        I =     |       0      |  gx*ex + is  |
                -------------------------------

        Where internal terms are defined as:
        - Ga: n-by-n matrix of intracellular axial conductance
        - Ga/cm: Ga matrix column-normalized by capacitance at each time step
        - Gp: n-by-n matrix of extracellular axial conductance
        - Gx: n-by-n sparse matrix of transverse extracellular conductance
        - Cx: n-by-n sparse matrix of transverse extracellular capacitance

        Note also that the gx conductance is connected in series with another transverse
        conductance of value 1e9, to mimick NEURON's default 2nd extracellular layer.
    '''

    def __init__(self, seclist, connections, has_ext_layer, is_dynamic_cm=False, verbose=False,
                 use_explicit_iax=False):
        ''' Initialization.

            :param seclist: list of sections.
            :param connections: list of index pairs (tuples) indicating connections across sections
            :param has_ext_layer: boolean indicating whether to implement an extracellular layer
            :param is_dynamic_cm: boolean indicating whether membrane capacitance is time-varying
            :param verbose: boolean stating whether to log details about the network
            :param use_explicit_iax: boolean stating whether to use explicit axial currents (as
                opposed to embedding them in the global conductance matrix)
        '''
        # Assign attributes
        self.seclist = seclist
        self.connections = connections
        self.has_ext_layer = has_ext_layer
        self.is_dynamic_cm = is_dynamic_cm
        self.use_explicit_iax = use_explicit_iax
        self.setGlobalComponents()
        self.setBaseLayer()
        if self.has_ext_layer:
            if self.use_explicit_iax:
                raise ValueError('explicit axial currents not defined for extracellular layers')
            self.setExtracellularLayer()
        if verbose:
            self.log(details=True)

    def __repr__(self):
        cm = {False: 'static', True: 'dynamic'}[self.is_dynamic_cm] + ' cm'
        dims = f'{self.nsec} sections, {self.nlayers} layers'
        s = f'{self.__class__.__name__}({dims}, {cm}'
        if self.use_explicit_iax:
            s = f'{s}, explicit iax'
        return f'{s})'

    @property
    def seclist(self):
        return self._seclist

    @seclist.setter
    def seclist(self, value):
        self._seclist = value
    
    @property
    def secnames(self):
        return [str(x).split('.')[-1] for x in self.seclist]

    @property
    def connections(self):
        return self._connections

    @connections.setter
    def connections(self, value):
        self._connections = value

    @property
    def has_ext_layer(self):
        return self._has_ext_layer

    @has_ext_layer.setter
    def has_ext_layer(self, value):
        self._has_ext_layer = value

    @property
    def nsec(self):
        ''' Number of sections in the network. '''
        return len(self.seclist)

    @property
    def nlayers(self):
        ''' Number of layers in the network. '''
        return 1 if not self.has_ext_layer else 2

    @property
    def size(self):
        ''' Overall size if the network (i.e. number of nodes). '''
        return self.nsec * self.nlayers

    def getVector(self, k):
        ''' Get a vector of values of a given parameter for each section in the list.

            :param k: parameter name
            :return: 1D numpy array with parameter values
        '''
        return np.array([getattr(sec, k) for sec in self.seclist])

    def setGlobalComponents(self):
        ''' Set the network's global components used by NEURON's linear Mechanism. '''
        self.C = SquareMatrix(self.size)  # capacitance matrix (mF/cm2)
        self.G = SquareMatrix(self.size)  # conductance matrix (S/cm2)
        self.y = h.Vector(self.size)      # charge density / extracellular voltage vector (mV)
        self.I = h.Vector(self.size)      # current vector (mA/cm2)

    def setBaseLayer(self):
        ''' Set components used to define the network's base layer. '''
        # Get required vectors
        self.Am = self.getVector('Am')   # membrane area (cm2)
        self.cm = self.getVector('Cm0')  # resting membrane capacitance (uF/cm2)
        self.ga = self.getVector('Ga')   # intracellular axial conductance (S)
        if self.use_explicit_iax:
            # Define qm and iax vectors
            self.qm = h.Vector(self.nsec)
            self.iax = PointerVector(self.nsec)
            self.iax.addRef(self.I, 0)

        # Define Gacm matrix and point it towards G top-left (if no explicit iax)
        self.Gacm = pointerMatrix(NormalizedConductanceMatrix)(
            self.ga, links=self.connections, rnorm=self.Am, cnorm=self.cm)
        if not self.use_explicit_iax:
            self.Gacm.addRef(self.G, 0, 0)

    def setExtracellularLayer(self):
        ''' Set components used to define the network's extracellular layer. '''
        # Get required vectors
        self.cx = self.getVector('cx')  # uF/cm2
        self.gx = self.getVector('gx')  # S/cm2
        self.gp = self.getVector('Gp')  # S

        # Define additional matrices and point them towards G
        self.Cx = pointerMatrix(DiagonalMatrix)(self.cx * UF_CM2_TO_MF_CM2)
        self.Ga = pointerMatrix(NormalizedConductanceMatrix)(
            self.ga, links=self.connections, rnorm=self.Am)
        self.Gp = pointerMatrix(NormalizedConductanceMatrix)(
            self.gp, links=self.connections, rnorm=self.Am)
        self.Gx = pointerMatrix(DiagonalMatrix)(self.gx)

        # Add references to global matrices
        self.Cx.addRef(self.C, self.nsec, self.nsec)  # Bottom-right: Cx * dvx/dt
        self.Gacm.addRef(self.G, self.nsec, 0)        # Bottom-left: Ga/cm * Qm
        self.Ga.addRef(self.G, 0, self.nsec)          # Top-right: Ga * vx
        self.Gp.addRef(self.G, self.nsec, self.nsec)  # Bottom-right: Gp * vx
        self.Ga.addRef(self.G, self.nsec, self.nsec)  # Bottom-right: Ga * vx
        self.Gx.addRef(self.G, self.nsec, self.nsec)  # Bottom-right: Gx * vx

        # Define pointing vectors toward extracellular voltage and currents
        self.istim = PointerVector(self.nsec)
        self.istim.addRef(self.I, self.nsec)
        self.iex = PointerVector(self.nsec)
        self.iex.addRef(self.I, self.nsec)
        self.vx = PointerVector(self.nsec)
        self.vx.addRef(self.y, self.nsec)

    def index(self, sec):
        return self.seclist.index(sec)

    def getVextRef(self, sec):
        ''' Get reference to a section's extracellular voltage variable. '''
        if not self.has_ext_layer:
            raise ValueError('Network does not have an extracellular layer')
        return self.y._ref_x[self.index(sec) + self.nsec]

    def startLM(self):
        ''' Feed network into a LinearMechanism object. '''
        # Define section list (using explicit append statements to ensure retro-compatibility)
        self.sl = h.SectionList()
        for sec in self.seclist:
            self.sl.append(sec=sec)
        # Define vector of relative positions in the sections (always at mid-point)
        self.relx = h.Vector([0.5] * self.nsec)
        # Set initial conditions vector
        self.y0 = h.Vector(self.size)
        # Define linear mechanism arguments
        lm_args = [self.C, self.G, self.y, self.y0, self.I, self.sl, self.relx]
        # Add update callback for dynamic cm
        if self.is_dynamic_cm:
            lm_args = [self.updateCmTerms] + lm_args
        # Add update callback for explicit iax
        elif self.use_explicit_iax:
            lm_args = [self.updateIax] + lm_args
        # Create LinearMechanism object
        self.lm = h.LinearMechanism(*lm_args)

    def clear(self):
        ''' Delete the network's LinearMechanism object. '''
        # Remove base layer components
        self.Am = None
        self.cm = None
        self.ga = None
        self.qm = None
        self.iax = None
        self.Gacm = None
        self.istim = None
        # Remove extracellular layer components
        self.cx = None
        self.gx = None
        self.gp = None
        self.Ga = None
        self.Gp = None
        self.Gx = None
        self.vx = None
        self.iex = None
        # Remove global components
        self.C = None
        self.G = None
        self.I = None
        self.y = None
        self.y0 = None
        self.sl = None
        self.relx = None
        self.lm = None
        self.seclist = None

    def updateCmTerms(self):
        ''' Update capacitance-dependent network components. '''
        # Update membrane capacitance vector
        self.cm = np.array([sec.getCm(x=0.5) for sec in self.seclist])

        # Modify Gacm matrix accordingly
        self.Gacm.setCNorm(self.cm)
        # Update iax vector if needed
        if self.use_explicit_iax:
            self.updateIax()

    def computeIax(self, Vi):
        ''' Compute axial current density vector from a internal potential vector.

            :param Vi: internal potential vector (mV)
            :return: axial current density vector (mA/cm2)
        '''
        Ga = self.Gacm.scaled(rnorm=np.ones(self.nsec))  # S/cm2
        return np.array(Ga.mulv(h.Vector(Vi)).to_python())

    def updateIax(self):
        ''' Update axial currents as an explicit current vector contribution. '''
        # Get membrane charge density vector
        self.qm.copy(self.y, 0, self.nsec - 1)
        # Compute axial current vector
        iax = -self.Gacm.mulv(self.qm)
        # Assign to I vector
        for i in range(self.nsec):
            self.iax.setVal(i, iax.get(i))

    def setIstim(self, sec, I):
        ''' Set the stimulation current of a given section to a new value. '''
        self.istim.setVal(self.index(sec), I / MA_TO_NA / sec.Am)

    def setEx(self, sec, old_ex, new_ex):
        ''' Set the imposed extracellular potential of a given section to a new value
            (propagating the discontinuity into the vx vector).
        '''
        i = self.index(sec)
        self.vx.addVal(i, new_ex - old_ex)
        self.iex.setVal(i, self.gx[i] * new_ex)

    def mat_to_table(self, k):
        ''' Convert matrix to table '''
        if not hasattr(self, k):
            raise ValueError(f'"{k}" matrix not found in {self}')
        mat = getattr(self, k)
        if not isinstance(mat, SquareMatrix):
            raise ValueError(f'"{k}" is not a square matrix')
        arr = mat.to_array()
        if arr.shape[0] == self.nsec:
            idxs = self.secnames
        elif arr.shape[0] == 2 * self.nsec:
            idxs = self.secnames + self.secnames
        else:
            raise ValueError(
                f'cannot convert {arr.shape} matrix to table: dimensions do not match network size') 
        return pd.DataFrame(data=arr, index=idxs, columns=idxs).style.format('{:.3g}')

    def logMat(self, k, unit):
        print(f'{k} ({unit}):\n{self.mat_to_table(k)}')
        # getattr(self, k).printf('%-8g' if self.size <= 22 else '%-5g')

    def logCmat(self, k):
        self.logMat(k, 'mF/cm2')

    def logGmat(self, k):
        self.logMat(k, 'S/cm2')

    def log(self, details=False):
        ''' Print network components. '''
        if details and self.has_ext_layer:
            self.logCmat('Cx')
        self.logCmat('C')
        if details:
            self.logGmat('Gacm')
            if self.has_ext_layer:
                self.logGmat('Ga')
                self.logGmat('Gp')
                self.logGmat('Gx')
        self.logGmat('G')
