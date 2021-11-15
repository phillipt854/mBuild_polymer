#!/usr/bin/env python
# coding: utf-8


from __future__ import print_function
import mbuild as mb
import builtins as __builtin__
import numpy as np
from foyer import Forcefield
#All the old import from compound
import collections
from collections import OrderedDict, defaultdict
from copy import deepcopy
import itertools
import os
import sys
import tempfile
from warnings import warn
##
from mbuild.utils.io import run_from_ipython, import_

# Define backbone (BB) classes
class BB(mb.Compound):
    def __init__(self, type_name='X', r0BB = 1, r0BBHB = 0.37, HB_orient = np.array([1,0,0])):
        super(BB, self).__init__(pos=[0,0,0],name = 'bb_'+type_name)#Initizlize an instance of abstract class
        full_name = 'BB_'+type_name
        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name=full_name)
        self.add(bead)

        port = mb.Port(anchor=list(self.particles_by_name(full_name))[0], orientation=[0, 1, 0], separation=r0BB/2)
        self.add(port, 'up')
        port = mb.Port(anchor=list(self.particles_by_name(full_name))[0],
            orientation=HB_orient, separation=r0BBHB/2)
        self.add(port,'toHB')
        port = mb.Port(anchor=list(self.particles_by_name(full_name))[0], orientation=[0, -1, 0], separation=r0BB/2)
        self.add(port, 'down')



 
# Define hydrogen bonding (HB) classes    
class HB(mb.Compound):
    def __init__(self, type_name='X', r0BBHB = 0.37, HB_orient = np.array([1,0,0])):
        full_name = 'HB_'+type_name
        super(HB, self).__init__(pos=[0.0, 0.0, 0.0], name='hb_'+type_name)

        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name=full_name)
        self.add(bead)

        port = mb.Port(anchor=list(self.particles_by_name(full_name))[0],
            orientation=-HB_orient, separation=r0BBHB/2)
        self.add(port, 'toBB')

# Define monomer (M) classes - combine backbone and hydrogen bonding beads
class M(mb.Compound):
    def __init__(self,type_name='X', r0BB = 1, r0BBHB = 0.37, HB_orient = np.array([1,0,0])):
        HB_orient = np.array(HB_orient)
        super(M,self).__init__()
        bb = BB(type_name, r0BB, r0BBHB, HB_orient)
        hb = HB(type_name, r0BBHB, HB_orient)
        self.add((bb,hb))
        #Move
        mb.force_overlap(move_this=hb, from_positions=hb['toBB'],to_positions=bb['toHB'])

class general_polymer:
    bead_dict = {}
    def __init__(self, beadtype_info):
        if isinstance(beadtype_info, dict):
            self.bead_dict = beadtype_info
        if isinstance(beadtype_info, str):
            from yaml import load
            from yaml import CLoader as Loader
            stream = open(beadtype_info, 'r')
            self.bead_dict = load(stream, Loader=Loader)
    
    def _get_M(self,beadtype_name):
        beadtype = self.bead_dict[beadtype_name]
        return M(beadtype_name, beadtype['r0BB'], beadtype['r0BBHB'], beadtype['HB_orient'])
                 
                 

    def gen_chain(self,sequence=[]):
        chain = mb.Compound()
        if len(sequence) > 0:
            last_M = self._get_M(sequence[0])
            chain.add(last_M)
        
            for i in range(1,len(sequence)):
                new_M = self._get_M(sequence[i])
                chain.add(new_M)
                mb.force_overlap(move_this = new_M, from_positions=(new_M.all_ports())[0],to_positions = (last_M.all_ports())[-1])
                last_M = new_M
        return chain
    def visualize(self, obj, show_ports=False):
        """Visualize the Compound using py3Dmol.
        Allows for visualization of a Compound within a Jupyter Notebook.
        Parameters
        ----------
        show_ports : bool, optional, default=False
            Visualize Ports in addition to Particles
        color_scheme : dict, optional
            Specify coloring for non-elemental particles
            keys are strings of the particle names
            values are strings of the colors
            i.e. {'_CGBEAD': 'blue'}
        Returns
        ------
        view : py3Dmol.view
        """
        py3Dmol = import_('py3Dmol')
        remove_digits = lambda x: ''.join(i for i in x if not i.isdigit()
                                              or i == '_')


        for particle in obj.particles():
            if not particle.name:
                particle.name = 'UNK'
                
        view = py3Dmol.view()
        
        
        for p in obj.particles(include_ports=False):
            if p.name[:2] == 'BB':
                 col = self.bead_dict[p.name[3:]]['BB_color']
                 rad = self.bead_dict[p.name[3:]]['r0BB']
            elif p.name[:2] == 'HB':
                 col = self.bead_dict[p.name[3:]]['HB_color']
                 rad = self.bead_dict[p.name[3:]]['r0HB']
            else:
                 col = 'black'
            view.addSphere({
                'center': {'x':p.pos[0], 'y':p.pos[1], 'z':p.pos[2]},
                'radius': rad,
                'color': col,
                'alpha': 0.9})
        view.zoomTo()
        view.show()

        return view



		








