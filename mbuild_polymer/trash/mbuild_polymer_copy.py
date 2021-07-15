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

r0BB_pvpy = 1.0 
r0BB_pvph = 1.134 
r0BBHB_pvpy = 0.37
r0BBHB_pvph = 0.42

# Define backbone (BB) classes
class BB(mb.Compound):
    def __init__(self):
        super(BB, self).__init__(pos=[0,0,0],name = 'bb')#Initizlize an instance of abstract class

        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='BB')
        self.add(bead)

        port = mb.Port(anchor=list(self.particles_by_name('BB'))[0], orientation=[0, 1, 0], separation=r0BB/2)
        self.add(port, 'up')
        port = mb.Port(anchor=list(self.particles_by_name('BB'))[0], orientation=[1, 0, 0], separation=r0BBHB/2)
        self.add(port,'toHB')
        port = mb.Port(anchor=list(self.particles_by_name('BB'))[0], orientation=[0, -1, 0], separation=r0BB/2)
        self.add(port, 'down')

class BB_pvpy(BB):
    def __init__(self):
        super(BB_pvpy,self).__init__(r0BB=r0BB_pvpy)
        for par in self.particles():
            par.name = '_BB_pvpy'
            #print(par.name)

        # Name the entire compound (particle+its ports) 
        self.name = '_bb_pvpy'

class BB_pvpy():
    def __init__(self):
        super(BB_pvpy, self).__init__(pos=[0,0,0],name = 'bb_pvpy')#Initizlize an instance of abstract class
        
        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='BB_pvpy')
        self.add(bead)
        
        port = mb.Port(anchor=list(self.particles_by_name('BB_pvpy'))[0], orientation=[0, 1, 0], separation=r0BB_pvpy/2)
        self.add(port, 'up')
        port = mb.Port(anchor=list(self.particles_by_name('BB_pvpy'))[0], orientation=[1, 0, 0], separation=r0BBHB_pvpy/2)
        self.add(port,'toHB')
        port = mb.Port(anchor=list(self.particles_by_name('BB_pvpy'))[0], orientation=[0, -1, 0], separation=r0BB_pvpy/2)
        self.add(port, 'down')
        
class BB_pvph(mb.Compound):
    def __init__(self):
        super(BB_pvph, self).__init__(pos=[0,0,0],name = 'bb_pvph')#Initizlize an instance of abstract class

        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='BB_pvph')
        self.add(bead)

        port = mb.Port(anchor=list(self.particles_by_name('BB_pvph'))[0], orientation=[0, 1, 0], separation=r0BB_pvph/2)
        self.add(port, 'up')
        port = mb.Port(anchor=list(self.particles_by_name('BB_pvph'))[0], orientation=[1, 0, 0], separation=r0BBHB_pvph/2)
        self.add(port,'toHB')
        port = mb.Port(anchor=list(self.particles_by_name('BB_pvph'))[0], orientation=[0, -1, 0], separation=r0BB_pvph/2)
        self.add(port, 'down')


# Define hydrogen bonding (HB) classes    
class HB_pvpy(mb.Compound):
    def __init__(self):
        super(HB_pvpy, self).__init__(pos=[0.0, 0.0, 0.0], name='hb_pvpy')
        
        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='HB_pvpy')
        self.add(bead)
        
        port = mb.Port(anchor=list(self.particles_by_name('HB_pvpy'))[0], orientation=[-1, 0, 0], separation=r0BBHB_pvpy/2)
        self.add(port, 'toBB')      

class HB_pvph(mb.Compound):
    def __init__(self):
        super(HB_pvph, self).__init__(pos=[0.0, 0.0, 0.0], name='hb_pvph')

        bead = mb.Particle(pos=[0.0, 0.0, 0.0], name='HB_pvph')
        self.add(bead)

        port = mb.Port(anchor=list(self.particles_by_name('HB_pvph'))[0], orientation=[-1, 0, 0], separation=r0BBHB_pvph/2)
        self.add(port, 'toBB')


# Define monomer (M) classes - combine backbone and hydrogen bonding beads
class M_pvpy(mb.Compound):
    def __init__(self):
        super(M_pvpy,self).__init__()
        bb = BB_pvpy()
        hb = HB_pvpy()
        self.add((bb,hb))
        #Move
        mb.force_overlap(move_this=hb, from_positions=hb['toBB'],to_positions=bb['toHB'])

class M_pvph(mb.Compound):
    def __init__(self):
        super(M_pvph,self).__init__()
        bb = BB_pvph()
        hb = HB_pvph()
        self.add((bb,hb))
        #Move
        mb.force_overlap(move_this=hb, from_positions=hb['toBB'],to_positions=bb['toHB'])

class polymer_pvpy(mb.Compound):
    chain_len = None
    def __init__(self,length=None):
        self.chain_len = length
        super(polymer_pvpy,self).__init__()
        if length != 0:
            last_M = M_pvpy()
            self.add(last_M)
        
            for i in range(length-1):
                new_M = M_pvpy()
                self.add(new_M)
                mb.force_overlap(move_this = new_M, from_positions=(new_M.all_ports())[0],to_positions = (last_M.all_ports())[-1])
                last_M = new_M

    def polymer_visualize_py3dmol(self, show_ports=False):
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


        for particle in self.particles():
            if not particle.name:
                particle.name = 'UNK'
                
        view = py3Dmol.view()
        rad = {'BB_pvpy':.5,'HB_pvpy':0.22}
        col = {'BB_pvpy':'#ffff14','HB_pvpy':'#000000'}
        
        for p in self.particles(include_ports=False):
            view.addSphere({
                'center': {'x':p.pos[0], 'y':p.pos[1], 'z':p.pos[2]},
                'radius' :rad[p.name],
                'color': col[p.name],
                'alpha': 0.9})
        view.zoomTo()
        view.show()

        return view

class polymer_pvph(mb.Compound):
    chain_len = None
    def __init__(self,length=None):
        self.chain_len = length
        super(polymer_pvph,self).__init__()
        if length != 0:
            last_M = M_pvph()
            self.add(last_M)

            for i in range(length-1):
                new_M = M_pvph()
                self.add(new_M)
                mb.force_overlap(move_this = new_M, from_positions=(new_M.all_ports())[0],to_positions = (last_M.all_ports())[-1])
                last_M = new_M

    def polymer_visualize_py3dmol(self, show_ports=False):
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


        for particle in self.particles():
            if not particle.name:
                particle.name = 'UNK'

        view = py3Dmol.view()
        rad = {'BB_pvph':.567,'HB_pvph':0.24}
        col = {'BB_pvph':'#15b01a','HB_pvph':'#e50000'}

        for p in self.particles(include_ports=False):
            view.addSphere({
                'center': {'x':p.pos[0], 'y':p.pos[1], 'z':p.pos[2]},
                'radius' :rad[p.name],
                'color': col[p.name],
                'alpha': 0.9})
        view.zoomTo()
        view.show()

        return view

class polymer_box(mb.Compound):
    def __init__(self,chain_len_pvpy = [],chain_len_pvph = [],dim = [0,0,0]):
        super(polymer_box,self).__init__()
        self.chain_len_pvpy = chain_len_pvpy
        self.chain_len_pvph = chain_len_pvph
        if (len(chain_len_pvpy) + len(chain_len_pvph)) > 0 and dim[0]*dim[1]*dim[2] != (len(chain_len_pvpy) + len(chain_len_pvph)):
            dim = [len(chain_len_pvpy) + len(chain_len_pvph),1,1]
        chain_num_pvpy = 0
        chain_num_pvph = 0
        total_chain_num = 0
        for i in range(dim[0]):
            for j in range(dim[1]):
                for k in range(dim[2]):
                    if (len(chain_len_pvpy) > 0) | (len(chain_len_pvph) > 0):
                        if (total_chain_num < len(chain_len_pvpy)):
                            new_polymer = polymer_pvpy(chain_len_pvpy[chain_num_pvpy])
                            new_polymer.translate_to([i*3.,j*3.,k*3.])
                            self.add(new_polymer)
                            chain_num_pvpy += 1          
                            total_chain_num += 1
                        elif ((total_chain_num - len(chain_len_pvpy)) < len(chain_len_pvph)):
                            new_polymer = polymer_pvph(chain_len_pvph[chain_num_pvph])
                            new_polymer.translate_to([i*3.,j*3.,k*3.])
                            self.add(new_polymer)
                            chain_num_pvph += 1
                            total_chain_num += 1

         
    def visualize(self):
        py3Dmol = import_('py3Dmol')
        view = py3Dmol.view()
        if (len(self.chain_len_pvpy) > 0) and (len(self.chain_len_pvph) > 0):
            rad = {'BB_pvpy':.5,'BB_pvph':.567,'HB_pvpy':0.22,'HB_pvph':0.24}
            col = {'BB_pvph':'#15b01a','BB_pvpy':'#ffff14','HB_pvph':'#e50000','HB_pvpy':'#000000'}
        elif (len(self.chain_len_pvpy) == 0):
            rad = {'BB_pvph':.567,'HB_pvph':0.24}
            col = {'BB_pvph':'#15b01a','HB_pvph':'#e50000'}
        elif (len(self.chain_len_pvph) == 0):
            rad = {'BB_pvpy':.5,'HB_pvpy':0.22}
            col = {'BB_pvpy':'#ffff14','HB_pvpy':'#000000'}
        
        remove_digits = lambda x: ''.join(i for i in x if not i.isdigit()
                                              or i == '_')

        #modified_color_scheme = {}
        
        for chain in self.children:
            for particle in chain.particles():
                #particle.name = remove_digits(particle.name).upper()
                if not particle.name:
                    particle.name = 'UNK'
                
        
        
                for p in chain.particles(include_ports=False):
                    view.addSphere({
                        'center': {'x':p.pos[0], 'y':p.pos[1], 'z':p.pos[2]},
                        'radius' :rad[p.name],
                        'color': col[p.name],
                        'alpha': 0.9})
        view.zoomTo()
        view.show()
        return view

		








