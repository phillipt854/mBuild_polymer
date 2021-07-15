import mbuild as mb
import numpy as np
import pytest
from mbuild_polymer import *

def test_chain_length_pvpy():
    chain = polymer_pvpy(20)
    assert len(chain.children) == 20 and chain._n_particles() == 40


def test_sc_particles():
    chain = polymer_pvph(1)
    name_list = []
    for aa in chain.children:
        for p in aa.children:
            name_list.append(p.name)
    assert name_list == ['bb_pvph', 'hb_pvph']


def test_HBBB_bond():
    # Check there is only 1 bond that has the only 2 particles in the system 
    # connected
    chain = polymer_pvpy(1)
    particle_set = set([p for p in chain.particles()])
    bonds = [b for b in chain.bonds()]
    firstBond = set(bonds[0])
    assert len(bonds) == 1 and firstBond == particle_set


