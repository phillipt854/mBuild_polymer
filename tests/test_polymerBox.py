import mbuild as mb
import numpy as np
import pytest
from mbuild_polymer import *

def test_boxChains():
    testBox = CLP_box([5],[10], dim=[2, 1, 1])
    chains = [c.chain_length_pvpy for c in testBox.children]
    testBox.save('polymerbox.gsd')
    assert chains == ['POG']
