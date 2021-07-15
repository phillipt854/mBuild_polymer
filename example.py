import mbuild as mb
import py3Dmol
from mbuild_polymer.mbuild_polymer import polymer_box

test_box = polymer_box([10],[1,1,2])
#test_box.visualize()
test_box.write_lammps('test.lammps')
test_box.create_lammps_input_script(sim_name='test.lammps')

