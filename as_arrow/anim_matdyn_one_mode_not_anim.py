from vpython import *
import sys
import numpy as np
import ph_atom_motion_function as func



class display():
    def __init__(self):
     self.COLORS=[color.red , color.yellow , color.green, color.purple , color.blue 	, color.cyan 	, color.orange 	, color.magenta ,color.orange, color.black	 ]
     self.crystal_lattice=[]
     self.arrows=[]
     self.atomic_balls=[]
     self.scene=()
     self.A=10 #amplitude. The displacement is multiplied by A
    def draw_lattice(self,crystal,crystal_primitive):
        self.crystal_lattice=func.draw_lattice(crystal,crystal_primitive)    
    def set_scene(self,crystal):
        self.scene=func.set_scene(crystal)
    def set_coord_system(self,alat):
        func.set_coord_system(alat)
    def draw_equilibrium_atoms(self,equil_atoms):
        self.atomic_balls=func.draw_equilibrium_atoms(equil_atoms,self.COLORS)   
    def init_arrows(self,equil_atoms,vib):
        self.arrows=func.init_arrows\
                    (equil_atoms,vib,self.A,self.COLORS)     
    def draw_displacement_arrows(self,atoms,vib,no_of_modes):
        func.draw_displacement_arrows\
                    (self.scene,self.arrows,atoms,vib,self.A,no_of_modes)    
    def legend(self,atoms):
        func.legend(atoms[:],self.COLORS)

class inputs():
 def __init__(self):
        self.file_scf_out='scf.out'
        self.file_matdyn_in='matdyn.in'
        self.file_matdyn_modes='matdyn.modes'

class system(inputs):
 def __init__(self):
        inputs.__init__(self)
        self.atoms=[]
        self.crystal=[]
        self.crystal_conv=[]
        self.crystal_primitive=[]
       # self.SYMM_OP=[]
        self.alat=0
 def read_crystal_info(self):
  self.atoms,self.crystal_conv,self.crystal_primitive,self.alat=\
       func.read_crystal_info(self.file_scf_out)
# def set_symm(self):
#  self.SYMM_OP= func.set_sym_bl(self.crystal_primitive)
 def move_atoms_to_cell(self):
  self.atoms=func.move_atoms_to_cell(self.atoms,self.crystal)
 def add_atoms_by_symmetry(self):
  self.atoms=func.add_atoms_by_symmetry(self.atoms,\
             self.crystal,self.crystal_primitive)

class motion(inputs):
 def __init__(self):
       inputs.__init__(self)
       self.FREQ=[]
       self.DISPL=[]
       self.Q=[]    
 def read_freqs_and_displacements(self):
  self.DISPL,self.FREQ,self.Q,self.no_of_modes=\
       func.read_freqs_and_displacements(self.file_matdyn_modes)

class chosen_motion(motion):
 def __init__(self):
       motion.__init__(self)
       self.freq=[]
       self.vib=[]
 def ask_which_q(self):
  self.vib,self.freq=func.ask_which_q(self.Q,self.DISPL,self.FREQ,
       self.no_of_modes)


##########DRAWING ATOMS AND CRYSTAL LATTICE
disp=display()
crystal_system=system()
crystal_system.read_crystal_info()
crystal_system.crystal=crystal_system.crystal_primitive
obj=chosen_motion()
obj.read_freqs_and_displacements()
obj.ask_which_q()
crystal_system.move_atoms_to_cell()
crystal_system.add_atoms_by_symmetry()
disp.set_scene(crystal_system.crystal)
disp.set_coord_system(crystal_system.alat)
disp.draw_lattice(crystal_system.crystal,crystal_system.crystal_primitive)
disp.draw_equilibrium_atoms(crystal_system.atoms)
disp.init_arrows(crystal_system.atoms,obj.vib[0])
disp.draw_displacement_arrows(crystal_system.atoms,obj.vib,obj.no_of_modes)

disp.legend(crystal_system.atoms)

crystal_system_conv=system()
crystal_system_conv.read_crystal_info()
crystal_system_conv.crystal=crystal_system.crystal_conv
crystal_system_conv.move_atoms_to_cell()
crystal_system_conv.add_atoms_by_symmetry()
disp.set_scene(crystal_system_conv.crystal)
disp.set_coord_system(crystal_system_conv.alat)
disp.draw_lattice(crystal_system_conv.crystal,crystal_system_conv.crystal_primitive)
disp.draw_equilibrium_atoms(crystal_system_conv.atoms)
disp.init_arrows(crystal_system_conv.atoms,obj.vib[0])
disp.draw_displacement_arrows(crystal_system_conv.atoms,obj.vib,obj.no_of_modes)

