from vpython import *
import sys
import numpy as np
import ph_atom_motion_function as func

class display():
    def __init__(self):
     self.COLORS=[color.red , color.yellow , color.green, color.purple , color.blue 	, color.cyan 	, color.orange 	, color.magenta ,color.orange, color.black,\
color.gray(0.9),color.gray(0.8),color.gray(0.7),color.gray(0.6)	 ]
     self.crystal_lattice=[]
     self.arrows=[]
     self.atomic_balls=[]
     self.moving_atoms=[]
     self.scene=()
     self.A=10 #amplitude. The displacement is multiplied by A
     self.coord_system=[]
    def set_scene(self,crystal):
        self.scene=func.set_scene(crystal[:])
    def set_coord_system(self,alat):
         self.coord_system=func.set_coord_system(alat,self.scene)
    def draw_lattice(self,crystal,crystal_primitive):
        self.crystal_lattice=func.draw_lattice(crystal[:],crystal_primitive[:])    
    def draw_equilibrium_atoms(self,equil_atoms):
        self.atomic_balls=func.draw_equilibrium_atoms(equil_atoms[:],self.COLORS[:])   
    def init_arrows(self,equil_atoms,vib):
        self.arrows,self.moving_atoms=func.init_arrows\
                                      (equil_atoms[:],vib,self.A,self.COLORS[:])     
    def draw_displacement_arrows_and_balls(self,atoms,vib,freq,no_of_modes):
        func.draw_displacement_arrows\
                    (self.scene,self.arrows[:],self.moving_atoms[:],atoms[:],vib[:],freq,\
                     self.A,no_of_modes)    
    def draw_moving_atoms(self,equil_atoms,vib,freq):
        func.draw_moving_atoms\
                                  (equil_atoms[:],\
                                   vib[:],freq[:],self.COLORS[:])   
    def draw_atomic_bondings(self,equil_atoms):
        self.atomic_bondings=func.draw_atomic_bondings(equil_atoms[:])
    def legend(self,atoms):
        func.legend(atoms[:],self.COLORS[:])
    def choose_color(self,irr_atoms,atoms):
        func.choose_color(irr_atoms,atoms,\
                  self.atomic_balls,self.moving_atoms,self.arrows,\
                  self.scene,self.COLORS[:])
    def if_display_tetrahedrons(self,atoms):
        func.if_display_tetrahedrons(atoms,self.COLORS,self.scene)

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
        self.crystal_primitive=[]
        self.alat=0
        self.celldm=[]
        self.ibrav=0
 def read_crystal_info(self):
  self.atoms,self.crystal_primitive,self.alat,self.celldm,self.ibrav=\
       func.read_crystal_info(self.file_scf_out)
# def set_symm(self):
#  self.SYMM_OP= func.set_sym_bl(self.crystal_primitive)
 def move_atoms_to_cell(self):
  self.atoms=func.move_atoms_to_cell(self.atoms[:],self.crystal[:])
 def add_atoms_by_symmetry(self):
  self.atoms=func.add_atoms_by_symmetry(self.atoms[:],\
             self.crystal[:],self.crystal_primitive[:])
 def make_conventional_cell(self):
  self.crystal=func.make_conv_cell(self.ibrav,self.crystal_primitive[:],self.celldm[:])



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
  self.vib,self.freq=func.ask_which_q(self.Q,self.DISPL[:],self.FREQ[:],
       self.no_of_modes)


##########DRAWING ATOMS AND CRYSTAL LATTICE
#scene
disp=display()
#primitive cell
crystal_system=system()
crystal_system.read_crystal_info()
irr_atoms=crystal_system.atoms
crystal_system.crystal=crystal_system.crystal_primitive
obj=chosen_motion()
obj.read_freqs_and_displacements()
obj.ask_which_q()
crystal_system.move_atoms_to_cell()
crystal_system.add_atoms_by_symmetry()
disp.set_scene(crystal_system.crystal)
disp.draw_atomic_bondings(crystal_system.atoms)
disp.set_coord_system(crystal_system.alat)
disp.draw_lattice(crystal_system.crystal,crystal_system.crystal_primitive)
disp.draw_equilibrium_atoms(crystal_system.atoms)
disp.init_arrows(crystal_system.atoms,obj.vib[0])
disp.draw_displacement_arrows_and_balls(crystal_system.atoms,obj.vib,obj.freq,obj.no_of_modes)

#legend
disp.legend(crystal_system.atoms)
#conventional cell
disp2=display()
crystal_system_conv=system()
crystal_system_conv.read_crystal_info()
crystal_system_conv.make_conventional_cell()
crystal_system_conv.move_atoms_to_cell()
crystal_system_conv.add_atoms_by_symmetry()
disp2.set_scene(crystal_system_conv.crystal)
disp2.draw_atomic_bondings(crystal_system_conv.atoms)
disp2.set_coord_system(crystal_system_conv.alat)
disp2.draw_lattice(crystal_system_conv.crystal,crystal_system_conv.crystal_primitive)
disp2.draw_equilibrium_atoms(crystal_system_conv.atoms)
disp2.init_arrows(crystal_system_conv.atoms,obj.vib[0])
disp2.if_display_tetrahedrons(crystal_system_conv.atoms)
disp2.draw_displacement_arrows_and_balls(crystal_system_conv.atoms,obj.vib,obj.freq,obj.no_of_modes)
disp2.choose_color(irr_atoms,crystal_system_conv.atoms)


