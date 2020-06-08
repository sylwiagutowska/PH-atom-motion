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

class inputs():
 def __init__(self):
        self.file_scf_out='scf.out'
        self.file_matdyn_in='matdyn.in'
        self.file_matdyn_modes='matdyn.modes'
        self.if_conv_cell=0
 def ask_if_conv_cell(self):
    self.if_conv_cell=func.ask_if_conv_cell()

class system(inputs):
 def __init__(self):
        inputs.__init__(self)
        self.atoms=[]
        self.crystal=[]
        self.crystal_primitive=[]
       # self.SYMM_OP=[]
        self.alat=0
 def read_crystal_info(self):
  self.atoms,self.crystal,self.crystal_primitive,self.alat=\
       func.read_crystal_info(self.file_scf_out,self.if_conv_cell)
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
crystal_system.ask_if_conv_cell()
obj=chosen_motion()
obj.read_freqs_and_displacements()
obj.ask_which_q()
obj.if_conv_cell=0
crystal_system.read_crystal_info()
#obj.if_conv_cell=crystal_system.if_conv_cell
#crystal_system.set_symm()
crystal_system.move_atoms_to_cell()
crystal_system.add_atoms_by_symmetry()
disp.set_scene(crystal_system.crystal)
disp.set_coord_system(crystal_system.alat)
disp.draw_lattice(crystal_system.crystal,crystal_system.crystal_primitive)
disp.draw_equilibrium_atoms(crystal_system.atoms)
disp.init_arrows(crystal_system.atoms,obj.vib[0])
disp.draw_displacement_arrows(crystal_system.atoms,obj.vib,obj.no_of_modes)




######BUTTONS
'''
def f(s):
 print (posi[int(s)])
 return 2*posi[int(s)]
def keyInput(evt):
    s =evt.key.split()[0]
    try: scene.camera.center(f(s))
    except: print(s+'dds')
    scene.capture('a') 
scene.bind('keydown', keyInput)
'''
'''
def camera_center(evt):
    loc = evt.pos
    scene.camera.axis=loc #+vector(0,0,0.5*alat)
scene.bind('click', camera_center)
'''
'''
def B(b):
    scene.camera.pos=scene.camera.pos-vector(10,0,0)
def C(b):
    scene.camera.pos=scene.camera.pos+vector(10,0,0)
def D(b):
    scene.camera.pos=scene.camera.pos-vector(0,10,0)
def E(b):
    scene.camera.pos=scene.camera.pos+vector(0,10,0)
def F(b):
    if but.value==0: 
     labele.visible=False
     but.value=1
    else:
     labele.visible=True
     but.value=0

def G(b):
    scene.up=vector(1,0,0)
    scene.camera.rotate(angle=90., axis=vector(0,0,1), origin=vector(0,0,0))
    scene.center=0.5*(crystal[0]+crystal[1]+crystal[2])
    
def H(b):
    scene.up=vector(-1,0,0)
    scene.camera.rotate(angle=-90, axis=vector(0,0,1), origin=vector(0,0,0))
    scene.center=0.5*(crystal[0]+crystal[1]+crystal[2])
def I(b):
    scene.up=vector(0,1,0)
    scene.camera.rotate(angle=90., axis=vector(0,0,1), origin=vector(0,0,0))
    scene.center=0.5*(crystal[0]+crystal[1]+crystal[2])
def J(b):
    scene.up=vector(0,-1,0)
    scene.camera.rotate(angle=-90., axis=vector(0,0,1),  origin=vector(0,0,0))
    scene.center=0.5*(crystal[0]+crystal[1]+crystal[2])
def K(b):
    scene.up=vector(0,0,0.1)
    scene.camera.rotate(angle=90., axis=vector(0,0,0.1),  origin=vector(0,0,0))
    scene.center=0.5*(crystal[0]+crystal[1]+crystal[2])
def L(b):
    scene.up=vector(0,0,-0.1)
    scene.camera.rotate(angle=-90., axis=vector(0,0,1),  origin=vector(0,0,0))
    scene.center=0.5*(crystal[0]+crystal[1]+crystal[2])

def M(b):
    scene.camera.pos=scene.camera.pos*0.5
    scene.center=0.5*(crystal[0]+crystal[1]+crystal[2])
def N(b):
    scene.camera.pos=scene.camera.pos*2
    scene.center=0.5*(crystal[0]+crystal[1]+crystal[2])

scene.append_to_caption('       ')
but=button( bind=D, text='TOP',height=100)
scene.append_to_caption('\n')
but=button( bind=B, text='LEFT',height=100)
but=button( bind=C, text='RIGHT',height=100)
scene.append_to_caption('\n    ')
but=button( bind=E, text='DOWN',height=100)
scene.append_to_caption('\n    ')
but=button( bind=F, text='DISABLE/ENABLE LABELS',height=100,value=0)
scene.append_to_caption('\n    ')

but=button( bind=G, text='ROTATE X+',height=100,value=0)
but=button( bind=H, text='ROTATE X-',height=100,value=0)
scene.append_to_caption('\n    ')
but=button( bind=I, text='ROTATE Y+',height=100,value=0)
but=button( bind=J, text='ROTATE Y-',height=100,value=0)
scene.append_to_caption('\n    ')
but=button( bind=K, text='ROTATE Z+',height=100,value=0)
but=button( bind=L, text='ROTATE Z-',height=100,value=0)
scene.append_to_caption('\n    ')
but=button( bind=N, text='ZOOM+',height=100,value=0)
but=button( bind=M, text='ZOOM-',height=100,value=0)
'''

