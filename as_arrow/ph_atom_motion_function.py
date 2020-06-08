from vpython import *
import sys
import numpy as np

eps1=1e-4
PRECIS=4

def round_vec(v):
 return np.array([ round(i,PRECIS) for i in v])

def read_crystal_info(file_scf_out, if_conv_cell):
  crystal,crystal_primitive,atoms=[],[],[]
  try: h=open(file_scf_out,'r')
  except: 
   print('Couldnt find scf.out file. Copy it here')
   sys.exit(1)
  tmp=h.readlines()
  h.close()
  for i in range(len(tmp)):
   if 'bravais-lattice index' in tmp[i]:
    ibrav=int(tmp[i].split()[-1])
   if ' number of atoms/cell' in tmp[i]: 
    nat=int(tmp[i].split()[-1])
   elif 'lattice parameter (alat)' in tmp[i]:
    alat=round(float(tmp[i].split()[-2]),PRECIS)
   elif 'site n.     atom                  positions (alat units)' in tmp[i]:
    for m in range(nat):
     i=i+1
     tmpi=tmp[i].split()         
     atoms.append([tmpi[1],\
                   round_vec([alat*float(tmpi[-4]),\
                             alat*float(tmpi[-3]),
                             alat*float(tmpi[-2])]),
                   m ])
   elif 'crystal axes: (cart. coord.' in tmp[i]:
    for j in range(3):
     i=i+1
     tmpi=tmp[i].split()
     crystal_primitive.append(\
     round_vec([alat*float(tmpi[3]),alat*float(tmpi[4]),alat*float(tmpi[5])]))
    if if_conv_cell==1 and ibrav in [1,2,3]:
     crystal=np.array([alat*round_vec(m) for m in [[1,0,0],[0,1,0],[0,0,1]]])
    else: crystal=crystal_primitive
  return atoms,crystal,crystal_primitive,alat

def move_one_atom_to_cell(atom_pos,crystal):
 cry=np.linalg.inv(np.transpose(crystal))
 v2=sum([atom_pos[m]*cry[:,m] for m in range(3)])
 for j in range(3):
   while v2[j]>1+eps1:
    v2[j]-=1
   while v2[j]<0-eps1:
    v2[j]+=1
 v2=round_vec(sum([v2[m]*np.transpose(crystal)[:,m] for m in range(3)]))
 return v2


def move_atoms_to_cell(atoms2,crystal):
 for i in atoms2:
  i[1]= move_one_atom_to_cell(i[1],crystal)
 return atoms2


def add_atoms_by_symmetry(atoms,crystal,crystal_primitive):
 #add atom at the faces
 pm=[-1.,0.,1.]
 cr=[round_vec(h*crystal_primitive[0]+k*crystal_primitive[1]+l*crystal_primitive[2]) for h in pm for k in pm for l in pm if not (h==0 and k==0 and l==0)]
 for i in atoms[:len(atoms)]:
  for e in cr:
   j= move_one_atom_to_cell(i[1]+e,crystal)
   sign=0
   for i2 in atoms:
    if abs(j[0]-i2[1][0])<eps1 and  abs(j[1]-i2[1][1])<eps1 and  abs(j[2]-i2[1][2])<eps1:
     sign=1
     break
   if sign==0:
    k=i[:]
    k[1]=j
    atoms.append(k)
 return atoms

def ask_if_conv_cell():
 try: if_conv_cell=int(raw_input('Display primitive [0] or conventional cell [1]?'))
 except: if_conv_cell=int(input('Display primitive [0] or conventional cell [1]?'))
 return if_conv_cell

def read_freqs_and_displacements(file_matdyn_modes):
 DISPL=[]
 FREQ=[]
 Q=[]
 h=open(file_matdyn_modes,'r')
 tmp=h.readlines()
 h.close() 
 i=2
 while i<len(tmp):
  if 'q' in tmp[i].split(): 
   DISPL.append([])
   FREQ.append([])
   Q.append(tmp[i])
   i=i+2
   while '***' not in tmp[i]:
    if 'freq' in tmp[i]: 
     FREQ[-1].append(float(tmp[i].split()[4]))
     DISPL[-1].append([])
    else:
     tmpi=tmp[i].split()
     DISPL[-1][-1].append([\
                   float(tmpi[1]), float(tmpi[3]),float(tmpi[5])])    
    i=i+1
  i=i+3
 for i in range(len(Q)):
  print( i,Q[i],),
 no_of_modes=len(FREQ[-1])
 return DISPL,FREQ,Q,no_of_modes

def ask_which_q(Q,DISPL,FREQ,no_of_modes):
 nq=input( 'I listed qs calculated by matdyn. which one you want to visualize? Write  the number 0 -'+str((len(Q)-1))+':    ')
 nq=int(nq)
 #nmode=input( 'Which mode you want to visualize? Write the number 0 -'+str((no_of_modes-1))+':    ')
 #nmode=int(nmode)
 vib=DISPL[nq] #[nmode]
 freq=FREQ[nq]
 return vib,freq 


def draw_lattice(crystal,crystal2):
 C=[np.array([0,0,0]), crystal[0],crystal[0]+crystal[1],\
   crystal[0]+crystal[1]+crystal[2],crystal[0]+crystal[2],\
   crystal[2],\
   np.array([0,0,0]),crystal[1],crystal[1]+crystal[0],\
   crystal[0]+crystal[1]+crystal[2],crystal[1]+crystal[2],\
   crystal[2],\
   crystal[2]+crystal[0],\
   crystal[0],np.array([0,0,0]),crystal[1],crystal[1]+crystal[2]]
 C2=[np.array([0,0,0]), crystal2[0],crystal2[0]+crystal2[1],\
   crystal2[0]+crystal2[1]+crystal2[2],crystal2[0]+crystal2[2],\
   crystal2[2],\
   np.array([0,0,0]),crystal2[1],crystal2[1]+crystal2[0],\
   crystal2[0]+crystal2[1]+crystal2[2],crystal2[1]+crystal2[2],\
   crystal2[2],\
   crystal2[2]+crystal2[0], crystal2[0],np.array([0,0,0]),crystal2[1],crystal2[1]+crystal2[2]]

 crystal_lattice=curve(color=color.black,radius=0.1)
 for i in C:
  crystal_lattice.append(vector(i[0],i[1],i[2]))
 return C,C2



def set_scene(crystal):
 crystal_vec=[ vector(i[0],i[1],i[2]) for i in crystal]
 scene=canvas(width=900,height=450,background=color.white)
 scene.center=0.5*(crystal_vec[0]+crystal_vec[1]+crystal_vec[2])
 #scene.parallel_projection = True
 scene.autocenter=0.5*(crystal_vec[0]+crystal_vec[1]+crystal_vec[2])
 scene.stereo = 'active'
 scene.stereodepth = 0
 default_camera_pos=scene.camera.pos
 def B1(b):
    scene.camera.pos-=0.1*(scene.camera.pos-scene.center)
 def B2(b):
    scene.camera.pos+=0.1*(scene.camera.pos-scene.center)
 scene.append_to_caption('       ')
 def C1(b):
    scene.camera.rotate(angle=-.1, axis=scene.up)
 def C2(b):
    scene.camera.rotate(angle=.1, axis=scene.up)
 def C3(b):
    scene.camera.pos-=scene.up
 def C4(b):
    scene.camera.pos+=scene.up

 but=button( bind=B2, text='ZOOM-',height=100,value=0)
 but=button( bind=B1, text='ZOOM+',height=100,value=0)
 scene.append_to_caption('\n')
 scene.append_to_caption('       ')
 but=button( bind=C3, text='UP',height=100)
 scene.append_to_caption('\n')
 but=button( bind=C1, text='LEFT',height=100,value=0)
 but=button( bind=C2, text='RIGHT',height=100,value=0)
 scene.append_to_caption('\n')
 scene.append_to_caption('    ')
 but=button( bind=C4, text='DOWN',height=100)
 return scene
 #for i in range(no_of_modes):
 # but=button( bind=M(i),
 
# [distant_light(direction=vector( 0.22,  0.44,  0.88),       color=color.gray(0.8)),
# distant_light(direction=vector(-0.88, -0.22, -0.44),       color=color.gray(0.3)),
#distant_light(direction=vector( -0.22,  -0.44,  -0.88),       color=color.gray(0.8)),
# distant_light(direction=vector(0.88, 0.22, 0.44),       color=color.gray(0.3))]



def set_coord_system(alat):
 coord_sys=[\
arrow(pos=-vector(4*alat,0,0), axis=vector(alat,0,0),color=color.green),\
arrow(pos=-vector(4*alat,0,0), axis=vector(0,alat,0),color=color.red),\
arrow(pos=-vector(4*alat,0,0), axis=vector(0,0,alat),color=color.blue)]
 coord_sys_lab=[\
label(pos=coord_sys[0].pos+coord_sys[0].axis,text='<b>x</b>',color=color.green,box=False),\
label(pos=coord_sys[1].pos+coord_sys[1].axis,text='<b>y</b>',color=color.red,box=False),\
label(pos=coord_sys[2].pos+coord_sys[2].axis,text='<b>z</b>',color=color.blue,box=False)]


def draw_equilibrium_atoms(atoms,COLORS):
 at1_equil=[ sphere(pos=vector(i[1][0],i[1][1],i[1][2]),\
                    radius=0.7,color=COLORS[i[2]]) \
        for i in (atoms)]

def init_arrows(atoms,vib,A,COLORS):
 arrows=[]
 for i in atoms:
   arrows.append(arrow(pos=vector(i[1][0],i[1][1],i[1][2]),\
                 axis=vector(A*vib[i[2]][0],A*vib[i[2]][1],A*vib[i[2]][2]),\
                 fixedwidth=True, shaftwidth=0.4,\
                 color=COLORS[atoms[i[2]][2]],shininess=0.1))
 return arrows

def draw_displacement_arrows(scene,arrows,atoms,vib,A,no_of_modes):
 scene.append_to_caption('\nChooose no of mode:\n')
 but=[]
 def F(b):
  m=int(b.text)
  for numi,i in enumerate(atoms):
   arrows[numi].axis=vector(A*vib[m][i[2]][0],A*vib[m][i[2]][1],A*vib[m][i[2]][2])
 for k in [i for i in range(no_of_modes)]:
  but.append(button( bind=F , text=str(k),height=100))

