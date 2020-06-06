from vpython import *
import sys
import numpy as np

eps1=1e-4
PRECIS=4

def read_crystal_info(file_scf_out, if_conv_cell):
  crystal,crystal2,atoms=[],[],[]
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
    alat=float(tmp[i].split()[-2])
   elif 'site n.     atom                  positions (alat units)' in tmp[i]:
    for m in range(nat):
     i=i+1
     tmpi=tmp[i].split()         
     atoms.append([tmpi[1],\
                   np.array([alat*float(tmpi[-4]),\
                             alat*float(tmpi[-3]),
                             alat*float(tmpi[-2])]),
                   m ])
   elif 'crystal axes: (cart. coord.' in tmp[i]:
    for j in range(3):
     i=i+1
     tmpi=tmp[i].split()
     crystal2.append(\
     np.array([alat*float(tmpi[3]),alat*float(tmpi[4]),alat*float(tmpi[5])]))
    if if_conv_cell==1 and ibrav in [1,2,3]:
     crystal=np.array([alat,alat,alat])
    else: crystal=crystal2
  return atoms,crystal,crystal2,alat

def move_one_atom_to_cell(atom_pos,crystal):
 cry=np.linalg.inv(np.transpose(crystal))
 v2=sum([atom_pos[m]*cry[:,m] for m in range(3)])
 for j in range(3):
   while v2[j]>1:
    v2[j]-=1
   while v2[j]<0:
    v2[j]+=1
 v2=sum([v2[m]*np.transpose(crystal)[:,m] for m in range(3)])
 return v2


def move_atoms_to_cell(atoms2,crystal):
 for i in atoms2:
  i[1]= move_one_atom_to_cell(i[1],crystal)
 return atoms2


def ruotaijk(s,k):
  return [s[0][0]*k[0]+s[0][1]*k[1]+s[0][2]*k[2],
  s[1][0]*k[0]+s[1][1]*k[1]+s[1][2]*k[2],
  s[2][0]*k[0]+s[2][1]*k[1]+s[2][2]*k[2]]

def add_atoms_by_symmetry(atoms,crystal,SYMM_OP):
 
 pm=[-1.,0.,1.]
 cr=[h*crystal[0]+k*crystal[1]+l*crystal[2] for h in pm for k in pm for l in pm if not (h==0 and k==0 and l==0)]
 for i in atoms:
  for e in cr:
   j= move_one_atom_to_cell(i[1]+e,crystal)
   sign=0
   for i2 in atoms:
    if (j[0]-i2[1][0])<eps1 and  (j[1]-i2[1][1])<eps1 and  (j[2]-i2[1][2])<eps1:
     sign=1
     break
   if sign==0:
    k=i[:]
    k[1]=j
    atoms.append(k)

 for i in atoms:
  for s in SYMM_OP:
   j=ruotaijk(s,i[1])
   j=move_one_atom_to_cell(j,crystal)
   sign=0
   for i2 in atoms:
    if (j[0]-i2[1][0])<eps1 and  (j[1]-i2[1][1])<eps1 and  (j[2]-i2[1][2])<eps1:
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

def ask_which_mode(Q,DISPL,FREQ,no_of_modes):
 nq=input( 'I listed qs calculated by matdyn. which one you want to visualize? Write  the number 0 -'+str((len(Q)-1))+':    ')
 nq=int(nq)
 nmode=input( 'Which mode you want to visualize? Write the number 0 -'+str((no_of_modes-1))+':    ')
 nmode=int(nmode)
 vib=DISPL[nq][nmode]
 freq=FREQ[nq]
 print('vib:'),
 for (numi,i) in enumerate(vib):
  print(numi)
 print('freq:', freq)
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

def add_atom(numi,i,atoms,atoms2,vv,vib2,vib):
    k=i[1]+vv
    atoms2.append([i[0],k,i[2]])
    vib2.append(vib[numi])
def if_zero(x):
 if abs(x)<0.02: return 1
 else: return 0


def set_scene(crystal):
 crystal_vec=[ vector(i[0],i[1],i[2]) for i in crystal]
 scene=canvas(width=900,height=450,background=color.white)
 scene.center=0.5*(crystal_vec[0]+crystal_vec[1]+crystal_vec[2])
 scene.parallel_projection = True
 [distant_light(direction=vector( 0.22,  0.44,  0.88),       color=color.gray(0.8)),
 distant_light(direction=vector(-0.88, -0.22, -0.44),       color=color.gray(0.3)),
distant_light(direction=vector( -0.22,  -0.44,  -0.88),       color=color.gray(0.8)),
 distant_light(direction=vector(0.88, 0.22, 0.44),       color=color.gray(0.3))]



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

def draw_displacement_arrows(atoms,vib,A,COLORS):
 arrows=[]
 for i in atoms:
   arrows.append(arrow(pos=vector(i[1][0],i[1][1],i[1][2]),\
                 axis=vector(A*vib[i[2]][0],A*vib[i[2]][1],A*vib[i[2]][2]),\
                 fixedwidth=True, shaftwidth=0.4,\
                 color=COLORS[atoms[i[2]][2]],shininess=0.1))


def set_sym_bl(a_vec):
#  ! Provides symmetry operations for all bravais lattices
#  ! Tests first the 24 proper rotations for the cubic lattice;
#  ! then the 8 rotations specific for the hexagonal axis (special axis c);
#  ! then inversion is added
 sin3 = 0.866025403784438597
 cos3 = 0.5
 msin3 =-0.866025403784438597
 mcos3 = -0.5
 # ! s0: the s matrices in cartesian axis
 # ! overlap: inverse overlap matrix between direct lattice
 # ! rat: the rotated of a direct vector ( cartesian )
 # ! rot: the rotated of a direct vector ( crystal axis )
 # ! value: component of the s matrix in axis basis
 # INTEGER :: jpol, kpol, mpol, irot, imat(32)
 # ! counters over the polarizations and the rotations

 S0= [[[ 1.,  0.,  0.],[  0.,  1.,  0.],[  0.,  0.,  1.]], 
          [[-1.,  0.,  0.],[  0., -1.,  0.],[  0.,  0.,  1.]],
          [[-1.,  0.,  0.],[  0.,  1.,  0.],[  0.,  0., -1.]],
          [[ 1.,  0.,  0.],[  0., -1.,  0.],[  0.,  0., -1.]],
          [[ 0.,  1.,  0.],[  1.,  0.,  0.],[  0.,  0., -1.]],
          [[ 0., -1.,  0.],[ -1.,  0.,  0.],[  0.,  0., -1.]],
          [[ 0., -1.,  0.],[  1.,  0.,  0.],[  0.,  0.,  1.]],
          [[ 0.,  1.,  0.],[ -1.,  0.,  0.],[  0.,  0.,  1.]],
          [[ 0.,  0.,  1.],[  0., -1.,  0.],[  1.,  0.,  0.]],
          [[ 0.,  0., -1.],[  0., -1.,  0.],[ -1.,  0.,  0.]],
          [[ 0.,  0., -1.],[  0.,  1.,  0.],[  1.,  0.,  0.]],
          [[ 0.,  0.,  1.],[  0.,  1.,  0.],[ -1.,  0.,  0.]],
          [[-1.,  0.,  0.],[  0.,  0.,  1.],[  0.,  1.,  0.]],
          [[-1.,  0.,  0.],[  0.,  0., -1.],[  0., -1.,  0.]],
          [[ 1.,  0.,  0.],[  0.,  0., -1.],[  0.,  1.,  0.]],
          [[ 1.,  0.,  0.],[  0.,  0.,  1.],[  0., -1.,  0.]],
          [[ 0.,  0.,  1.],[  1.,  0.,  0.],[  0.,  1.,  0.]],
          [[ 0.,  0., -1.],[ -1.,  0.,  0.],[  0.,  1.,  0.]],
          [[ 0.,  0., -1.],[  1.,  0.,  0.],[  0., -1.,  0.]],
          [[ 0.,  0.,  1.],[ -1.,  0.,  0.],[  0., -1.,  0.]],
          [[ 0.,  1.,  0.],[  0.,  0.,  1.],[  1.,  0.,  0.]],
          [[ 0., -1.,  0.],[  0.,  0., -1.],[  1.,  0.,  0.]],
          [[ 0., -1.,  0.],[  0.,  0.,  1.],[ -1.,  0.,  0.]],
          [[ 0.,  1.,  0.],[  0.,  0., -1.],[ -1.,  0.,  0.]],
          [[ cos3,  sin3, 0.],[ msin3,  cos3, 0.],[ 0., 0.,  1.]],
          [[ cos3, msin3, 0.],[  sin3,  cos3, 0.],[ 0., 0.,  1.]],
          [[mcos3,  sin3, 0.],[ msin3, mcos3, 0.],[ 0., 0.,  1.]],
          [[mcos3, msin3, 0.],[  sin3, mcos3, 0.],[ 0., 0.,  1.]],
          [[ cos3, msin3, 0.],[ msin3, mcos3, 0.],[ 0., 0., -1.]],
          [[ cos3,  sin3, 0.],[  sin3, mcos3, 0.],[ 0., 0., -1.]],
          [[mcos3, msin3, 0.],[ msin3,  cos3, 0.],[ 0., 0., -1.]],
          [[mcos3,  sin3, 0.],[  sin3,  cos3, 0.],[ 0., 0., -1.]]]
 
 S0NAME=['identity                                     ',
                '180 deg rotation - cart. axis [0,0,1]        ',
                '180 deg rotation - cart. axis [0,1,0]        ',
                '180 deg rotation - cart. axis [1,0,0]        ',
                '180 deg rotation - cart. axis [1,1,0]        ',
                '180 deg rotation - cart. axis [1,-1,0]       ',
                ' 90 deg rotation - cart. axis [0,0,-1]       ',
                ' 90 deg rotation - cart. axis [0,0,1]        ',
                '180 deg rotation - cart. axis [1,0,1]        ',
                '180 deg rotation - cart. axis [-1,0,1]       ',
                ' 90 deg rotation - cart. axis [0,1,0]        ',
                ' 90 deg rotation - cart. axis [0,-1,0]       ',
                '180 deg rotation - cart. axis [0,1,1]        ',
                '180 deg rotation - cart. axis [0,1,-1]       ',
                ' 90 deg rotation - cart. axis [-1,0,0]       ',
                ' 90 deg rotation - cart. axis [1,0,0]        ',
                '120 deg rotation - cart. axis [-1,-1,-1]     ',
                '120 deg rotation - cart. axis [-1,1,1]       ',
                '120 deg rotation - cart. axis [1,1,-1]       ',
                '120 deg rotation - cart. axis [1,-1,1]       ',
                '120 deg rotation - cart. axis [1,1,1]        ',
                '120 deg rotation - cart. axis [-1,1,-1]      ',
                '120 deg rotation - cart. axis [1,-1,-1]      ',
                '120 deg rotation - cart. axis [-1,-1,1]      ',
                ' 60 deg rotation - cryst. axis [0,0,1]       ',
                ' 60 deg rotation - cryst. axis [0,0,-1]      ',
                '120 deg rotation - cryst. axis [0,0,1]       ',
                '120 deg rotation - cryst. axis [0,0,-1]      ',
                '180 deg rotation - cryst. axis [1,-1,0]      ',
                '180 deg rotation - cryst. axis [2,1,0]       ',
                '180 deg rotation - cryst. axis [0,1,0]       ',
                '180 deg rotation - cryst. axis [1,1,0]       ',
                'inversion                                    ',
                'inv. 180 deg rotation - cart. axis [0,0,1]   ',
                'inv. 180 deg rotation - cart. axis [0,1,0]   ',
                'inv. 180 deg rotation - cart. axis [1,0,0]   ',
                'inv. 180 deg rotation - cart. axis [1,1,0]   ',
                'inv. 180 deg rotation - cart. axis [1,-1,0]  ',
                'inv.  90 deg rotation - cart. axis [0,0,-1]  ',
                'inv.  90 deg rotation - cart. axis [0,0,1]   ',
                'inv. 180 deg rotation - cart. axis [1,0,1]   ',
                'inv. 180 deg rotation - cart. axis [-1,0,1]  ',
                'inv.  90 deg rotation - cart. axis [0,1,0]   ',
                'inv.  90 deg rotation - cart. axis [0,-1,0]  ',
                'inv. 180 deg rotation - cart. axis [0,1,1]   ',
                'inv. 180 deg rotation - cart. axis [0,1,-1]  ',
                'inv.  90 deg rotation - cart. axis [-1,0,0]  ',
                'inv.  90 deg rotation - cart. axis [1,0,0]   ',
                'inv. 120 deg rotation - cart. axis [-1,-1,-1]',
                'inv. 120 deg rotation - cart. axis [-1,1,1]  ',
                'inv. 120 deg rotation - cart. axis [1,1,-1]  ',
                'inv. 120 deg rotation - cart. axis [1,-1,1]  ',
                'inv. 120 deg rotation - cart. axis [1,1,1]   ',
                'inv. 120 deg rotation - cart. axis [-1,1,-1] ',
                'inv. 120 deg rotation - cart. axis [1,-1,-1] ',
                'inv. 120 deg rotation - cart. axis [-1,-1,1] ',
                'inv.  60 deg rotation - cryst. axis [0,0,1]  ',
                'inv.  60 deg rotation - cryst. axis [0,0,-1] ',
                'inv. 120 deg rotation - cryst. axis [0,0,1]  ',
                'inv. 120 deg rotation - cryst. axis [0,0,-1] ',
                'inv. 180 deg rotation - cryst. axis [1,-1,0] ',
                'inv. 180 deg rotation - cryst. axis [2,1,0]  ',
                'inv. 180 deg rotation - cryst. axis [0,1,0]  ',
                'inv. 180 deg rotation - cryst. axis [1,1,0]  ' ]

#####finding the symmetries 
#    compute the overlap matrix for crystal axis
 ROT=np.array([[sum( [a_vec[kpol][i]*a_vec[jpol][i] for i in range(3)]) for kpol in range(3)] for jpol in range(3)])
 OVERLAP=np.linalg.inv(ROT)
 S=np.zeros((48,3,3))
 SNAME=[ [] for i in range(48)]
 IMAT=[ 0 for i in range(32)]
 nrot=0
 for irot in range(32):
  isign=0
  for jpol in range(3):
   #compute, in cartesian coordinates the rotated vector
   RAT=[ sum([S0[irot][i][mpol]*a_vec[jpol][i] for i in range(3)]) for mpol in range(3)]
   ROT[jpol]=[ sum([a_vec[kpol][i]*RAT[i] for i in range(3)]) for kpol in range(3)]
  # and the inverse of the overlap matrix is applied
  for jpol in range(3):
   for kpol in range(3):
    value=round(sum([OVERLAP[i][jpol]*ROT[kpol][i] for i in range(3)]),PRECIS)
    if abs((int(value))-value)>eps1:
              # if a noninteger is obtained, this implies that this operation
              # is not a symmetry operation for the given lattice
     isign=1
     break
    else: S[nrot][jpol][kpol]=(int(value))
   if isign==1: break
  if isign==1: continue
  SNAME[nrot]=S0NAME[irot]
  IMAT[nrot]=irot
  nrot=nrot+1
#####end of finding the symmetries

 #check number of symmetries
 temp_numbers=[1,2,4,6,8,12,24]
 print (nrot)
 if nrot not in temp_numbers:
  print("NOTICE: Bravais lattice has wrong number  of symmetries - symmetries are disabled")
  nrot = 0

 #set the inversion symmetry ( Bravais lattices have always inversion
 # !     symmetry )
 for irot in range(nrot):
  SNAME[irot+nrot]=S0NAME[IMAT[irot]+32]
  for kpol in range(3):
   for jpol in range(3):
    S[irot+nrot][jpol][kpol]=-S[irot][jpol][kpol]
 nrot=2*nrot

 print("Found "+str(nrot)+" symmetry operations")
 return  [np.transpose(a) for a in S[:nrot]]



#### END OF FUNCTIONS


