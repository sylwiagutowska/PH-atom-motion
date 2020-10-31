from vpython import *
import sys
import numpy as np

eps1=1e-4
PRECIS=4
eps2=1e-2

def round_vec(v):
 return np.array([ round(i,PRECIS) for i in v])

def make_vector(v):
 return vector(v[0],v[1],v[2])

def read_crystal_info(file_scf_out):
  crystal_conv,crystal_primitive,atoms=[],[],[]
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
   elif 'celldm(1)=' in tmp[i]:
    c_tmp=tmp[i].split()+tmp[i+1].split()
    celldm=[c_tmp[1],c_tmp[3],c_tmp[5],c_tmp[7], c_tmp[9],c_tmp[11]]
    celldm=[ round(float(m),PRECIS) for m in celldm]
    for numi in range(1,3):
     celldm[numi]=celldm[numi]*celldm[0]   
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
    crystal_primitive=np.array(crystal_primitive)
  return atoms,crystal_primitive,alat,celldm,ibrav

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
 pm=[-1.,0.,1.,-2.,2.,-3.,3.]
 cr=[round_vec(h*crystal_primitive[0]+k*crystal_primitive[1]+l*crystal_primitive[2]) for h in pm for k in pm for l in pm if not (h==0 and k==0 and l==0)]
 for i in atoms[:]:
  for e in cr:
   j= move_one_atom_to_cell(i[1]+e,crystal)
   sign=0
   for i2 in atoms:
     if abs(j[0]-i2[1][0])<eps2 and  abs(j[1]-i2[1][1])<eps2 and  abs(j[2]-i2[1][2])<eps2:
       sign=1
       break
   if sign==0:
    k=i[:]
    k[1]=j
    atoms.append(k)
 return atoms

def make_conv_cell(ibrav,crystal_primitive,celldm):
 cart=[round_vec(m) for m in [[1,0,0],[0,1,0],[0,0,1]]]
 if ibrav in [1,2,3]: crystal_conv=np.array([round_vec(celldm[0]*m) for m in cart])
 elif ibrav in [4,6,8]: crystal_conv=crystal_primitive[:]
 elif ibrav in [9]:  crystal_conv=np.array([round_vec(celldm[m]*cart[m]) for m in range(3)])
 elif ibrav==5: 
  M=np.array([[-1.,1.,0.],[1.,0.,-1.],[1.,1.,1.]])
  crystal_conv=np.dot(M,np.array(crystal_primitive))
 return crystal_conv



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
     DISPL[-1][-1].append(round_vec([\
                   float(tmpi[1]), float(tmpi[3]),float(tmpi[5])]))    
    i=i+1
  i=i+3
 for i in range(len(Q)):
  print( i,Q[i][:-1]),
 no_of_modes=len(FREQ[-1])
 return DISPL,FREQ,Q,no_of_modes

def ask_which_q(Q,DISPL,FREQ,no_of_modes):
 nq=input( 'I listed qs calculated by matdyn. which one you want to visualize? Write  the number 0 -'+str((len(Q)-1))+':    ')
 nq=int(nq)
 q=Q[nq]
 vib=DISPL[nq] #[nmode]
 freq=FREQ[nq]
 return vib,freq,q 



def set_scene(crystal,disp):
 crystal_vec=[ make_vector(i) for i in crystal]
 scene=canvas(width=900,height=450,background=color.white)
 scene.center=0.5*(crystal_vec[0]+crystal_vec[1]+crystal_vec[2])
 #scene.parallel_projection = True
 scene.autocenter=0.5*(crystal_vec[0]+crystal_vec[1]+crystal_vec[2])
 scene.stereo = 'active'
 scene.stereodepth = 0
 scene.fov=0.001
 default_camera_pos=scene.camera.pos
 def B1(b):
    scene.camera.pos-=0.1*(scene.camera.pos-scene.center)
 def B2(b):
    scene.camera.pos+=0.1*(scene.camera.pos-scene.center)
 scene.append_to_caption('       ')
 def C1(b):
    rotate(disp,np.pi/4,vec(0,0,1),scene.center)
 def C2(b):
    rotate(disp,np.pi/4,vec(1,0,0),scene.center)

 but=button( bind=B2, text='ZOOM-',height=100,value=0)
 but=button( bind=B1, text='ZOOM+',height=100,value=0)
 scene.append_to_caption('\n')
 but=button( bind=C1, text='ROTATE XY',height=100,value=0)
 but=button( bind=C2, text='ROTATE YZ',height=100,value=0)
 scene.append_to_caption('\n')
 return scene
 #for i in range(no_of_modes):
 # but=button( bind=M(i),
 
# [distant_light(direction=vector( 0.22,  0.44,  0.88),       color=color.gray(0.8)),
# distant_light(direction=vector(-0.88, -0.22, -0.44),       color=color.gray(0.3)),
#distant_light(direction=vector( -0.22,  -0.44,  -0.88),       color=color.gray(0.8)),
# distant_light(direction=vector(0.88, 0.22, 0.44),       color=color.gray(0.3))]


def set_coord_system(alat,scene):
 coord_sys=[\
arrow(pos=-vector(4*alat,0,0), axis=vector(alat,0,0),color=color.green),\
arrow(pos=-vector(4*alat,0,0), axis=vector(0,alat,0),color=color.red),\
arrow(pos=-vector(4*alat,0,0), axis=vector(0,0,alat),color=color.blue)]
 coord_sys_lab=[\
label(pos=coord_sys[0].pos+coord_sys[0].axis,text='<b>x</b>',color=color.green,box=False),\
label(pos=coord_sys[1].pos+coord_sys[1].axis,text='<b>y</b>',color=color.red,box=False),\
label(pos=coord_sys[2].pos+coord_sys[2].axis,text='<b>z</b>',color=color.blue,box=False)]
 coord_system=[coord_sys,coord_sys_lab]
 def D(b):
    for i in coord_system:
     if i[0].visible==True:
      for j in i:
       j.visible=False
     else:  
      for j in i:
       j.visible=True
 but=button( bind=D, text='Coord.sys.',height=100)


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
  crystal_lattice.append(make_vector(i))
 return crystal_lattice



def draw_equilibrium_atoms(atoms,COLORS):
 at_equil=[ sphere(pos=make_vector(i[1]),\
                    radius=0.7,color=COLORS[i[2]]) \
        for i in (atoms)]
 return at_equil

def init_arrows(atoms,vib,A,COLORS):
 arrows=[]
 moving_atoms=[]
 for i in atoms:
   arrows.append(arrow(pos=make_vector(i[1]),\
                 axis=make_vector(A*vib[i[2]]),\
                 fixedwidth=True, shaftwidth=0.4,\
                 color=COLORS[atoms[i[2]][2]],shininess=0.1))
   moving_atoms=[0,[ sphere(pos=make_vector(i[1]),\
                    radius=0.7,color=COLORS[i[2]], make_trail=True,visible=False) \
        for i in (atoms)]] #first element stands for no of mode
 return arrows,moving_atoms


def draw_displacement_arrows(scene,arrows,moving_atoms,atoms,vib,freq,A,no_of_modes):
 mode_buttons=[]
 #moving atoms
 def G(b):
   t,dt=0,0.05
   if b.checked==False:
    for i in moving_atoms[1]: 
     i.clear_trail() 
     i.visible=False
   else:
    for i in moving_atoms[1]: i.visible=True  
   while b.checked==True:
    rate(10)
    for numi,i in enumerate(arrows):
     moving_atoms[1][numi].pos=(i.pos+(i.axis*sin(freq[moving_atoms[0]]*t)))
    t+=dt  

 moving_atoms_button=radio(bind=G,text="moving atoms on/off")
 #arrows
 def F(b):
  m=int(b.text.split('\n')[0])-1
  #dont know why, but it HAS TO be done twice, otherwise not all atoms change their arrows :(
  for numi,i in enumerate(atoms):
   arrows[numi].length=make_vector(A*vib[m][i[2]]).mag
   arrows[numi].axis=make_vector(A*vib[m][i[2]])
   moving_atoms[0]=m
  for numi,i in enumerate(atoms):
   arrows[numi].length=make_vector(A*vib[m][i[2]]).mag
   arrows[numi].axis=make_vector(A*vib[m][i[2]])
   moving_atoms[0]=m
 scene.append_to_caption('\nChoose mode:\n')
 for k in range(no_of_modes):
  mode_buttons.append(button( bind=F , text=str(k+1)+'\n'+str(round(freq[k],2)),height=100))


def draw_atomic_bondings(atoms): 
 atomic_bondings=[]
 mini=[]
 for numi,i in enumerate(atoms):
  dist=[]
  for numj,j in enumerate(atoms):
   dist.append([(sum([m**2 for m in i[1]-j[1]]))**0.5,i,j])
  min_dist=min([j[0] for j in dist if j[0]>1e-1])
  for numj, j in enumerate(dist):
   if abs(j[0]-min_dist)<1e-1: mini.append(j)
 atomic_bondings=[curve(make_vector(i[1][1]),make_vector(i[2][1]),radius=.05) for i in mini]
 def E(b):
  if atomic_bondings[0].visible==False: 
   for i in atomic_bondings:
    i.visible = True
  else: 
   for i in atomic_bondings:
    i.visible = False
 but=button( bind=E, text='bondings on/off',height=100)
 return atomic_bondings


def legend(atoms,COLORS): 
 scene=canvas(width=900,height=60,background=color.white)
 scene.parallel_projection = True
 scene.center=vector(0.5,0,0)
 scene.fov=0.001
 scene.camera.pos=vector(0.5,0,-.25)
 scene.stereo = 'active'
 scene.stereodepth = 0
 scene.append_to_caption('\n')
 label(pos=vector(-1,0,0),text='Legend',box=False) 
 legend_atoms=[]
 for i in atoms:
  sign=0
  for j in legend_atoms:
   if i[2]==j[2]:
    sign=1
    break
  if sign==0: legend_atoms.append(i)
 scene.center=vector(int(len(legend_atoms)/2)*0.25-0.1,0,0)
 scene.camera.pos=vector(int(len(legend_atoms)/2)*0.25-0.1,0,-2)
 atom_bals=[ sphere(pos=vector(numi*0.25-0.1,0,0), radius=0.04,color=COLORS[i[2]]) \
        for numi,i in enumerate(legend_atoms)]
 atom_labels=[ label(pos=vector(numi*0.25-.1,-.1,0),\
                    text=i[0],box=False) \
        for numi,i in enumerate(legend_atoms)]

def choose_color(atoms,all_atoms,balls,moving_atoms,arrows,scene,COLORS):
 color_buttons=[]
 names=[]
 scene.append_to_caption('\nChoose color of atoms :\n')
 for i in COLORS:
   if i.x!=0 and i.y==0 and i.z==0: name='red'+str(i.x)
   elif i.x!=0 and i.x==i.y and i.z==0: name='yellow'+str(i.y)
   elif i==vector(0,0,0): name='black'
   elif i.x==0 and i.y!=0 and i.z==0: name='green'+str(i.y)
   elif i==vector(1,0.6,0): name='orange'
   elif i==vector(1,1,1): name='white'
   elif i.x==0 and i.y==0 and i.z!=0: name='blue'+str(i.z) 
   elif i==vector(0,1,1): name='cyan' 	 
   elif i==vector(0.4,0.2,0.6): name='purple'
   elif i==vector(1,0,1): name='magenta'
   elif i.x==i.y and i.y==i.z: name='gray'+str(i.x)
   else: name=str(i.x)+' '+str(i.y)+' '+str(i.z)
   names.append(name)
 def F(b):
  m=int(b.text)
  for numi,i in enumerate(all_atoms):
   if i[2]==m: 
    balls[numi].color=COLORS[b.index]
    moving_atoms[1][numi].color=COLORS[b.index]
    arrows[numi].color=COLORS[b.index]
 for k in range(len(atoms)):
  color_buttons.append(menu( bind=F, text=str(k), height=100,\
                choices=names, selected=names[k]))


def make_tetrahedrons(atoms_balls,COLORS,scene,tetrahedrons,maxbonding):
 scene.select()
 at_pos=[ at.pos for at in atoms_balls]
 tetra=[]
 for i in tetrahedrons:
  for j in i:  
   j.visible=False
   del j
 for numi,i in enumerate(atoms_balls):
  dist=[]
  mini=[]
  for numj,j in enumerate(atoms_balls[numi+1:]):
   dist.append([mag(i.pos-j.pos),i,j])
  if len(dist)==0: continue
  min_dist=min([j[0] for j in dist if j[0]>1e-1])
  if min_dist>maxbonding: continue
  for numj, j in enumerate(dist):
   if abs(j[0]-min_dist)<1e-1: mini.append(j)
  for numj, j in enumerate(mini):
   for numk, k in enumerate(mini[numj+1:]):
    if abs(mag(j[2].pos-k[2].pos)-min_dist)<1e-1: 
     for numl, l in enumerate(mini[numk+1:]):
      if abs(mag(l[2].pos-k[2].pos)-min_dist)<1e-1 \
         and abs(mag(l[2].pos-j[2].pos)-min_dist)<1e-1: 
       tetra.append([j[1],j[2],k[2],l[2]])
 print(len(tetra))
 for i in tetra:
  if len(i)<4: continue
  ats=[vertex(pos=j.pos,opacity=.5,color=j.color,canvas=scene) for j in i[:4]]
  tetrahedrons.append([triangle(v0=ats[0],v1=ats[1],v2=ats[2]),\
                       triangle(v0=ats[0],v1=ats[2],v2=ats[3]),\
                       triangle(v0=ats[1],v1=ats[2],v2=ats[3]),curve(canvas=scene)])
  for k in ats: tetrahedrons[-1][-1].append(k.pos)
  tetrahedrons[-1][-1].append(ats[0].pos)


def tetrahedrons_menu(atoms_balls,COLORS,scene,tetrahedrons):
 maxbonding=10.0
 make_tetrahedrons(atoms_balls,COLORS,scene,tetrahedrons,maxbonding)
 def G(b): 
  if b.text=='all on':
   for i in tetrahedrons: 
     for j in i: j.visible=True 
  elif b.text=='all off':
   for i in tetrahedrons: 
     for j in i: j.visible=False
  else:
   m=int(b.text)-1
   if tetrahedrons[m][0].visible==True:
    for i in tetrahedrons[m]: i.visible=False
   else:
    for i in tetrahedrons[m]: i.visible=True

 def F(b):
  maxbonding=float(b.text)  
  make_tetrahedrons(atoms_balls,COLORS,scene,tetrahedrons,maxbonding)

 scene.append_to_caption('\nChoose maxbonding and tetrahedron :')
 winput(text=maxbonding, prompt='Type maxbonding',bind=F,pos=scene.caption_anchor)
 tetra_buttons=[\
     button( bind=G , text='all on',height=100,pos=scene.caption_anchor),\
     button( bind=G , text='all off',height=100,pos=scene.caption_anchor)]
 for k in range(len(tetrahedrons)):
  tetra_buttons.append(button( bind=G , text=str(k+1),height=100,pos=scene.caption_anchor))

def if_display_tetrahedrons(atoms_balls,COLORS,scene,tetrahedrons):
 def F(b):
  tetrahedrons_menu(atoms_balls,COLORS,scene,tetrahedrons)
 button(text='tetrahedrons menu', bind=F)
 
def add_plane(alat):
 try:
  h=open('plane.in')
  tmp=h.readlines()
  h.close()
  line=tmp[0]
  filename=tmp[1].split()[0]
  if '1 1 0' in line: #przekatna
   rt = shapes.rectangle(width=alat*(2**0.5), height=alat)
   v=vector(.5*alat,.5*alat,.5*alat)
   extrusion(shape=rt, path=[v,v+vector(0.1,0.,0.1)],texture=filename)
  elif '1 1 1' in line: #[111]
   rt = shapes.triangle(length=alat*(2**0.5),rotate=pi/2-pi/6)
   v=vector(.33*alat,2/3.*alat,.33*alat)
   extrusion(shape=rt, path=[v,v+vector(0.1,-0.1,0.1)], texture=filename)
  elif '1 0 0 center' in line: #sciana w srodku
   rt = shapes.rectangle(width=alat, height=alat)
   v=vector(.5*alat,.5*alat,.5*alat)
   extrusion(shape=rt, path=[v,v+vector(0.01,0.,0.)], texture='tot_pot.png')
  elif '1 0 0 face' in line: #sciana
   rt = shapes.rectangle(width=alat, height=alat)
   v=vector(0,.5*alat,.5*alat)
   extrusion(shape=rt, path=[v,v+vector(0.01,0.,0.)], texture='tot_pot.png')
 except: 1

def rotate_one_obj(obj,ang,ax,orig): 
 obj.rotate(angle=ang,  axis=ax,  origin=orig)


def rotate(disp,ang,ax,orig):
 for v in vars(disp): 
  try: 
   rotate_one_obj(vars(disp)[v],ang,ax,orig)
  except:
   if type(vars(disp)[v]) is not list: continue 
   for v1 in vars(disp)[v]:
    try:
     rotate_one_obj(v1,ang,ax,orig)
    except:
     if type(v1) is not list: continue
     for v2 in v1:
      try:
       rotate_one_obj(v2,ang,ax,orig)
      except: 
       if type(v2) is not list: continue
       for v3 in v2:
        try:
         rotate_one_obj(v3,ang,ax,orig)       
        except: 1

def plot_dispersion(Q,FREQ,chosen_q,no_of_modes):
 Q2=[[float(m) for m in i.split()[2:]] for i in Q]
 for i in range(len(Q2)):
   if i==0: dist=0
   else: 
    dq=[Q2[i][m]-Q2[i-1][m] for m in range(3)]
    dist=dist+(dq[0]**2+dq[1]**2+dq[2]**2)**0.5
   if chosen_q==Q[i]: 
    f2=gcurve(color=color.red)
    f2.plot([dist,0],[dist, FREQ[i][-1]])
 plot_modes=[]
 for j in range(no_of_modes):
  plot_modes.append(gcurve())
  for i in range(len(Q2)):
   if i==0: dist=0
   else: 
    dq=[Q2[i][m]-Q2[i-1][m] for m in range(3)]
    dist=dist+(dq[0]**2+dq[1]**2+dq[2]**2)**0.5
   plot_modes[-1].plot(dist,  FREQ[i][j])

'''
def move_to_center(disp,crystal):
 o=0.5*(crystal[0]+crystal[1]+crystal[2])
 for v in vars(disp): 
  try: 
   rotate_one_obj(vars(disp)[v],ang,ax)
  except:
   if type(vars(disp)[v]) is not list: continue 
   for v1 in vars(disp)[v]:
    try:
     v1.pos=v1.pos+o
    except:
     if type(v1) is not list: continue
     for v2 in v1:
      try:
       v2.pos=v2.pos+o
      except: 
       if type(v2) is not list: continue
       for v3 in v2:
        try:
         v3.pos=v3.pos+o
        except: 1
'''
