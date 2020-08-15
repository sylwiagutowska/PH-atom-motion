Here you can find files produced by PHonon - part of Quantum Espresso.
The include the informations about vibrations of atoms.
If you want to visualize these vibrations, just 
1) type: python3 ../as_balls_and_arrows/run.py
2) type any of number (every number represent different wave vector of phonons) 
3) you will see vpython window with two windows, which represents primitive and conventional crystal cell.
   At the bottom of each window you have to choose one of numbers, which represent the number of phonon mode 
   (there are 3*N modes, N=number of atoms in unit cell. First 3 modes are always acoustic, rest are the optical ones)
   You have various of options: you can remove bonding of atoms, switch on moving atoms
