#!/usr/bin/python3
import numpy as np 
import math 
import os
import itertools
import sys

#Classes definitions

#Class Atom collect objects atoms with attributes as in Crystal output file 
class Atom():
    
    def __init__(self,atom_number,atom_status,atom_Z,atom_label,atom_coord):
        self.atom_number=atom_number
        self.atom_status=atom_status
        self.atom_label=atom_label
        self.atom_Z=atom_Z
        self.atom_coord=np.asarray(atom_coord)
        
#Class CrystalFiles includes methods to handle crystal output and input files
class CrystalFiles():

    #Read data from crystal output file, both transformation matrix and geometries at each optimization step are read
    def read_crys_out(self,ifile):
        with open(ifile, 'r') as f:
            # Read line

            self.geom_list = []

            for line in f:

                # Find and read slab structural parameters
                if " LATTICE PARAMETERS (ANGSTROMS AND DEGREES) - BOHR = 0.5291772083 ANGSTROM" in line:

                    self.slab_parameters = []
                    new_line = f.readline()
                    new_line = f.readline()
                    new_line = f.readline()
                    new_line = new_line.split(" ");

                    while ("" in new_line):
                        new_line.remove("")

                    new_line = [float(x) for x in new_line]

                    #Definition of transformation matrix from data read
                    self.slab_parameters = np.array([[new_line[0], new_line[1] * math.cos(math.radians(new_line[5])), 0],[0, new_line[1] * math.sin(math.radians(new_line[5])), 0], [0, 0, 1]])

                #Find and read geometry block as expressed in fractional coordinates
                elif "ATOMS IN THE ASYMMETRIC UNIT" in line:
                    n_atoms = ([int(s) for s in line.split() if s.isdigit()])

                    #n_atoms_primitive is the number of atoms with label T while n_atoms_withsymm include even atoms labelled as F

                    n_atoms_primitive = n_atoms[0]
                    self.n_atoms_withsymm = n_atoms[1]

                    new_line = f.readline()
                    new_line = f.readline()

                    atom_list = []

                    #Read attributes for each atom and an object atom is inizialized

                    for i in range(self.n_atoms_withsymm):

                        new_line = f.readline()
                        new_line = new_line.split(" ");

                        # Remove blanks from the line
                        while ("" in new_line):
                            new_line.remove("")

                        atom_list.append(Atom(new_line[0], new_line[1], new_line[2], new_line[3],
                                              atom_coord=[float(new_line[4]), float(new_line[5]), float(new_line[6])]))



                    self.geom_list.append(atom_list)
                  #  print(self.slab_parameters)

            return self.geom_list, self.slab_parameters

    #This method subtract atom coordinates to calculate displacements at the end of the optimization
    def subtract_geoms(self):

        self.displacement = []  # difference between first and last geom

        for i in range(self.n_atoms_withsymm):

          #  print('Starting geometry:', self.geom_list[0][i].atom_coord)
          #  print('Last geometry:', self.geom_list[-1][i].atom_coord)

            self.displacement.append(Atom(self.geom_list[0][i].atom_number, self.geom_list[0][i].atom_status, self.geom_list[0][i].atom_Z,
                     self.geom_list[0][i].atom_label, atom_coord=(np.subtract(self.geom_list[-1][i].atom_coord, self.geom_list[0][i].atom_coord))))

            #Here start the procedure to check if some atom has been substituted by the corresponding image and displacements need to be rescaled
            chk_x_incr = 1-abs(float(self.geom_list[0][i].atom_coord[0]))
            chk_x_decr = 1-abs(float(self.geom_list[0][i].atom_coord[0]))
            chk_y_incr = 1-abs(float(self.geom_list[0][i].atom_coord[1]))
            chk_y_decr = 1-abs(float(self.geom_list[0][i].atom_coord[1]))
            tmp_x = float(self.displacement[i].atom_coord[0])
            tmp_y = float(self.displacement[i].atom_coord[1])
         #   print(self.geom_list[0][i].atom_coord)
#            print(self.geom_list[0][0].atom_coord,self.geom_list[-1][0].atom_coord)

            if (tmp_x >= 0) and (tmp_x >= abs(chk_x_incr)):
 #               print('Changing X coordinate for:', i)

                self.displacement[i].atom_coord[0]= -self.geom_list[-1][i].atom_coord[0]-self.geom_list[0][i].atom_coord[0]

            elif (tmp_x <= 0) and (abs(tmp_x) > chk_x_decr):
  #              print('Changing X coordinate for:', i)

                self.displacement[i].atom_coord[0] = 1 - abs(self.geom_list[-1][i].atom_coord[0]) - self.geom_list[0][i].atom_coord[0]

            if (tmp_y >= 0) and (tmp_y >= abs(chk_y_incr)):
   #             print('Changing Y coordinate for:', i)

                self.displacement[i].atom_coord[1]= -1+self.geom_list[-1][i].atom_coord[1]-self.geom_list[0][i].atom_coord[1]

            elif (tmp_y <= 0) and ((abs(float(self.displacement[i].atom_coord[1])))>abs(chk_y_decr)):
    #            print('Changing Y coordinate for:', i)

                self.displacement[i].atom_coord[1] = 1 + tmp_y

            #Conversion of displacements from fractional to cartesian
            self.displacement[i].atom_coord=np.dot(self.slab_parameters,self.displacement[i].atom_coord)
         #   print(self.displacement[i].atom_coord)


     #       print(self.geom_list[0][0].atom_coord, self.geom_list[-1][0].atom_coord)

        return self.displacement

    #This method write a new Crystal input including section of ATOMDISP and GEOMTEST
    #A file New_Cry is generated and can be copied and pasted in your new input
    def write_input(self, ofile):

        with open(ofile, 'a') as ofile:
            ofile.write('ATOMDISP\n')
            ofile.write(str(self.n_atoms_withsymm) + '\n')
            for i in range(self.n_atoms_withsymm):
                ofile.write(
                    str(self.displacement[i].atom_number) + " " + str(self.displacement[i].atom_coord)[1:-1] + "\n")
            ofile.write('TESTGEOM\n' + 'END\n' + '\n' + '\n')

#Crystal output file to process taken from command line
CRYout_file=sys.argv[-1]

class_instance= CrystalFiles()
class_instance.read_crys_out(CRYout_file)
class_instance.subtract_geoms()
class_instance.write_input('New_Cry.d12')








