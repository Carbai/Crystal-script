#!/usr/bin/python3
import sys
import os
import re
import numpy as np

#Definition of an Atom class, each atom has attributes required for basis set in CRYSTAL
class Atom():

    def __init__(self,AtSym=None,NEle=None,NSHL=None,LAT=None,NG=None,CHG=None,exponent=None):

        self.AtSym = AtSym
        self.NEle = NEle
        self.LAT = LAT
        self.NSHL = NSHL
        self.NG = NG
        self.CHG = CHG
        self.exponent = exponent

    #Simple function to calculate number of electrons once the atomic symbol is known
    def At_Sym_to_NEle(self,line):

        #Atom Flag is a flag used to identify that an atom has been found -- atom found than AtomFlag = True
        self.AtomFlag = False
        self.AtSym_list  = ['H', 'He',
                      'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                      'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
                      'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
                      'Br',
                      'Kr',
                      'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
                      'I',
                      'Xe',
                      'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
                      'Rn',
                      'Fr', 'Ra']

        for i in range(0, len(self.AtSym_list)):
            if self.AtSym_list[i] in line:
                self.AtSym = self.AtSym_list[i]
                self.NEle = i + 1
                self.AtomFlag = True

        return self.NEle, self.AtomFlag, self.AtSym

    #Obtain electronic configuration from the total number of electrons
    def get_electronic_configuration(self,NEle):

        self.orbital_list = ['1s',
                             '2s','2p',
                             '3s','3p',
                             '4s','3d','4p',
                             '5s','4d','5p',
                             '6s','4f','5d','6p',
                             '7s','5f','6d','7p']

        max_occ_s = 2
        max_occ_p = 6
        max_occ_d = 10
        max_occ_f = 14

        self.ele_conf = []

        tot_occ = 0
        s_electrons = 0
        p_electrons = 0
        d_electrons = 0
        f_electrons = 0

        for k in range(len(self.orbital_list)):
            s_occ = 0
            p_occ = 0
            d_occ = 0
            f_occ = 0

            if 's' in self.orbital_list[k]:

                for j in range(max_occ_s):

                    if s_occ < 2 and tot_occ < NEle:

                        s_electrons += 1
                        s_occ += 1
                        tot_occ += 1

            elif 'p' in self.orbital_list[k]:

                for j in range(max_occ_p):

                    if p_occ < 6 and tot_occ < NEle:

                        p_electrons += 1
                        p_occ += 1
                        tot_occ += 1

            elif 'd' in self.orbital_list[k]:

                for j in range(max_occ_d):

                    if d_occ < 10 and tot_occ < NEle:

                        d_electrons += 1
                        d_occ += 1
                        tot_occ += 1

            elif 'f' in self.orbital_list[k]:

                for j in range(max_occ_f):

                    if f_occ < 14 and tot_occ < NEle:

                        f_electrons += 1
                        f_occ += 1
                        tot_occ += 1

        self.ele_conf.append(s_electrons)
        self.ele_conf.append(p_electrons)
        self.ele_conf.append(d_electrons)
        self.ele_conf.append(f_electrons)

        return self.ele_conf

    #Convert the shell type S, SP, P or D to the corresponding LAT value in Crystal
    def shell_to_LAT(self, shell_type):
        item = shell_type

        if item == 'S':
            self.LAT = 0
        elif item == 'SP':
            self.LAT = 1
        elif item == 'P':
            self.LAT = 2
        elif item == 'D':
            self.LAT = 3
        elif item == 'F':
            self.LAT = 4

        return self.LAT

    #Distribute electrons in each shell according to LAT values. CHG keywords is defined
    def ele_per_shell_to_CHG(self, LAT_list_per_atom, n_ele_per_shell_per_atom):

        LAT_s = 0
        LAT_sp = 0
        LAT_p = 0
        LAT_d = 0
        LAT_f = 0
        tot_LAT_list = []
        self.CHG=[]

        for i in range(len(LAT_list_per_atom)):

            if LAT_list_per_atom[i] == 0:
                LAT_s += 1
            elif LAT_list_per_atom[i] == 1:
                LAT_sp += 1
            elif LAT_list_per_atom[i] == 2:
                LAT_p += 1
            elif LAT_list_per_atom[i] == 3:
                LAT_d += 1
            elif LAT_list_per_atom[i] == 4:
                LAT_f += 1

        if LAT_s != 0:
            tot_LAT_list.append(LAT_s)
        if LAT_sp != 0:
            tot_LAT_list.append(LAT_sp)
        if LAT_p != 0:
            tot_LAT_list.append(LAT_p)
        if LAT_d != 0:
            tot_LAT_list.append(LAT_d)
        if LAT_f != 0:
            tot_LAT_list.append(LAT_f)

        for j in range(len(tot_LAT_list)):
            res = 0

            if n_ele_per_shell_per_atom[0][j] == 0:
                for orbital in range(tot_LAT_list[j]):
                    self.CHG.append(0)
                occ_LAT = 0
            elif n_ele_per_shell_per_atom[0][j] != 0 and j == 0 and LAT_sp == 0:
                occ_LAT = (int(n_ele_per_shell_per_atom[0][j] / 2))
                for orbital in range(occ_LAT):
                    self.CHG.append(2)
                res = (n_ele_per_shell_per_atom[0][j] - (occ_LAT * 2))
                tot_orb = tot_LAT_list[j]
            elif n_ele_per_shell_per_atom[0][j] != 0 and j == 0 and LAT_sp != 0:
                self.CHG.append(2)
                n_ele_sp = n_ele_per_shell_per_atom[0][0] + n_ele_per_shell_per_atom[0][1] - 2
                occ_LAT = int(n_ele_sp / 8)
                if n_ele_sp <= 8:
                    self.CHG.append(n_ele_sp)
                    occ_LAT = 1
                    res = 0
                else:
                    for orbital in range(occ_LAT):
                        self.CHG.append(8)
                    res = (n_ele_sp - (occ_LAT * 8))
                occ_LAT += 1
                tot_orb = tot_LAT_list[j] + tot_LAT_list [j + 1]
            elif n_ele_per_shell_per_atom[0][j] != 0 and j == 1 and LAT_sp == 0:
                occ_LAT = (int(n_ele_per_shell_per_atom[0][j] / 6))
                tot_orb = tot_LAT_list[j]
                if n_ele_per_shell_per_atom[0][j] <= 6:
                    self.CHG.append(n_ele_per_shell_per_atom[0][j])
                    occ_LAT = 1
                    res = 0
                else:
                    for orbital in range(occ_LAT):
                        self.CHG.append(6)
                    res = (n_ele_per_shell_per_atom[0][j] - (occ_LAT * 6))
            elif n_ele_per_shell_per_atom[0][j] != 0 and j == 2:
                tot_orb = tot_LAT_list[j]
                occ_LAT = (int(n_ele_per_shell_per_atom[0][j] / 10))
                if n_ele_per_shell_per_atom[0][j] <= 10:
                    self.CHG.append(n_ele_per_shell_per_atom[0][j])
                    occ_LAT = 1
                    res = 0
                else:
                    for orbital in range(occ_LAT):
                        self.CHG.append(10)
                    res = (n_ele_per_shell_per_atom[0][1] - (occ_LAT * 10))
            elif n_ele_per_shell_per_atom[0][j] != 0 and j == 3:
                tot_orb = tot_LAT_list[j]
                occ_LAT = (int(n_ele_per_shell_per_atom[0][j] / 14))
                if n_ele_per_shell_per_atom[0][j] <= 14:
                    self.CHG.append(n_ele_per_shell_per_atom[0][j])
                    occ_LAT = 1
                    res = 0
                else:
                    for orbital in range(occ_LAT):
                        self.CHG.append(14)
                    res = (n_ele_per_shell_per_atom[0][1] - (orbital * 14))

            if res != 0:
                self.CHG.append(res)
                occ_LAT += 1
                res_orb = tot_orb - occ_LAT
            else:
                res_orb = tot_orb - occ_LAT

                if res_orb != 0:
                    for i in range(res_orb):
                        self.CHG.append(0)

        return self.CHG

    #Read exponent for shell in numpy arrays
    def exponent_array(self,ifile,n_line_counter,NG):

        with open(ifile, 'r') as f:

            lines = f.readlines()[n_line_counter+1:n_line_counter+NG+1]
            self.exponent_arr = np.asarray(lines)

        return self.exponent_arr

#Method class for parsing gaussian basis file
class basis_set_files_handling(Atom):

    def atoms_in_gaus_basis_file2(self, ifile):

        tmp = Atom()
        self.atoms_in_basis = []
        self.AtomBlockEndFlag = False
        n_line_counter_1 = 0
        self.NSHL1 = []

        # Values in second line crystal basis file
        self.LAT1 = []
        self.NG1 = []
        self.CHG1 = []
        self.ele_per_shell = []

        self.exponent1 = []

        with open(ifile, 'r') as f:
            for line in f:
                if line.startswith('!') == False and line.startswith('\n') == False:
                    line = line.rstrip('\n')
                    line = line.split()
                    if line[0].lower().islower() is True:
                        if self.AtomBlockEndFlag == False:
                            self.NEle1, self.AtomBlockEndFlag, self.AtSym1 = tmp.At_Sym_to_NEle(line)
                            self.ele_per_shell.append(tmp.get_electronic_configuration(self.NEle1))
                        else:
                            if 'D+' or 'D-' in line[0]:
                                tmp_value = line[0].replace('D', 'E')
                            try:
                                float(tmp_value)
                                res = True
                            except:
                                res = False
                            if res == False:
                                self.LAT1.append(tmp.shell_to_LAT(line[0]))
                                self.NG1.append(line[1])
                                self.exponent1.append(tmp.exponent_array(ifile, n_line_counter_1, int(line[1])))
                    elif self.AtomBlockEndFlag == True and '****' in line:
                        self.AtomBlockEndFlag = False
                        self.NSHL1 = len(self.LAT1)
                        self.CHG1.append(tmp.ele_per_shell_to_CHG(self.LAT1, self.ele_per_shell))
                        self.atoms_in_basis.append(Atom(self.AtSym1,self.NEle1,self.NSHL1,self.LAT1,self.NG1,self.CHG1,self.exponent1))
                        self.NSHL1 = []
                        self.LAT1 = []
                        self.NG1 = []
                        self.CHG1 = []
                        self.ele_per_shell = []
                        self.exponent1 = []

                n_line_counter_1 += 1

            return self.atoms_in_basis

#Use all the values to write a CRY_BAS file with all the basis in CRYSTAL format
    def crystal_basis_set_files(self,ofile):

        with open(ofile, 'a') as f:

            ITYP = '0'
            SCAL = '1.000'
            n_atoms = len(self.atoms_in_basis)

            for i in range(n_atoms):
                f.write(str(self.atoms_in_basis[i].NEle) + ' ' + str(self.atoms_in_basis[i].NSHL) + '\n')
                for j in range(len(self.atoms_in_basis[i].LAT)):
                    f.write('{0} {1} {2} {3} {4}\n'.format(ITYP, str(self.atoms_in_basis[i].LAT[j]),
                                                          str(self.atoms_in_basis[i].NG[j]),
                                                          self.atoms_in_basis[i].CHG[0][j], SCAL))
                    for item in self.atoms_in_basis[i].exponent[j]:
                        f.write(item)
        return


GAUSSIANBASIS_file=sys.argv[-1]
classinstance=basis_set_files_handling()
classinstance.atoms_in_gaus_basis_file2(GAUSSIANBASIS_file)
classinstance.crystal_basis_set_files('CRYSTAL_BAS')
