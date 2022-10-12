# If you use NNScore in your research, please cite the following reference:
# NNScore: A Neural-Network-Based Scoring Function for the Characterization of
# Protein-Ligand Complexes. Jacob D. Durrant, J. Andrew McCammon. Journal of
# Chemical Information and Modeling, 2010, 50 (10), pp 1865-1871.

import __future__

import sys
import math
import time
import os
import random

######################################## Variables to Modify if AutoDock Output Files are to be Used ##############################################################################
mglenv = '' # If you need to use source to set environmental variables, put the filename of the file to source here.
prepare_ligand4_location = '/net/linux/pkg/autodocktools-1.5.1/MGLTools-1.5.0/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py' # the location of prepare_ligand4.py
###################################################################################################################################################################################

program_info = "NNScore is based in part on ffnet, coded by Marek Wojciechowski. \nIt is distributed under the GNU General Public License, version 3. \nSee gpl-3.0.txt for more details."
version = '1.1'

'''Basic Feed Forward Neural Network functionality.'''
class FFNet:
    '''* Network weights '''
    weights = []
    '''* Connection array '''
    conec = [] # [][]
    '''* List of input nodes '''
    inno = []
    '''* List of output nodes '''
    outno = []
    '''* Input normalization parameters '''
    eni = [] # [][]
    '''* Output denormalization parameters '''
    deo = [] # [][]
    '''* Network size '''
    unitCount = 0

    def findMax2d(self, table):
        max = -sys.maxsize - 1
        for row in table:
            for value in row:
                if value > max:
                    max = value
        return max


    '''Sigmoid transfer function'''
    def sigmoid(self, x):
        return 1.0 / (1.0 +  math.exp(-x));

    def file_to_double_array(self, filename):
        weights = []
        conec = []
        inno = []
        outno = []
        eni = []
        deo = []
        text_file = open(filename, "r")
        lines = text_file.readlines()
        for line in lines:
            if line[:8] == "WEIGHTS:": weights.append(float(line[8:]))
            if line[:6] == "CONEC:":
                pair = line[6:]
                pair = pair.split(" ")
                pair[0] = int(pair[0])
                pair[1] = int(pair[1])
                conec.append(pair)
            if line[:5] == "INNO:": inno.append(int(line[5:]))
            if line[:6] == "OUTNO:": outno.append(int(line[6:]))
            if line[:4] == "ENI:":
                pair = line[4:]
                pair = pair.split(" ")
                pair[0] = float(pair[0])
                pair[1] = float(pair[1])
                eni.append(pair)
            if line[:4] == "DEO:":
                pair = line[4:]
                pair = pair.split(" ")
                pair[0] = float(pair[0])
                pair[1] = float(pair[1])
                deo.append(pair)

        text_file.close()

        output = [weights, conec, inno, outno, eni, deo]

        return output

    '''def file_to_int_array(self, filename):
        ints = []
        text_file = open(filename, "r")
        lines = text_file.readlines()
        for line in lines:
            ints.append(int(line))
        text_file.close()

        return ints'''

    # def __init__(self, weights, conec, inno, outno, eni, deo):
    def __init__(self, filename):

        # first, load the weights into a double array
        arrays = self.file_to_double_array(filename)

        weights = arrays[0]
        conec = arrays[1]
        inno = arrays[2]
        outno = arrays[3]
        eni = arrays[4]
        deo = arrays[5]

        #

        prev_target = conec[0][1]
        # Check that connections are properly ordered
        for connection in conec:
            target = connection[1]
            prev_target = target

        self.weights = weights
        self.conec = conec
        self.inno = inno
        self.outno = outno
        self.eni = eni
        self.deo = deo

        self.unitCount = self.findMax2d(conec) + 1

    '''*
     * Normalizes inputs and sets network status for propagation.
     *
     * @param input
     * @return units - network state before propagation
     '''
    def setInput(self, input):
        # double[] units = new double[unitCount];
        units = [0 for i in range(self.unitCount)]
        k = 0
        while k < len(self.inno):
            units[self.inno[k]] = self.eni[k][0] * input[k] + self.eni[k][1]
            k = k + 1

        return units

    '''*
     * Gets network state with input already set and calculates all activations.
     * Identity input and sigmoid activation function for other units is assumed.
     '''
    def prop(self, units):

        ''' Connections are arranged so, that inputs to one node are computed together '''
        ctrg = self.conec[0][1]
        i = 0
        while i < len(self.conec):
            src = self.conec[i][0]
            trg = self.conec[i][1]

            # If next target, apply sigmoid
            if trg != ctrg:
                units[ctrg] = self.sigmoid(units[ctrg])
                ctrg = trg
                units[ctrg] = 0

            if src == -1: # bias
                units[ctrg] = units[ctrg] + self.weights[i]
            else:
                units[ctrg] = units[ctrg] + units[src] * self.weights[i]

            i = i + 1

        units[ctrg] = self.sigmoid(units[ctrg])


    '''*
     * Gets output from network state and denormalizes it
     '''
    def getOutput(self, units):
        # double[] output = new double[outno.length];
        output = [0 for i in range(len(self.outno))]

        k = 0
        while k < len(self.outno):
            output[k] = self.deo[k][0] * units[self.outno[k]] + self.deo[k][1];
            k = k + 1

        return output;

    '''*
     * Calls the network
     *
     * @param input Input vector
     * @return Output vector
     '''
    def call(self, input):

        units = self.setInput(input);
        self.prop(units)
        return self.getOutput(units)

start_time = time.time()

def therange(start, end, step):
    arange = []
    delta = end - start
    numsteps = delta / step + 1
    for i in range(0, int(numsteps)):
        if start + i * step<= end:
            arange.append(start + i * step)

    return arange

def format_num(num):
    astr = "%.3f" % num
    if len(astr) < 7:
        astr = ' ' + astr
    return astr

def convert_atomname_to_element(atomname):
    element = atomname
    element = element.replace("0","")
    element = element.replace("1","")
    element = element.replace("2","")
    element = element.replace("3","")
    element = element.replace("4","")
    element = element.replace("5","")
    element = element.replace("6","")
    element = element.replace("7","")
    element = element.replace("8","")
    element = element.replace("9","")
    element = element[:1]
    return element

def get_vdw(element):
        ans = 1.0 # the default
        if element == "H":
                ans = 1.2
        if element == "C":
                ans = 1.7
        if element == "N":
                ans = 1.55
        if element == "O":
                ans = 1.52
        if element == "F":
                ans = 1.47
        if element == "P":
                ans = 1.8
        if element == "S":
                ans = 1.8
        return ans

class point:
    x = 99999.0
    y = 99999.0
    z = 99999.0

    def __init__ (self, x, y ,z):
        self.x = x
        self.y = y
        self.z = z

    def print_coors(self):
        print(str(self.x)+"\t"+str(self.y)+"\t"+str(self.z))

    def snap(self,reso): # snap the point to a grid
        self.x = round(self.x / reso) * reso
        self.y = round(self.y / reso) * reso
        self.z = round(self.z / reso) * reso

    def dist_to(self,apoint):
        return math.sqrt(math.pow(self.x - apoint.x,2) + math.pow(self.y - apoint.y,2) + math.pow(self.z - apoint.z,2))

    def description(self):
        return str(self.x) + " " + str(self.y) + " " + str(self.z)

class region:

    def __init__(self):
        self.center = point(9999.9, 9999.9, 9999.9)
        self.radius = 9999.9 # in case the region is a sphere
        self.box_dimen = point(9999.9, 9999.9, 9999.9) # in case the region is a box

        self.region_type = "SPHERE" # could also be BOX

    def volume(self):
        vol = 0.0
        if self.region_type == "SPHERE":
            vol = (4.0 / 3.0) * math.pi * math.pow(self.radius,3)
        elif self.region_type == "BOX":
            vol = self.box_dimen.x * self.box_dimen.y * self.box_dimen.z
        return vol

    def point_in_region(self, point, padding):
        response = False
        if self.region_type == "SPHERE":
            dist = math.sqrt(math.pow(point.x-self.center.x,2) + math.pow(point.y-self.center.y,2) + math.pow(point.z-self.center.z,2))
            if dist<= self.radius + padding:
                response = True
        elif self.region_type == "BOX":
            x_min = self.center.x - self.box_dimen.x/2 - padding
            x_max = self.center.x + self.box_dimen.x/2 + padding
            y_min = self.center.y - self.box_dimen.y/2 - padding
            y_max = self.center.y + self.box_dimen.y/2 + padding
            z_min = self.center.z - self.box_dimen.z/2 - padding
            z_max = self.center.z + self.box_dimen.z/2 + padding

            if point.x>= x_min and point.x<= x_max and point.y>= y_min and point.y<= y_max and point.z>= z_min and point.z<= z_max:
                response = True

        return response

    def print_out(self):
        self.center.print_coors()
        print(self.radius)
        self.box_dimen.print_coors()
        print(self.region_type)

    def points_set(self, reso):
        points = []
        if self.region_type == "BOX":
            # self.center = point(9999.9, 9999.9, 9999.9)
            # self.box_dimen
            for x in therange(self.center.x - self.box_dimen.x / 2,self.center.x + self.box_dimen.x / 2, reso):
                for y in therange(self.center.y - self.box_dimen.y / 2, self.center.y + self.box_dimen.y / 2, reso):
                    for z in therange(self.center.z - self.box_dimen.z / 2,self.center.z + self.box_dimen.z / 2, reso):
                        thepoint = point(x,y,z)
                        thepoint.snap(reso)
                        # thepoint.print_coors()
                        points.append(thepoint)
        elif self.region_type == "SPHERE":
            for x in therange(self.center.x - self.radius,self.center.x + self.radius, reso):
                for y in therange(self.center.y - self.radius, self.center.y + self.radius, reso):
                    for z in therange(self.center.z - self.radius,self.center.z + self.radius, reso):
                        thepoint = point(x,y,z)
                        thepoint.snap(reso)
                        if self.point_in_region(thepoint,0): # padding always 0 for inclusion objects
                            points.append(thepoint)
        return points


class atom:

    def __init__ (self):
        self.atomname = ""
        self.residue = ""
        self.coordinates = point(99999, 99999, 99999)
        self.undo_coordinates = point(99999, 99999, 99999)
        self.element = ""
        self.PDBIndex = ""
        self.line = ""
        self.atomtype = ""
        self.IndeciesOfAtomsConnecting = []
        self.charge = 0

    def ReadPDBLine(self, Line):
        self.line = Line
        self.atomname = Line[11:16].strip()

        if len(self.atomname) == 1:
            self.atomname = self.atomname + "  "
        elif len(self.atomname) == 2:
            self.atomname = self.atomname + " "
        elif len(self.atomname) == 3:
            self.atomname = self.atomname + " " # This line is necessary to work, though many PDBs in the PDB would have this line commented out

        self.coordinates = point(float(Line[30:38]), float(Line[38:46]), float(Line[46:54]))

        # now atom type (for pdbqt)
        self.atomtype = Line[77:79].strip().upper()

        self.charge = float(Line[69:76])

        # if len(Line)>= 79:
        # self.element = Line[76:79].strip().upper() # element specified explicitly at end of life
        # el
        if self.element == "": # try to guess at element from name
            two_letters = self.atomname[0:2].strip().upper()
            if two_letters == 'BR':
                self.element = 'BR'
            elif two_letters == 'CL':
                self.element = 'CL'
            elif two_letters == 'BI':
                self.element = 'BI'
            elif two_letters == 'AS':
                self.element = 'AS'
            elif two_letters == 'AG':
                self.element = 'AG'
            elif two_letters == 'LI':
                self.element = 'LI'
            # elif two_letters == 'HG':
            # self.element = 'HG'
            elif two_letters == 'MG':
                self.element = 'MG'
            elif two_letters == 'MN':
                self.element = 'MN'
            elif two_letters == 'RH':
                self.element = 'RH'
            elif two_letters == 'ZN':
                self.element = 'ZN'
            elif two_letters == 'FE':
                self.element = 'FE'
            else: # So, just assume it's the first letter.
                # Any number needs to be removed from the element name
                self.element = self.atomname
                self.element = self.element.replace('0','')
                self.element = self.element.replace('1','')
                self.element = self.element.replace('2','')
                self.element = self.element.replace('3','')
                self.element = self.element.replace('4','')
                self.element = self.element.replace('5','')
                self.element = self.element.replace('6','')
                self.element = self.element.replace('7','')
                self.element = self.element.replace('8','')
                self.element = self.element.replace('9','')
                self.element = self.element.replace('@','')

                self.element = self.element[0:1].strip().upper()


        self.PDBIndex = Line[6:12].strip()
        self.residue = Line[16:20]
        if self.residue.strip() == "": self.residue = " MOL"

    def CreatePDBLine(self):

        # if len(self.atomname) > 1: self.atomname = self.atomname[:1].upper() + self.atomname[1:].lower()

        output = "ATOM "
        # output = output + str(index).rjust(6) + self.atomname.rjust(5) + self.residue.rjust(4)
        output = output + self.PDBIndex.rjust(6) + self.atomname.rjust(5) + self.residue.rjust(4)
        output = output + ("%.3f" % self.coordinates.x).rjust(18)
        output = output + ("%.3f" % self.coordinates.y).rjust(8)
        output = output + ("%.3f" % self.coordinates.z).rjust(8)
        output = output + self.element.rjust(24) # + "   " + str(uniqueID) # This last part must be removed
        return output


class PDB:

    def __init__ (self):
        self.AllAtoms = {}

    def LoadPDB(self, FileName):

        self.__init__()

        # Now load the file into a list
        file = open(FileName,"r")
        lines = file.readlines()
        file.close()

        self.LoadPDB_from_list(lines)


    def print_out_info(self):
        for index in self.AllAtoms:
            print(self.AllAtoms[index].CreatePDBLine())

    def LoadPDB_from_list(self, lines):

        autoindex = 1
        self.entropy_count = 0

        '''for line in file.readlines():
            if "between atoms" in line and " A " in line:
                entropy_count = entropy_count + 1
        file.close()'''


        for t in range(0, len(lines)):
            line = lines[t]
            if len(line)>= 7:

                if "between atoms" in line and " A " in line:
                    self.entropy_count = self.entropy_count + 1

                if line[0:5] == "ATOM " or line[0:7] == "HETATM ": # Load atom data (coordinates, etc.)
                    TempAtom = atom()
                    TempAtom.ReadPDBLine(line)

                    self.AllAtoms[autoindex] = TempAtom # because points files have no indecies
                    autoindex = autoindex + 1

class Complex:
    def not_trained(self,type, key, ligand_atom,receptor_atom, ligand_filename):
            toreturn = "The program can't deal with the " + type + " calculation between ligand,receptor atom types: " + key + ", ligand = " + ligand_filename + "\n"
            toreturn = toreturn + "LIGAND:   " + ligand_atom.line.strip() + "\n"
            toreturn = toreturn + "RECEPTOR: " + receptor_atom.line.strip()
            return toreturn
            # sys.exit(0)

    def make_key(self,string1, string2):
        # compare 2 strings cmp(string1, string2)
        #   if string_1 < string_2 then return -1
        #   elif string_1 == string_2   return zero
        #   elif string_1 > string_2    retrun +1
        #  Example "OB" > "OA" > "NA" > "MA"
        # This is to replace the cmp() function which was removed in python 3
        # replace cmp() with (a > b) - (a < b) as an equivalent substitution.


        compare_value = (string1 > string2) - (string1 < string2)

        if compare_value == -1:
            key = string1 + "_" + string2
        else:
            key = string2 + "_" + string1
        return key


    def __init__ (self, ligand, receptor): # ligand_name, receptor_name):

        # self.IC50 = IC50 * math.pow(10,15) # so IC fifty is given in M, converted to fM.

        '''if IC50<= math.pow(10,-8): # less than 10 nM
            self.IC50 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        if IC50<= math.pow(10,-7): # less than 100 nM
            self.IC50 = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
        elif IC50<= math.pow(10,-6): # less than 1 uM
            self.IC50 = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
        elif IC50<= math.pow(10,-5): # less than 10 uM
            self.IC50 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        elif IC50<= math.pow(10,-4): # less than 100 uM
            self.IC50 = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        else: # > 100 uM, doesn't bind
            self.IC50 = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]'''

        # self.IC50_orig = IC50

        # if IC50<= 2.5 * math.pow(10,-5): # less than 25 uM
        # self.IC50 = [1.0, 0.0]
        # else:
        # self.IC50 = [0.0, 1.0]

        # picklename = "./pickles/"+os.path.basename(ligand_name)+"_"+os.path.basename(receptor_name)+".pickle"


        # else:

        # First, load in the ligand
        # ligand = PDB()
        # ligand.LoadPDB(ligand_name)

        # receptor = PDB()
        # receptor.LoadPDB(receptor_name)

        # get info
        # all atom types, too big
        # charge_type_combos = ["FE_N", "A_BR", "MG_NA", "CL_SA", "A_A", "A_C", "A_CL", "A_F", "A_FE", "A_HD", "A_I", "A_N", "A_NA", "A_OA", "A_P", "A_S", "A_SA", "A_ZN", "BR_C", "BR_HD", "BR_N", "BR_OA", "C_C", "C_CL", "C_F", "C_FE", "C_HD", "C_I", "CL_HD", "CL_N", "CL_OA", "C_MG", "C_MN", "C_N", "C_NA", "C_OA", "C_P", "C_S", "C_SA", "C_ZN", "FE_HD", "FE_OA", "F_HD", "F_N", "F_OA", "F_SA", "HD_HD", "HD_I", "HD_MG", "HD_MN", "HD_N", "HD_NA", "HD_OA", "HD_P", "HD_S", "HD_SA", "HD_ZN", "I_N", "I_OA", "MG_OA", "MG_P", "MN_N", "MN_OA", "MN_P", "NA_OA", "NA_SA", "NA_ZN", "N_N", "N_NA", "N_OA", "N_P", "N_S", "N_SA", "N_ZN", "OA_OA", "OA_P", "OA_S", "OA_SA", "OA_ZN", "P_ZN", "SA_SA", "SA_ZN", "S_ZN"]
        charge_type_combos = ["A_A", "A_BR", "A_C", "A_CL", "A_F", "A_FE", "A_HD", "A_I", "A_N", "A_NA", "A_OA", "A_P", "A_S", "A_SA", "A_ZN", "BR_C", "BR_HD", "BR_N", "BR_OA", "C_C", "C_CL", "C_F", "C_FE", "C_HD", "C_I", "CL_HD", "CL_N", "CL_OA", "CL_SA", "C_MG", "C_MN", "C_N", "C_NA", "C_OA", "C_P", "C_S", "C_SA", "C_ZN", "FE_HD", "FE_N", "FE_OA", "F_HD", "F_N", "F_OA", "F_SA", "HD_HD", "HD_I", "HD_MG", "HD_MN", "HD_N", "HD_NA", "HD_OA", "HD_P", "HD_S", "HD_SA", "HD_ZN", "I_N", "I_OA", "MG_NA", "MG_OA", "MG_P", "MN_N", "MN_OA", "MN_P", "NA_OA", "NA_SA", "NA_ZN", "N_N", "N_NA", "N_OA", "N_P", "N_S", "N_SA", "N_ZN", "OA_OA", "OA_P", "OA_S", "OA_SA", "OA_ZN", "P_ZN", "SA_SA", "SA_ZN", "S_ZN"]
        # proximity_2_type_combos = ["FE_HD", "A_HD", "C_HD", "C_OA", "C_SA", "HD_HD", "HD_MG", "HD_N", "HD_NA", "HD_OA", "HD_ZN", "MG_OA", "NA_ZN", "OA_ZN"]
        proximity_2_type_combos = ["A_HD", "C_HD", "C_OA", "C_SA", "FE_HD", "HD_HD", "HD_MG", "HD_N", "HD_NA", "HD_OA", "HD_ZN", "MG_OA", "NA_ZN", "OA_ZN"]
        # proximity_4_type_combos = ["FE_N", "A_BR", "MG_NA", "CL_SA", "A_A", "A_C", "A_CL", "A_F", "A_FE", "A_HD", "A_I", "A_N", "A_NA", "A_OA", "A_P", "A_S", "A_SA", "A_ZN", "BR_C", "BR_HD", "BR_N", "BR_OA", "C_C", "C_CL", "C_F", "C_FE", "C_HD", "C_I", "CL_HD", "CL_N", "CL_OA", "C_MG", "C_MN", "C_N", "C_NA", "C_OA", "C_P", "C_S", "C_SA", "C_ZN", "FE_HD", "FE_OA", "F_HD", "F_N", "F_OA", "F_SA", "HD_HD", "HD_I", "HD_MG", "HD_MN", "HD_N", "HD_NA", "HD_OA", "HD_P", "HD_S", "HD_SA", "HD_ZN", "I_N", "I_OA", "MG_OA", "MG_P", "MN_N", "MN_OA", "MN_P", "NA_OA", "NA_SA", "NA_ZN", "N_N", "N_NA", "N_OA", "N_P", "N_S", "N_SA", "N_ZN", "OA_OA", "OA_P", "OA_S", "OA_SA", "OA_ZN", "P_ZN", "SA_SA", "SA_ZN", "S_ZN"]
        proximity_4_type_combos = ["A_A", "A_BR", "A_C", "A_CL", "A_F", "A_FE", "A_HD", "A_I", "A_N", "A_NA", "A_OA", "A_P", "A_S", "A_SA", "A_ZN", "BR_C", "BR_HD", "BR_N", "BR_OA", "C_C", "C_CL", "C_F", "C_FE", "C_HD", "C_I", "CL_HD", "CL_N", "CL_OA", "CL_SA", "C_MG", "C_MN", "C_N", "C_NA", "C_OA", "C_P", "C_S", "C_SA", "C_ZN", "FE_HD", "FE_N", "FE_OA", "F_HD", "F_N", "F_OA", "F_SA", "HD_HD", "HD_I", "HD_MG", "HD_MN", "HD_N", "HD_NA", "HD_OA", "HD_P", "HD_S", "HD_SA", "HD_ZN", "I_N", "I_OA", "MG_NA", "MG_OA", "MG_P", "MN_N", "MN_OA", "MN_P", "NA_OA", "NA_SA", "NA_ZN", "N_N", "N_NA", "N_OA", "N_P", "N_S", "N_SA", "N_ZN", "OA_OA", "OA_P", "OA_S", "OA_SA", "OA_ZN", "P_ZN", "SA_SA", "SA_ZN", "S_ZN"]
        lig_types_combos = ["A", "BR", "C", "CL", "F", "HD", "I", "N", "NA", "OA", "P", "S", "SA"]

        # autodock_atom_types = ["HD", "HS", "C", "A", "N", "NA", "OA", "F", "P", "SA", "S", "CL", "BR", "I"]
        # autodock_atom_types = ["H", "C", "N", "O", "P", "S", "CL"]

        # Autogenerate keys
        # ligand_atom_types = ["H", "C", "N", "O", "F", "P", "S", "CL", "BR", "I"]
        # receptors_atom_types = ["C", "FE", "H", "MG", "MN", "N", "O", "P", "S", "ZN"]
        # ligand_atom_types = ["H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS", "F", "MG", "P", "SA", "S", "CL", "CA", "MN", "FE", "ZN", "BR", "I"]
        # receptors_atom_types = ["H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS", "F", "MG", "P", "SA", "S", "CL", "CA", "MN", "FE", "ZN", "BR", "I"]
        # type_combos = []
        # for lig_type in ligand_atom_types:
        # for recep_type in receptors_atom_types:
        # type_combos.append(self.make_key(lig_type, recep_type))
        # charge_type_combos = type_combos[:]
        # proximity_2_type_combos = type_combos[:]
        # proximity_4_type_combos = type_combos[:]

        # Enhancement: sort alphabetically. You actually have double the needed here.
        # charge_type_combos = ["BR_C", "BR_H", "BR_N", "BR_O", "C_C", "C_CL", "C_F", "C_FE", "C_H", "C_I", "CL_H", "CL_N", "CL_O", "C_MG", "C_MN", "C_N", "C_O", "C_P", "C_S", "C_ZN", "FE_H", "FE_O", "FE_S", "F_H", "F_N", "F_O", "F_S", "F_ZN", "H_H", "H_I", "H_MG", "H_MN", "H_N", "H_O", "H_P", "H_S", "H_ZN", "I_N", "I_O", "MG_N", "MG_O", "MG_P", "MN_N", "MN_O", "MN_P", "N_N", "N_O", "N_P", "N_S", "N_ZN", "O_O", "O_P", "O_S", "O_ZN", "P_S", "P_ZN", "S_S", "S_ZN"]
        # proximity_2_type_combos = ["C_H", "C_N", "C_O", "C_S", "F_H", "H_H", "H_MG", "H_N", "H_O", "H_ZN", "MG_O", "MN_O", "N_ZN", "O_O", "O_ZN"]
        # proximity_4_type_combos = ["BR_C", "BR_H", "BR_N", "BR_O", "C_C", "C_CL", "C_F", "C_FE", "C_H", "C_I", "CL_H", "CL_N", "CL_O", "C_MG", "C_MN", "C_N", "C_O", "C_P", "C_S", "C_ZN", "FE_H", "FE_O", "FE_S", "F_H", "F_N", "F_O", "F_S", "F_ZN", "H_H", "H_I", "H_MG", "H_MN", "H_N", "H_O", "H_P", "H_S", "H_ZN", "I_N", "I_O", "MG_N", "MG_O", "MG_P", "MN_N", "MN_O", "MN_P", "N_N", "N_O", "N_P", "N_S", "N_ZN", "O_O", "O_P", "O_S", "O_ZN", "P_S", "P_ZN", "S_S", "S_ZN"]

        self.coulomb_energy = {} # where to store energy info
        self.proximity_2 = {} # where to store proximity info (within 2 A)
        self.proximity_4 = {} # where to store proximity info (within 3.5 A)
        self.lig_types = {}

        for key in charge_type_combos:
            self.coulomb_energy[key] = 0
        for key in proximity_2_type_combos:
            self.proximity_2[key] = 0
        for key in proximity_4_type_combos:
            self.proximity_4[key] = 0
        for key in lig_types_combos:
            self.lig_types[key] = 0

        # Here you could remove receptor atoms that are far from the ligand... not currently implemented, but would speed things up a bit.
        self.bad_training = ""

        for ligand_index in ligand.AllAtoms:
            ligand_atom = ligand.AllAtoms[ligand_index]

            # also keep track of the different ligand types
            lig_type = ligand_atom.atomtype

            if lig_type not in self.lig_types:
                self.bad_training = self.bad_training + "\n\n" + "The program can't deal with ligands that have the atom type: " + lig_type
            else:
                self.lig_types[lig_type] = self.lig_types[lig_type] + 1

            for receptor_index in receptor.AllAtoms:
                receptor_atom = receptor.AllAtoms[receptor_index]
                dist = ligand_atom.coordinates.dist_to(receptor_atom.coordinates)
                if dist < 0.5: print("There may be steric clashes between " +ligand_name + ", " + receptor_name)

                if dist < 4.0:

                    # get the key
                    # lig_type = ligand_atom.element
                    # recep_type = receptor_atom.element
                    recep_type = receptor_atom.atomtype

                    key = self.make_key(lig_type, recep_type)
                    key_ordered = lig_type + " " + recep_type

                    if key in charge_type_combos:
                        # calculate charges
                        lig_charge = ligand_atom.charge
                        recep_charge = receptor_atom.charge
                        coulomb = lig_charge * recep_charge / dist # ignore all constants. Just let the neural net take care of that
                        self.coulomb_energy[key] = self.coulomb_energy[key] + coulomb

                    else:
                        self.bad_training = self.bad_training + "\n\n" + self.not_trained("CHARGE", key_ordered, ligand_atom,receptor_atom,ligand_name)

                    if dist < 2.0: # so it's a close contact
                        if key in proximity_2_type_combos:
                            self.proximity_2[key] = self.proximity_2[key] + 1

                        else:
                            self.bad_training = self.bad_training + "\n\n" + self.not_trained("PROXIMITY_2", key_ordered, ligand_atom,receptor_atom,ligand_name)
                    else: # so it's between 2.0 and 4.0
                        if key in proximity_4_type_combos:
                            self.proximity_4[key] = self.proximity_4[key] + 1

                        else:
                            self.bad_training = self.bad_training + "\n\n" + self.not_trained("PROXIMITY_4", key_ordered, ligand_atom,receptor_atom, ligand_name)

        # now we need to find out how many active torsions there are in the ligand (so the network can account for some entropy effects).

        entropy_count = ligand.entropy_count

        '''file = open(ligand_name, "r")
        entropy_count = 0
        for line in file.readlines():
            if "between atoms" in line and " A " in line:
                entropy_count = entropy_count + 1
        file.close()'''

        self.nn_input = []
        for key in charge_type_combos:
            self.nn_input.append(self.coulomb_energy[key])
        for key in proximity_2_type_combos:
            self.nn_input.append(self.proximity_2[key])
        for key in proximity_4_type_combos:
            self.nn_input.append(self.proximity_4[key])
        for key in lig_types_combos:
            self.nn_input.append(self.lig_types[key])

        self.nn_input.append(entropy_count)

networks = []
receptor_name = ""
ligand_name = ""
vina_output = ""
autodock_output = ""

for index in range(1,len(sys.argv)):
    var = sys.argv[index].strip()
    if var.upper() == "-RECEPTOR": receptor_name = sys.argv[index+1]
    if var.upper() == "-LIGAND": ligand_name = sys.argv[index+1]
    if var.upper() == "-VINA_OUTPUT": vina_output = sys.argv[index+1]
    if var.upper() == "-AUTODOCK_OUTPUT": autodock_output = sys.argv[index+1]  ###### ADD TO USER MANUAL #####
    if var.upper() == "-NETWORK": networks.append(sys.argv[index+1])
    if var.upper() == "-NETWORKS_DIR": # a directory containing only networks
        sys.argv[index+1] = sys.argv[index+1].replace("\\","/")
        if sys.argv[index+1][-1:] != "/": sys.argv[index+1] = sys.argv[index+1] + "/"
        for filename in os.listdir(sys.argv[index+1]):
            if os.path.isfile(sys.argv[index+1] + filename): networks.append(sys.argv[index+1] + filename)

print("")
print("NNScore " + version)
print("")
print(program_info)
print("")
print("If you use NNScore in your research, please cite the following reference:")
print("  NNScore: A Neural-Network-Based Scoring Function for the Characterization")
print("  of Protein-Ligand Complexes. Jacob D. Durrant, J. Andrew McCammon. Journal")
print("  of Chemical Information and Modeling, 2010, 50 (10), pp 1865-1871.")
print("")

error = False
if receptor_name == "":
    error = True
    extra_message = "Error! Required parameters were not passed to NNScore!"
if len(networks) == 0:
    error = True
    extra_message = "Error! Required parameters were not passed to NNScore!"
if ligand_name == "" and vina_output == "" and autodock_output == "":
    error = True
    extra_message = "Error! Required parameters were not passed to NNScore!"
if ligand_name != "" and vina_output != "":
    error = True
    extra_message = "Error! You cannot use both -ligand and -vina_output!"
if ligand_name != "" and autodock_output != "":
    error = True
    extra_message = "Error! You cannot use both -ligand and -autodock_output!"
if autodock_output != "" and vina_output != "":
    error = True
    extra_message = "Error! You cannot use both -autodock_output and -vina_output!"

if error is True:
    print(extra_message)
    print("\nParameters are:")
    print("\t-receptor <pdbqt filename>")
    print("\t-ligand <pdbqt filename>")
    print("\t-vina_output <vina output filename>")
    print("\t-autodock_output <autodock output filename>")
    print("\t-network <network filename>")
    print("\t-networks_dir <directory>")
    print("")
    print("Note: It is best to use multiple neural networks to judge ligand binding by")
    print("consensus. Commandline parameters can be used to add neural-network files")
    print("to the list of those that will be used. To add a single neural network to")
    print("the list, use the -network parameter to specify a single network file. To")
    print("add mutliple networks to the list, create a directory containing only")
    print("network files and specify the path to that directory using the -networks_dir")
    print("parameter.")
    print("")
    print("Note: Only pdbqt files of the receptor and ligand are accepted. Scripts to")
    print("convert from pdb to pdbqt are included in the AutoDockTools package:")
    print("http://autodock.scripps.edu/resources/adt")
    print("")
    print("Note: If an AutoDock Vina or AutoDock output file is specified, NNScore will")
    print("evaluate all docked poses and return the best NNScore calculated. To score an")
    print("AutoDock output file, modify the mglenv and prepare_ligand4_location variables")
    print("at the begining of this python file.")
    print("")
    print("Examples:")
    print("\t python NNScore.py -receptor neuraminidase.pdbqt -ligand oseltamivir.pdbqt -network ./networks/top_3_networks/12.net")
    print("\t python NNScore.py -receptor integrase.pdbqt -ligand raltegravir.pdbqt -networks_dir ./networks/top_3_networks/")
    print("\t python NNScore.py -receptor integrase.pdbqt -vina_output docked_poses.vina.out -networks_dir ./networks/top_3_networks/")
    print("\t python NNScore.py -receptor integrase.pdbqt -autodock_output docked_poses.autodock.out -networks_dir ./networks/top_3_networks/")
    print("\t python NNScore.py -receptor protease.pdbqt -ligand tipranavir.pdbqt -networks_dir ./networks/top_24_networks/ -network ./networks/top_3_networks/16.net")
    print("")
    sys.exit()


print("Receptor: "+ receptor_name)
print("")

def process_ligand(ligand_name, ligand, receptor):
    global networks
    print("Ligand: "+ ligand_name)
    print(" = " * len("Ligand: "+ ligand_name))
    if len(networks) == 1:
            print("\tNetwork: ", networks[0])
    else:
            print("\tNetworks: ")
            for net in networks:
                    print("\t\t", net)

    print("")

    acomplex = Complex(ligand, receptor)
    scores = []

    # describe binding
    print("\tAtom types (one ligand, one receptor) within 2 angstroms of each other:")
    first = 1
    something = 1
    for item in acomplex.proximity_2:
            value = acomplex.proximity_2[item]
            if value != 0:
                    something = 0
                    addin = ""
                    if first != 1: addin = "; "
                    first = 0

                    print(addin + "(" + str(item.replace("_",', ')) + "), " + str(value))
                    if value == 1:
                            sys.stdout.write(" time")
                    else:
                            sys.stdout.write(" times")
    if something == 1:
        print("(None)")

    print("")
    print("")
    print("\tAtom types (one ligand, one receptor) within 4 angstroms of each other:")
    first = 1
    something = 1
    for item in acomplex.proximity_4:
            value = acomplex.proximity_4[item]
            if value != 0:
                    something = 0
                    addin = ""
                    if first != 1: addin = "; "
                    first = 0

                    print(addin + "(" + str(item.replace("_",', ')) + "), " + str(value))
                    if value == 1:
                            sys.stdout.write(" time")
                    else:
                            sys.stdout.write(" times")
    if something == 1: print("(None)")
    print("")
    print("")
    print("\tRelative coulombic energy between atom types (one ligand, one receptor) within 4 angtroms of each other:")
    first = 1
    something = 1
    for item in acomplex.coulomb_energy:
            value = round(acomplex.coulomb_energy[item],3)
            if value != 0:
                    something = 0
                    addin = ""
                    if first != 1: addin = "; "
                    first = 0

                    print(addin + "(" + str(item.replace("_",', ')) + "), " + str(value))
                    if value == 1:
                            sys.stdout.write(" unit")
                    else:
                            sys.stdout.write(" units")
    if something == 1: print("(None)")
    print("")
    print("")
    print("\tAtom types in the ligand:")
    first = 1
    something = 1
    for item in acomplex.lig_types:
            value = acomplex.lig_types[item]
            if value != 0:
                    something = 0
                    addin = ""
                    if first != 1: addin = "; "
                    first = 0

                    print(addin + str(item.replace("_",', ')) + ", " + str(value))
                    if value == 1:
                            sys.stdout.write(" time")
                    else:
                            sys.stdout.write(" times")
    if something == 1: print("(None)")
    print("")
    print("")



    # self.coulomb_energy = {} # where to store energy info
    # self.proximity_2 = {} # where to store proximity info (within 2 A)
    # self.proximity_4 = {} # where to store proximity info (within 3.5 A)
    # self.lig_types = {}


    for net in networks:
            print("\tUsing network " + net + " to predict binding: ")
            net_name = net
            net = FFNet(net_name)
            result = net.call(acomplex.nn_input)
            score = result[0] - result[1]
            print("\t", score)
            scores.append(score)
            if score < 0:
                    print("(bad binder)")
            else:
                    print("(good binder)")

    # compute average score
    total = 0
    for score in scores: total = total + score
    average_score = total/len(scores)

    if acomplex.bad_training != "":
        print(acomplex.bad_training)
        average_score = -999999.9

    return average_score

if ligand_name != "": # so a single pdbqt ligand was provided
    ligand = PDB()
    ligand.LoadPDB(ligand_name)

    receptor = PDB()
    receptor.LoadPDB(receptor_name)

    average_score = process_ligand(ligand_name, ligand, receptor)

    print("")
    print("Average score: ", average_score)

    if average_score < 0:
            print("(bad binder)")
    else:
            print("(good binder)")

    print("")

elif vina_output != "": # so a vina output file has been passed
    receptor = PDB()
    receptor.LoadPDB(receptor_name)

    file = open(vina_output,"r")
    lines = file.readlines()
    file.close()

    thelines = []
    label = ""

    best_binder = -10000000.0
    best_binder_name = ""

    for line in lines:
        if "MODEL " in line:
            if len(thelines) != 0:

                # so create a complex of this frame and the receptor, get the score
                ligand = PDB()
                ligand.LoadPDB_from_list(thelines)
                average_score = process_ligand(vina_output + ", " +label, ligand, receptor)
                if best_binder < average_score:
                    best_binder = average_score
                    best_binder_name = vina_output + ", " +label
                print("\n\tAverage score: ", average_score)
                if average_score < 0: print("(bad binder)")
                else: print("(good binder)")
                print("")
                thelines = [line.strip()]

            else: thelines = [line.strip()]

            label = line.strip()

        else: thelines.append(line.strip())

    # so create a complex of this frame and the receptor, get the score
    ligand = PDB()
    ligand.LoadPDB_from_list(thelines)
    average_score = process_ligand(vina_output + ", " +label, ligand, receptor)
    if best_binder < average_score:
        best_binder = average_score
        best_binder_name = vina_output + ", " +label
    print("\n\tAverage score: ", average_score)
    if average_score < 0: print("(bad binder)")
    else: print("(good binder)")
    print("")
    print("Best score:", best_binder,"(" + best_binder_name + ")")
    print("")
elif autodock_output != "":
    receptor = PDB()
    receptor.LoadPDB(receptor_name)

    file = open(autodock_output,"r")
    lines = file.readlines()
    file.close()

    for t in range(len(lines)):
        if "Keeping original residue number (specified in the input PDBQ file) for outputting." in lines[t]: break

    pdbs = []
    thelines = []
    for i in range(t+2,len(lines)):
        line = lines[i]
        if " l " in line:
            print("Atom name \"l\" being interpreted as \"Cl.\"")
            line = line.replace(' l ','Cl ')
        if " r " in line:
            print("Atom name \"r\" being interpreted as \"Br.\"")
            line = line.replace(' r ','Br ')
    if " A " in line:
        print("Atom name \"A\" being interpreted as \"C\" (i.e. an aromatic carbon).")
        line = line.replace(' A ',' C ')
        if "ENDMDL" in line:
            thelines.append(line)
            pdbs.append(thelines)
            thelines = []
        else: thelines.append(line)

    # now, printOut PBD files in temp directory
    if not os.path.exists(sys.path[0] + os.sep + "tmp"):
        os.mkdir(sys.path[0] + os.sep + "tmp")
        print("Created directory " + sys.path[0] + os.sep + 'tmp' + os.sep)

    # now, pick a random id number
    id = random.randrange(0, 1000000)

    filenames = []

    for pdb in pdbs:
        name = pdb[0]
        name = name.replace("\t"," ")
        while "  " in name: name = name.replace("  "," ")
        name = name.replace(" ", "_").strip()
        filename = sys.path[0] + os.sep + 'tmp' + os.sep + name +"."+ str(id) + ".pdb"
        filenames.append(filename)

        f = open(filename,'w')
        f.writelines(pdb)
        f.close()

    # now convert the pdbs into pdbqts
    runstring = ''
    if mglenv.strip() != "": runstring = 'source ' + mglenv + '; '
    runstring = runstring + prepare_ligand4_location
    print('Converting frames from AutoDock output file into individual pdbqt files...')
    for filename in filenames:
        torun = runstring + " -l " + filename + " -o " + filename + "qt"
        print("\t"+torun)
        os.system(torun)

    best_binder = -10000000.0
    best_binder_name = ""

    for filename in filenames:
        filename = filename+"qt"

        ligand = PDB()
        ligand.LoadPDB(filename)

        tmp = os.path.basename(filename)
        tmp = tmp.split(".")
        tmp = tmp[0]
        ligand_name = autodock_output + ", " + tmp.replace("_",' ')

        average_score = process_ligand(ligand_name, ligand, receptor)
        if best_binder < average_score:
            best_binder = average_score
            best_binder_name = ligand_name
        print("\n\tAverage score: ", average_score)
        if average_score < 0: print("(bad binder)")
        else: print("(good binder)")
        print("")
    print("Best score:", best_binder,"(" + best_binder_name + ")")
    print("")

    # now delete files
    for filename in filenames:
        os.remove(filename)
        os.remove(filename + "qt")
