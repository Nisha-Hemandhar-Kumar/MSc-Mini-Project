#!/usr/bin/env python3

"""A code to calulate the atom-atom distane, average B-factor and to extract the ligand information from the pdb file"""

#argparse module parse the pdb file and extract a particular function from the command line
import argparse
# module to perform math functions
import math
# module to perform statistics
import statistics

class Pdb():

 """Pdb class can deal with pdb files to perform three functions get_Distane(), get_Bfactor(), and get_Ligand()"""

 def __init__(self):

  """Created a constructor for the class Pdb which assigns a file object for the pdb file which is inputed from the terminal"""

  self.pdb_file = args.filename

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
 def get_Distance(self): 

  """Returns all the distance between the C and CA atoms of a pdb file"""
  
  #readlines() method to the pdb file to return a list containing lines
  all_lines = self.pdb_file.readlines()

  # empty list to hold the X, Y, Z coordinates, atom number, residue name, chain and residue number of the Ca atom
  atom_1 = list()

  # empty list to hold the X, Y, Z coordinates, atom number, residue name, chain and residue number of the C atom
  atom_2 = list()
  
  # to get first atom input from the user
  atom1 = input("Input the first atom: ")

  # to get second atom input from the user
  atom2 = input("Input the second atom: ")

  # for loop to iterate over the lines in the pdb file
  for i in all_lines:

    #split() method to return a list of words splitted by space
    j = i.split()
     
    #condition to specify the row starting with ATOM and the row where the alpha carbon atom is present
    if j[0]=='ATOM'and j[2]== atom1:

        #appends the atom_number, residue name, chain id, residue number,X,Y,Z coordinates of Ca to an empty list
        atom_1.append((j[1],j[3],j[4],j[5],float(j[6]), float(j[7]),float(j[8])))

    #condition to specify the row starting with ATOM and the row where the carbon atom is present
    if j[0]=='ATOM'and j[2]== atom2:

        #appends the atom_number, residue name, chain id, residue number,X,Y,Z coordinates of C atom to an empty list
        atom_2.append((j[1],j[3],j[4],j[5],float(j[6]), float(j[7]),float(j[8])))

  #calculation of the euclidean distance between the Ca and C atom 
  #
  #zip function is mainly used to combining data of two iterable elements ca_atom and c_atom 

  dist_calc = [( a[0],a[1],a[2],a[3],(a[4]-b[4])**2, (a[5] - b[5])**2, (a[6] - b[6])**2) for a, b in zip(atom_1, atom_2) ]
  ca_c_dist = [(atm_no,res,chain,res_num,round(math.sqrt( x+y+z),3)) for atm_no, res,chain,res_num,x,y,z in dist_calc ]

  #output file to write the calculated distances between Ca and C atom
  with open(atom1+'_'+atom2+".txt","w") as dist_file: 

   # writing the headings to the distance.txt file 
   dist_file.write("Distance between all the "+atom1+ " and "+atom2+" are: "+"\n"+"Residue"+"\t "+"Chain"+"\t"+"Residue_num"+"\t"+"Distance"+"\n")  

    #writing the calculated distance to a file name distance.txt 
   for atm_no, res, chain,res_num,distance  in ca_c_dist:

    #center() method will center align the string, with the specified length
    dist_file.write(res.center(9)+chain.center(8)+ res_num.center(8)+(str(distance)).center(12)+"\n" )
  
  print("The atom-atom distance calculations are written to "+atom1+"_"+atom2+".txt file")  
 
  return 

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
 def get_Bfactor(self):

   """Returns the average overall B-factor, the backbone average B-factors and average alpha-carbon B-factors"""

   # readlines() method to the pdb file to return a list containing lines
   all_lines = self.pdb_file.readlines() 

   # empty list to hold the bfactor values of all atoms 
   avg_bfactor = list() 

   # empty list to hold the bfactor values of main chain atoms
   mainchain_b_factor = list()

   # empty list to hold the bfactor values of Ca atom
   ca_bfactor = list()

   # for loop to iterate over the lines in the pdb file
   for i in all_lines:

    #split() method to return a list of words splitted by space 
    j=i.split()

    #condition to specify the row where there is atom 
    if j[0]=='ATOM':

        # extracts the bfactor values of all atoms
        b_fac = (float(j[10])) 

        # appends the bfactor values of all atoms to an empty list
        avg_bfactor.append(b_fac) 

    #condition to specify the row which starts with ATOM and rows where there are backbone atoms
    if j[0]=='ATOM'and (j[2]=='CA')|(j[2]=='C')|(j[2]=='N')|(j[2]=='O'):

        # extracts the bfactor values of main chain atoms
        b_factor_mainchain = (float(j[10]))

        # appends the bfactor values of main chain atoms to an empty list
        mainchain_b_factor.append(b_factor_mainchain)

    #condition to specify the row starting with ATOM and the row where the alpha carbon atom is present
    if j[0]=='ATOM' and (j[2]=='CA'):

        # extracts the bfactor values of Ca atom
        b_fac_ca = (float(j[10])) 

        # appends the bfactor values of Ca atom to an empty list
        ca_bfactor.append(b_fac_ca) 

   # printing the calculated B-factor values to the terminal 
   print (" Average overall B-factor is: "+ str(round(statistics.mean(avg_bfactor),3)))
   print (" Backbone averaged B-factor is: "+str(round(statistics.mean(mainchain_b_factor),3)))
   print (" Average alpha carbon B-factor is: "+str(round(statistics.mean(ca_bfactor),3)))

   return
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
 def get_Sequence(self):

  """Returns a file with the fasta sequence for the chain id specified by the user as an input"""

  # readlines() method to the pdb file to return a list containing lines
  all_lines = self.pdb_file.readlines() 

  # dictionary with key as three code of amino acid and values as a single letter code
  amino_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',}
  # empty list to hold SEQRES section of the pdb file
  seqres = list()

  # creating an input statement from the user to obtain the chain id 
  chain = input('Input the chain id for which for the seqeunce :')
  
  # for loop to iterate over the lines in the pdb file
  for i in all_lines:

        #condition to specify the row which starts with SEQRES     
        if i.startswith('SEQRES'):
        
            #split() method to return a list of words splitted by space    
            lines = i.split()

            #condition to specify rows for the input chain ID
            if lines[2]== chain : 
 
              # appends the SEQRES for the input chain to an empty list  
              seqres.append(lines[4:])
              
             
             
  #converting seqres list of list to a single list
  flat_list = [residue for sublist in seqres for residue in sublist]

  #converting list to a string
  amino_3_letter = ''.join(flat_list)

  #converting the amino acid three letter to a single letter 
  amino_3_to_1 = "".join(amino_dict[amino_3_letter[x:x+3]] for x in range(0, len(amino_3_letter), 3))

  #extracting the PDB three letter code by using the input file
  PDB_ID = self.pdb_file.name[:4]
  
  #writing the seqeunce for the input chain to chain.fasta file 
  with open(chain+".fasta","w") as seq_file:
   seq_file.write(">"+PDB_ID+":"+chain+"|"+"PDBID|CHAIN|SEQUENCE"+"\n"+amino_3_to_1)

  print("The fasta sequence is written to "+chain+".fasta file")
  
  return
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
 def get_Ligand(self): 

  """Returns a file with the ligand information of the pdb file with the chain"""  
  # readlines() method to the pdb file to return a list containing lines
  all_lines = self.pdb_file.readlines()

  # creating an input statement from the user to obtain the chain id
  chain = chain = input('Input the chain id for the ligand information :')  

  # output file to write the ligand information from the pdb file
  with open(chain+".pdb","w") as lig_file:

   # for loop to iterate over the lines in the pdb file
   for i in all_lines:

    #split() method to return a list of words splitted by space
    j=i.split()

    # condition to specify the row starting with HETATM and the row where there no HOH (water) present
    if j[0]=='HETATM' and j[3] != 'HOH' and j[4]==chain:

     # center() method will center align the string, with the specified length, rjust() method to return a string right justified, ljust() method to return a string left justified    
     lig_file.write(( j[0].ljust(8)+j[1].ljust(5)+j[2].center(4)+j[3].ljust(5))+j[4].rjust(1)+j[5].rjust(4)+ str('%8.3f' % (float(j[6]))).rjust(8)+ str('%8.3f' % (float(j[7]))).rjust(8)+str('%8.3f' % (float(j[8]))).rjust(8)+str('%6.2f'%(float(j[9]))).rjust(6)+str('%6.2f'%(float(j[10]))).ljust(6)+j[11].rjust(12)+'\n')

  print("The ligand from the pdb file for input chain is written to "+chain+".pdb file") 

  return 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#parsing over the pdb file and functions in Pdb class
parser = argparse.ArgumentParser(description = "Calculating atom-atom distance, Average B-factor, extacting the fasta seqeunce of a chain and extracting the      ligand section of a pdb file")
parser.add_argument('filename',type=argparse.FileType('r'),help="input a pdb file")
parser.add_argument('--distance', help='distance calculates the atom-atom distances', action = 'store_true')
parser.add_argument('--bfactor', help='bfactor calculates the average B-factor for all atoms, main chain atoms and CA atom', action = 'store_true')
parser.add_argument('--sequence', help='sequence extracts the fasta sequence', action = 'store_true')
parser.add_argument('--ligand', help='ligand extract the ligand section of the pdb file', action = 'store_true')
args = parser.parse_args()

#conditional statements to check the input in the command line and execute the function
if args.ligand:     
    l = Pdb()
    l.get_Ligand()
if args.distance:
    l = Pdb()
    l.get_Distance()
if args.sequence:
    l = Pdb()
    l.get_Sequence()
if args.bfactor:
    l = Pdb()
    l.get_Bfactor()





