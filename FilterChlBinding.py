import Bio.PDB
import sys
import numpy as np
import pandas as pd

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def get_plane_angle(CLA, pivot, atom):
   v1 = -CLA['MG'].get_coord() + CLA['NC'].get_coord()
   v2 = -CLA['MG'].get_coord() + CLA['NB'].get_coord()
   v3 = -np.cross(v1,v2)
   v4 = atom.get_coord() - pivot.get_coord()
   angle = angle_between(v4, v3)
   return 90 - np.degrees(angle)
   #if angle<90: return 90-angle
   #if angle>90: return angle-90

def predict_chl(fname, pivot, radius, angleT=90, DN=None, AN=None, exclude_alpha_helix=True):

    """This script identify the chl binding sites that have hydrogen bond donors within
       a certain (radius) around the (pivot) atom. DN is the list of atoms that are considered
       hydrogen bond donor and AN is the list of hydrogen bond acceptor. There is an option to
       exclude the the backbone nitrogen in alpha helcies (exclude_alpha_helix) because they 
       make hydrogen with the oxygens. With the angleT option you may exclude hydrogen
       bonds that are on the top of the chl plan, which is defined by the Mg and the N atoms"""  

    N = ['NB','NA','NC','ND']
    if DN is None:
        D = ['OMB','N','NE','NH1','NH2','NE2','ND2','ND1','NZ','OG','OG1','NE1','OH','SG']
    else:
        D = DN
    if AN is None:
        A = ['O1A','O2A','O1D', 'OBD','OD1','OD2','OE1','OE2','NZ','OG','OG1','NA','NB','NC','ND','O']
    else:
        A = AN
    chl = ['CLA','CL7','F6C','CHL']
    parser = Bio.PDB.PDBParser(QUIET=True) # QUIET=True avoids comments on errors in the pdb.
    structures = parser.get_structure('1prot', fname)
    structure = structures[0] # 'structures' may contain several proteins in this case only one.
    target = []
    for chain in structure:
        for res in chain:
            if res.get_resname() in chl:
                target.append((res,res[pivot]))
          
    atoms  = Bio.PDB.Selection.unfold_entities(structure, 'A')
    ns = Bio.PDB.NeighborSearch(atoms)
    df = pd.DataFrame(columns=['CLA','Num','HBcountD','HBcountA','subtract','Pred','DistD','DistA',
                               'Dnames','Dresname', 'Anames','Aresname'])
    for c,i in enumerate(target):
        close_res = ns.search(i[1].coord, radius, level='A')
        countA = 0
        countD = 0
        dA = 0
        dD = 0
        Dnames =[]
        Dresname=[]
        Anames =[]
        Aresname=[]
        for atom in close_res:
            if atom.element=="O" or atom.element=="N" or atom.element=="S":
                if atom.get_parent()!=i[0] and atom.get_name() not in N:
                    if atom.get_parent().get_resname()=="HOH":
                        countD+=1
                        dD+=(atom-i[1])            
                        Dnames.append(atom)
                        Dresname.append(atom.get_parent().get_resname())
                    elif atom.get_parent().get_resname()!="HIS" or atom.get_name()=="N":
                        angle = get_plane_angle(i[0], i[1], atom)
                        if abs(angle) < angleT or True:
                            if atom.get_name() in D:
                            #print(i[0], atom.get_name(), get_plane_angle(i[0], i[1], atom))
                                countD+=1
                                dD+=(atom-i[1])
                                Dnames.append(atom)
                                Dresname.append(atom.get_parent().get_resname())
               	            elif atom.get_name() in A:
                                countA+=1
                                dA+=(atom-i[1])
                                Anames.append(atom)
                                Aresname.append(atom.get_parent().get_resname()) 


        if exclude_alpha_helix:
            di, ai = None, None 
            for l,datom in enumerate(Dnames):
                if datom.get_name()=="N":
                   for k,aatom in enumerate(Anames):
                      if aatom.get_name()=="O" and aatom.get_parent().get_resname()!="HOH" and datom-aatom<4.:
                          di, ai = l, k
                          dD -= (datom-i[1])
                          dA -= (aatom-i[1]) 
                          break
            if di is not None and ai is not None:
               Dnames.pop(di);Dresname.pop(di)
               Anames.pop(ai);Aresname.pop(ai)
               countD-=1; countA-=1;

        if countA>0: dA=dA/countA
        if countD>0: dD=dD/countD
        subtract = countD-countA
        if subtract>0: pred = "!CLA"
        elif subtract<0: pred = "CLA"
        else: pred = "!CLA/CLA"    
        df.loc[c] = [i[0].get_resname(),
                     i[0].get_id()[1],
                     countD, countA, subtract, pred, dD, dA, Dnames, Dresname, Anames, Aresname]
    return df

if __name__ == "__main__":
    #df = predict_chl(fname=sys.argv[1],pivot=sys.argv[2],radius=float(sys.argv[3]))
    #print(df[df["CLA"]=="F6C"])
    D = ['OMB','N','NE','NH1','NH2','NE2','ND2','ND1','NZ','OG','OG1','NE1','OH','SG','CE1']
    df = predict_chl(fname="FRL-PSII_monomer_10-1-2022_FINAL.pdb",pivot="CMB",radius=4., DN=D)
    print(df)
    #df.to_csv(sys.argv[1].split(".")[0]+"_HB.csv")
