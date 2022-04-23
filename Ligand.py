import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Draw
import pysmiles


#m = Chem.MolFromSmiles('N1C=Cc2ccccc12')
#m.GetSubstructMatches(Chem.MolFromSmarts('c'))
#m
""" This function takes input as list of smiles and preprocesses it further"""
def input_smile(lig, smiles=Chem.MolFromSmiles(pd_1)):
    print("Enter the smile comma separated ")

    x = input()
    pc = []
    pc.append(x)
    print(pc)
    for s in pc:

        sub_list = s.split(",")

        print(sub_list)


        for pd in sub_list:

            pd= Chem.MolFromSmiles(pd)
            print(pd)
            show_implicit_h(pd)
            minimization(pd)
#s1 = 'OCCn2c(=N)n(CCOc1ccc(Cl)cc1Cl)c3ccccc23'  # aromatic
#s2 = 'OCCN2C(=N)N(CCOC1=CC=C(Cl)C=C1Cl)C3=CC=CC=C23'  # kekulized

""" This function implicitly adds hydrogen atoms """
def show_implicit_h(smiles):
    m = Chem.MolFromSmiles(smiles)
    for atom in m.GetAtoms():
        atom.SetProp('atomLabel', str(atom.GetIdx()))
    m = Chem.AddHs(m)

    return (m)

#show_implicit_h(s1)


#input_smile('N1C=Cc2ccccc12')
""" This function minimizes the small molecule using a forcefield providing mdp and sh files to run gromacs"""
def minimization(smiles):
    for ligand in pd:
                        path_save = ligand()
                        #ck=chain.get_id()

                        #print(ck)
                        fullpath= path_save+'\\'
                        print(fullpath)
                        #save_c=chain.get_id()
                    #c_p=path_sav+save_c
                        #print(save_c)
                        name_of_f= (path_save+"_")
                        c_name= os.path.join(fullpath,name_of_f+".mdp")

                        md=open(c_name,"w")

                        md.write("integrator  = steep \n emtol       = 1000.0 \n emstep      = 0.01\n nsteps      = 50000 \n nstxout = 10\n"
                         "nstlist =1 \n cutoff-scheme   = Verlet \n ns_type         = grid \n coulombtype     = PME \n rcoulomb   = 1.0 \n"
                         "rvdw  = 1.0 \n pbc = xyz")

                        md.close()
                        print(md)
                        min_md=os.path.join(fullpath,name_of_f+".sh")
                        min_m= open(min_md, "w")
                        #id_chain=chain.get_id()
                        md_p=ligand+"_"
                        print(md_p)
                        min_m.write("gmx pdb2gmx -f " +ligand+".pdb -o protein.gro -ignh -ff oplsaa -water spce\n "
                        "gmx editconf -f protein.gro -o box.gro -bt cubic -d 1.0 \n"
                        "gmx grompp -f "+md_p+".mdp -c solv.gro -p topol.top -o ions.tpr\n"
                        "gmx genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral\n"
                        "gmx grompp -f minim.mdp -c ions.gro -p topol.top -o em.tpr -maxwarn 1\n"
                        "gmx mdrun -v -deffnm em\n"
                        "gmx editconf -f em.gro -o protein_minimized.pdb")
                        min_m.close()
                        print(min_m)

    return True
