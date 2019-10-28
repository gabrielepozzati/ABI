from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from codes import nnlab as nl
import pickle
import networkx as nx
import multiprocessing as mp
import pandas as pd


allowed_atoms = ['N','CA','C','O','OXT','CB',
                 'CZ','CZ1','CZ2','CE','CE1','CE2',
                 'CD','CD1','CD2','CG','CG1','CG2','CG3',
                 'CH2','OE','OE1','OE2','OD','OD1','OD2',
                 'OH','OG','OG1','OG2','NZ','NE','NE1','NE2',
                 'NH','NH1','NH2','ND','ND1','ND2',
                 'SG','SD','SE']

def chainfind(code, hotspots):

    try: structure = parser.get_structure(code, 'hotspot_pdb/'+code+'.pdb')
    except:
        print (code, ' no structure!')
        return 'none'
    model_chains=[]
    for chain in structure[0]: model_chains.append(chain.get_id())

    main_chains=[]
    for key in hotspots:
        if hotspots[key][1] not in main_chains and hotspots[key][1] in model_chains: 
            main_chains.append(hotspots[key][1])

    if main_chains == []: 
        residue_list = []
        for chain in structure[0]:
            for residue in chain:
                if not is_aa(residue) or residue.get_id()[0] != ' ': continue
                residue_list.append(residue.get_resname()+str(residue.get_id()[0]))
            for key in hotspots: 
                if key in residue_list: main_chains.append(chain.get_id())
    
    if main_chains == []:
        print (code, ' no hotspot chain!')
        return 'none'

    chain_couples = {}
    for chain_id in main_chains:
        for hot_residue in structure[0][chain_id]:
            hot_residue_name = hot_residue.get_resname()+str(hot_residue.get_id()[1])
            if hot_residue_name in hotspots:
                for chain in structure[0]:
                    atoms = []
                    if chain.get_id() == chain_id: continue
                    for residue in chain:
                        if not is_aa(residue) or residue.get_id()[0] != ' ': continue 
                        for atom in residue: atoms.append(atom)

                    if atoms == []: continue
                    if NeighborSearch(atoms).search(hot_residue['CA'].get_coord(), 12, level='R') != []:
                        chain_couples[chain_id+'_'+chain.get_id()] = chain_couples.get(chain_id+'_'+chain.get_id(), {})
                        chain_couples[chain_id+'_'+chain.get_id()][hot_residue_name] = hotspots[hot_residue_name]
                        del hotspots[hot_residue_name]
                    break

    if len(chain_couples) == 0: 
        print (code, ' no pair!')
        return 'none'

    return chain_couples


def chaingraph(chain, thr):

    G = nx.Graph()
    for residue1 in chain:
        r1 = residue1.get_resname()+str(residue1).split('=')[2].split()[0]
        G.add_node(r1)
        atom_list=[]
        for atom in residue1: atom_list.append(atom)
        for residue2 in chain:
            r2 = residue2.get_resname()+str(residue2).split('=')[2].split()[0]
            contact = 'n'
            if residue1 == residue2: continue
            for atom in residue2: 
                if NeighborSearch(atom_list).search(atom.get_coord(), thr, level='R') != []:
                    contact = 'y'
                    break
            if contact == 'y':
                dist_r1_r2 = thr+0.1
                for atom1 in atom_list:
                    for atom2 in residue2:
                        dist_a1_a2 = nl.distance_find(atom1.get_coord(), atom2.get_coord())
                        dist_r1_r2 = min(dist_r1_r2, dist_a1_a2)
                G.add_edge(r1, r2, weight = dist_r1_r2)
    return G
            

def worker(args, thr=6):

    code = args[0].rstrip()
    c1 = args[1]
    c2 = args[2]

    G = nx.Graph()
    broken = 'n'

    try: structure = parser.get_structure(code, 'hotspot_pdb/'+code+'.pdb')
    except: broken = 'y'
    chain1 = structure[0][c1]
    chain2 = structure[0][c2]
    print (code, chain1, chain2)

    residues = []
    chain_hotspots = {}
    if chain1 == chain2 and model1 == model2: broken = 'y'
    if chain1 == '' or chain2 == '': broken = 'y'


    if broken == 'n':
        interface1 = []
        interface2 = []
        interface_contacts = {}
        for residue1 in chain1:
            atom_list1 = []
            if not is_aa(residue1): continue
            for atom in residue1: atom_list1.append(atom)

            for residue2 in chain2:

                distance = []
                if not is_aa(residue2): continue
                for atom in residue2: 
                    distance += NeighborSearch(atom_list1).search(atom.get_coord(), thr, level='A')

                if distance:
                    interface1.append(residue1)
                    interface2.append(residue2) 
                    dist_r1_r2 = thr+0.1
                    for atom1 in distance: 
                        for atom2 in residue2: 
                            dist_a1_a2 = nl.distance_find(atom1.get_coord(), atom2.get_coord())
                            dist_r1_r2 = min(dist_r1_r2, dist_a1_a2)
                            if dist_r1_r2 == dist_a1_a2: midpoint = nl.midpoint_find(atom1.get_coord(), atom2.get_coord()) 

                    r1 = residue1.get_resname()+str(residue1).split('=')[2].split()[0]
                    r2 = residue2.get_resname()+str(residue2).split('=')[2].split()[0]
                    interface_contacts[r1+'_'+r2] = [dist_r1_r2, midpoint]
                    G.add_node(r1+'_'+r2)
    
        for key1 in interface_contacts:
            k1c1 = key1.split('_')[0]
            k1c2 = key1.split('_')[1]
    
            for key2 in interface_contacts:
                k2c1 = key2.split('_')[0]
                k2c2 = key2.split('_')[1]
    
                if key1 == key2: continue
                elif k1c1 == k2c1: 
                    dk1_k2 = nl.distance_find(interface_contacts[key1][1], interface_contacts[key2][1])
                    if interface_contacts[key1][0] <= interface_contacts[key2][0]: G.add_edge(key1, key2, weight = dk1_k2)
                    else: G.add_edge(key2, key1, weight = dk1_k2)
                elif k1c2 == k2c2: 
                    dk1_k2 = nl.distance_find(interface_contacts[key1][1], interface_contacts[key2][1])
                    if interface_contacts[key1][0] <= interface_contacts[key2][0]: G.add_edge(key1, key2, weight = dk1_k2)
                    else: G.add_edge(key2, key1, weight = dk1_k2)

        G1 = chaingraph(interface1, thr)
        G2 = chaingraph(interface2, thr)

        return [code+'_'+c1+'_'+c2, G, G1, G2, args[3]]

    else: 
        print (code.rstrip()+' failed!')
        return 'none'

if __name__ == '__main__':
    parser = PDBParser()
    joblist = []

    df = pd.read_csv('all.csv')
    hot_list = []
    for row in range(df.shape[0]):
        pdb = str(df['pdb'][row]).lower()
        chain = str(df['chain'][row])
        try: delta = float(df['ddG'][row])
        except: continue
        mutation = df['mutant'][row]
        if mutation[-1] == 'A': mutation = mutation[:-1]
        try: mutation = nl.one2three[mutation[0]]+mutation[1:] 
        except: mutation = nl.one2three[mutation[-1]]+mutation[:-1]
        if [pdb, chain, mutation, delta] not in hot_list: 
            hot_list.append([pdb, chain, mutation, delta])

    fail_list = []
    fail_count = 0
    for code in open('pdb_list','r'):
        count = 0
        code_spots = {}
        for row in range(len(hot_list)):
            if code.rstrip() == hot_list[row][0]: code_spots[hot_list[row][2]] = [hot_list[row][3], hot_list[row][1]]
        print (code.rstrip()+'...')
        chain_couples = chainfind(code.rstrip(), code_spots)
        if chain_couples == 'none':
            fail_list.append(code)
            fail_count += 1
            continue

        for couple in chain_couples: 
            joblist.append([code.rstrip(), couple.split('_')[0], couple.split('_')[1], chain_couples[couple]])
        print (code, chain_couples)

    pool = mp.Pool(processes=7)
    results = pool.map(worker, joblist)

    for r in results:
        if r == 'none': continue
        with open('graphs/'+r[0], 'wb') as f:
            print (r[0], r[4])
            pickle.dump({'1_2':r[1], '1':r[2], '2':r[3], 'H':r[4]},  f)
    with open('fail_list','w') as out:
        for code in fail_list: out.write(code)
    print ('failed: ', fail_count)
