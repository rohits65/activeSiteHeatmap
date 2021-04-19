from Bio.PDB import *
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


# Get active site residues
def activeSite(modelName, modelFile, drug, drugFile, bindingAffinity, dist):
    parser = PDBParser()
    structure = parser.get_structure(modelName, modelFile)
    otherStructure = parser.get_structure(drug, drugFile)

    proteinChain = None
    drugResidue = None



    for c in structure[0].get_chains():
        # Finds last chain in protein
        if c.id != ' ':
            print('Chain ' + c.id + ' detected.')
            proteinChain = c
    for c in otherStructure[0].get_chains():
        # Finds last element not associated with a chain
        if c.id == ' ' or c.id == '':
            print('Added element detected.')
            iterator = c.get_residues()
            for r in iterator:
                drugResidue = r

    # Get atoms out of protein chain
    proteinAtoms = []
    for a in proteinChain.get_atoms():
        proteinAtoms.append(a)
    neighbor = NeighborSearch(proteinAtoms)

    # Gets atoms from residue
    resAtoms = []
    for a in drugResidue.get_atoms():
        resAtoms.append(a)
    
    residueMap = {}
    anotherMap = {}

    for resAtom in resAtoms:
        if resAtom.name[0] != 'H':
            residueMap.update({resAtom.name: list(map(lambda x: x.resname, neighbor.search(resAtom.get_coord(), dist, level="R")))})
            anotherMap.update({resAtom.name: list(map(lambda x: x.element, neighbor.search(resAtom.get_coord(), 4.0, level="A")))})


    residueMap.update({'BindingAffinity': bindingAffinity})
    anotherMap.update({'BindingAffinity': bindingAffinity})
    
    return [residueMap, anotherMap]

def getResidueMap(mapArray, anotherArray, targets, scores, dist):
    
    for i in range(0, len(targets)):
        for j in range(1, len(scores[i])+1):
            if scores[i][j-1] != 99:
                output = activeSite(str(targets[i]), str(targets[i]) + "0.pdb", "UNK", targets[i] + "3." + str(j) + ".pdb", scores[i][j-1], dist)
                mapArray.append(output[0])
                anotherArray.append(output[1])

    print("MAP ARRAY")
    print(mapArray)
    print("ANOTHER ARRAY")
    print(anotherArray)

    aminoAcids = ['Phe','Met','Leu','Ile','Val','Pro','Tyr','Trp','Cys','Ala','Gly','Ser','Thr','His','Glu','Gln','Asp','Asn','Lys','Arg']
    polarAtoms = ['N', 'O']
    data = {'AA': aminoAcids}

    cols = ['AA']
    # data = {}
    # cols = []
    for a in mapArray[0].keys():
        print(a)
        results = [0] * 20
        for i in range(0, len(mapArray)):
            site = mapArray[i]
            arr = np.array(scores).flatten()
            # factor = site['BindingAffinity'] - np.delete(arr, np.argwhere(arr == 99)).mean()
            factor = site['BindingAffinity'] - 7.5
            # factor = 1

            
            try:
                for i in range(0,20):
                    for j in range(0, len(site[a])):
                        if site[a][j] == aminoAcids[i].upper():
                            
                            if any(atom in anotherArray[i][a][j] for atom in polarAtoms) and a[0] in polarAtoms:
                                factor = -1
                            else:
                                factor = 1
                            results[i] += factor
            except:
                continue
        data.update({a: results})
        cols.append(a)
    return data, cols

# main
if __name__ == "__main__":
    mapArray = []
    anotherArray = []
    # targets = ["thpk", "pim1", "prk1", "mapk5", "mek1", "mapk", "rac1", "tiam1", "mth1", "mtor", "pi3k", 
    #     "rho", "rac1", "prl3", "mmp2", "ck2", "hsp90", "ahr", "titin", "mapk12", "mapk2",
    #     "mapk9", "mapk10", "mapk2-4", "dsp6"]
    targets = ["thpk", "pim1", "prk1", "mapk5", "mek1", "titin", "mapk12", "mapk2",
        "mapk9", "mapk10"]
    # targets = ["mapk", "rac1", "tiam1", "mth1", "mtor", "pi3k", 
    #     "rho", "rac1", "prl3", "mmp2", "ck2", "hsp90", "ahr", "mapk2-4", "dsp6","dsp7"]
    # targets = ["mapk2-4", "dsp6","dsp7"]
    scores = [
        [9.7, 9.5, 9.3, 9.1, 8.9, 8.7, 8.2, 8.1, 8.0],
        [9.2, 9.2, 9.1, 9.0, 9.0, 8.4, 8.3, 8.1, 7.9],
        [8.5, 8.4, 7.9, 7.9, 7.8, 7.8, 7.4, 7.4, 7.3],
        [8.3, 8.2, 7.1, 7.0, 7.0, 7.0, 6.9, 6.7, 6.4],
        [8.8, 8.6, 8.5, 8.3, 7.9, 7.9, 7.8, 7.8, 7.5],
        [6.5, 6.3, 6.1, 6.0, 5.8, 5.7, 5.6, 5.6, 5.6],
        [7.4, 7.3, 7.2, 7.0, 6.9, 6.8, 6.7, 6.4, 6.4],
        [5.9, 5.8, 5.7, 5.7, 5.3, 5.3, 5.2, 5.2, 5.2],
        [7.8, 7.7, 7.7, 7.2, 7.1, 6.6, 6.4, 6.3, 6.3],
        [8.0, 7.7, 7.6, 7.5, 7.5, 7.5, 7.5, 7.4, 7.3],
        [9.1, 9.1, 8.4, 8.2, 8.1, 7.8, 7.8, 7.7, 7.7],
        [6.5, 6.5, 6.3, 6.2, 6.2, 6.2, 6.1, 6.0, 6.0],
        [7.2, 7.1, 6.6, 6.2, 6.0, 6.0, 5.9, 5.8, 5.8],
        [7.4, 7.3, 7.2, 7.0, 6.9, 6.8, 6.7, 6.4, 6.4],
        [6.9, 6.8, 6.7, 6.6, 6.5, 6.5, 6.5, 6.4, 6.4],
        [5.9, 5.8, 5.7, 5.7, 5.3, 5.3, 5.2, 5.2, 5.2],
        [6.6, 6.5, 6.5, 6.4, 6.4, 6.3, 6.3, 6.2, 6.2],
        [7.8, 7.5, 7.5, 7.4, 7.1, 7.1, 6.9, 6.5, 6.3],
        [8.3, 8.2, 7.1, 7.0, 7.0, 7.0, 6.9, 6.7, 6.4],
        [8.2, 7.9, 7.5, 7.3, 7.3, 7.3, 7.2, 7.1, 7.1],
        [9.5, 9.2, 8.5, 8.2, 8.2, 7.8, 7.3, 7.2, 6.9],
        [9.3, 8.7, 8.1, 8.1, 7.9, 7.8, 7.7, 7.4, 7.0],
        [9.7, 9.6, 9.5, 9.5, 9.4, 9.4, 7.5, 7.1, 7.0],
        [9.6, 9.1, 8.9, 8.9, 8.8, 8.1, 8.1, 8.0, 7.9],
        [7.8, 7.6, 7.6, 7.4, 7.4, 7.1, 7.1, 7.1, 6.9],
        [6.6, 6.5, 6.3, 6.3, 6.2, 6.2, 6.1, 6.0, 5.9]
    ]

    data, cols = getResidueMap(mapArray, anotherArray, targets, scores, 2.75)

    finaldf = pd.DataFrame(data, columns=cols).set_index('AA')
    finaldf = finaldf.reindex(columns= ['C7', 'C15', 'C9', 'C12', 'C8', 'O1', 'C2', 'C1', 'C6', 'O3', 'C5', 'O2', 'C4', 'C3', 'C10', 'O4', 'C14', 'C13', 'O5', 'C11'])

    # for i in np.arange(1.5, 4.5, 0.5):
    #     data, cols = getResidueMap(mapArray, targets, scores, i)
    #     finaldf += pd.DataFrame(data, columns=cols).set_index('AA') * ((4-i)/4)

    print(finaldf)

    # plt.figure(figsize=(16,9))
    sns.heatmap(finaldf, center=0, vmax=10, vmin=-6)
    #ax = sns.heatmap(arr)
    plt.show()  

    aminoAcids = ['Phe','Met','Leu','Ile','Val','Pro','Tyr','Trp','Cys','Ala','Gly','Ser','Thr','His','Glu','Gln','Asp','Asn','Lys','Arg']
    data = {'AA': aminoAcids}

    cols = ['AA'] 

    for i in range(0, len(mapArray)):
        final = [0]*20

        results = []
        for k in mapArray[i].keys():
            try:
                for e in mapArray[i][k]:
                    results.append(e)
            except:
                continue
        
        for j in range(0,20):
            for aaF in results:
                if aaF == aminoAcids[j].upper():
                    final[j] += 1
        
        data.update({str(i): final})
        cols.append(str(i))
        print(final)
    
    print(data)
    df = pd.DataFrame(data, columns=cols).set_index('AA')

    print(df)

    # plt.figure(figsize=(16,9))
    sns.heatmap(df)
    #ax = sns.heatmap(arr)
    plt.show()  
