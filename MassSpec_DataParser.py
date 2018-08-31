# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 12:17:55 2018

@author: darne
"""

from rdkit import Chem
import os
import numpy as np
import pickle
import cirpy

MS_Records = []

MassBank_Folder = os.listdir("C:/Users/darne/MassBank-data")

all_canon_SMILES = []

for subfolder in MassBank_Folder:
    if "." in subfolder:
        continue
    else:
        subfolder_path = "C:/Users/darne/MassBank-data" + "/" + subfolder
        subfolder_files = os.listdir(subfolder_path) #list of filename strings

        EXTENSION = "txt"

        for file in subfolder_files:
            if not file.endswith(EXTENSION):
                continue
            else:
                textfilepath = subfolder_path + "/" + file
                with open(textfilepath) as f:
                    text = f.read()
             
                    mz_list = []
                    
                    spectrum_data = text.split("PK$PEAK:")[-1]    
                    spectrum_data = spectrum_data.splitlines()
                    spectrum_data.remove(' m/z int. rel.int.')
                    spectrum_data.remove('//')
                    split_spectrum_data = []  
                    for line in spectrum_data:
                        line = line.split()
                        mz_entry = (line[0], line[1], line[2])
                        mz_list.append(mz_entry)
                        
                    mz_array = np.array(mz_list)
                    
                    text_lines = text.splitlines()
                    detail_lines = []
    
                    for line in text_lines:
        
                        if "CH$SMILES:" in line:
                            SMILES_text = line[11:]

                        elif "CH$IUPAC:" in line:
                            InChI_text = line[10:]
                            
                        elif "CH$NAME:" in line:
                            Mol_name = line[9:]
                           
                        elif (line[0:3] == "AC$") or (line[0:3] == "MS$"):
                           detail_lines.append(line[3:])
                    Details = detail_lines
                  
                    
                    if SMILES_text != "N/A":
                        mol = Chem.MolFromSmiles(SMILES_text)
                        if not mol:
                            continue
                        else:
                            canonical_SMILES = Chem.MolToSmiles(mol)
                            canonical_InChI = Chem.MolToInchi(mol)
                        
                    elif InChI_text != "N/A":
                        mol = Chem.MolFromInchi(InChI_text)
                        if not mol:
                            continue
                        else:
                            canonical_SMILES = Chem.MolToSmiles(mol)
                            canonical_InChI = Chem.MolToInchi(mol)
                    
                    else: 
                        SMILES_text = cirpy.resolve(Mol_name, "smiles")
                        if SMILES_text == None:
                            
                            InChI_text = cirpy.resolve(Mol_name, "stdinchi")
                            
                            if InChI_text == None:
                                continue
                                
                            else:
                                mol = Chem.MolFromInchi(InChI_text)
                                if not mol:
                                    continue
                                    
                                else:
                                    canonical_SMILES = Chem.MolToSmiles(mol)
                                    canonical_InChI = Chem.MolToInchi(mol)                                        

                        else:
                            mol = Chem.MolFromSmiles(SMILES_text)
                            if not mol:
                                continue
                            else:
                                canonical_SMILES = Chem.MolToSmiles(mol)
                                canonical_InChI = Chem.MolToInchi(mol)


            molecule_record = {"Name":Mol_name,
                               "SMILES":SMILES_text,
                               "InChI":InChI_text,
                               "Canonical SMILES":canonical_SMILES,
                               "Canonical InChI":canonical_InChI,
                               "m/z Array":mz_array,
                               "Details":Details
                               }
            
            MS_Records.append(molecule_record)

MS_Records_pickle = open("MS_Records.pickle", "wb")
pickle.dump(MS_Records, MS_Records_pickle)
MS_Records_pickle.close()