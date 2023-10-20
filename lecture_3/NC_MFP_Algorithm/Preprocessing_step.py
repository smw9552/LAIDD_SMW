from rdkit import Chem

class Preprocessing:


    def get_Query_Mol_single(self):

        #Query molecule name: Vismiaguianin A
        qMol = Chem.MolFromSmarts('[#8]1-[#6](=[#6](-[#6](=[#8])-[#8]2)-[#6](-[#1])=[#6](-[#1])-[#6]1(-[#6](-[#1])(-[#1])-[#1])-[#6](-[#1])(-[#1])-[#1])-[c]3:[c]2:[c]4:[c](:[c](-[#1]):[c](:[c](-[#1]):[c]4-[#8]-[#1])-[#8]-[#6](-[#1])(-[#1])-[#1]):[c](-[#1]):[c]3-[#6](-[#1])(-[#1])-[#1]')

        return qMol


    def get_Query_Mols(self, FilePath):

        f = open(FilePath, 'r')
        qMols = []
        while True:
            line = f.readline()
            if not line:break
            new_line = line.split("\t")
            qMols.append(Chem.MolFromSmarts(new_line[1].strip()))

        return qMols


    def get_Query_Mol_single_Smarts(self):

        #Query molecule name: Vismiaguianin A
        qMol_Smarts = '[#8]1-[#6](=[#6](-[#6](=[#8])-[#8]2)-[#6](-[#1])=[#6](-[#1])-[#6]1(-[#6](-[#1])(-[#1])-[#1])-[#6](-[#1])(-[#1])-[#1])-[c]3:[c]2:[c]4:[c](:[c](-[#1]):[c](:[c](-[#1]):[c]4-[#8]-[#1])-[#8]-[#6](-[#1])(-[#1])-[#1]):[c](-[#1]):[c]3-[#6](-[#1])(-[#1])-[#1]'

        return qMol_Smarts

    def get_Query_Mols_Smarts(self, FilePath):

        f = open(FilePath, 'r')
        qMols_Smarts = []
        while True:
            line = f.readline()
            if not line:break
            new_line = line.split("\t")
            qMols_Smarts.append(new_line[1].strip())

        return qMols_Smarts