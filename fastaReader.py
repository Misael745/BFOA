import numpy
from multiprocessing import Manager

import numpy

class fastaReader():
    def __init__(self, path="multifasta.fasta"):
        self.path = path 
        self.seqs = []
        self.names = []
        self.read()

    def read(self):
        with open(self.path, "r") as f:
            lines = f.readlines()

        seq = ""
        for line in lines:
            if line.startswith(">"):
                self.names.append(line[1:].strip())
                if seq:
                    self.seqs.append(seq)
                seq = ""
            else:
                seq += line.strip()
        if seq:
            self.seqs.append(seq)

    
    
    