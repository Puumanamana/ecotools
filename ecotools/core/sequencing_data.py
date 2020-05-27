from tempfile import mkdtemp

from Bio import SeqIO, AlignIO
from ecotools.core.biotable import BioData
    
class SequencingData(BioData):

    def __init__(self, **kwargs):
        BioData.__init__(self, **kwargs)

    def __repr__(self):
        return "OTUs: {} sequences".format(len(self.data))

    def __len__(self):
        return len(self.data)

    def subset_otus(self, otus):
        self.data = {x: self.data[x] for x in otus}

    def trim_aln(self):
        print('Not implemented')

    def to_fasta(self, output='sequences.fasta'):
        SeqIO.write(self.data.values(), output, 'fasta')

    def to_phylip(self, output='sequences.phylip'):
        tmp = mkdtemp()

        fasta_tmp_file = f'{tmp}/sequences.fasta'
        self.to_fasta(fasta_tmp_file)
        AlignIO.convert(fasta_tmp_file, 'fasta', output, 'phylip-relaxed')

        return output
        
    def _to_h5(self, output):
        print('to_h5() not implemented for {} object. Using to_fasta() instead'
              .format(self.__class__.__name__))
        self.to_fasta(output.replace('.h5', '.fasta'))
