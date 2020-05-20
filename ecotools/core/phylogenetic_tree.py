import sys
import subprocess
from tempfile import mkdtemp

from Bio.Phylo.Applications import PhymlCommandline

from ecotools.core.biotable import BioData

class PhylogeneticTree(BioData):
    def __init__(self, **kwargs):
        BioData.__init__(self, **kwargs)

    def root_tree(self):
        self.tree = self.tree.root_at_midpoint()

    def to_nwk(self, output):
        self.tree.write(output)

    def update_tree(self, fasta, app='FastTree'):

        tmp = mkdtemp()
        phylip_path = fasta.to_phylip(f'{tmp}/sequences.phylip')

        if app == 'phyml' and len(fasta) < 2000:
            phyml_cmd = PhymlCommandline(input=phylip_path)
            phyml_cmd()
            self.tree_path = f'{self.outdir}/sequences_phyml_tree.txt'
            
        elif app == 'FastTree':
            if app == 'phyml':
                print('Fallback to FastTree (too many OTUs)')

            self.tree_path = f'{self.outdir}/tree.nwk'

            with open(self.tree_path, 'w') as file_handle:
                proc = subprocess.Popen(['FastTree', '-nt', phylip_path],
                                        stdout=file_handle)
            proc.wait()

        else:
            sys.exit('Unknown application {}'.format(app))

        self.root_tree()
        
        
