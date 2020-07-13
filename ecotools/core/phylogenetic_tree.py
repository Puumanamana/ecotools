import sys
import subprocess
from tempfile import mkdtemp

from Bio.Phylo.Applications import PhymlCommandline
# from ete3 import Tree
# from Bio import Phylo
from skbio import TreeNode

from ecotools.core.biotable import BioData
from ecotools.decorators import timer

class PhylogeneticTree(BioData):
    def __init__(self, **kwargs):
        BioData.__init__(self, **kwargs)
        self.tree_path = None

    def load(self):
        self.data = TreeNode.read(self.tree_path, 'newick')

    def root(self):
        self.data = self.data.root_at_midpoint()
        # midpoint = self.data.get_midpoint_outgroup()
        # self.data.set_outgroup(midpoint)

    @timer
    def prune(self, otus):
        self.data = self.data.common_ancestor(otus)
        # self.data.prune(otus, preserve_branch_length=True)

    def to_nwk(self, output):
        if self.data is not None:
            self.data.write(output)
        else:
            print('No sequence data. Skipping')

    def _to_h5(self, output):
        print('to_h5() not implemented for {} object. Using to_nwk() instead'
              .format(self.__class__.__name__))
        self.to_nwk(output.replace('.h5', '.nwk'))

    @timer
    def compute_tree(self, sequencing_data, app='FastTree'):

        tmp = mkdtemp()
        phylip_path = sequencing_data.to_phylip(f'{tmp}/sequences.phylip')

        if app == 'phyml' and len(sequencing_data) < 2000:
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
