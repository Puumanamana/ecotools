from skbio.stats.distance import DistanceMatrix
from skbio.stats.distance import permanova, permdisp

from ecotools.decorators import pairwise

@pairwise
def permanova_test(metagenome, factor, permutations=999, pairwise=False):

    distances = DistanceMatrix(metagenome.distance_matrix)
    result = permanova(distances, factor)

    return result

@pairwise
def betadisper_test(metagenome, factor, permutations=999, pairwise=False):

    distances = DistanceMatrix(metagenome.distance_matrix)
    result = permdisp(distances, factor)

    return result
