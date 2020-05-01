from skbio.stats.distance import DistanceMatrix
from skbio.stats.distance import permanova, permdisp

def permanova_test(metagenome, column, permutations=999):
    
    result = permanova(DistanceMatrix(metagenome.distance_matrix),
                       metagenome.medata.data[column])

    print(result)


def betadisp(metagenome, column, permutations=999):
    
    result = permdisp(DistanceMatrix(metagenome.distance_matrix),
                      metagenome.medata.data[column])

    print(result)

    
