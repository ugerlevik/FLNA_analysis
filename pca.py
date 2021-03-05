import sys
from prody import *
from matplotlib.pylab import *

import matplotlib.pyplot as plt

id1, id2 = sys.argv[1:]

str1 = parsePDB("pca_"+id1+"_cMD100ns_CA.pdb")
ens1 = parseDCD("pca_"+id1+"_cMD100ns_CA.dcd")

str2 = parsePDB("pca_"+id2+"_cMD100ns_CA.pdb")
ens2 = parseDCD("pca_"+id2+"_cMD100ns_CA.dcd")

ens1.setCoords(str1)
ens1.setAtoms(str1.calpha)
ens1.superpose()

ens2.setCoords(str2)
ens2.setAtoms(str2.calpha)
ens2.superpose()

eda_ens1 = EDA(id1+' Ensemble')
eda_ens1.buildCovariance( ens1 )
eda_ens1.calcModes()

eda_ens2 = EDA(id2+' Ensemble')
eda_ens2.buildCovariance( ens2 )
eda_ens2.calcModes()

# showProjection(ens1, eda_ens1[:3], color = 'black', marker = '.', label = 'Wild-type', alpha = 0.5)
# showProjection(ens2, eda_ens2[:3], color = 'red', marker = '.', label = id2, alpha = 0.5)
# plt.legend(loc='upper left')
# plt.savefig('PCA_'+id1+'_'+id2+'_3D.pdf')
# plt.close()

showProjection(ens1, eda_ens1[:2], rmsd = False, color = 'black', marker = '.', label = 'Wild-type', alpha = 0.5)
showProjection(ens2, eda_ens2[:2], rmsd = False, color = 'red', marker = '.', label = id2, alpha = 0.5)
plt.legend(loc='upper right')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.xlim(-85,60)
plt.ylim(-60,75)
plt.savefig('PCA_'+id1+'_'+id2+'_2D.pdf', bbox_inches = 'tight', pad_inches = 0)
