import numpy as np
import matplotlib.pyplot as plt

##t-maps

X = ['visLV5','visRV5','tacLS1', 'tacLV5', 'tacRV5']
visVertVsHor = [0.16,
0.66,
0.58,
0.33,
0.58]
tacVertVsHor = [0.41,
0.58,
0.58,
0.33,
0.50]

X_axis = np.arange(len(X))
plt.bar(X_axis - 0.2, visVertVsHor, 0.4, label = 'visVertVsHor')
plt.bar(X_axis + 0.2, tacVertVsHor, 0.4, label = 'tacVertVsHor')
plt.plot(1.5*X_axis-0.5, [0.5,0.5,0.5,0.5,0.5], linestyle = 'dashed', color ="black")
plt.xticks(X_axis, X)
plt.xlabel("ROIs")
plt.ylabel("Decoding Accuracy")
plt.xlim([-0.5,5])
plt.ylim([0,1.])
plt.title("Within-Modality Decoding\nSub-002, Features = 363, smoothing = 0, t-maps")
plt.legend()
#plt.show()
plt.savefig('/Users/shahzad/Files/fMRI/TacMotionAnalysis/derivatives/cosmoMvpa/2mmvoxel/Sub-002_WithinModalityDecoding_Features363_smoothing0_tmaps.png')