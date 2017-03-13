import numpy as np
import random
from matplotlib import pyplot as plt

import pickle
with open('off_target_scores.p', 'rb') as inp:
  off_target_scores = pickle.load(inp)
vals = off_target_scores.values()
bins = np.arange(0,570,1)
plt.xlim([0, 50])
plt.hist(vals, bins=bins, alpha=0.5)
plt.title('G-score distribution')
plt.xlabel('G-score')
plt.ylabel('count')
plt.yscale('log', nonposy='clip')
plt.show()
