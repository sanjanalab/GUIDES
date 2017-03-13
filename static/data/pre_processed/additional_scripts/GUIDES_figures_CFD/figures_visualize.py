# Visualize results

# import packages
import pickle
from matplotlib import pyplot as plt

# Read in results
with open('all_scores.p', 'rb') as all_scores_f:
  all_scores = pickle.load(all_scores_f)

with open('selected_scores.p', 'rb') as selected_scores_f:
  selected_scores = pickle.load(selected_scores_f)


print 'visualizing all_scores'
vals = all_scores
max_bin = int(max(vals)) + 13
bins = range(0,max_bin,1)
plt.xlim([0, max_bin])
plt.hist(vals, bins=bins, alpha=0.5)
plt.title('G-score distribution')
plt.xlabel('G-score')
plt.ylabel('count')
plt.yscale('log', nonposy='clip')
plt.savefig('all_scores_log.pdf', bbox_inches='tight')

plt.clf()

print 'visualizing all_scores'
vals = all_scores
min_bin = 1
max_bin = int(max(vals)) + 13
bins = range(min_bin,max_bin,1)
plt.xlim([min_bin, max_bin])
plt.hist(vals, bins=bins, alpha=0.5)
plt.title('G-score distribution')
plt.xlabel('G-score')
plt.ylabel('count')
plt.savefig('all_scores_lin.pdf', bbox_inches='tight')

plt.clf()

print 'visualizing selected_scores'
vals = selected_scores
max_bin = int(max(vals)) + 13
bins = range(0,max_bin,1)
plt.xlim([0, max_bin])
plt.hist(vals, bins=bins, alpha=0.5)
plt.title('G-score distribution')
plt.xlabel('G-score')
plt.ylabel('count')
plt.yscale('log', nonposy='clip')
plt.savefig('selected_scores_log.pdf', bbox_inches='tight')

plt.clf()

print 'visualizing selected_scores'
vals = selected_scores
min_bin = 1
max_bin = int(max(vals)) + 13
bins = range(min_bin,max_bin,1)
plt.xlim([min_bin, max_bin])
plt.hist(vals, bins=bins, alpha=0.5)
plt.title('G-score distribution')
plt.xlabel('G-score')
plt.ylabel('count')
plt.savefig('selected_scores_lin.pdf', bbox_inches='tight')
