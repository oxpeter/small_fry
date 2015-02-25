#!/usr/bin/env python

#!/usr/bin/env python

import sys

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# analyse fluorescence signals:
file_name = sys.argv[1]
results = pd.read_csv( file_name,sep='\t')
# remove all the blank lines:
clean = results.dropna(how='all')
# transpose the table:
clean = clean.transpose()
# remove some empty rows
cleaned_t = clean.dropna(how='all')
# change the header names to be the well:
cleaned_t.columns = cleaned_t.iloc[0]
# remove the non-text rows:
final = cleaned_t.drop(cleaned_t.index[[0,1]])
# convert to numbers:
final = final.convert_objects(convert_numeric=True)

# calculate the signal to noise ratio:
max = final[36:].max()
background = final[:36].mean()
background_std = final[:36].std()
signal2noise = (max-background)/background_std


# identify the treatment type from a tab-separated template:
template_name = sys.argv[2]
layout = pd.read_csv( template_name,sep='\t')
# create dictionary { well:treatment }
layout_d = { (layout[[0]].loc[idx].values[0]+col):row.values[0] for col in layout.columns for idx, row in layout[[col]].iterrows()}

# convert signal to noise list into dataframe:
s2n = pd.DataFrame(signal2noise)
# add a new column to the signal 2 noise which contains the treatment type:
s2n['treatment'] = s2n.index.map(lambda x: layout_d[x])
# create list of all unique treatment types:
treatments = set(s2n.treatment.dropna())
# sort all inotocin treatment types by concentration:
itc_sorted = sorted([ val for val in set(s2n.treatment.dropna()) if val[:8]=='inotocin'], key=(lambda x: int(x[8:]) if x[8:].isdigit() else x[8:]))


# calculate mean signal to noise ratio +/- standard deviation for each treatment:
treat_means = {}
treat_std = {}

for treatment in sorted(treatments):
    t_mean = s2n[s2n.treatment == treatment].mean()
    t_std = s2n[s2n.treatment == treatment].std()
    t_sem = s2n[s2n.treatment == treatment].std()/len(s2n[s2n.treatment == treatment])
    treat_means[treatment] = t_mean
    treat_std[treatment] = t_sem
    print "%-14s %.2f +/- %.2f" % (treatment, t_mean, t_sem)

# plot inotocin signal to noise:
itc_pos = np.arange(0.5, len(itc_sorted)*0.8+0.5, 0.8  )
itc_means = [ float(treat_means[val]) for val in itc_sorted ]
itc_sem = [ treat_std[val]*2 for val in itc_sorted ]


plt.bar(itc_pos, itc_means, yerr=itc_sem, ecolor='red' )
plt.title("Inotocin average signal to noise")
plt.ylim(0,200)
plt.xticks(itc_pos+0.4, [name[8:] for name in itc_sorted])
plt.show()