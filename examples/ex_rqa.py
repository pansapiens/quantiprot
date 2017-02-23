import os
import sys
sys.path.insert(0, os.path.abspath('..'))

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.feature import FeatureSet
from quantiprot.metrics.aaindex import get_aa2mj
from quantiprot.metrics.rqa import RQAFeatureSet
from quantiprot.metrics.basic import average

from matplotlib import pyplot as plt

# Load the HET-E1 sequence with WD40 repeats:
hete1_seq = load_fasta_file("data/HETE1_PODAS.fasta")

# Prepare FeatureSet for conversion from aa to Miyazawa-Jernigan hydrophobicity:
mj_fs = FeatureSet("mj")
mj_fs.add(get_aa2mj())

# Prepare specialized FeatureSet with basic RQA parameters calculated
# over 100aa window, then smoothed over the 10aa window:
rqa_fs = RQAFeatureSet("rqa", features=['recurrence','determinism'],
                       window=100, metric='taxi', radius=4, dim=4, det_len=8)
rqa_fs.then_all(average, window=10)
print rqa_fs

# From fasta to MJ hydrophobicity to RQA
hete1_rqa_seq = rqa_fs(mj_fs(hete1_seq))

plt.plot(hete1_rqa_seq.columns(feature="recurrence>average"),label="REC")
plt.plot(hete1_rqa_seq.columns(feature="determinism>average"),label="DET")
plt.xlabel('position')
plt.legend()
plt.show()

# Now, try to find coordinates of the repeats region based on determinism:
repeats_pos_fs = RQAFeatureSet("repeat indices", features=['determinism'],
                               window=100, metric='taxi', radius=4, dim=4, det_len=8)
repeats_pos_fs.then_all(average, window=10)

# The following function returns indices at positions where value > 'threshold':
def ind_if_gt(data, **params):
    threshold = params['threshold'] if 'threshold' in params else 0
    return [i for i in range(len(data)) if data[i] > threshold]

# Now add the function to the feature chain with the DET threshold of 0.5:
repeats_pos_fs.then_all(ind_if_gt, threshold=0.5)

# From fasta to MJ hydrophobicity to positions where DET > 0.5
hete1_repeats_pos_seq= repeats_pos_fs(mj_fs(hete1_seq))

for seq in hete1_repeats_pos_seq:
    print min(seq.data), max(seq.data)+99
