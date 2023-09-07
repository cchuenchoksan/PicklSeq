import pickle
import matplotlib.pyplot as plt
import numpy as np

with open("output.pkl", "rb") as file:
    loaded_list = pickle.load(file)

data = []
for seq in loaded_list:
    data.append(list(seq[2]))

data_trimmed = [row[0:3000000] for row in data]

transdict = {"A":0, "C": 1, "G":2, "T":3}
alignments_ints = np.vectorize(transdict.get)(np.array(data_trimmed[0:200]))
plt.imshow(alignments_ints, cmap='jet')
plt.show()