from Bio import SeqIO

#handle= open("BE MAFB5.ab1", "rb")
#for record in SeqIO.parse(handle, "abi"):
#    print(record)


record = SeqIO.read("BE MAFB5.ab1", "abi")
print(list(record.annotations.keys()),"\n")
print(list(record.annotations["abif_raw"].keys()),"\n")

from collections import defaultdict
channels = ['DATA1', 'DATA2', 'DATA3', 'DATA4', 'DATA5', 'DATA6', 'DATA7', 'DATA8', 'DATA9', 'DATA10', 'DATA11', 'DATA12']
trace = defaultdict(list)
for c in channels:
    trace[c] = record.annotations["abif_raw"][c]

import matplotlib


plt.plot(trace["DATA9"], color="blue")
plt.plot(trace["DATA10"], color="red")
plt.plot(trace["DATA11"], color="green")
plt.plot(trace["DATA12"], color="yellow")
plt.show()
