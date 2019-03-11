import matplotlib.pyplot as plt
import numpy


def plot_data(seq):
       fig=plt.figure(figsize=(10, 8), dpi=800)
       cryo_data=numpy.genfromtxt('vol_vs_lvl_cryo.csv', dtype=float, delimiter=',', names=True)
       lvl_cryo = cryo_data['lvl']
       vol_cryo = cryo_data['vol']
       plt.scatter(lvl_cryo,vol_cryo,s=3, c='black')
       plt.title('Volume at corresponding sensor level from Cryogenics Documentation')
       plt.ylabel('Volume (liters)')
       plt.xlabel('Sensor Level (cm)')
       filename="Volume_vs_levels.png"
       fig.savefig(filename)

