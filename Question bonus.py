# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 15:05:41 2024

@author: lolo
"""

from scipy.stats import mode
import numpy as np

norm_v = [] # liste des toutes les normes de v pour une seule particule
all_norm = [] # liste des toutes les normes de v de l'ensemble des particules
instant = [] # liste de la position de chaque particule à un instant T

for i in range(Natoms):
    TrackParticule(i)
    v = valeurs_f[i]
    for vector in v:
        norm_v.append(mag(vector))  # Calcul des normes pour chaque vecteur de vitesse
    all_norm.append(norm_v)
    

for part in all_norm:
    instant.append(part[4]) # Parcours de la norme de chaque particule à un instant précis
    

fig, axs = plt.subplots(2, figsize=(8, 12)) 

# Histogramme de la norme de l'ensemble des particules en fonction de t
axs[0].hist(norm_v, bins=50, color='skyblue', edgecolor='black')
axs[0].set_title("Fig. 1 - Distribution de la norme des vitesses de l'ensemble des particules.")
axs[0].set_xlabel('Norme des vitesses [m/s]')
axs[0].set_ylabel('Fréquence [-]')
axs[0].legend()


# Histogramme de la norme de v en fonction de la particule
axs[1].hist(instant, bins=50, color='lightgreen', edgecolor='black')
axs[1].set_title("Fig. 1 - Distribution de la norme des vitesses de l'ensemble des particules.")
axs[1].set_xlabel('Norme des vitesses [m/s]')
axs[1].set_ylabel('Fréquence [-]')
axs[1].legend()

plt.tight_layout()

plt.show()