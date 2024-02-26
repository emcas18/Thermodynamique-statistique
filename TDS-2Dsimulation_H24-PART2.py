"""
Created on Fri Feb  5 15:40:27 2021

#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen
"""


import matplotlib.pyplot as plt
from vpython import *
import numpy as np
from scipy.stats import maxwell as maxwell

#Variables
Natom = 200  # Change this to have more or fewer free electrons
Nion = 36   # Nombre de coeurs
dt = 1E-7    # pas d'incrémentation temporel
mass = 9.1E-31 # Masse de electron
Ratom = 0.01 # wildly exaggerated size of an electron
Rion = 0.03 # Rayon des ions
k = 1.4E-23   # Boltzmann constant
T = 300     # around room temperature
e = 1.6E-19 # charge d'un électron

#### CANEVAS DE FOND ####
L = 1  # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container and spheres below
animation = canvas(width=750, height=500)
animation.range = L

#### ARÊTES DE BOÎTE 2D ####
d = L/2+Ratom
r = 0.005
cadre = curve(color=gray, radius=r)
cadre.append([vector(-d, -d, 0), vector(d, -d, 0), vector(d, d, 0), vector(-d, d, 0), vector(-d, -d, 0)])


#### POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES ÉLECTRONS ####
Atom = [] # objet qui contiendra les sphères d'électrons pour l'animation
p = [] # quantité de mouvement des electrons
apos = [] # position des electrons
pavg = sqrt(2 * mass * 1.5 * k * T) #Principe de l'équipartition de l'énergie en thermodynamique statistique classique

for i in range(Natom):
    x = L * random() - L / 2 # position aléatoire qui tient compte que l'origine est au centre de la boîte
    y = L * random() - L / 2

    Atom.append(simple_sphere(pos=vector(x, y, 0), radius=Ratom, color=gray))
    apos.append(vec(x, y, 0))  # Liste de la position initiale des electrons

    phi = 2 * pi * random()  # Direction aléatoire pour la quantité de mouvement
    px = pavg * cos(phi) # Quantité de mouvement initiale
    py = pavg * sin(phi)
    p.append(vector(px, py, 0))  # Liste de la quantité de mouvement initiale de toutes les sphères d'électrons
    
# CREATION COEURS IONIQUEs 
x_ion, step_x = np.linspace(-d + 4 * Rion, d - 4 * Rion, int(np.sqrt(Nion)), retstep=True)
y_ion = np.linspace(-d + 2 * Rion, d - 2 * Rion, int(np.sqrt(Nion)))

ipos = []  # Liste des positions de tous les ions

# Nombre de lignes et de colonnes dans la grille
num_rows = len(x_ion)
num_cols = len(y_ion)

# Espace entre chaque ion
spacing_x = step_x / 2
spacing_y = 2 * Rion

for i in range(num_rows):
    for j in range(num_cols):
        # Position en x
        x_pos = x_ion[i]

        # Position en y
        y_pos = y_ion[j]

        # Créer la sphère à la position calculée
        simple_sphere(pos=vector(x_pos, y_pos, 0), radius=Rion, color=color.blue)
        ipos.append(vec(x_pos, y_pos, 0))
        
def champElec(valeur_champ):
    x = -e*valeur_champ*dt
    dp = vector(x,0,0)
    return dp
    

def checkCollisions():
    """
    Fonction pour identifier les collisions entre les électrons et les ions.
    :return:
        hitlist (list) : Liste contenant une paire électron-ion qui ont subi une collision
    """
    hitlist = []

    # Distance critique où les 2 sphères entrent en contact à la limite de leur rayon
    r2 = Ratom + Rion  # Distance critique où les 2 sphères (électron et ion) entrent en contact à la limite de leur rayon
    r2 *= r2  # Produit scalaire pour éviter une comparaison vectorielle ci-dessous

    for i in range(Natom):
        ai = apos[i]
        for j in range(Nion):
            ionj = ipos[j]

            # Calcul de la distance vectorielle entre un électron et les ions
            dr = ai - ionj
            if mag2(dr) < r2:
                hitlist.append([i, j])

    return hitlist

#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
temps = 0   
p_atom = [] #Qte mvt electron           
p_moy = [] #Qte mvt moyenne pour tout les electron
x_moy = []
y_moy = []         
temps_tot = [] #Temps total 
valeur_E = -10     
for itr in range(1000):
    rate(300)  # Limite la vitesse de calcul de la simulation pour que l'animation soit visible à l'œil humain!

    #### DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
    vitesse = []  # Vitesse instantanée de chaque électron
    deltax = []  # Pas de position de chaque sphère correspondant à l'incrément de temps dt
    dp = vector(0,0,0)
    for i in range(Natom):
        dp = champElec(valeur_E)
        dv = dp / mass
        vitesse.append((p[i] / mass) + dv)  # Par définition de la quantité de mouvement pour chaque électron
        deltax.append(vitesse[i] * dt)  # Différence avant pour calculer l'incrément de position
        apos[i] = apos[i] + deltax[i]
        loc = apos[i]
        if abs(loc.x) > L / 2:
            if loc.x < 0:
                p[i].x = abs(p[i].x)  # Renverse la composante x au mur de gauche
                apos[i].x = -(L / 2) - (apos[i].x + (L / 2))
            else:
                p[i].x = -abs(p[i].x)  # Renverse la composante x au mur de droite
                apos[i].x = (L / 2) - (apos[i].x - (L / 2))
        if abs(loc.y) > L / 2:
            if loc.y < 0:
                p[i].y = abs(p[i].y)  # Renverse la composante y au mur du bas
                apos[i].y = -(L / 2) - (apos[i].y + (L / 2))
            else:
                p[i].y = -abs(p[i].y)  # Renverse la composante y au mur du haut
                apos[i].y = (L / 2) - (apos[i].y - (L / 2))
        Atom[i].pos = apos[i]  # Nouvelle position de l'atome après l'incrément de temps dt

    #### LET'S FIND THESE COLLISIONS!!! ####
    hitlist = checkCollisions()

   # Boucle pour traiter les collisions identifiées dans hitlist
    for ij in hitlist:
        i = ij[0]               
        j = ij[1]               
        posi = apos[i]         
        posj = ipos[j]         
        rrel = posi - posj      
        vrel = -p[i] / mass    

    # Vérifier si la vitesse relative est nulle, si c'est le cas, passer à la prochaine collision
        if vrel.mag2 == 0: continue

    # Vérifier si la distance entre l'électron et l'ion est supérieure à la somme de leurs rayons, si c'est le cas, passer à la prochaine collision
        if rrel.mag > Rion + Ratom: continue

    # Calcul de la composante x du déplacement relatif entre l'électron et l'ion
        dx = dot(rrel, vrel.hat)
    # Calcul de la composante y du déplacement relatif entre l'électron et l'ion
        dy = cross(rrel, vrel.hat).mag
    # Calcul de l'angle alpha entre la composante y et la somme des rayons de l'électron et de l'ion
        alpha = asin(dy / (Rion + Ratom))
    # Calcul de la distance de chevauchement entre l'électron et l'ion
        d = (Rion + Ratom) * cos(alpha) - dx
    # Calcul du temps nécessaire pour atteindre le point de contact entre l'électron et l'ion
        deltat = d / vrel.mag

    # Mettre à jour la position de l'électron après la collision
        posi = posi + vrel * deltat

    # Générer une nouvelle quantité de mouvement pour l'électron après la collision
        p_norm = mass * maxwell.rvs(scale=sqrt(k * T / mass))
        phi = 2 * pi * random()
        px = p_norm * cos(phi)
        py = p_norm * sin(phi)
        p[i] = vector(px, py, 0)

    # Mettre à jour la position de l'électron après la collision et son nouveau déplacement
        apos[i] = posi + (p[i] / mass) * deltat

# Calcul de la température moyenne du système après les collisions
    T = np.mean([mag2(vec) for vec in p]) / (3 * k * mass)
    temps += dt

# Enregistrer le temps écoulé depuis le début de la simulation
    temps_tot.append(temps)

# Enregistrer la norme de la quantité de mouvement du premier électron
    p_atom.append(p[0].mag)

# Enregistrer la quantité de mouvement moyenne de tous les électrons
    p_moy.append(np.mean([mag(vec) for vec in p]))

# Enregistrer la posititon moyenne de tous les électrons
    x_moy.append(np.mean([horiz.x for horiz in apos]))
    y_moy.append(np.mean([vert.y for vert in apos]))
