"""
Created on Fri Feb  5 15:40:27 2021

#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen
"""

import matplotlib.pyplot as plt
from vpython import *
from scipy.stats import maxwell as maxwell
import numpy as np

import math


# win = 500 # peut aider à définir la taille d'un autre objet visuel comme un histogramme proportionnellement à la taille du canevas.

# Déclaration de variables influençant le temps d'exécution de la simulation
Natoms = 200                 # change this to have more or fewer atoms
num_ions = 6                 # Nombre d'ions
dt = 1E-7                    # Pas d'incrémentation temporelle
blue = vector(0, 0, 1)

# Déclaration de variables physiques "Typical values"
mass = 9.10E-31    # Masse electron
Ratom = 0.01    # Taille electron
Rions = 0.03                # Taille des ions
k = 1.4E-23    # Constante de Boltzmann
T = 300                  # Température

#### CANEVAS DE FOND ####
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container and spheres below
animation = canvas( width=750, height=500) # , align='left')
animation.range = L
# animation.title = 'Cinétique des gaz parfaits'
# s = """  Simulation de particules modélisées en sphères dures pour représenter leur trajectoire ballistique avec collisions. Une sphère est colorée et grossie seulement pour l’effet visuel permettant de suivre sa trajectoire plus facilement dans l'animation, sa cinétique est identique à toutes les autres particules.

# """
# animation.caption = s

#### ARÊTES DE BOÎTE 2D ####
d = L/2+Ratom
r = 0.005
cadre = curve(color=gray, radius=r)
cadre.append([vector(-d,-d,0), vector(d,-d,0), vector(d,d,0), vector(-d,d,0), vector(-d,-d,0)])


#### POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES SPHÈRES ####
Atoms = [] # Objet qui contiendra les sphères pour l'animation
p = [] # quantité de mouvement des sphères
apos = [] # position des sphères
pavg = sqrt(2*mass*1.5*k*T) #Principe de l'équipartition de l'énergie en thermodynamique statistique classique

for i in range(Natoms):
    x = L*random()-L/2 # position aléatoire qui tient compte que l'origine est au centre de la boîte
    y = L*random()-L/2
    z = 0
    if i == 0:  # garde une sphère plus grosse et colorée parmis toutes les grises
        Atoms.append(simple_sphere(pos=vector(x,y,z), radius=0.03, color=color.magenta))
    else: Atoms.append(simple_sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z)) # liste de la position initiale de toutes les sphères
#    theta = pi*random() # direction de coordonnées sphériques, superflue en 2D
    phi = 2*pi*random() # direction aléatoire pour la quantité de mouvement
    px = pavg*cos(phi)  # quantité de mvt initiale selon l'équipartition
    py = pavg*sin(phi)
    pz = 0
    p.append(vector(px,py,pz)) # liste de la quantité de mvt initiale de toutes les sphères

#### POSITION DES IONS ET CREATION####
ions = []
ion_positions = []  # Liste pour stocker les vecteurs de positions des ions
spacing = L / num_ions

for i in range(num_ions):
    for j in range(num_ions):
        x = -L / 2 + (i + 0.5) * spacing
        y = -L / 2 + (j + 0.5) * spacing
        z = 0
        ions.append(simple_sphere(pos=vector(x,y,z), radius = Rions, color=blue))
        ion_positions.append(vector(x, y, z))  # Création du vecteur de position de l'ion


#### FONCTION POUR IDENTIFIER LES COLLISIONS ####
def checkCollisions():
    hitlist = []

    r2 = (Ratom + Rions)
    r2 *= r2  # produit scalaire pour éviter une comparaison vectorielle ci-dessous

    for i in range(Natoms):
        ai = apos[i]
        for j in range(num_ions**2):
            ipos = ion_positions[j]

            # Calcul de la distace vectorielle un électron les ions
            dr = ai - ipos
            if mag2(dr) < r2:
                hitlist.append([i, j])

    return hitlist
paire_collision = []
#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
temps = 0
p_atom = []  # Quantité de mouvement des atomes (électrons)
p_moyen = []  # Quantité de mouvement moyenne
temps_ecoule = []

for itr in range(20000):
    rate(300)  # limite la vitesse de calcul de la simulation pour que l'animation soit visible à l'oeil humain!
    #### DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
    vitesse = []   # vitesse instantanée de chaque sphère
    deltax = []  # pas de position de chaque sphère correspondant à l'incrément de temps dt
    for i in range(Natoms):
        vitesse.append(p[i]/mass)   # par définition de la quantité de nouvement pour chaque sphère
        deltax.append(vitesse[i] * dt)   # différence avant pour calculer l'incrément de position
        Atoms[i].pos = apos[i] = apos[i] + deltax[i]  # nouvelle position de l'atome après l'incrément de temps dt
     #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS AVEC LES PAROIS DE LA BOÎTE ####
    for i in range(Natoms):
        loc = apos[i]
        if abs(loc.x) > L/2:
            if loc.x < 0: p[i].x =  abs(p[i].x)  # renverse composante x à la paroi de gauche
            else: p[i].x =  -abs(p[i].x)   # renverse composante x à la paroi de droite
        if abs(loc.y) > L/2:
            if loc.y < 0: p[i].y = abs(p[i].y)  # renverse composante y à la paroi du bas
            else: p[i].y =  -abs(p[i].y)  # renverse composante y à la paroi du haut

    #### IDENTIFICATION ET GESTION DES COLLISIONS ####
    hitlist = checkCollisions()
    if hitlist:
        paire_collision.extend(hitlist)
    for ij in hitlist:
        posi = apos[i]              # position de l'électron
        posj = ion_positions[j]              # position de l'ion
        rrel = posi - posj 
        vrel = -p[i]/mass
        
        if vrel.mag2 == 0: continue  
        if rrel.mag > Rions + Ratom: continue
            
        dx = dot(rrel, vrel.hat)                # rrel.mag*cos(theta) où theta is the angle between vrel and rrel:
        dy = cross(rrel, vrel.hat).mag          # rrel.mag*sin(theta)
        alpha = asin(dy/(Rions + Ratom))  
        d = d = (Rions + Ratom)*cos(alpha) - dx
        deltat = d/vrel.mag 
        #### CHANGE L'INTERPÉNÉTRATION DES SPHÈRES PAR LA CINÉTIQUE DE COLLISION ####
        posi = posi + vrel * deltat  

        # Calcul de la qte de mvt des électrons en collision
        # selon le modèle de Drude
        p_norm = mass * maxwell.rvs(scale=sqrt(k * T / mass))           
        phi = 2 * pi * random()                                         
        px = p_norm * cos(phi)                                        
        py = p_norm * sin(phi)
        p[i] = vector(px, py, 0)

        apos[i] = posi + (p[i]/mass)*deltat         



    T = np.mean([mag2(vec) for vec in p])/(2*k*mass)
    temps += dt

    temps_ecoule.append(temps)
    p_atom.append((p[0].mag)) 
    p_moyen.append(np.mean([mag(vec) for vec in p]))

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    #for ij in hitlist:
        

        # définition de nouvelles variables pour chaque paire de sphères en collision
        #i = ij[0]  # extraction du numéro des 2 sphères impliquées à cette itération
        #j = ij[1]
        #ptot = p[i]   # quantité de mouvement totale des 2 sphères
        #mtot = masse_reduite    # masse totale des 2 sphères
        ##Vcom = ptot/mtot   # vitesse du référentiel barycentrique/center-of-momentum (com) frame
        ##posi = apos[i]   # position de chacune des 2 sphères
        ##posj = apos[j]
        ##vi = p[i]/mass   # vitesse de chacune des 2 sphères
        ##vj = p[j]/mass
        ##rrel = posi-posj  # vecteur pour la distance entre les centres des 2 sphères
        ##vrel = vj-vi   # vecteur pour la différence de vitesse entre les 2 sphères

        # exclusion de cas où il n'y a pas de changements à faire
        ##if vrel.mag2 == 0: continue  # exactly same velocities si et seulement si le vecteur vrel devient nul, la trajectoire des 2 sphères continue alors côte à côte
        ##if rrel.mag > Ratom: continue  # one atom went all the way through another, la collision a été "manquée" à l'intérieur du pas deltax

        # calcule la distance et temps d'interpénétration des sphères dures qui ne doit pas se produire dans ce modèle
        ##dx = dot(rrel, vrel.hat)       # rrel.mag*cos(theta) où theta is the angle between vrel and rrel:
        ##dy = cross(rrel, vrel.hat).mag # rrel.mag*sin(theta)
        ##alpha = asin(dy/(2*Ratom))  # alpha is the angle of the triangle composed of rrel, path of atom j, and a line from the center of atom i to the center of atom j where atome j hits atom i
        ##d = (2*Ratom)*cos(alpha)-dx # distance traveled into the atom from first contact
        ##deltat = d/vrel.mag         # time spent moving from first contact to position inside atom

        #### CHANGE L'INTERPÉNÉTRATION DES SPHÈRES PAR LA CINÉTIQUE DE COLLISION ####
        #posi = posi-vi*deltat   # back up to contact configuration
        #posj = posj-vj*deltat
        #pcomi = p[i]-mass*Vcom  # transform momenta to center-of-momentum (com) frame
        #pcomj = p[j]-mass*Vcom
        #rrel = hat(rrel)    # vecteur unitaire aligné avec rrel
        #pcomi = pcomi-2*dot(pcomi,rrel)*rrel # bounce in center-of-momentum (com) frame
        #pcomj = pcomj-2*dot(pcomj,rrel)*rrel
        #p[i] = pcomi+mass*Vcom # transform momenta back to lab frame
        #p[j] = pcomj+mass*Vcom
        #apos[i] = posi+(p[i]/mass)*deltat # move forward deltat in time, ramenant au même temps où sont rendues les autres sphères dans l'itération
        #apos[j] = posj+(p[j]/mass)*deltat