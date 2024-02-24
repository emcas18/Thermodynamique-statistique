#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 15:40:27 2021

#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen
"""
from vpython import *
import numpy as np
import math
import matplotlib.pyplot as plt

# win = 500 # peut aider à définir la taille d'un autre objet visuel comme un histogramme proportionnellement à la taille du canevas.

# Déclaration de variables influençant le temps d'exécution de la simulation
Natoms =   # change this to have more or fewer atoms
dt = 1E-6  # pas d'incrémentation temporel

# Déclaration de variables physiques "Typical values"
mass = 4E-3/6E23 # helium mass
Ratom = 0.02 # wildly exaggerated size of an atom
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature

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



# Ajouter un maillage de sphères fixes
Coeurs = []
coeur_radius = 0.04  # Rayon des sphères du maillage
coeur_pos = []
pcoeur = []

# Définir les positions du maillage
num_coeurs = 5
coeurs_spacing = L / (num_coeurs + 1)
for i in range(1, num_coeurs + 1):
    for j in range(1, num_coeurs + 1):
        x = i * coeurs_spacing - L / 2
        y = j * coeurs_spacing - L / 2
        z = 0
        coeur_pos.append(vec(x,y,z))
        Coeurs.append(sphere(pos=vector(x,y,z), radius=coeur_radius, color=color.red))
        pcoeur_x = 0  # quantité de mvt initiale selon l'équipartition
        pcoeur_y = 0
        pcoeur_z = 0
        pcoeur.append(vector(pcoeur_x,pcoeur_y,pcoeur_z))

# Ajouter les sphères du maillage au scénario
for coeur_sphere in Coeurs:
    coeur_sphere.visible = True



#### POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES SPHÈRES ####
Atoms = [] # Objet qui contiendra les sphères pour l'animation
p = [] # quantité de mouvement des sphères
apos = [] # position des sphères
pavg = sqrt(2*mass*1.5*k*T) #Principe de l'équipartition de l'énergie en thermodynamique statistique classique

for i in range(Natoms):
    x = L*random()-L/2 # position aléatoire qui tient compte que l'origine est au centre de la boîte
    y = L*random()-L/2
    z = 0
    apos.append(vec(x,y,z)) # liste de la position initiale de toutes les sphères
    for j in range(num_coeurs):
        while (apos[i]-coeur_pos[j]).mag < (coeur_radius+Ratom):  # garde une sphère plus grosse et colorée parmis toutes les grises
            x += coeur_radius
            apos[i] = vec(x,y,z)
            
    phi = 2*pi*random() # direction aléatoire pour la quantité de mouvement
    px = pavg*cos(phi)  # quantité de mvt initiale selon l'équipartition
    py = pavg*sin(phi)
    pz = 0
    p.append(vector(px,py,pz)) # liste de la quantité de mvt initiale de toutes les sphères
    Atoms.append(simple_sphere(pos=vector(x,y,z), radius=Ratom, color=color.green))
    
    
def randomDirection(p: list):
    angle = random.uniform(0, 2 * math.pi) # Génère un angle aléatoire entre 0 et 2*pi
    p[i].x = p[i].mag * math.cos(angle)
    p[i].y = p[i].mag * math.sin(angle)     
    return p


#### FONCTION POUR IDENTIFIER LES COLLISIONS, I.E. LORSQUE LA DISTANCE ENTRE LES CENTRES DE 2 SPHÈRES EST À LA LIMITE DE S'INTERPÉNÉTRER ####
def checkCollisions():
    hitlist = []   # initialisation
    r2 = coeur_radius+Ratom   # distance critique où les 2 sphères entre en contact à la limite de leur rayon
    r2 *= r2   # produit scalaire pour éviter une comparaison vectorielle ci-dessous
    for i in range(Natoms):
        ai = apos[i]
        for j in range(num_coeurs) :
            aj = coeur_pos[j]
            dr = ai - aj   # la boucle dans une boucle itère pour calculer cette distance vectorielle dr entre chaque paire de sphère
            if mag2(dr) < r2:   # test de collision où mag2(dr) qui retourne la norme élevée au carré de la distance intersphère dr
                hitlist.append([i,j]) # liste numérotant toutes les paires de sphères en collision
    return hitlist


#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
## ATTENTION : la boucle laisse aller l'animation aussi longtemps que souhaité, assurez-vous de savoir comment interrompre vous-même correctement (souvent `ctrl+c`, mais peut varier)
## ALTERNATIVE : vous pouvez bien sûr remplacer la boucle "while" par une boucle "for" avec un nombre d'itérations suffisant pour obtenir une bonne distribution statistique à l'équilibre

while True:
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

    #### LET'S FIND THESE COLLISIONS!!! ####
    hitlist = checkCollisions()

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    for ij in hitlist:

        i = ij[0]  # extraction du numéro des 2 sphères impliquées à cette itération
        j = ij[1]
        
        # Calcule la direction de rebond après la collision
        normal_vector = apos[i] - coeur_pos[j] # Vecteur normal au point de collision sur le cercle
        length = normal_vector.mag # Normalisation du vecteur
        normal_vector = normal_vector / length

        new_direction = randomDirection()
        
        # Vérification si la nouvelle trajectoire contourne le cercle
        new_pos = apos[i] + deltax[i]  # Nouvelle position après la collision
        vector_to_center = new_pos - coeur_pos[j]  # Vecteur allant du centre du cercle à la nouvelle position
        while acos(dot(vector_to_center.normalize(), -normal_vector)) > (math.pi / 4):
            new_direction = randomDirection()

