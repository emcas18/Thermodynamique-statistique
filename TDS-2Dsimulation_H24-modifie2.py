"""
Created on Fri Feb  5 15:40:27 2021

#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen
"""
import numpy.random as random
from vpython import *
from vpython import vector
import numpy as np
import math
import matplotlib.pyplot as plt

# win = 500 # peut aider à définir la taille d'un autre objet visuel comme un histogramme proportionnellement à la taille du canevas.
blue = vector(0, 0, 1)
# Déclaration de variables influençant le temps d'exécution de la simulation
Natoms = 200  # change this to have more or fewer atoms
dt = 1E-5  # pas d'incrémentation temporel

# Déclaration de variables physiques "Typical values"
mass = 4E-3/6E23 # helium mass
Ratom = 0.01 # wildly exaggerated size of an atom
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature
ions_mass = 1.67E-27
Rions = 0.04

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

# Ajout des cœurs ioniques
ions = []
num_ions = 6
spacing = L / num_ions


for i in range(num_ions):
    for j in range(num_ions):
        x = -L / 2 + (i + 0.5) * spacing
        y = -L / 2 + (j + 0.5) * spacing
        z = 0
        ions.append(vector(x, y, z))

for pos in ions:
    sphere(pos=pos, radius=Rions, color=blue)
    
##FONCTION POUR UPDATE LA QTE DE MVT

def update_momentum(T):
    # Calcul de la nouvelle norme de la quantité de mouvement à partir de la distribution de Maxwell-Boltzmann
    p_norm = np.sqrt(2 * mass * k * T)

    # Génération aléatoire d'un angle theta pour le plan xy
    theta = np.random.uniform(0, 2 * np.pi)  

    # Calcul des composantes de la nouvelle quantité de mouvement en 2D
    px = p_norm * np.cos(theta)
    py = p_norm * np.sin(theta)

    # Retourne un objet vector VPython avec les composantes calculées
    return vector(px, py, 0)  # La composante z est 0 puisque nous sommes en 2D

##FONCTION POUR TRACK UNE PARTICULE

def TrackParticule(particule):
    global valeurs_f
    tau = []
    # Initialisation des listes nécessaires
    temps_final = []
    v = []
    t = []
    d = []

    # Calcul des quantités de mouvement de la particule
    qte_mvt_part = np.array([sous_liste[particule - 1] for sous_liste in qte_mvt if sous_liste])

    # Parcours des collisions
    for index, paire in enumerate(paire_collision):
        if particule in paire:
            t.append(itr_collisions[index])
            # Calcul des temps finaux entre les collisions
            temps_final = [0] + [(t[i] - t[i-1]) * dt for i in range(1, len(t))]

    # Calcul des vitesses aux instants des collisions
    for i in t:
        v.append(qte_mvt_part[i] / mass)

    # Calcul des distances parcourues entre les collisions
    for i in range(len(v)):
        d.append(mag(v[i]) * temps_final[i])
    return v, d, temps_final, qte_mvt_part
    # Affectation du résultat à la variable globale
     #valeurs_f = (v, d, temps_final)



    
#### FONCTION POUR IDENTIFIER LES COLLISIONS, I.E. LORSQUE LA DISTANCE ENTRE LES CENTRES DE 2 SPHÈRES EST À LA LIMITE DE S'INTERPÉNÉTRER ####
def checkCollisions():
    hitlist = []   # initialisation
    r2 = (Ratom + Rions)   # distance critique où les 2 sphères entrent en contact à la limite de leur rayon
    r2 *= r2   # produit scalaire pour éviter une comparaison vectorielle ci-dessous
    for i in range(Natoms):
        # Vérifier si la sphère i est un électron mobile
        if Atoms[i].radius == Ratom:
            ai = apos[i]
            for j in range(len(ions)):
                # Vérifier si la sphère j est un ion
                aj = ions[j]
                dr = ai - aj   # calcul de la distance entre l'électron et l'ion
                if mag2(dr) < r2:   # test de collision
                    hitlist.append([i,j]) # ajouter la paire en collision à la liste
    return hitlist

 

#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
## ATTENTION : la boucle laisse aller l'animation aussi longtemps que souhaité, assurez-vous de savoir comment interrompre vous-même correctement (souvent `ctrl+c`, mais peut varier)
## ALTERNATIVE : vous pouvez bien sûr remplacer la boucle "while" par une boucle "for" avec un nombre d'itérations suffisant pour obtenir une bonne distribution statistique à l'équilibre
itr_collisions = []
paire_collision = []
qte_mvt = []
qte_mvt_moyenne = []
Temp = []
valeurs_f = np.array([]).reshape(3, 0)
for itr in range(10):

    p_iter = []
    temps = itr*dt
    rate(300)  # limite la vitesse de calcul de la simulation pour que l'animation soit visible à l'oeil humain!
    #### DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
    vitesse = []   # vitesse instantanée de chaque sphère
    deltax = []  # pas de position de chaque sphère correspondant à l'incrément de temps dt
    for i in range(Natoms):
        vitesse.append(p[i]/mass)   # par définition de la quantité de nouvement pour chaque sphère
        deltax.append(vitesse[i] * dt)   # différence avant pour calculer l'incrément de position
        Atoms[i].pos = apos[i] = apos[i] + deltax[i]  # nouvelle position de l'atome après l'incrément de temps dt
        p_iter.append(p[i])
    qte_mvt.append(p_iter)
    for iteration in qte_mvt:
        for qte in iteration:
            qte_mvt_moyenne.append(np.mean(mag(qte)))
    for i in qte_mvt_moyenne:
        Temp.append(i**2 / (2 * mass * k))

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
    if hitlist:
        paire_collision.extend(hitlist) 
    if hitlist:
        # Parcours de chaque paire de collision dans la hitlist
        for collision_pair in hitlist:
            # Ajoute la valeur de itr au moment de la collision à la liste des itr_collisions
            itr_collisions.append(itr)
            
    if hitlist:
    # Parcours de chaque paire de collision dans la hitlist
        for collision_pair in hitlist:
        # Mise à jour de la quantité de mouvement des électrons après la collision inélastique
            electron_index = collision_pair[0] 
            temperature_at_collision = Temp[itr_collisions[hitlist.index(collision_pair)]]# Indice de l'électron impliqué dans la collision
            new_momentum = update_momentum(temperature_at_collision)  # Génère une nouvelle quantité de mouvement pour l'électron
            p[electron_index] = new_momentum  # Met à jour la quantité de mouvement de l'électron
            apos[electron_index] += new_momentum * dt / mass
    

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