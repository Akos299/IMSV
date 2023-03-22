
source("ens_rules.R", chdir = TRUE)


#################################### Définition de la grille ##################################

size_x = 100     #distance maximale à partir du centre en abscisse
size_y = 100     #distance maximale à partir du centre en ordonnée

## center coordonne
center_x = size_x + 1
center_y = size_y + 1

## Taille totale de la lattice
n_x = 2*size_x + 1
n_y = 2*size_y + 1
n_tot = n_x * n_y

##############################################################################################

################################## Paramètres ################################################

## Parameters of automata
k = sample(1:8,size=1)   # capacité proliferative utilisé dans l'article peut varier entre 1 et 8 dans l'article fixé à 3



## Initialisation

init_x = 5 
init_y = size_y

## Nombre d'itération (valeur observee tous les n_inner)
n_outer = 5
n_inner = 10

## Initialisation des variables et conteneurs

my_lattice = matrix(data = 0,n_x , n_y) # matrice de 0 et de 1 décrivant l'état des agents
profils = matrix(data = NA,n_x, n_outer+1) # nombre de agent sur chaque ligne
population = rep(NA,n_outer+1) # stocke les tailles de population toutes les n_outer iterations
ii = rep(NA, n_tot)  #abscisses de toutes les cellules de la population
jj = rep(NA, n_tot) #ordonnees de toutes les cellules de la population

################################################################################################

############################### Modele ########################################################

##Initialisation uniforme des agents actives

n_agents = 0
while(n_agents < 4*init_x*init_y)
{
    i = sample(1:(2*init_x),1, replace = FALSE)
    j = sample(1:(2*init_y),1, replace = FALSE)
    n_agents = n_agents + 1

    my_lattice[i,j] = 1
    ii[n_agents] = i
    jj[n_agents] = j
}


# Population à initiale
population[1] = n_agents

# Nombre d'agents sur chaque ligne de la grille à l'instant init

profils[,1] = rowSums(my_lattice)

# Visualisation de la configuration initiale de la grille
plot(ii,jj,xlim = c(1,n_x), ylim = c(1, n_y), col = "red4", bg ="red",
pch = 22, cex = 40/n_x, xlab = "x", ylab = "y", main = "Initial field")

# Iterations

for (i_outer in 1:n_outer){
    out = ens(my_lattice, n_x, n_y, n_agents, ii, jj, k, n_inner)
    my_lattice = out$field
    n_agents = out$nbr
    ii       = out$i_coord
    jj       = out$j_coord
    profils[, i_outer+1] = rowSums(my_lattice) # nombre de cellules sur chaque ligne de la matrice
    population[i_outer+1] = n_agents # nombre de cellules total
    plot(ii,jj,xlim = c(1, n_x), ylim = c(1, n_y), col = "red4", bg = "red",
    pch = 22, cex = 40/n_x, xlab = "x", ylab = "y",
    main = paste(toString((i_outer)*n_inner), " iterations"))
}

# Distribution de la densité au fil du temps
plot (profils[,1], type = "l", col = "black", xlab = "Position",
ylab = "Concentration", main = "Profils", lwd = 2)

grid()

for( i in 1:n_outer) {
    lines(profils[,i+1], col = i+1, lwd =2)
}

plot(population) # nombre de cellules (constant si pas de proliferation)

grid()

