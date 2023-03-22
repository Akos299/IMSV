
source("ens_rules.R", chdir = TRUE)


#################################### Définition de la grille ##################################

size_x = 25     #distance maximale à partir du centre en abscisse
size_y = 8     #distance maximale à partir du centre en ordonnée

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
#k = sample(1:8,size=1)   # capacité proliferative utilisé dans l'article peut varier entre 1 et 8 dans l'article fixé à 3

k = 3


## Initialisation

init_x =5 
init_y = size_y

## Nombre d'itération (valeur observee tous les n_inner)
n_outer = 10
n_inner = 20

## Initialisation des variables et conteneurs

my_lattice = matrix(data = 0,n_x , n_y) # matrice de 0 et de 1 décrivant l'état des agents
profils = matrix(data = NA,n_x, n_outer+1) # nombre de agent sur chaque ligne
population = rep(NA,n_outer+1) # stocke les tailles de population toutes les n_outer iterations
ii = rep(NA, n_tot)  #abscisses de toutes les cellules de la population
jj = rep(NA, n_tot) #ordonnees de toutes les cellules de la population
cel_colors= rep(NA, n_tot) #couleurs des cellules
cel_freq_tabs = matrix(data = "", n_tot,n_outer+1) # matrice des coleurs toutes n_outers iterations
################################################################################################

############################### Modele ########################################################

##Initialisation uniforme des agents actives et atribution des couleurs aux cellules souches

electives_colors = c("black", "green", "yellow4", "red4", "cyan3", "coral","mediumblue", "darkviolet", "lightsalmon4","khaki1")
n_agents = 0
#while(n_agents < 4*init_x*init_y)
##

i =  rep(1,10)   # on décide pour l'initialisation de prendre uniquement les agents de la première ligne
j = sample(1:n_y,10,replace = FALSE) # ordonnées des agents à initialiser
my_lattice[i,j] = rep(1,10);         # mise à jour de la grille
ii[1:10] =  i                      # mise à jour du tableau des abscisses
jj[1:10] = j                       # mise à jour du tableau des ordonnées

n_agents = 10
## Pour changer la configuration d'initialisation 
pos1 = sample(1:10,10,replace = FALSE)
electives_colors = electives_colors[pos1] 
cel_colors[1:10] = electives_colors
#col_freq 
# Population à la date initiale
population[1] <- n_agents

# Couleurs des agents actifs à la date initiale

cel_freq_tabs[,1] <- cel_colors

# Nombre d'agents sur chaque ligne de la grille à l'instant init

profils[,1] <- rowSums(my_lattice)

# Visualisation de la configuration initiale de la grille

plot(ii,jj,xlim = c(1,n_x), ylim = c(1, n_y), col = cel_colors,
     pch = 19, cex = 40/n_x, xlab = "x", ylab = "y", main = " t = 0 iteration")

# Iterations
agents_evol = matrix(data=NA,10,n_outer+1)     # Evolution des agents toutes i_outer iteration
tmp_data1 = data.frame(table(cel_freq_tabs[,1]))
tmp_df1   = tmp_data1[order(tmp_data1[,1]),]
agents_evol[,1] = tmp_df1$Freq
for (i_outer in 1:n_outer){
  out = ens(my_lattice, n_x, n_y, n_agents, ii, jj, cel_colors, k, n_inner)
  my_lattice = out$field
  n_agents = out$nbr
  ii       = out$i_coord
  jj       = out$j_coord
  cel_colors = out$my_colors
  cel_freq_tabs[,i_outer+1] = cel_colors;
  profils[, i_outer+1] = rowSums(my_lattice) # nombre de cellules sur chaque ligne de la matrice
  population[i_outer+1] = n_agents # nombre de cellules total
 
  plot(ii,jj,xlim = c(1, n_x), ylim = c(1, n_y), col = cel_colors,
       pch = 19, cex = 40/n_x, xlab = "x", ylab = "y",
       main = paste(toString((i_outer)*n_inner), " iterations"))
  
  tmp_data1 = data.frame(table(cel_freq_tabs[,i_outer+1]))
  tmp_df1   = tmp_data1[order(tmp_data1[,1]),]
  agents_evol[,i_outer+1] = tmp_df1$Freq
}

row.names(agents_evol) = c("black","coral","cyan3","darkviolet","green","khaki1","lightsalmon4","mediumblue","red4","yellow4")

# # Distribution de la densité au fil du temps
 plot (profils[,1], type = "l", col = "black", xlab = "Position",
       ylab = "Concentration", main = "Profils", lwd = 2)

grid()

for( i in 1:n_outer) {
  lines(profils[,i+1], col = i+1, lwd =2)
}

plot(population) # nombre de cellules (constant si pas de proliferation)

grid()

