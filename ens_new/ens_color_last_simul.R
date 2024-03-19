
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

electives_colors = c("black", "green", "yellow4", "red4", "cyan3", "coral","mediumblue", "darkviolet", "lightsalmon4","khaki1")

last_agents = matrix(data = NA, 10,100)   #Une matrice pour collecter les populations finales 
row.names(last_agents) = c("black","coral","cyan3","darkviolet","green","khaki1","lightsalmon4","mediumblue","red4","yellow4")
for(iter in 1:100)
{
## Initialisation des variables et conteneurs

my_lattice = matrix(data = 0,n_x , n_y) # matrice de 0 et de 1 décrivant l'état des agents
ii = rep(NA, n_tot)  #abscisses de toutes les cellules de la population
jj = rep(NA, n_tot) #ordonnees de toutes les cellules de la population
cel_colors= rep(NA, n_tot) #couleurs des cellules
cel_freq_tabs = matrix(data = "", n_tot,n_outer+1) # matrice des coleurs toutes n_outers iterations
################################################################################################

############################### Modele ########################################################

##Initialisation uniforme des agents actives et attribution des couleurs aux cellules souches

n_agents = 0

##


i =  rep(1,10)   # on décide pour l'initialisation de prendre uniquement les agents de la première ligne
j = sample(1:n_y,10,replace = FALSE) # ordonnées des agents à initialiser
my_lattice[i,j] = rep(1,10);         # mise à jour de la grille
ii[1:10] =  i                      # mise à jour du tableau des abscisses
jj[1:10] = j                       # mise à jour du tableau des ordonnées

# pos1 = sample(1:10,10,replace = FALSE)
# electives_colors = electives_colors[pos1] 
n_agents = 10



  cel_colors[1:10] = electives_colors
 
  # Population à la date initiale
  #population[1] <- n_agents
  
  # Couleurs des agents actifs à la date initiale
  
  cel_freq_tabs[,1] <- cel_colors
  
  # Visualisation de la configuration initiale de la grille
  
  plot(ii,jj,xlim = c(1,n_x), ylim = c(1, n_y), col = cel_colors,
       pch = 19, cex = 40/n_x, xlab = "x", ylab = "y", main = " t = 0 iteration")
  
  # Iterations

  for (i_outer in 1:n_outer){
    out = ens(my_lattice, n_x, n_y, n_agents, ii, jj, cel_colors, k, n_inner)
    my_lattice = out$field
    n_agents = out$nbr
    ii       = out$i_coord
    jj       = out$j_coord
    cel_colors = out$my_colors
    cel_freq_tabs[,i_outer+1] = cel_colors;
    
    plot(ii,jj,xlim = c(1, n_x), ylim = c(1, n_y), col = cel_colors,
         pch = 19, cex = 40/n_x, xlab = "x", ylab = "y",
         main = paste(toString((i_outer)*n_inner), " iterations"))
  }
  
  tmp_data1 = data.frame(table(cel_freq_tabs[,ncol(cel_freq_tabs)]))
  tmp_df1   = tmp_data1[order(tmp_data1[,1]),]
  last_agents[,iter] = tmp_df1$Freq
  
  rm(tmp_df1,tmp_data1,cel_colors,cel_freq_tabs,ii,jj,my_lattice,n_agents,i,j)
  
  
}


