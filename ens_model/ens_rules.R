
##################################################
# Les differents Etats d'une agents ---> {1 : agent vivant , 0 : agent mort}
# ------------------@@@@@@@@@@@@@@---------------
# Deux architectures pour les déplacements :
#       1 . Von Neumann : 4 directions
#             1 ---> NORD (x,y+1)
#             2 ---> SUD  (x,y-1)
#             3 ---> EST  (x-1,y)
#             4 ---> OUEST (x+1,y)
#
#       2 . Moore   : 8 directions 
#            11 ---> NORD (x,y+1)
#            13 ---> NORD-EST (x-1,y+1)
#            14 ---> NORD-OUEST (x+1,y+1)
#            21 ---> SUD  (x,y-1)
#            23 ---> SUD-EST (x-1,y-1)
#            14 ---> SUD-OUEST (x+1,y-1)
#            31 ---> EST  (x-1,y)
#            41 ---> OUEST (x+1,y)
# -------------------@@@@@@@@@@@@@@@--------------
# Proliferation de type "Isotropic mitosis" (prob de migration egale)
#--------------------@@@@@@@@@@@@@@---------------
# Migration de type "Mediated by Contact Inhibition" (+ dense ---> - dense) avec décès si region cible occupée
# -------------------@@@@@@@@@@@@@-----------------
# Un agent peut :
#       1. Proliférer
#       2. mouvoir    
#       3. rester immobile
#--------------------@@@@@@@@@@@@------------------
#Les cellules filles se déplace en direction opposées 
#
#
#
####################################################

ens <-function(my_lattice,n_x,n_y,n_agents,ii,jj,k,n_inner)
{
    for (inner in 1:n_inner)
    {
        ##Mobilité
        num_agents = sample(1:n_agents,n_agents, replace=TRUE) #Tirage avec remises des cellules affectés

        for(k in 1:n_agents) #Pour chaque agents
        {
            id_agent = num_agents [k]   # récupérer le k-ième agent
            i_loc    = ii[id_agent]     # l'absisse actuelle de l'agent
            j_loc    = jj[id_agent]     # l'ordonnée actuelle de l'agent

            cible_i = 0;
            cible_j = 0;

            #Choisir une cellule cible parmis les quatres cellules voisine en config Von Nuemann
            pos_to_move_to = sample(1:4,size=1,Replace=TRUE)  # Choix de l'une des 4 diretion avec probabilité uniforme
            #1:droite; 2:gauche; 3:haut; 4:bas

            ## calculer la position de la cellule ciblée pour son mouvement
            if(pos_to_move_to == 1)
            {
                cible_i = ((i_loc) %% n_x) + 1
                cible_j = j_loc
            }

            else if (pos_to_move_to == 2) {
               cible_i = ((i_loc - 2) %% n_x) +1
               cible_j = j_loc
            }

            else if (pos_to_move_to == 3) {
                cible_i = i_loc
                cible_j = ((j_loc) %% n_y) + 1
            }

            else {
                cible_i = i_loc
                cible_j = ((j_loc - 2) %% n_y) +1
            }

            ## Récupérer les cordonnées des huits cellules voisines de la cellules cible en config Moore
            #Nord 
            N_j = ((cible_j) %% n_y) + 1
            N_i = cible_i
            #Sud
            S_j = ((cible_j - 2) %% n_y) + 1
            S_i = cible_i
            #Est
            E_j = cible_j
            E_i = ((cible_i - 2) %% n_x) + 1
            #Ouest
            O_j = cible_j
            O_i = ((cible_i) %% n_x) + 1
            #Nord_Est
            NE_j = ((cible_j) %% n_y) + 1
            NE_i = ((cible_i - 2) %% n_x) + 1
            #Nord_Ouest
            NO_j = ((cible_j) %% n_y) + 1
            NO_i = ((cible_i) %% n_x) + 1
            #Sud_Est
            SE_j = ((cible_j - 2) %% n_y) + 1
            SE_i = ((cible_i - 2) %% n_x) + 1
            #Sud_Ouest
            SO_j = ((cible_j - 2) %% n_y) + 1
            SO_i = ((cible_i) %% n_x) + 1

            ## Calculer M pour la cellule ciblé
            M_m = my_lattice[N_i,N_j] + my_lattice[S_i,S_j] + my_lattice[E_i,E_j] + my_lattice[O_i,O_j] +
             my_lattice[NE_i,NE_j] + my_lattice[NO_i,NO_j] + my_lattice[SE_i,SE_j] + my_lattice[SO_i,SO_j]
            
            ## faire le test pour voir si mouvements
            if(M_m <= k)
            {
                #Vérifier si bien que le mouvement est autoriser la cellule cible est occupé ou non
                if(my_lattice[cible_i,cible_j] == 0)
                {
                    i_nb = i_cible
                    j_nb = j_cible
                    my_lattice[i_nb,j_cible] = 1; my_lattice[i_loc,j_loc]=0;
                    ii[id_agent] = i_nb; jj[id_agent]=j_nb;
                }
            }
        }


        ## Proliferation
        num_agents = sample(1:n_agents,n_agents, replace=TRUE) #Tirage avec remises des cellules affectés
        for(k in 1:n_agents)
        {
            id_agent = num_agents [k]   # récupérer le k-ième agent
            i_loc    = ii[id_agent]     # l'absisse actuelle de l'agent
            j_loc    = jj[id_agent]     # l'ordonnée actuelle de l'agent

            ## Pour initier une mitose il faut que M >= k
            ## Récupérer les cordonnées des huits cellules voisines de la cellules cible en config Moore
            #Nord 
            N_j = ((j_loc) %% n_y) + 1
            N_i = i_loc
            #Sud
            S_j = ((j_loc - 2) %% n_y) + 1
            S_i = i_loc
            #Est
            E_j = j_loc
            E_i = ((i_loc) - 2) %% n_x) + 1
            #Ouest
            O_j = j_loc
            O_i = ((i_loc) %% n_x) + 1
            #Nord_Est
            NE_j = ((j_loc) %% n_y) + 1
            NE_i = ((i_loc - 2) %% n_x) + 1
            #Nord_Ouest
            NO_j = ((j_loc) %% n_y) + 1
            NO_i = ((i_loc) %% n_x) + 1
            #Sud_Est
            SE_j = ((j_loc - 2) %% n_y) + 1
            SE_i = ((i_loc - 2) %% n_x) + 1
            #Sud_Ouest
            SO_j = ((j_loc - 2) %% n_y) + 1
            SO_i = ((i_loc) %% n_x) + 1

            ## Calculer M pour la cellule ciblé
            M_p = my_lattice[N_i,N_j] + my_lattice[S_i,S_j] + my_lattice[E_i,E_j] + my_lattice[O_i,O_j] +
             my_lattice[NE_i,NE_j] + my_lattice[NO_i,NO_j] + my_lattice[SE_i,SE_j] + my_lattice[SO_i,SO_j]

            #Faire le test pour voir si on peut ou non initier la mitose
            if(M_p <= p) 
            {
                # Direction du mouvement des cellules filles
                # 1:Nord et Sud ; 2:Est et Ouest ; 3:Nord-Est et Sud-Ouest 4:Nord-Ouest et Sud-Est
                dir_to_move_to = sample(1:4,size=1,Replace=TRUE,prob=1/4) # Choix de l'une des 4 diretion avec probabilité uniforme

                fille1_i = i_loc; fille1_j = j_loc; fille2_i = i_loc; fille2_j = j_loc
                if(dir_to_move_to == 1)
                {
                    fille1_i = N_i; fille1_j = N_j
                    fille2_i = S_i; fille2_j = S_j
                }
                else if (dir_to_move_to == 2) {
                   fille1_i = E_i; fille1_j = E_j;fille2_i = O_i; fille2_j = O_j
                }
                else if(dir_to_move_to == 3) {
                    fille1_i = NE_i; fille1_j = NE_j; fille2_i = SO_i; fille2_j = SO_j
                }
                else {
                    fille1_i = NO_i; fille1_j = NO_j; fille2_i = SE_i; fille2_j = SE_j
                }

                ## déterminer si les cellules ciblées par les filles sont libres ou non
                if((my_lattice[fille1_i,fille1_j]==0) and (my_lattice[fille2_i,fille2_j]==0))
                {
                    my_lattice[i_loc,j_loc] = 1
                    my_lattice[fille1_i,fille1_j]=1; n_agents = n_agents + 1; ii[n_agents]=fille1_i; jj[n_agents]=fille1_j
                    my_lattice[fille2_i,fille2_j]=1; n_agents = n_agents + 1; ii[n_agents]=fille2_i; jj[n_agents]=fille2_j
                }
            }
        }
    }

    return (list(field = my_lattice, nbr = n_agents, i_coord = ii, j_coord = jj))
}