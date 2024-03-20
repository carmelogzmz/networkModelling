# network_builder function

# INFO
# Author: Carmelo Gomez-Martinez
# Last update: 21/03/2024
# contact: carmelogzmz@gmail.com
# github: https://github.com/carmelogzmz/networkModelling/

# DESCRIPTION
# This function build a pair of networks with a defined structure according to the following parameters:
# sp_num: number of species in each trophic level (if 100, each network will have 100 rows and 100 columns)
# sp_unique: number of species unique of each network. if sp_num = 100 and sp_unique = 50, 150 species will be shared
#            between network, and 50 will be unique of each network.
# sp_inis: number of species that occur in each community not interacting in its network but in the other. If sp_inis = 10,
#          in community at time 1 will be 10 species that interact in network 2 but not in network 1.
# interactions: number of interactions (qualitative) in each network. If 1000, a thousand cells will be filled with "1".
# perc_conserved: the percentage of interactions that will be conserved (shared) between networks. If 50, the 50% of
#                 the interactions in network 1 will be also in network 2.
# perc_rewiring: from the percentage of interactions changes (100-perc_conserved), perc_rewiring is the percentage of
#                interactions changes due to rewiring. If perc_conserved = 50 and perc_rewiring = 80, the 80% of the
#                50% of the total interactions will be unique of each network but involving shared species. The rest
#                (20% in this example), will be interactions involving unique species of each network.

# OUTPUT
# The function return a 4-items list, with items 1 and 2 being network at time 1 and at time 2, respectively, and
# items 3 and 4 the list of the species in the community at time 1 and at time 2, respectively.

# CODE
network_builder <- function(sp_num = 100, sp_unique = 50, sp_inis = 10, interactions = 1000, perc_conserved = 50, perc_rewiring = 80) {

  build_comm <- function(sp_num = 100, trophic, start = 1) {
    if(sp_num < 2 || sp_num > 999) {
      stop("sp_num must be a number between 2 and 999")
    }
    else {
      species <- NULL
      for(i in 1:sp_num) {
        if(i+start-1 < 10) {
          species <- c(species, paste0(trophic,"_00",i+start-1))
        }
        else if(i+start-1 >= 10 && i+start-1 < 100) {
          species <- c(species, paste0(trophic,"_0",i+start-1))
        }
        else {
          species <- c(species, paste0(trophic,"_",i+start-1))
        }
      }
      return(species)
    }
  }

  con <- round((interactions*perc_conserved)/100) # conserved interactions
  cha <- interactions-con # interaction changes
  rew <- round((cha*perc_rewiring)/100) # interactions changes due to rewiring
  tur <- cha-rew # interaction changes due to turnover
  sha <- con+rew # interactions in the shared subnetwork (conserved and rewired)
  sp_shared <- sp_num - round(sp_unique/2, 0)
  cells_shared <- sp_shared*sp_shared
  cells_rew <- (cells_shared - con)/2

  if(tur < round(sp_unique/2, 0)) {
    stop("Not enough interactions due to species turnover. Decrease the number of unique species or increase the number of interactions")
  } else if(cells_rew < rew) {
    stop("Not enough cells for rewiring efficiently. Decrease sp_unique or interactions")
  } else {
    print("Correct parameters. Networks can be built")
  }

  inturnover <- round(sp_unique/2, 0)
  plturnover <- sp_unique-inturnover

  plants_n1 <- build_comm(sp_num, "plant", 1)
  insects_n1 <- build_comm(sp_num, "insect", 1)

  plants_n2 <- plants_n1
  plants_n2 <- plants_n2[1:(length(plants_n2) - plturnover)]
  shared_pl <- plants_n2
  plants_n2 <- c(plants_n2, build_comm(plturnover, "plant", sp_num+1))

  insects_n2 <- insects_n1
  insects_n2 <- insects_n2[1:(length(insects_n2) - inturnover)]
  shared_in <- insects_n2
  insects_n2 <- c(insects_n2, build_comm(inturnover, "insect", sp_num+1))

  # Creating the shared subnetwork 1 and fill it with the number of interactions designated by "con" (% of conserved interactions)
  network_1 <- matrix(0, nrow = length(shared_pl), ncol = length(shared_in), dimnames = list(shared_pl,shared_in))
  x <- con
  while(x != 0) {
    i <- sample(seq_len(nrow(network_1)), 1)
    j <- sample(seq_len(ncol(network_1)), 1)
    if(network_1[i,j] == 0) {
      network_1[i,j] <- 1
      x <- x - 1
    }
  }

  # Copying network_1 to network_2 (conserved interactions)
  network_2 <- network_1

  # Adding to network_1 the interactions designated by "rew" (% of interactions due to rewiring)
  continue <- "NO"
  while(continue == "NO") {
    net <- network_1
    x <- rew
    for(i in 1:nrow(net)) {
      if(rowSums(net)[i] == 0) {
        j <- sample(seq_len(ncol(net)), 1)
        net[i,j] <- 1
        x <- x - 1
      }
    }

    for(j in 1:ncol(net)) {
      if(colSums(net)[j] == 0) {
        i <- sample(seq_len(nrow(net)), 1)
        net[i,j] <- 1
        x <- x - 1
      }
    }

    while(x > 0) {
      i <- sample(seq_len(nrow(net)), 1)
      j <- sample(seq_len(ncol(net)), 1)
      if(net[i,j] == 0) {
        net[i,j] <- 1
        x <- x - 1
      }
    }

    if(!0 %in% colSums(net) && !0 %in% rowSums(net) && sum(net) == sha) {
      continue <- "YES"
      network_1 <- net
    }
  }

  # Adding to network_2 the interactions designated by "rew" (% of interactions due to rewiring)
  # First fill interactions in the species not interacting in network_1, after that fill them randomly.
  continue <- "NO"
  while(continue == "NO") {
    net <- network_2
    x <- rew
    for(i in 1:nrow(net)) {
      if(rowSums(net)[i] == 0) {
        next_i <- "NO"
        while(next_i == "NO") {
          j <- sample(seq_len(ncol(net)), 1)
          if(net[i,j] == 0 && network_1[i,j] == 0) {
            net[i,j] <- 1
            x <- x - 1
            next_i <- "YES"
          }
        }
      }
    }

    for(j in 1:ncol(net)) {
      if(colSums(net)[j] == 0) {
        next_j <- "NO"
        while(next_j == "NO") {
          i <- sample(seq_len(nrow(net)), 1)
          if(net[i,j] == 0 && network_1[i,j] == 0) {
            net[i,j] <- 1
            x <- x - 1
            next_j <- "YES"
          }
        }
      }
    }

    while(x != 0) {
      i <- sample(seq_len(nrow(net)), 1)
      j <- sample(seq_len(ncol(net)), 1)
      if(net[i,j] == 0 && network_1[i,j] == 0) {
        net[i,j] <- 1
        x <- x - 1
      }
    }

    if(!0 %in% colSums(net) && !0 %in% rowSums(net) && sum(net) == sha) {
      continue <- "YES"
      network_2 <- net
    }
  }

  shared_rew <- sum(abs(network_2-network_1))
  shared_total <- sum(sum(network_1), sum(network_2))
  if(shared_rew == (sha+sha-con-con)) { print("Well built networks") }
  # extraer las especies únicas de cada red y nivel trofico
  unique_pl1 <- plants_n1[(length(plants_n1)-plturnover+1):length(plants_n1)]
  unique_in1 <- insects_n1[(length(insects_n1)-inturnover+1):length(insects_n1)]
  unique_pl2 <- plants_n2[(length(plants_n2)-plturnover+1):length(plants_n2)]
  unique_in2 <- insects_n2[(length(insects_n2)-inturnover+1):length(insects_n2)]

  # Adding the unique species of network_1 and the interactions designated by "tur" (percentage of interactions due to species turnover)
  if(length(unique_pl1) != 0) {
    x <- rep(0, ncol(network_1))
    for(i in seq_along(unique_pl1)) {
      network_1 <- rbind(network_1, x)
    }
    filas <- rownames(network_1)
    p <- length(unique_pl1)-1
    for(i in seq_along(unique_pl1)) {
      filas[length(filas)-p] <- unique_pl1[i]
      p <- p - 1
    }
    rownames(network_1) <- filas
  }

  if(length(unique_in1) != 0) {
    x <- rep(0, nrow(network_1))
    for(i in seq_along(unique_in1)) {
      network_1 <- cbind(network_1, x)
    }
    columnas <- colnames(network_1)
    p <- length(unique_in1)-1
    for(i in seq_along(unique_in1)) {
      columnas[length(columnas)-p] <- unique_in1[i]
      p <- p - 1
    }
    colnames(network_1) <- columnas
  }

  res <- tur
  size <- nrow(network_1)-length(unique_pl1)
  if(size < nrow(network_1)) {
    size1 <- ((size+1):nrow(network_1))
    for(k in size1) {
      i <- k
      continue <- "NO"
      while(continue == "NO") {
        j <- sample((ncol(network_1)-length(unique_in1)+1):ncol(network_1), 1)
        if(sum(network_1[,j]) == 0) {
          network_1[i,j] <- 1
          continue <- "YES"
          res <- res - 1
        }
      }
    }
  } else { size1 <- NULL }

  size <- ncol(network_1)-length(unique_in1)
  colsum <- colSums(network_1)[(size+1):100]
  pos <- NULL
  if(0 %in% colsum) {
    pos <- which(colsum == 0)
  }
  if(size < ncol(network_1)) {
    size2 <- ((size+1):ncol(network_1))
    if(!is.null(pos)) {
      continue <- "NO"
      while(res > 0 || continue == "NO") {
        for(k in pos) {
          j <- k
          l <- 1
          while(l == 1) {
            i <- sample(seq_len(nrow(network_1)), 1)
            if(network_1[i,j] == 0){
              network_1[i,j] <- 1
              l <- 0
              res <- res - 1
            }
          }
          if(res == 0) {
            break
          }
        }
        continue <- "YES"
      }
    }
    continue <- "NO"
    while(res > 0 || continue == "NO") {
      for(k in size2) {
        j <- sample(size2, 1)
        l <- 1
        while(l == 1) {
          i <- sample(seq_len(nrow(network_1)), 1)
          if(network_1[i,j] == 0){
            network_1[i,j] <- 1
            l <- 0
            res <- res - 1
          }
        }
        if(res == 0) {
          break
        }
        continue <- "YES"
      }
    }
  } else { size2 <- NULL }

  # Adding the unique species of network_2 and the interactions designated by "tur" (percentage of interactions due to species turnover)
  if(length(unique_pl2) != 0) {
    x <- rep(0, ncol(network_2))
    for(i in seq_along(unique_pl2)) {
      network_2 <- rbind(network_2, x)
    }
    filas <- rownames(network_2)
    p <- length(unique_pl2)-1
    for(i in seq_along(unique_pl2)) {
      filas[length(filas)-p] <- unique_pl2[i]
      p <- p - 1
    }
    rownames(network_2) <- filas
  }

  if(length(unique_in2) != 0) {
    x <- rep(0, nrow(network_2))
    for(i in seq_along(unique_in2)) {
      network_2 <- cbind(network_2, x)
    }
    columnas <- colnames(network_2)
    p <- length(unique_in2)-1
    for(i in seq_along(unique_in2)) {
      columnas[length(columnas)-p] <- unique_in2[i]
      p <- p - 1
    }
    colnames(network_2) <- columnas
  }

  res <- tur
  size <- nrow(network_2)-length(unique_pl2)
  if(size < nrow(network_2)) {
    size1 <- ((size+1):nrow(network_2))
    for(k in size1) {
      i <- k
      continue <- "NO"
      while(continue == "NO") {
        j <- sample((ncol(network_2)-length(unique_in2)+1):ncol(network_2), 1)
        if(sum(network_2[,j]) == 0) {
          network_2[i,j] <- 1
          continue <- "YES"
          res <- res - 1
        }
      }
    }
  } else { size1 <- NULL }

  size <- ncol(network_2)-length(unique_in2)
  colsum <- colSums(network_2)[(size+1):100]
  pos <- NULL
  if(0 %in% colsum) {
    pos <- which(colsum == 0)
  }
  if(size < ncol(network_2)) {
    size2 <- ((size+1):ncol(network_2))
    if(!is.null(pos)) {
      continue <- "NO"
      while(res > 0 || continue == "NO") {
        for(k in pos) {
          j <- k
          l <- 1
          while(l == 1) {
            i <- sample(seq_len(nrow(network_2)), 1)
            if(network_2[i,j] == 0){
              network_2[i,j] <- 1
              l <- 0
              res <- res - 1
            }
          }
          if(res == 0) {
            break
          }
        }
        continue <- "YES"
      }
    }
    continue <- "NO"
    while(res > 0 || continue == "NO") {
      for(k in size2) {
        j <- sample(size2, 1)
        l <- 1
        while(l == 1) {
          i <- sample(seq_len(nrow(network_2)), 1)
          if(network_2[i,j] == 0){
            network_2[i,j] <- 1
            l <- 0
            res <- res - 1
          }
        }
        if(res == 0) {
          break
        }
        continue <- "YES"
      }
    }
  } else { size2 <- NULL }
  dim(network_2); sum(network_2)

  # Creating lists of species in the community at time 1 and 2 (c1 and c2). c1 will be composed by the species in
  # network_1 plus a number of the unique species of network_2 defined by sp_inis.
  if(sp_inis == 0) {
    c1 <- list(rownames(network_1), colnames(network_1))
    c2 <- list(rownames(network_2), colnames(network_2))
  } else {
    pl1 <- rownames(network_1)
    po1 <- colnames(network_1)

    pl2 <- rownames(network_2)
    po2 <- colnames(network_2)

    plu1 <- pl1[!pl1 %in% pl2]
    plu2 <- pl2[!pl2 %in% pl1]
    pl1 <- c(pl1, plu2[1:sp_inis])
    pl2 <- c(pl2, plu1[1:sp_inis])

    pou1 <- po1[!po1 %in% po2]
    pou2 <- po2[!po2 %in% po1]
    po1 <- c(po1, pou2[1:sp_inis])
    po2 <- c(po2, pou1[1:sp_inis])

    c1 <- list(pl1, po1)
    c2 <- list(pl2, po2)
  }

  networks <- list(network_1, network_2, c1, c2)

  corrected_shared <- sha+sha
  corrected_rew <- rew+rew
  corrected_total <- sha+sha+tur+tur

  print(paste0("Total number of interaccions in both shared networks should be: ", corrected_shared, " interaccions in total"))
  print(paste0("Total number of interactions in both shared networks is: ", shared_total, " interactions in total"))
  print(paste0("Total number of rewired interactions should be: ", corrected_rew, " interactions in total"))
  print(paste0("Total number of rewired interactions is: ", shared_rew, " interactions in total"))
  print(paste0("Total number of interaccions in both total networks should be: ", corrected_total, " interactions in total"))
  print(paste0("Total number of interaccions in both total networks is: ", sum(sum(network_1),sum(network_2)), " interactions in total"))
  print(paste0("Total number of interactions in network 1: ", sum(network_1)))
  print(paste0("Total number of interactions in network 2: ", sum(network_2)))
  print(paste0("Total number of species in 1: ", nrow(network_1), " plant species and ", ncol(network_1), " insect species."))
  print(paste0("Total number of species 2: ", nrow(network_2), " plant species and ", ncol(network_2), " insect species."))
  return(networks)
}

# END OF SCRIPT
