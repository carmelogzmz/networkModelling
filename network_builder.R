network_modelling_comp <- function(lower = NULL, upper = NULL, sp_num = 100, interactions = 240, perc_conserved = 40, perc_rewiring = 80, sp_turnover = 10) {

  con <- round((interactions*perc_conserved)/100) # interacciones conservadas
  cha <- interactions-con # itneraccioens que cambian
  rew <- round((cha*perc_rewiring)/100) # interacciones que cambian por rewiring
  tur <- cha-rew # interacciones que cambian por turnover
  sha <- con+rew # interacciones en la red compartida (conservadas y rewiring)

  # comprobando que la red se puede construir
  int_tur <- sp_turnover * sp_turnover # número de celdillas que comparten las especies no compartidas entre ellas
  tot <- 2*(sp_turnover * sp_num) # número de celdillas que involucran a las especies no compartidas y el resto de la red
  resto <- tot-int_tur # número de celdillas que involucran a las especies no compartidas y el resto de la red, corregido

  if(resto < tur) {
    stop("Hay más interacciones que celdas posibles en la red no compartida. Disminuye el porcentaje de interaccioens que se deben a species turnover o aumenta el número de sp_turnover")
  }
  if(tur < sp_turnover*2) {
    stop("El porcentaje de interaction turnover debido a species turnover es demasiado bajo de acuerdo al número de interacciones total de la red, debes aumentar el numero de interacciones o disminuir el numero de sp_turnover")
  }

  # aleatoriza el número de especies que cambiará de insectos y plantas de acuerdo al total de especies a cambiar definido por "sp_turnover"
  inturnover <- sample(0:sp_turnover, 1)
  if(inturnover == 0) { print("NOTE: No insect species turnover") }
  plturnover <- sp_turnover - inturnover
  if(plturnover == 0) { print("NOTE: No plant species turnover") }

  # crea la lista de especies de plantas e insectos que habrá en la red 1 a partir de las comunidades incluidas en "lower" y "upper"
  plants_n1 <- sort(sample(lower, sp_num, replace = FALSE)) # extraigo 20 especies al azar de la comunidad
  insects_n1 <- sort(sample(upper, sp_num, replace = FALSE)) # extraigo 20 especies al azar de la comunidad

  # crea la lista de especies de plantas que habrá en la red 2 a partir de la lista de plantas de la red 1, eliminando aleatoriamente el numero de plantas definido por plturnover
  plants_n2 <- plants_n1
  if(plturnover != 0) {
    plextract <- sample(seq_along(plants_n2), plturnover, replace = FALSE)
    plants_n2 <- plants_n2[-plextract]
  } else {
    print("plturnover es 0")
  }
  while(plturnover != 0) {
    sp <- sample(lower, 1)
    if(!sp %in% plants_n2 && !sp %in% plants_n1) {
      plants_n2 <- c(plants_n2, sp)
      plturnover <- plturnover - 1
    }
  }

  # crea la lista de especies de insectos que habrá en la red 2 a partir de la lista de insectos de la red 1, eliminando aleatoriamente el numero de insects definido por inturnover
  insects_n2 <- insects_n1
  if(inturnover != 0) {
    inextract <- sample(seq_along(insects_n2), inturnover, replace = FALSE)
    insects_n2 <- insects_n2[-inextract]
  } else {
    print("inturnover es 0")
  }
  while(inturnover != 0) {
    sp <- sample(upper, 1)
    if(!sp %in% insects_n2 && !sp %in% insects_n1) {
      insects_n2 <- c(insects_n2, sp)
      inturnover <- inturnover - 1
    }
  }

  # Define las especies compartidas en ambas redes
  shared_pl <- plants_n1[plants_n1 %in% plants_n2]
  shared_in <- insects_n1[insects_n1 %in% insects_n2]

  # Crea la red 1 compartida y la rellena con el número de interacciones designado en "con" (que es el que corresponde el porcentaje determinado en "int_conserved")
  network_1 <- matrix(0, nrow = length(shared_pl), ncol = length(shared_in), dimnames = list(shared_pl,shared_in)) # monto la matriz
  x <- con
  while(x != 0) {
    i <- sample(seq_len(nrow(network_1)), 1)
    j <- sample(seq_len(ncol(network_1)), 1)
    if(network_1[i,j] == 0) {
      network_1[i,j] <- 1
      x <- x - 1
    }
  }
  sum(network_1)
  # Copio en la red 2 para que las compartan
  network_2 <- network_1

  # Crea la red 1 compartida y la rellena con el número de interacciones definidas por "rew"

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
    print(x)
    for(j in 1:ncol(net)) {
      if(colSums(net)[j] == 0) {
        i <- sample(seq_len(nrow(net)), 1)
        net[i,j] <- 1
        x <- x - 1
      }
    }
    print(x)
    while(x != 0) {
      i <- sample(seq_len(nrow(net)), 1)
      j <- sample(seq_len(ncol(net)), 1)
      if(net[i,j] == 0) {
        net[i,j] <- 1
        x <- x - 1
      }
    }
    print(x)
    if(!0 %in% colSums(net) && !0 %in% rowSums(net) && sum(net) == sha) {
      continue <- "YES"
      network_1 <- net
      print(paste0("Interacciones de la red 1 compartida: ", sum(network_1)))
    }
  }

  # Crea la red 2 compartida y la rellena con el número de interacciones definidas por "rew". Primero rellena
  # interacciones en las especies que no interactuaban en la red 1, y después coloca las demás de forma aleatoria.
  # Estos dos bucles for rellenan interacciones en las especies que no interactuaban en la red 1 (primer bucle plantas, segundo bucle insectos)

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
    print(x)
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
    print(x)
    while(x != 0) {
      i <- sample(seq_len(nrow(net)), 1)
      j <- sample(seq_len(ncol(net)), 1)
      if(net[i,j] == 0 && network_1[i,j] == 0) {
        net[i,j] <- 1
        x <- x - 1
      }
    }
    print(x)
    if(!0 %in% colSums(net) && !0 %in% rowSums(net) && sum(net) == sha) {
      continue <- "YES"
      network_2 <- net
      print(paste0("Interacciones de la red 2 compartida: ", sum(network_2)))
    }
  }

  shared_rew <- sum(abs(network_2-network_1))
  shared_total <- sum(sum(network_1), sum(network_2))
  if(shared_rew == (sha+sha-con-con)) { print("Redes compartidas bien construidas") }
  # extraer las especies únicas de cada red y nivel trofico
  unique_pl1 <- plants_n1[!plants_n1 %in% plants_n2]
  unique_in1 <- insects_n1[!insects_n1 %in% insects_n2]
  unique_pl2 <- plants_n2[!plants_n2 %in% plants_n1]
  unique_in2 <- insects_n2[!insects_n2 %in% insects_n1]

  dim(network_1); sum(network_1)
  dim(network_2); sum(network_2)

  # Añado las especies nuevas en la red 1 y relleno tantas interacciones como marca el valor "tur" (int_turnover)
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

  size <- nrow(network_1)-length(unique_pl1)
  if(size < nrow(network_1)) {
    size1 <- ((size+1):nrow(network_1))
    for(k in size1) {
      i <- k
      j <- sample(seq_len(ncol(network_1)), 1)
      network_1[i,j] <- 1
    }
  } else { size1 <- NULL }

  size <- ncol(network_1)-length(unique_in1)
  if(size < ncol(network_1)) {
    size2 <- ((size+1):ncol(network_1))
    for(k in size2) {
      j <- k
      l <- 1
      while(l == 1) {
        i <- sample(seq_len(nrow(network_1)), 1)
        if(network_1[i,j] == 0){
          network_1[i,j] <- 1
          l <- 0
        }
      }
    }
  } else { size2 <- NULL }
  dim(network_1); sum(network_1)

  x <- tur-length(size1)-length(size2)
  while(x != 0) {
    sel <- sample(c("i","j"), 1)
    i <- NULL
    j <- NULL
    if(sel == "i" && length(size1) > 0) {
      w <- nrow(network_1)-length(unique_pl1)
      if(w < nrow(network_1)) { w <- w+1 }
      i <- as.numeric(sample(as.character(c((w):nrow(network_1))), 1))
      #z <- ncol(network_1)-length(unique_in1)
      j <- sample(seq_len(ncol(network_1)), 1)
    } else if(sel == "j" && length(size2) > 0) {
      #w <- nrow(network_1)-length(unique_pl1)
      i <- sample(seq_len(nrow(network_1)), 1)
      z <- ncol(network_1)-length(unique_in1)
      if(z < ncol(network_1)) { z <- z+1 }
      j <- as.numeric(sample(as.character(c((z):ncol(network_1))), 1))
    }
    if(!is.null(i) && !is.null(j)) {
      if(network_1[i,j] == 0) {
        network_1[i,j] <- 1
        x <- x - 1
      }
    }
  }
  dim(network_1); sum(network_1)

  # Añado las especies nuevas en la red 2 y relleno tantas interacciones como marca el valor "tur" (int_turnover)
  if(!is.null(unique_pl2)) {
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

  if(!is.null(unique_in2)) {
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

  rownames(network_1) == rownames(network_2)
  colnames(network_1) == colnames(network_2)

  size <- nrow(network_2)-length(unique_pl2)
  if(size < nrow(network_2)) {
    size1 <- ((size+1):nrow(network_2))
    for(k in size1) {
      i <- k
      j <- sample(seq_len(ncol(network_2)), 1)
      network_2[i,j] <- 1
    }
  } else { size1 <- NULL }

  size <- ncol(network_2)-length(unique_in2)
  if(size < ncol(network_2)) {
    size2 <- ((size+1):ncol(network_2))
    for(k in size2) {
      j <- k
      l <- 1
      while(l == 1) {
        i <- sample(seq_len(nrow(network_2)), 1)
        if(network_2[i,j] == 0){
          network_2[i,j] <- 1
          l <- 0
        }
      }
    }
  } else { size2 <- NULL }
  dim(network_2); sum(network_2)

  x <- tur-length(size1)-length(size2)
  while(x != 0) {
    sel <- sample(c("i","j"), 1)
    i <- NULL
    j <- NULL
    if(sel == "i" && length(size1) > 0) {
      w <- nrow(network_2)-length(unique_pl2)
      if(w < nrow(network_2)) { w <- w+1 }
      i <- as.numeric(sample(as.character(c((w):nrow(network_2))), 1))
      #z <- ncol(network_2)-length(unique_in2)
      j <- sample(seq_len(ncol(network_2)), 1)
    } else if(sel == "j" && length(size2) > 0) {
      #w <- nrow(network_2)-length(unique_pl2)
      i <- sample(seq_len(nrow(network_2)), 1)
      z <- ncol(network_2)-length(unique_in2)
      if(z < ncol(network_2)) { z <- z+1 }
      j <- as.numeric(sample(as.character(c((z):ncol(network_2))), 1))
    }
    if(!is.null(i) && !is.null(j)) {
      if(network_2[i,j] == 0) {
        network_2[i,j] <- 1
        x <- x - 1
      }
    }
  }
  dim(network_2); sum(network_2)

  tot_net <- sum(abs(network_2-network_1))
  c1 <- list(rownames(network_1), colnames(network_1))
  c2 <- list(rownames(network_2), colnames(network_2))
  networks <- list(network_1, network_2, c1, c2)

  corrected_shared <- sha+sha
  corrected_rew <- rew+rew
  corrected_total <- sha+sha+tur+tur

  print(paste0("Total number of interaccions in both shared networks should be: ", corrected_shared, " interaccions in total"))
  print(paste0("Total number of interactions in both shared networks is: ", shared_total, " interactions in total"))
  print(paste0("Total number of rewired interactions should be: ", corrected_rew, " interactions in total"))
  print(paste0("Total number of rewired interactions should is: ", shared_rew, " interactions in total"))
  print(paste0("Total number of interaccions in both total networks should be: ", corrected_total, " interactions in total"))
  print(paste0("Total number of interaccions in both total networks is: ", sum(sum(network_1),sum(network_2)), " interactions in total"))
  print(paste0("Interacciones totales de la red 1: ", sum(network_1)))
  print(paste0("Interacciones totales de la red 2: ", sum(network_2)))
  print(paste0("Especies totales de la red 1: ", nrow(network_1), " especies de plantas y ", ncol(network_1), " especies de insectos."))
  print(paste0("Especies totales de la red 2: ", nrow(network_2), " especies de plantas y ", ncol(network_2), " especies de insectos."))
  return(networks)
}
