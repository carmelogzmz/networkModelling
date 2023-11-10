subsampling <- function(d1, d2) {

  require(bipartite)
  detectability <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

  totLevels <- list()
  for(i in 1:10) {
    det <- detectability[[i]]
    rep_level <- list()
    for(j in 1:100) {
      f1 <- sample(seq_len(nrow(d1)), (det*nrow(d1))/100, replace = FALSE)
      f1 <- sort(f1)
      dx1 <- d1[f1,]
      dx1$net <- "1"
      f2 <- sample(seq_len(nrow(d2)), (det*nrow(d2))/100, replace = FALSE)
      f2 <- sort(f2)
      dx2 <- d2[f2,]
      dx2$net <- "2"
      dx <- rbind(dx1, dx2)
      redes <- frame2webs(dx, varnames = c("lower", "upper", "net", "interaction"))
      rep_level[[j]] <- redes
    }
    totLevels[[i]] <- rep_level
  }
  return(totLevels)
}