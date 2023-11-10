webs2frames <- function(n1, n2) {
  row1 <- nrow(n1)
  col1 <- ncol(n1)
  row2 <- nrow(n2)
  col2 <- ncol(n2)
  l1 <- NULL; u1 <- NULL
  l2 <- NULL; u2 <- NULL
  for(i in 1:row1) {
    for(j in 1:col1) {
      if(n1[i,j] == 1) {
        l1 <- c(l1, rownames(n1)[i])
        u1 <- c(u1, colnames(n1)[j])
      }
    }
  }
  print(paste0(length(l1),"-",length(u1)))
  for(i in 1:row2) {
    for(j in 1:col2) {
      if(n2[i,j] == 1) {
        l2 <- c(l2, rownames(n2)[i])
        u2 <- c(u2, colnames(n2)[j])
      }
    }
  }
  print(paste0(length(l2),"-",length(u2)))
  d1 <- data.frame("lower" = l1, "upper" = u1, "interaction" = 1)
  d2 <- data.frame("lower" = l2, "upper" = u2, "interaction" = 1)
  d <- list(d1, d2)
  names(d) <- c("n1", "n2")
  return(d)
}