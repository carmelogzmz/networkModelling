build_comm <- function(sp_num = 200, start = 1) {
  if(sp_num < 2 || sp_num > 999) {
    stop("sp_num must be a number between 2 and 999")
  }
  else {
    plants <- NULL
    insects <- NULL
    for(i in 1:sp_num) {
      if(i+start-1 < 10) {
        plants <- c(plants, paste0("plant_00",i+start-1))
        insects <- c(insects, paste0("insect_00",i+start-1))
      }
      else if(i+start-1 >= 10 && i+start-1 < 100) {
        plants <- c(plants, paste0("plant_0",i+start-1))
        insects <- c(insects, paste0("insect_0",i+start-1))
      }
      else {
        plants <- c(plants, paste0("plant_",i+start-1))
        insects <- c(insects, paste0("insect_",i+start-1))
      }
    }
    community <- list(plants, insects)
    return(community)
  }
}