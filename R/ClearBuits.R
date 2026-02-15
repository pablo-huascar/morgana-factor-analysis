

#########################################
###                                   ###
###             CLEARBUITS            ### 
###                                   ###
#########################################

ClearBuits <- function(B){
  B <- as.matrix(B, ncol=1)
  m <- nrow(B)
  h <- m 
  for(i in 1:m){
    if(B[i]== -1){
      h <- h-1
    }
  }
  Cc <- B[1:h, 1]
  return(Cc)
}
