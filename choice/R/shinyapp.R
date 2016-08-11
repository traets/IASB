#' nominal design matrix
#'
#' Transforms an effect coded design matrix into the design matrix containing the real attributelevels.
#' This design can be used to present to respondents.
#' @param design An effect coded design matrix.
#' @param lvl_names A list containing the values of each level of each attribute.
#' @param n_alts Number of alternatives per choice set.
#' @return A desing matrix with presentable attributelevels.
#' @export
present<-function (design, lvl_names, levels, n_alts){


  levels<-numeric(length(lvl_names))
  for(i in 1:length(levels)){levels[i]<-length(lvl_names[[i]])}
  n_att<-length(lvl_names)
  design<-as.matrix(design)
  mat<-matrix(nrow = nrow(design), ncol = n_att)
  jumps<-levels - 1
  end<-cumsum(jumps)
  begin<-1+end-jumps

  for (col in 1:n_att){
    x<-as.matrix(design[, seq(begin[col], end[col], 1)])
    mat[,col]<- attr(uniquecombs(x),"index")
   }

  fac_des<-as.matrix(mat)

  for (i in 1:ncol(fac_des)){
    fac_des[ ,i]<-mapvalues(fac_des[,i], from = sort(unique(fac_des[,i])),
                            to = lvl_names[[i]][1:max(as.numeric(fac_des[,i]))])
    fac_des<-as.data.frame(fac_des)
  }

 #setnr and alts
 set<-rep(seq(1, nrow(design)/n_alts, 1), each=n_alts)
 alt<-rep(seq(1, n_alts, 1), nrow(design)/n_alts )
 fac_des<-cbind(fac_des, set=set, alt=alt)

 return(fac_des)

}

#' Transform responses
#'
#' Transforms input responses to binary response vector
#' @param resp String vector containing input responses
#' @param resp_options String vector containing all possible responses.
#' The response options should be specified in increasing order, starting with (if included), the neutral response.
#' @param n_alts Number of alternatives per choice set.
#' @param neutral Logical value indicating whether a neutral option is provided or not. Default = TRUE.
#' @return A binary response vector
#' @export
map_resp<-function(resp, resp_options, n_alts, neutral=T){

  map<-match(resp, resp_options)
  l<-list()

  for(i in 1:length(map)){

    if (neutral){
    l[[i]] <- rep(0, n_alts)
    l[[i]][map[i]-1]<-1
    }else{
      l[[i]] <- rep(0, n_alts)
      l[[i]][map[i]]<-1
    }
  }
  v<-unlist(l)
  return(v)
}

#' Load from dropbox
#'
#' Load design from dropbox map
#'@param inputDir A dropbox directory that contains the design.
#'@return the file in that directory (or concatenated files).
#'@export
 loaddrop <- function(dir) {

  filesInfo <- drop_dir(dir)
  filePaths <- filesInfo$path
  data <- lapply(filePaths, drop_read_csv, stringsAsFactors = FALSE)
  # Concatenate all data together into one data frame
  data <- do.call(rbind, data)
  return(data)
}


#' Save to dropbox
#'
#' Save design and responses to dropbox map
#' @param d Matrix containing the design.
#' @param Y A response vector.
#' @export
savedrop <- function(d, Y, dir, filename) {

  data<-rbind(d, Y)
  fileName <- sprintf("%s.csv", filename )
  filePath <- file.path(tempdir(), fileName)
  write.csv(data, filePath, row.names = FALSE, quote = TRUE)

  drop_upload(filePath, dest = dir)
}





