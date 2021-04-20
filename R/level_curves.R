#' Create level curves for copulae and other models
#'
#' Take data on uniform margins and fit the Ledford and Tawn (1997) joint tail model. Also contains the method where additional information from values that are extreme in at most one variable is used.
#'
#' @param pxf Uniform values of the first margin with a mixed distribution (empirical below and gpd above a threshold)
#' @param pyf Uniform values of the second margin with a mixed distribution (empirical below and gpd above a threshold)
#' @param mar1 Values of the first margin
#' @param mar2 Values of the second margin
#' @param pos part of the curve to be modelled '\code{l}' for the left part, '\code{m}' for the moddle part and '\code{r}' for the right part
#' @param pobje probability of the level curve to be modeled
#' @param ng number of points to be interpolated
#' @param inter type of hazard interrelation '\code{comb}' for compound and '\code{casc}' for cascade,
#' @param coco a copula function from the following: GHcop,NORMcop,FGMcop,GLcop
#' @param c1 the parameter of the copula
#' @importFrom copBasic surfuncCOP
#' @export
#' @return Level curves for a given joint probability
#'

curve.funct<-function(pxf,pyf,mar1,mar2,pos,pobje,ng=100,inter="comb",coco,c1){
  if(pos=="l"){
    xmin=0
    xmax=0.9
    ymin=0.99
    ymax=1
  }
  if(pos=="m"){
    xmin=0.9
    xmax=0.995
    ymin=0.9
    ymax=.995
  }
  if(pos=="r"){
    xmin=0.99
    xmax=1
    ymin=0
    ymax=.9
  }
  godx<-spline(pxf,mar1, n = ng, method = "fmm",
               xmin = xmin, xmax = xmax, ties = mean)

  gody<-spline(pyf,mar2, n = ng, method = "fmm",
               xmin = ymin, xmax = ymax, ties = mean)

  coxi<-approx(godx$y,godx$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)

  coyi<-approx(gody$y,gody$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)
  plot(coxi$x,coxi$y)
  godx<-coxi$y
  gody<-coyi$y
  coxi<-coxi$x
  coyi<-coyi$x

  repeat{
    idx<-which(diff(coxi)<=0)
    if(length(idx)<1){
      break
    }
    coxi[idx+1]=coxi[idx]+0.001

  }
  repeat{
    idy<-which(diff(coyi)<=0)
    if(length(idy)<1){
      break
    }
    coyi[idy+1]=coyi[idy]+0.001

  }


  acp3<-matrix(NA, nrow = ng, ncol = ng)
  for (k in 1:length(gody)){
    for (j in 1:length(godx)){
      if (inter=="comb"){
        acp3[j,k]=surfuncCOP(godx[j], gody[k], cop=coco, para=c1)}
      if (inter=="casc"){
        acp3[j,k]= surfuncCOP(godx[j], gody[k], cop=coco, para=c1)/(1-godx[j])}

    }
  }

  cl2<-contourLines(coxi,coyi, acp3, levels = pobje)
  if(length(cl2)>0){
    cl2<-as.matrix(data.frame(cl2[[1]]$x,cl2[[1]]$y))} else{cl2<-NA}

}

#' Create level curves for copulae and other models
#'
#' Take data on uniform margins and fit the Ledford and Tawn (1997) joint tail model. Also contains the method where additional information from values that are extreme in at most one variable is used.
#'
#' @param px Uniform values of the first margin with a mixed distribution (empirical below and gpd above a threshold)
#' @param py Uniform values of the second margin with a mixed distribution (empirical below and gpd above a threshold)
#' @param mar1 Values of the first margin
#' @param mar2 Values of the second margin
#' @param pos part of the curve to be modelled '\code{l}' for the left part, '\code{m}' for the moddle part and '\code{r}' for the right part
#' @param pobje probability of the level curve to be modelled
#' @param ng number of points to be interpolated
#' @param inter type of hazard interrelation '\code{comb}' for compound and '\code{casc}' for cascade,
#' @param logm log tranformation of the margins '\code{T}'of '\code{F}'
#' @param coco a copula function from the following: GHcop,NORMcop,FGMcop,GLcop
#' @param c1 the parameter of the copula
#' @export
#' @importFrom copBasic surfuncCOP
#' @return Level curves for a given joint probability

curve.funct.a<-function(px,py,mar1,mar2,pos,pobje,ng=100,inter="comb",logm=FALSE,coco,c1){
  if(logm==TRUE) {
    mar1<-log(mar1)
    mar2<-log(mar2)
  }
  xmin=0
  xmax=1
  ymin=0
  ymax=1
  tieef<-mean
  ngx=100000


  godx<-approx(px,mar1, n = ngx, method = "linear",
               yleft = min(mar1), yright = max(mar1), ties = tieef)

  gody<-approx(py,mar2, n = ngx, method = "linear",
               yleft = min(mar1), yright = max(mar1), ties = tieef)


  coxi<-approx(godx$y,godx$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)

  coyi<-approx(gody$y,gody$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)

  godx<-coxi$y
  gody<-coyi$y
  coxi<-coxi$x
  coyi<-coyi$x
  repeat{
    idx<-which(diff(coxi)<=0)
    if(length(idx)<1){
      break
    }
    coxi[idx+1]=coxi[idx]+0.001

  }
  repeat{
    idy<-which(diff(coyi)<=0)
    if(length(idy)<1){
      break
    }
    coyi[idy+1]=coyi[idy]+0.001

  }

  acp3<-matrix(NA, nrow = ng, ncol = ng)
  for (k in 1:length(gody)){
    for (j in 1:length(godx)){
      if (inter=="comb"){
        acp3[j,k]=surfuncCOP(godx[j], gody[k], cop=coco, para=c1)}
      if (inter=="casc"){
        acp3[j,k]= surfuncCOP(godx[j], gody[k], cop=coco, para=c1)/(1-godx[j])}

    }
  }

  clf<-contourLines(coxi,coyi, acp3, levels = pobje)
  if(is.list(clf)|length(clf)>1){
    clf<-as.matrix(data.frame(clf[[1]]$x,clf[[1]]$y))} else{clf<-"error"}
  if(logm==TRUE) clf<-exp(clf)
  return(clf)

}

#' Create level curves for copulae and other models
#'
#' Take data on uniform margins and fit the Ledford and Tawn (1997) joint tail model. Also contains the method where additional information from values that are extreme in at most one variable is used.
#'
#' @param pxf Uniform values of the first margin with a mixed distribution (empirical below and gpd above a threshold)
#' @param pyf Uniform values of the second margin with a mixed distribution (empirical below and gpd above a threshold)
#' @param mar1 Values of the first margin
#' @param mar2 Values of the second margin
#' @param pos part of the curve to be modelled '\code{l}' for the left part, '\code{m}' for the middle part and '\code{r}' for the right part
#' @param pobje probability of the level curve to be modelled
#' @param ng number of points to be interpolated
#' @param coco a copula function from the following: GHcop,NORMcop,FGMcop,GLcop
#' @param c1 the parameter of the copula
#' @param inter type of hazard interrelation '\code{comb}' for compound and '\code{casc}' for cascade,
#' @export
#' @return Level curves for a given joint probability
#' @importFrom grDevices contourLines
#' @importFrom copBasic surfuncCOP
#' @importFrom graphics abline contour segments

curve.funct.b<-function(pxf,pyf,mar1,mar2,pos,pobje,ng=100,inter="comb",coco,c1){
  if(pos=="l"){
    xmin=0
    xmax=0.99
    ymin=0.99
    ymax=1
  }
  if(pos=="m"){
    xmin=0.98
    xmax=1
    ymin=0.98
    ymax=1
  }
  if(pos=="r"){
    xmin=0.99
    xmax=1
    ymin=0
    ymax=0.99
  }
  ngx=10000
  godx<-spline(pxf,mar1, n = ngx, method = "natural",
               xmin = xmin, xmax = xmax, ties = mean)

  gody<-spline(pyf,mar2, n = ngx, method = "natural",
               xmin = ymin, xmax = ymax, ties = mean)

  coxi<-approx(godx$y,godx$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)

  coyi<-approx(gody$y,gody$x, n = ng, method = "linear",
               yleft = 0, yright = 1, ties = mean)
  plot(coxi$x,coxi$y)
  godx<-coxi$y
  gody<-coyi$y
  coxi<-coxi$x
  coyi<-coyi$x

  repeat{
    idx<-which(diff(coxi)<=0)
    if(length(idx)<1){
      break
    }
    coxi[idx+1]=coxi[idx]+0.001

  }
  repeat{
    idy<-which(diff(coyi)<=0)
    if(length(idy)<1){
      break
    }
    coyi[idy+1]=coyi[idy]+0.001

  }


  acp3<-matrix(NA, nrow = ng, ncol = ng)
  for (k in 1:length(gody)){
    for (j in 1:length(godx)){
      if (inter=="comb"){
        acp3[j,k]=surfuncCOP(godx[j], gody[k], cop=coco, para=c1)}
      if (inter=="casc"){
        acp3[j,k]= surfuncCOP(godx[j], gody[k], cop=coco, para=c1)/(1-godx[j])}

    }
  }
  grid <- expand.grid(lon=godx, lat=gody)
  levelplot(acp3 ~ lon * lat, data=grid, cuts=20, pretty=TRUE)
  cl2<-contourLines(coxi,coyi, acp3, levels = pobje)
  if(length(cl2)>0){
    cl2<-as.matrix(data.frame(cl2[[1]]$x,cl2[[1]]$y))} else{cl2<-NA}

}
