#' prepro
#'
#' Data preprocessing function
#' @param x data
#' @export
prepro <- function(x) {
  # Remove missing values
  x_clean <- x[!is.na(x)]

  # Apply trimming if DescTools is available
  if ("Trim" %in% ls("package:DescTools")) {
    x_trimmed <- Trim(x_clean, trim = 0.01)
  } else {
    warning("Trim() function is not available. Using untrimmed data.")
    x_trimmed <- x_clean
  }

  return(x_trimmed)
}

#' Pickands
#'
#' Adaptive Pn selection function for Moment method
#' @param n sample size
#' @export
Pn_Mom <- function(n) {
  ifelse(n > 1850, 0.9,
         ifelse(n > 1110, 0.8,
                ifelse(n > 370, 0.6, 0.3)))
}

#' Pickands
#'
#' Adaptive Pn selection function for Pickands's method
#' @param n sample size
#' @export
Pn_Pick <- function(n) {
  ifelse(n > 1110, 0.9,
         ifelse(n > 740, 0.8,
                ifelse(n > 370, 0.6, 0.3)))
}

#' Pickands
#'
#' Pickands method
#' @param x data
#' @param pn the nominal coverage probability
#' @param alpha false alarm rate
#'
#' @returns A list containing the estimated tail index, control limits, the width of the control limits, sample size, and the nominal coverage probability.
#'   \item{Tail_Index_Upper}{uppler tail index}
#'   \item{Tail_Index_Lower}{lower tail index}
#'   \item{UCL}{upper control limit}
#'   \item{LCL}{lower control limit}
#'   \item{Band.Pickands}{width of control limits}
#'   \item{n}{sample size}
#'   \item{pn}{the nominal coverage probability}
#' @examples
#' data(ucisecom)
#' x1 = prepro(ucisecom$V2)[1:61]
#' n = length(x1)
#' Pickands(x1, pn=Pn_Pick(n))
#' @export
Pickands = function(x, pn, alpha = 0.0027){
  #  calculate tail index  Pickands (1975)
  x=unique(x)
  Rx=sort(x)
  n=length(x)
  range=seq(1,floor(n/4), by=1)
  r.p=function(k){
    log((Rx[n-1*k+1]-Rx[n-2*k+1])/(Rx[n-2*k+1]-Rx[n-4*k+1]))/log(2)
  }

  result <- dfphase1::rsp(r.p(range), plot = FALSE, alpha = 0.1)
  fitted_level <- result$fit[, "level"]
  fitted_scale <- result$fit[, "scale"]

  # -----level ----
  D.level=c()
  for (i in 1:length(fitted_level)-1) {D.level[i]=fitted_level[i+1]-fitted_level[i]}
  N.level=D.level[D.level!=0]
  S.level=c()
  for (i in 1:length(N.level)) {S.level[i]=sort(which(D.level==N.level[i]))}
  N.fitted.level=c()
  for (i in 1: length(N.level)){N.fitted.level[i]=fitted_level[S.level[i]]+N.level[i]}


  # -----Scale ----
  D.scale=c()
  for (i in 1:length(fitted_scale)-1) {D.scale[i]=fitted_scale[i+1]-fitted_scale[i]}
  N.scale=D.scale[D.scale!=0]
  S.scale=c()
  for (i in 1:length(N.scale)) {S.scale[i]=sort(which(D.scale==N.scale[i]))}
  N.fitted.scale=c()
  for (i in 1: length(N.scale)){N.fitted.scale[i]=fitted_scale[S.scale[i]]+N.scale[i]}

  # HJ TEST m
  s.level.scale=sort(c(S.level,S.scale))
  s.level.scale=unique(s.level.scale)

  if(is.null(s.level.scale)) {
    m=5
  } else if (length(s.level.scale)==1) {
    m=s.level.scale+1
  } else {
    D.S=c()
    for (i in 1:length(s.level.scale)) {
      if(is.na(s.level.scale[i+1])) {
        D.S[i]=floor(n/4)-s.level.scale[i]
      } else {D.S[i]=s.level.scale[i+1]-s.level.scale[i]}
    }
    final.max=sort(which(D.S==max(D.S)));final.max
    s.level.scale[final.max]
    m=s.level.scale[final.max]+1
  }

  m=m[1]

  k_u=m;m
  # ---------------------------------------------------------------------------------------------
  xm=max(range)
  ym1=min(r.p(range))
  ym2=max(r.p(range))


  SLC=s.level.scale+1;SLC
  SLC=unique(SLC)
  rnk=which(SLC==m);rnk


  if(is.null(s.level.scale)){SLC[rnk+1]=floor(n/4)
  } else {
    if(length(SLC)<=1 | length(SLC)==rnk) {SLC[rnk+1]=floor(n/4)} else {SLC[rnk+1]}}

  r.p.final=r.p(m);r.p.final.u=r.p.final


  #==== Step 5 : Exceed. Prob. UCL ===============================================

  var.x=function(r){
    (2^(2*r+1))*(r^2)/(2^r-1)^2
  }


  z=stats::qnorm((1+pn)/2)
  m=max(floor(n*alpha/2),1)

  UCL=Rx[n-m+1]+z*(Rx[n-m+1]-Rx[n-2*m+1])*(var.x(r.p.final)/(2*m))^0.5 # Theorem 3.1

  r_upper=r.p.final

  #center
  center=stats::median(x)

  # =============================================================================
  #  lower limit: LCL
  # =============================================================================

  y=-x

  ## 2nd order parameter lo estimation ---------------------------

  Rx=sort(y)
  n=length(y)

  #  calculate tail index  Pickands (1975) --------------------------------------

  Rx=sort(y)

  r.p=function(k){
    log((Rx[n-1*k+1]-Rx[n-2*k+1])/(Rx[n-2*k+1]-Rx[n-4*k+1]))/log(2)
  }

  range=seq(1,floor(n/4), by=1)
  r.p=function(k){
    log((Rx[n-1*k+1]-Rx[n-2*k+1])/(Rx[n-2*k+1]-Rx[n-4*k+1]))/log(2)
  }

  result <- dfphase1::rsp(r.p(range), plot = FALSE,alpha = 0.1)

  # Extract the fitted level estimates
  fitted_level <- result$fit[, "level"]
  fitted_scale <- result$fit[, "scale"]


  # -----level ----
  D.level=c()
  for (i in 1:length(fitted_level)-1) {D.level[i]=fitted_level[i+1]-fitted_level[i]}
  N.level=D.level[D.level!=0]
  S.level=c()
  for (i in 1:length(N.level)) {S.level[i]=sort(which(D.level==N.level[i]))}
  N.fitted.level=c()
  for (i in 1: length(N.level)){N.fitted.level[i]=fitted_level[S.level[i]]+N.level[i]}

  # -----Scale ----
  D.scale=c()
  for (i in 1:length(fitted_scale)-1) {D.scale[i]=fitted_scale[i+1]-fitted_scale[i]}
  N.scale=D.scale[D.scale!=0]
  S.scale=c()
  for (i in 1:length(N.scale)) {S.scale[i]=sort(which(D.scale==N.scale[i]))}
  N.fitted.scale=c()
  for (i in 1: length(N.scale)){N.fitted.scale[i]=fitted_scale[S.scale[i]]+N.scale[i]}

  # HJ TEST m
  s.level.scale=sort(c(S.level,S.scale))
  s.level.scale=unique(s.level.scale)

  if(is.null(s.level.scale)) {
    m=5
  } else if (length(s.level.scale)==1) {
    m=s.level.scale+1
  } else {
    D.S=c()
    for (i in 1:length(s.level.scale)) {
      if(is.na(s.level.scale[i+1])) {
        D.S[i]=floor(n/4)-s.level.scale[i]
      } else {D.S[i]=s.level.scale[i+1]-s.level.scale[i]}
    }
    final.max=sort(which(D.S==max(D.S)));final.max
    s.level.scale[final.max]
    m=s.level.scale[final.max]+1
  }

  m=m[1]

  k_l=m
  # ---------------------------------------------------------------------------------------------
  xm=max(range)
  ym1=min(r.p(range))
  ym2=max(r.p(range))

  SLC=s.level.scale+1
  SLC=unique(SLC)
  rnk=which(SLC==m)

  if(is.null(s.level.scale)){SLC[rnk+1]=floor(n/4)
  } else {
    if(length(SLC)<=1 | length(SLC)==rnk) {SLC[rnk+1]=floor(n/4)} else {SLC[rnk+1]}}

  r.p.final=r.p(m)
  r.p.final.l=r.p.final

  z=stats::qnorm((1+pn)/2)
  m=max(floor(n*alpha/2),1)
  LCL=-(Rx[n-m+1]+z*(Rx[n-m+1]-Rx[n-2*m+1])*(var.x(r.p.final)/(2*m))^0.5) # Theorem 3.1

  r_lower=r.p.final


  UCL.Pickands=UCL
  LCL.Pickands=LCL


  Band.Pickands=UCL-LCL

  # =========== OOC =====================================

  ooc=0
  for (i in 1: length(x)){
    if(x[i]<LCL.Pickands| x[i]>UCL.Pickands) ooc=ooc+1
  }
  ooc.Pickands=ooc
  ooc.Pickands.p=round(100*ooc.Pickands/length(x),2)


  return(list(Tail_Index_Upper=r_upper,Tail_Index_Lower=r_lower,
              UCL = UCL, LCL = LCL,
              Band.Pickands = Band.Pickands,
              n = n, pn=pn))
}

#' Moment
#'
#' Moment method
#' @param x data
#' @param pn the nominal coverage probability
#' @param alpha false alarm rate
#'
#' @returns A list containing the estimated tail index, control limits, the width of the control limits, sample size, and the nominal coverage probability.
#'   \item{Tail_Index_Upper}{uppler tail index}
#'   \item{Tail_Index_Lower}{lower tail index}
#'   \item{UCL}{upper control limit}
#'   \item{LCL}{lower control limit}
#'   \item{Band.Moment}{width of control limits}
#'   \item{n}{sample size}
#'   \item{pn}{the nominal coverage probability}
#' @examples
#' data(ucisecom)
#' x1 = prepro(ucisecom$V2)[1:61]
#' n = length(x1)
#' Moment(x1, pn=Pn_Mom(n))
#' @export
Moment = function(x, pn, alpha = 0.0027){
  Rx=sort(x)
  n=length(Rx)
  q=0.5
  xx=floor(n*q)+1
  Rxx=sum(Rx[seq(1:1)])


  Ry=Rx-Rxx
  n=length(Ry)
  range=seq(1,floor(n/4), by=1)

  ## Moment estimation

  Mn.1=function(k){
    apply(matrix(log(Ry[n-seq(1,k)+1])-log(Ry[n-k]),ncol=1),2,mean)
  }

  Mn.2=function(k){
    apply(matrix((log(Ry[n-seq(1,k)+1])-log(Ry[n-k]))^2,ncol=1),2,mean)
  }

  r.m=function(k){
    if(k==1){
      apply(matrix(log(Ry[n-seq(1,k)+1])-log(Ry[n-k]),ncol=1),2,mean)
    } else
      Mn.1(k)+1-0.5*((1-(Mn.1(k)^2)/Mn.2(k))^(-1))
  }

  rm=c()
  for (i in 1:floor(n/4)) {rm[i]=r.m(i)} # HJ test
  rm=stats::na.omit(rm)

  # result <- rsp(rm, plot = TRUE, alpha = 0.1)
  result <- dfphase1::rsp(rm, plot = FALSE, alpha = 0.1)

  # Extract the fitted level estimates
  fitted_level <- result$fit[, "level"]
  fitted_scale <- result$fit[, "scale"]

  # -----level ----
  D.level=c()
  for (i in 1:length(fitted_level)-1) {D.level[i]=fitted_level[i+1]-fitted_level[i]}
  N.level=D.level[D.level!=0]
  S.level=c()
  for (i in 1:length(N.level)) {S.level[i]=sort(which(D.level==N.level[i]))}
  N.fitted.level=c()
  for (i in 1: length(N.level)){N.fitted.level[i]=fitted_level[S.level[i]]+N.level[i]}


  # -----Scale ----
  D.scale=c()
  for (i in 1:length(fitted_scale)-1) {D.scale[i]=fitted_scale[i+1]-fitted_scale[i]}
  N.scale=D.scale[D.scale!=0]
  S.scale=c()
  for (i in 1:length(N.scale)) {S.scale[i]=sort(which(D.scale==N.scale[i]))}
  N.fitted.scale=c()
  for (i in 1: length(N.scale)){N.fitted.scale[i]=fitted_scale[S.scale[i]]+N.scale[i]}

  # HJ TEST m
  s.level.scale=sort(c(S.level,S.scale))
  s.level.scale=unique(s.level.scale)

  if(is.null(s.level.scale)) {
    m=5
  } else if (length(s.level.scale)==1) {
    m=s.level.scale+1
  } else {
    D.S=c()
    for (i in 1:length(s.level.scale)) {
      if(is.na(s.level.scale[i+1])) {
        D.S[i]=floor(n/4)-s.level.scale[i]
      } else {D.S[i]=s.level.scale[i+1]-s.level.scale[i]}
    }
    final.max=sort(which(D.S==max(D.S)));final.max
    s.level.scale[final.max]
    m=s.level.scale[final.max]+1
  }

  m=m[1]

  k_u=m
  # ---------------------------------------------------------------------------------------------
  rm=rm[is.finite(rm)]
  xm=max(range)
  ym1=min(rm)
  ym2=max(rm)

  SLC=s.level.scale+1
  SLC=unique(SLC)
  rnk=which(SLC==m)

  if(is.null(s.level.scale)){SLC[rnk+1]=floor(n/4)
  } else {
    if(length(SLC)<=1 | length(SLC)==rnk) {SLC[rnk+1]=floor(n/4)} else {SLC[rnk+1]}}

  rm1=function(k){rm[k]}



  # HJ revised below ###########
  rmx=r.m(m)
  rmx=if(is.na(rmx)){rm[m]} else rmx
  r.m.final=rmx
  r.m.final.u=r.m.final

  z=stats::qnorm((1+pn)/2)
  m=max(floor(n*alpha/2),1)
  m0=floor(n*alpha/2)
  UCL1=Rx[n-m0]+z*Rx[n-m0]*Mn.1(m)*abs(1-min(0,r.m.final))/sqrt(m) # Theorem 5.1

  UCL=Rxx+Ry[n-m0]+z*Ry[n-m0]*Mn.1(m)*abs(1-min(0,r.m.final))/sqrt(m) # HJ test


  r_upper=r.m.final

  #center
  center=stats::median(x)

  r_upper=r.m.final

  #center
  center=stats::median(x)

  # =============================================================================
  #  lower limit: LCL
  # =============================================================================
  y=-x

  ## 2nd order parameter lo estimation ---------------------------

  Rx=sort(y)
  n=length(y)

  ## Moment estimation
  Rx1=Rx[1]
  Rx2=Rx[2]
  Rxx=sum(Rx[seq(1:1)])
  Ry=Rx-Rxx
  n=length(Ry)
  range=seq(1,floor(n/4), by=1)

  ## Moment estimation

  Mn.1=function(k){
    apply(matrix(log(Ry[n-seq(1,k)+1])-log(Ry[n-k]),ncol=1),2,mean)
  }

  Mn.2=function(k){
    apply(matrix((log(Ry[n-seq(1,k)+1])-log(Ry[n-k]))^2,ncol=1),2,mean)
  }

  r.m=function(k){
    if(k==1){
      apply(matrix(log(Ry[n-seq(1,k)+1])-log(Ry[n-k]),ncol=1),2,mean)
    } else
      Mn.1(k)+1-0.5*((1-(Mn.1(k)^2)/Mn.2(k))^(-1))
  }

  rm=c()
  for (i in 1:floor(n/4)) {rm[i]=r.m(i)} # HJ test
  rm=stats::na.omit(rm)
  result <- dfphase1::rsp(rm, plot =FALSE, alpha = 0.1)

  # Extract the fitted level estimates
  fitted_level <- result$fit[, "level"]
  fitted_scale <- result$fit[, "scale"]

  # -----level ----
  D.level=c()
  for (i in 1:length(fitted_level)-1) {D.level[i]=fitted_level[i+1]-fitted_level[i]}
  N.level=D.level[D.level!=0]
  S.level=c()
  for (i in 1:length(N.level)) {S.level[i]=sort(which(D.level==N.level[i]))}
  N.fitted.level=c()
  for (i in 1: length(N.level)){N.fitted.level[i]=fitted_level[S.level[i]]+N.level[i]}


  # -----Scale ----
  D.scale=c()
  for (i in 1:length(fitted_scale)-1) {D.scale[i]=fitted_scale[i+1]-fitted_scale[i]}
  N.scale=D.scale[D.scale!=0]
  S.scale=c()
  for (i in 1:length(N.scale)) {S.scale[i]=sort(which(D.scale==N.scale[i]))}
  N.fitted.scale=c()
  for (i in 1: length(N.scale)){N.fitted.scale[i]=fitted_scale[S.scale[i]]+N.scale[i]}

  # HJ TEST m
  s.level.scale=sort(c(S.level,S.scale))
  s.level.scale=unique(s.level.scale)
  # s.level.scale
  if(is.null(s.level.scale)) {
    m=5
  } else if (length(s.level.scale)==1) {
    m=s.level.scale+1
  } else {
    D.S=c()
    for (i in 1:length(s.level.scale)) {
      if(is.na(s.level.scale[i+1])) {
        D.S[i]=floor(n/4)-s.level.scale[i]
      } else {D.S[i]=s.level.scale[i+1]-s.level.scale[i]}
    }
    final.max=sort(which(D.S==max(D.S)))
    s.level.scale[final.max]
    m=s.level.scale[final.max]+1
  }

  m=m[1]

  k_l=m
  # ---------------------------------------------------------------------------------------------
  rm=rm[is.finite(rm)]
  xm=max(range)
  ym1=min(rm)
  ym2=max(rm)
  SLC=s.level.scale+1
  SLC=unique(SLC)
  rnk=which(SLC==m)


  rm1=function(k){rm[k]}


  rmx=r.m(m)
  rmx=if(is.na(rmx)){rm[m]} else rmx

  r.m.final=rmx;r.m.final
  r.m.final.l=r.m.final

  #rsp(rm,alpha = 0.1, plot=TRUE)


  z=stats::qnorm((1+pn)/2)
  m=max(floor(n*alpha/2),1)
  m0=floor(n*alpha/2)
  LCL1=-(Rx[n-m0]+z*Rx[n-m0]*Mn.1(m)*abs(1-min(0,r.m.final))/sqrt(m)) # Theorem 5.1
  LCL=-(Rxx+Ry[n-m0]+z*Ry[n-m0]*Mn.1(m)*abs(1-min(0,r.m.final))/sqrt(m)) # HJ test


  r_lower=r.m.final

  Band.Moment=UCL-LCL

  UCL.Moment=UCL
  LCL.Moment=LCL


  # =========== OOC =====================================

  ooc=0
  for (i in 1: length(x)){
    if(x[i]<LCL.Moment| x[i]>UCL.Moment) ooc=ooc+1
  }
  ooc.Moment=ooc
  ooc.Moment.p=round(100*ooc.Moment/length(x),2)


  return(list(Tail_Index_Upper=r_upper,Tail_Index_Lower=r_lower,
              UCL = UCL, LCL = LCL,
              Band.Moment = Band.Moment,
              n = n, pn=pn))
}
