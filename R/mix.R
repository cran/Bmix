
"mix" <- function(Z, alpha, g0params, times=NULL,  rho=NULL,  cat=0, 
                  state=NULL, read=FALSE, print=FALSE, N=100, niter=0, kout=FALSE)
{
  if(is.null(state)){ state <- sample(seq(0,999), 3) }

  Z <- as.matrix(Z)
  dim <- ncol(Z)
  Total <- nrow(Z)
  if(cat > dim) cat <- 0
  dim <- dim-cat
  
    if(cat>0){ levels <- c()
               for(i in 1:cat){
                 levels <- c(levels,length(unique(Z[,dim+i])))
                 if(min(Z[,dim+i] != 0)){ stop("cat min starts at zero")  }
               }}
    else levels = 0
    
  if(N < 1 || niter < 0)
    { stop("nums parameter is invalid") }

  if(is.null(rho)){ rho <- 0 }
  if(alpha <= 0 || rho > 1 || rho < 0)
    { stop("alpha/rho is invalid") }

  if( length(times) != Total ){ times <- rep(0,Total) }

  k <- integer(N*(niter+1)*Total*kout)
  
  out <- .C("mixsample",
            rstate=as.integer(state),
            iopar=as.integer(c(read,print)),
            nums=as.integer(c(N,niter)),
            dims=as.integer(c(dim, cat, levels)),
            Total=as.integer(Total),
            times=as.integer(times),
            Z=as.double(t(Z)),
            params=as.double(c(alpha, rho, g0params)),
            margllhd=double(Total),
            m=integer(Total),
            k=k,
            PACKAGE="Bmix")


  
  if(kout){ k <- array(out$k, dim=c(Total,N,niter+1),
                       dimnames=list(paste("o", 1:Total, sep=""),
                         paste("p", 1:N, sep=""),
                         paste("i", 1:(niter+1), sep=""))) }
  
  invisible( list(state=state, read=read, print=print, Z=Z,
                  N=N, niter=niter, dim=dim, cat=cat, times=times,
                  levels=levels, Total=Total, g0params=g0params,
                  alpha = alpha, rho = rho, logprob = out$margllhd, m=out$m, k=k) )
                                 
}



"particle" <- function(i, mixobj, t, rho=0){
  pmat <- read.table(paste(".particle", i,".", t,".", rho,".txt", sep="")) 
  counts <- c()  
  if(mixobj$cat > 0){ for(b in 1:mixobj$cat){
    counts <- c(counts, paste("counts",b,1:mixobj$levels[b], sep=".")) }}
  names(pmat)  <- c("n", paste("mean", 1:mixobj$dim,sep="."), paste("S", 1:mixobj$dim^2,sep="."),
                        counts, "p", paste("a",1:mixobj$dim,sep="."), paste("B",1:mixobj$dim^2,sep="."),"c")
  row.names(pmat) <- c(0:(nrow(pmat)-1))
  return(pmat)
}
                       
