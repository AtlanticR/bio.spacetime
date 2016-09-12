
hdf5_timer = function( engine="hdf5", nr=2000, nc=2000, 
  cs=c(8, 16, 24, 32, 40, 48, 64, 78, 100, 128, 184, 256, 512, 1024), n=1, times=3 ) {
  #\\ find optimal chuncksize for hdf5
  #\\ also compare with bigmemory .. just a bit faster than hdf5 (1 to 2X) but no flexibility with data structures

  dd = matrix(runif(nr*nc), ncol=nc )

  if(engine=="hdf5") {
    require(h5)
    fn = tempfile(tmpdir=".", fileext=".h5")
    hdf5 = function(fn, nr, nc, dd, n, cs) {
      if (file.exists(fn)) file.remove(fn)
      b <- h5file(name =fn, mode = "a")
      b["dd", compression=0L, chunksize=cs] <- dd
      h5close(b)
      for (i in 1:n){
        b <- h5file(name =fn, mode = "a")
        b["dd"][1,] <- log(b["dd"][1,] + b["dd"][2,])
        h5close(b)
      }
      return(NULL)
    }
    #system.time( hdf5(nr,nc,dd, n, cs) )
    out = matrix( NA, ncol=2, nrow=length(cs) )
    out[,1] = cs
    for ( i in 1:length(cs) ) {
      mb = microbenchmark::microbenchmark( hdf5( fn=fn, nr=nr, nc=nc, dd=dd, cs=cs[i], n=n), 
          times=times, unit="s" )
      a = summary(mb)
      #a = system.time( hdf5(fn, nr, nc, dd, n, cs[i]) )[1]
      print( paste( cs[i], ":",  round(a$mean,5) , attr(a, "unit")))
      out[i,2]  = as.numeric(a$mean)
    }
     plot( out[,2] ~ log(out[,1]), xlab="chunksize", ylab="time", pch=20, xaxt="n", type="b" )
     axis(1, log(cs), as.character(cs) )
     if (file.exists(fn)) file.remove(fn)
     return(out)
  }
  
  if(engine=="bigmemory") {
    require(bigmemory)
    a = filebacked.big.matrix( nrow=nr, ncol=nc, type="double", dimnames=NULL, separated=FALSE,
         backingpath=".", backingfile="bm.test.bigmemory", descriptorfile="bm.test.desc" )
    a[] = dd 
    bm = function(nr, nc, dd, n) {
      for(i in 1:n) {
        a = bigmemory::attach.big.matrix("bm.test.desc", path="." )
        a[1,] = log(a[1,] +a[2,])
      }
      return(NULL)
    }
    #system.time( bm(nr,nc,dd, n) )
    s0 = c(10, 50, 100, 500, 1000, 2000, 4000)
    out = matrix( NA, ncol=2, nrow=length(s0) )
    out[,1] = s0
    for ( i in 1:length(s0) ) {
      mb = microbenchmark::microbenchmark( bm( nr=s0[i], nc=s0[i], dd=dd, n=n),
        times=times, unit="s" )
      ab = summary(mb)
      print( paste( s0[i], ":",  round(ab$mean,5) , attr(ab, "unit")))
      out[i,2] = as.numeric(ab$mean)
    }
    return(out)
  }
}
