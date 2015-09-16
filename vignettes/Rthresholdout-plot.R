

## defines a plot function for Rthresholdout performance
## styling using default Rcssplot style
plotHresults = function(ar, types=c("train.1","train.2","test"),
  k=20,
  col=c("#0000ff","#00ff00","#ff0000"),
  bgcol=c("#ccccff","#ccffcc","#ffcccc"), main="", Rcssclass=c()) {
  
  xlim = c(0, max(ar[,"k"]))
  ylim=c(0.4, 1)
  plot(xlim, ylim, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", type="n",
       xlab="", ylab="", main="", axes=F, frame=F)

  ## draw divider lines
  for (nowy in seq(0.6, 1, by=0.1)) {
    Rcsslines(xlim, rep(nowy, 2), Rcssclass="guides")
  }  
  Rcsslines(xlim, rep(0.5, 2), Rcssclass="random")
  
  ## draw polygons
  ari = seq(1, nrow(ar))
  arirev = seq(nrow(ar), 1)
  for (i in 1:length(types)) {
    ilab = types[i];
    Rcsspolygon(c(ar[ari, "k"], ar[arirev,"k"]),
            c(ar[ari, paste0(ilab,".low")], ar[arirev, paste0(ilab,".high")]),
            col=bgcol[i], border=NA)
    Rcsslines(ar[ari,"k"], ar[ari, paste0(ilab,".median")], col=col[i], lwd=2)
    Rcsspoints(ar[ari,"k"], ar[ari, paste0(ilab,".median")], col=col[i])
  }

  Rcsslines(xlim, rep(0.5, 2), Rcssclass="random")
  Rcsslines(rep(k,2), ylim, Rcssclass="random")
  
  Rcssaxis(1, Rcssclass="x")
  aty = seq(0.4, 1, by=0.1)
  Rcssaxis(2, lwd=0, Rcssclass="y")
  Rcssaxis(2, pos=0, at=aty, labels=rep("", length(aty)), Rcssclass="y")
  
  Rcssmtext("k", side=1, Rcssclass="x")
  Rcssmtext("Accuracy", side=2, Rcssclass="y")
  if (Rcssclass=="nopdf") {
    main = gsub(".\\(", "\n\\(", main);
  }
  Rcssmtext(main, side=3, Rcssclass=c("main", Rcssclass))

  ## draw legend
  ## (ad-hoc translation from codes to human-readable descriptors)
  types[types=="train"] = "Historical"
  types[types=="train.1"] = "Training"
  types[types=="train.2"] = "Holdout"
  types[types=="test"] = "Test"
  Rcsslegend("bottomright", legend=types, col=col, Rcssclass=Rcssclass)
  
}
