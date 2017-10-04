# -!- Stats -!-
# ConTest functions, version 1.0
# 20/07/2007
# By Julien Dutheil

#Copyright or ?? or Copr. Julien Dutheil
#
#Julien.Dutheil@univ-montp2.fr
#
#This software is a computer program whose purpose is to detect positions
#within a set of aligned sequence that are evolutionarily constrained.#
#
#This software is governed by the CeCILL license under French law and
#abiding by the rules of distribution of free software.  You can  use, 
#modify and/ or redistribute the software under the terms of the CeCILL
#license as circulated by CEA, CNRS and INRIA at the following URL
#"http://www.cecill.info". 
#
#As a counterpart to the access to the source code and  rights to copy,
#modify and redistribute granted by the license, users are provided only
#with a limited warranty  and the software's author,  the holder of the
#economic rights,  and the successive licensors  have only  limited
#liability. 
#
#In this respect, the user's attention is drawn to the risks associated
#with loading,  using,  modifying and/or developing or reproducing the
#software by the user in light of its specific status of free software,
#that may mean  that it is complicated to manipulate,  and  that  also
#therefore means  that it is reserved for developers  and  experienced
#professionals having in-depth computer knowledge. Users are therefore
#encouraged to load and test the software's suitability as regards their
#requirements in conditions enabling the security of their systems and/or 
#data to be ensured and,  more generally, to use and operate it in the 
#same conditions as regards security. 
#
#The fact that you are presently reading this means that you have had
#knowledge of the CeCILL license and that you accept its terms.
#

library(ade4)

# Generate p-values code
pval<-function(x)
{
  return(as.character(symnum(x, cutpoints=c(0,0.001,0.01,0.05,0.1,1), symbols=c("***","**","*",".","NS"))))
}

# Three built-in tag function for volume, polarity and charge.
# When detecting conserved positions, we vcan check which value of the property is conserved
# (for instance for volume conservation, small or large?).
# 'tag' function allow add the appropriate columns in the output.
# You can add your own...
volume.tag<-function(prop)
{
  tag<-character(length=length(prop))
  for(i in 1:length(prop))
  {
    #if(prop[i] <= 55) tag[i]<-"Small"
    #else if(prop[i] <= 100) tag[i]<-"Medium"
    #else tag[i]<-"Large"
    tag[i]<-ifelse(prop[i] < 83.5, "Small", "Large")
  }
  return(tag)
}

polarity.tag<-function(prop)
{
  tag<-character(length=length(prop))
  for(i in 1:length(prop))
  {
    #if(prop[i] <= 6) tag[i]<-"Non-polar"
    #else if(prop[i] <= 100) tag[i]<-"Medium"
    #else tag[i]<-"Polar"
    tag[i]<-ifelse(prop[i] < 8.35, "Non-polar", "Polar")
  }
  return(tag)
}

charge.tag<-function(prop)
{
  tag<-character(length=length(prop))
  for(i in 1:length(prop))
  {
    if(prop[i] <= -0.5) tag[i]<-"Negative"
    if(prop[i] <= 0.5) tag[i]<-"Neutral"
    else tag[i]<-"Positive"
  }
  return(tag)
}

# Main testing function:
test<-function(
    data,
    sim,
    property,
    level,
    exact=TRUE,
    FDR="fdr",
    alt="lower",
    do.plot="b",
    pch=list(simulation=1,observed=19,detected=19),
    col=list(simulation=gray(0.5),observed="blue",detected="red"),
    col.line="red",
    xlab=NULL,
    ylab=NULL,
    nsim=1000,                                                        #Number of simulated points to show
    label=FALSE,
    tag.function=NULL,
    ...)                                                              #Other parameters passed to plot
{
  p.norm<-paste("norm", property, sep=".")
  p.min <-paste("min" , property, sep=".")
  p.max <-paste("max" , property, sep=".")
  p.mean<-paste("mean", property, sep=".")
  d<-dudi.pca(sim[, c("norm", p.norm)], scannf=F, nf=2)
  dsup<-suprow(d, data[, c("norm", p.norm)])
  signif<-NA

  #We need to orient 2nd axis, so that negative values correspond to negative selection...
  mx<-max(sim$norm) #The fastest rate in the simulations
  tmp<-data.frame(x=mx, y=0) #Dummy point with highest speed and strongest conservation
  names(tmp)<-c("norm", p.norm)
  orig<-suprow(d, tmp)
  if(orig$lisup[,2] > 0)
  {
    d$li[,2]<- -d$li[,2]
    dsup$lisup[,2]<- -dsup$lisup[,2]
  }

  if(alt=="lower")
  {
    f<-function(x,level)
    {
      return((sum(d$li[,2] <= x) + 1) / (nrow(d$li) + 1) - level)
    }
    r<-uniroot(f, lower=-10, upper=10, level=level)
    signif<-dsup$lisup[,2] <= r$root
    pred<-data[signif, c("Group", "pr", "norm", p.norm, p.min, p.max, p.mean)]
    pred$Property<-rep(property, nrow(pred))
  }
  else
  {
    f<-function(x,level)
    {
      return((sum(d$li[,2] >= x) + 1) / (nrow(d$li) + 1) - level)
    }
    r<-uniroot(f, lower=-10, upper=10, level=level)
    signif<-dsup$lisup[,2] >= r$root
    pred<-data[signif, c("Group", "pr", "norm", p.norm, p.min, p.max, p.mean)]
    pred$Property<-rep(property, nrow(pred))
  }
  names(pred)<-c("Group", "pr", "norm", "norm.property", "min.property", "max.property", "mean.property", "property")
  pred[,"Axis1"]<-dsup$lisup[signif,1]
  pred[,"Axis2"]<-dsup$lisup[signif,2]

  #Compute exact p-values:
  if(exact & nrow(pred) > 0)
  {
    stats<-dsup$lisup[signif, 2]
    n<-nrow(d$li)
    pred$p.value<-NA
    for(i in 1:nrow(pred))
    {
      if(alt=="lower") pred[i,"p.value"]<-(sum(d$li[,2] <= stats[i]) + 1) / (n + 1)
      else             pred[i,"p.value"]<-(sum(d$li[,2] >= stats[i]) + 1) / (n + 1)
    }
    pred$code<-pval(pred$p.value)
  }
  if(exact & fdr != "none" & nrow(pred) > 0)
  {
    m<-nrow(data)
    pred$adjust.p.value<-p.adjust(pred$p.value, method = fdr, n = m)
    pred$adjust.p.value.code<-pval(pred$adjust.p.value)
  }
  #Correct for mutliple testing:
  signif2<-pred[,"adjust.p.value"] <= level
  signif[signif]<-signif2
 
  if(do.plot=="b")
  {
    close.screen(all=T)
    split.screen(c(1,2))
  }

  if(do.plot==1 | do.plot=="b")
  {
    if(do.plot=="b") screen(1)
    plot(dsup$lisup, type="n", xlab=ifelse(is.null(xlab),"Axis 1",xlab), ylab=ifelse(is.null(ylab),"Axis 2",ylab), ...)
    points(d$li[1:nsim,], col=col$simulation)
    points(dsup$lisup[!signif,], pch=pch$observed, col=col$observed)
    abline(v=0, col=col.line)
    abline(h=0, col=col.line)
    points(dsup$lisup[signif,], pch=pch$detected, col=col$detected)
    if(label) text(dsup$lisup[signif,], labels=data[signif,"Group"], col="black")
  }

  if(do.plot==2 | do.plot=="b")
  {
    if(do.plot=="b") screen(2)
    plot(data[, p.norm]~data[,"norm"], type="n", xlab=ifelse(is.null(xlab),"Norm",xlab), ylab=ifelse(is.null(ylab),property,ylab), ...)
    points(sim[1:nsim, p.norm]~sim[1:nsim, "norm"], col=col$simulation)
    points(data[!signif, p.norm]~data[!signif, "norm"], pch=pch$observed, col=col$observed)
    #abline(v=mean(sim[, "norm"]), col=col.line)
    #abline(h=mean(sim[, p.norm]), col=col.line)
    #Axis 1:
    p1<-d$cent - (as.matrix(d$co) %*% c(-100,0)) * d$norm
    p2<-d$cent - (as.matrix(d$co) %*% c(100,0)) * d$norm
    lines(c(p1[1,1],p2[1,1]),c(p1[2,1],p2[2,1]), col=col.line)
    #Axis 2:
    p1<-d$cent - (as.matrix(d$co) %*% c(0,-100)) * d$norm
    p2<-d$cent - (as.matrix(d$co) %*% c(0,100)) * d$norm
    lines(c(p1[1,1],p2[1,1]),c(p1[2,1],p2[2,1]), col=col.line)
    points(pred[signif2, c("norm", "norm.property")], pch=pch$detected, col=col$detected)
    if(label) text(pred[, c("norm", "norm.property")], labels=pred$Group, col="black")
  }

  if(!is.null(tag.function))
  {
    pred$label.property<-tag.function(pred$mean.property)
  }
  
  return(pred)
}
