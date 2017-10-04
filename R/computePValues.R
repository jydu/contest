#######################################################################################
# Options:

dataFile<-"../ConTest/TIM.csv"
simuFile<-"../ConTest/TIM.null.csv"
outputFile<-"TIM0.05.csv"
# Output file for positive selection detection.
# Set this option to "" if you don't want to perform positive selection detection.
outputFilePos<-"TIMpos0.05.csv"


# Correction for multiple testing:
fdr<-"fdr"

# Set the detection level (p-value threshold, final false discovery rate if multiple testing correction is performed. 
level<-0.05

# Properties names.
# Add as many properties as tested in ConTest.
# Each property must correspond to the column name in dataFile and outputFile.
# Built-in properties are "GVolume", "GPolarity" and "Charge".
properties<-c(
"Volume",
"Polarity",
"Charge"
)

#######################################################################################
# End of options. Do not edit this part...

data<-read.table(dataFile, header=T, sep="\t")
 sim<-read.table(simuFile, header=T, sep="\t")

source("ConTestFunctions.R")

tags<-list();
for(prop in properties) { tags[[prop]]<-NULL }
tags[["Volume"]]<-volume.tag
tags[["Polarity"]]<-polarity.tag
tags[["Charge"]]<-charge.tag
# You may add your own tag function here... see comments in ConTestFunctions.R.

predList<-list();
for(prop in properties)
{
  predList[[prop]]<-test(data, sim, prop, level, FDR=fdr, tag.function=tags[[prop]]);
}
pred<-NULL
if(length(properties) == 0) {
  cat("ERROR!!! No property specified!\n")
} else if(length(properties) == 1) {
  pred<-predList[[1]]
} else {
  pred<-merge(predList[[1]], predList[[2]], all=T)
  if(length(properties) > 2) {
    for(i in 3:length(properties)) {
      pred<-merge(pred, predList[[i]], all=T)
    }
  }
}
write.table(pred, outputFile, row.names=F, sep="\t", quot=F)

if(outputFilePos != "") {
  predList<-list();
  for(prop in properties) {
    predList[[prop]]<-test(data, sim, prop, level, FDR=fdr, alt="greater");
  }
  pred<-NULL
  if(length(properties) == 0) {
    cat("ERROR!!! No property specified!\n")
  } else if(length(properties) == 1) {
    pred<-predList[[1]]
  } else {
    pred<-merge(predList[[1]], predList[[2]], all=T)
    if(length(properties) > 2) {
      for(i in 3:length(properties)) {
        pred<-merge(pred, predList[[i]], all=T)
      }
    }
  }
  write.table(pred, outputFilePos, row.names=F, sep="\t",quot=F)
}

#######################################################################################
