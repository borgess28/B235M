### Analyze B-235 pct cover data in google drive ###
library(googlesheets)
library(plyr)
library(dplyr)
library(reshape2)
library(vegan)

gs_auth(new_user = TRUE)

pctCover <- as.data.frame(gs_read(gs_key('1usTHRIt4PfwAu_xJH7bxtRfmRdn48ZwSxtm6Fh3V_iA', 
                             lookup = TRUE), ws = "pct cover"))
round5 <- grep("round5|total|Total", names(pctCover))
pctCover1 <- pctCover[, -round5]
plotDat <- as.data.frame(gs_read(gs_key('1usTHRIt4PfwAu_xJH7bxtRfmRdn48ZwSxtm6Fh3V_iA', 
                                        lookup = TRUE), ws = "plotInfo"))

ptIntDat <- as.data.frame(gs_read(gs_key('1gnSDblUAbubdeaGwsXG-rMdzUf3m5UFzs-Rp-fCLtec', 
                                         lookup = TRUE) ))
mapping <- as.data.frame(gs_read(gs_key('1t17W_XX8yeIamxOp0LSlqDhvW9mt9QhOUQdGxaoqaBE',
                                        lookup=TRUE) ))
 
## DATA PREP ##

pctCover1$siteID <- ""
pctCover1$siteID[grep("CRE", pctCover1$plotID)] <- "Crescent"
pctCover1$siteID[grep("CAN", pctCover1$plotID)] <- "Canada"
pctCover1$siteID[grep("BOW", pctCover1$plotID)] <- "Bowles"
pctCover1$siteID[grep("PLA", pctCover1$plotID)] <- "Relict"
pctCover1$siteID[grep("MCK", pctCover1$plotID)] <- "McKnight"
pctCover1$siteID[grep("LOS", pctCover1$plotID)] <- "Lost Seal"
navals <- which(pctCover1$siteID=="" | is.na(pctCover1$siteID) )
if(length(navals) >0) {
  pctCover1 <- pctCover1[-navals, ]
}

pctCover1$plotNum <- gsub("([[:alpha:]]+)_P([1-3])", "\\2", pctCover1$plotID)
pctCover1$newPlotID <- paste0(pctCover1$siteID, "-", pctCover1$plotNum, "-", pctCover1$collectDate)

plotIDs <- pctCover1$newPlotID
pctCover2 <- pctCover1 %>% select(-plotID, -collectDate, -subplotID, -imageNumber, -siteID,
                                -remarks, -plotNum, -water, -newPlotID)


pctCoverMeans <- pctCover2[1, ]
for(i in 1:length(unique(plotIDs)) ) {
  tmpPlot <- unique(plotIDs)[i]
  tmp <- na.omit(pctCover2[unique(plotIDs)==tmpPlot,] )
  pctCoverMeans[i,] <- colSums(tmp)/nrow(tmp)* 0.01
  rownames(pctCoverMeans)[i] <- tmpPlot
}

pctCoverQ <- decostand(pctCoverMeans, method="total")    ### This is the final percent coverage data for Quadrat Survey data

rowSums(pctCoverQ)  # verify totals add up to 1 



## Point intercept data ##

ptIntDat$siteID <- ""
ptIntDat$siteID[grep("Bowl", ptIntDat$transectID)] <- "Bowles"
ptIntDat$siteID[grep("Pla", ptIntDat$transectID)] <- "Relict"
ptIntDat$siteID[grep("McKn", ptIntDat$transectID)] <- "McKnight"
ptIntDat$siteID[grep("Cres", ptIntDat$transectID)] <- "Crescent"
navals <- which(ptIntDat$siteID=="" | is.na(ptIntDat$siteID) )
if(length(navals) >0) {
  pctCover <- pctCover[-navals, ]
}

ptIntDat$collectDate <- gsub("([[:alpha:]]+)(01|02|03)_(201[0-9]{5})_T[0-9]{4}", "\\3", ptIntDat$transectID)
ptIntDat$plotNum <- gsub("([a-zA-Z])+(0)([1-3])_.*", "\\3", ptIntDat$transectID)
ptIntDat$plotID <- paste0(ptIntDat$siteID, "-", ptIntDat$plotNum, "-", ptIntDat$collectDate)

ptIntMelted <- melt(ptIntDat, id.vars = c("transectID", "siteID", "collectDate",
                                                      "plotNum", "plotID") )
summ <- as.data.frame(with(ptIntMelted, table(value, plotID) ) )

summWide <- spread(summ, value, Freq)
rownames(summWide) <- summWide$plotID
summWide <- summWide[,-1]

pctCoverPtInt <- decostand(summWide, method="total")    # Final percent cover data set for Line intercept
rowSums(pctCoverPtInt)    # verify totals add up to 1

# For Now: Export to map matching categories
write.csv(pctCoverQ, "/Users/lstanish/Github/devTOS/microbes/B235quadPctDat.csv",
          na="")

write.csv(pctCoverPtInt, "/Users/lstanish/Github/devTOS/microbes/B235ptIntPctDat.csv",
          na="")


### MANUAL EDITING OF MAPPED CATEGORIES: TOO TIRED To DO PROGRAMMATICALLY!

# Load edited data
quadFinal <- read.csv('/Users/lstanish/Github/devTOS/microbes/B235quadPctDatManEdited.csv',
                       stringsAsFactors=FALSE)
ptIntFinal <- read.csv('/Users/lstanish/Github/devTOS/microbes/B235ptIntPctDatManEdited.csv',
                        stringsAsFactors=FALSE)
quadFinal$method <- "quadrat"
ptIntFinal$method <- "pointIntercept"
#quadFinal <- select(quadFinal, -plotID)
#ptIntFinal <- select(ptIntFinal, -plotID)

### Compare data sets

datMerged <- bind_rows(quadFinal, ptIntFinal, na=0)
dupePlots <- datMerged$plotID[duplicated(datMerged$plotID)]
datMerged <- datMerged[datMerged$plotID %in% dupePlots, ]
write.csv(datMerged, '/Users/lstanish/Github/devTOS/microbes/B235mergedData.csv')
# manually edited the data sheet to remove NA values

# re-import manually edited data set
datNoNa <- read.csv('/Users/lstanish/Github/devTOS/microbes/B235mergedDataTotalSimplestCategories.csv')

for(i in 2:ncol(datNoNa)-2) {
  i=18
  tmp <- names(datNoNa)[i]
  mod <- aov(datNoNa[,i]~method, data=datNoNa)
  print(tmp)
  print(summary(mod) )
}

datforNMDS <- select(datNoNa, -method, -plotID)
mds <- metaMDS(datforNMDS, distance = 'bray', k=2, trymax = 50)
mds.sc <- data.frame(scores(mds, display="sites", choices=c(1:2)) )
mds.sp <- data.frame(scores(mds, display="species", choices=c(1:2)) )

plot(mds, type="none")
points(mds.sc[,1], mds.sc[,2], col=as.integer(datNoNa$method), 
       pch=as.integer(datNoNa$plotID), cex=1.5 )
text(mds.sp[,1], mds.sp[,2], pch=19,  cex=0.7, labels=names(datforNMDS) )

legend("topright", legend=c('line intercept', 'quadrat'), col=c(1,2), pch=19)


# convert to mats versus abiotic substrates



