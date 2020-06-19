# Auxiliary functions




# this function reads from .H5 files and returns average intensities
# for all the m/z channels included in the peak table. it uses
# functionalities of "TofDaqR" and "rhdf5" libraries

readH5File <- function(WorkDir, SubDirPattern){
  library("TofDaqR")
  library("rhdf5")
  DirNames <- list.files(path = WorkDir, pattern = SubDirPattern)
  for (i in 1:length(DirNames)) {
    setwd(paste0(WorkDir, DirNames[i]))
    FileNames <- list.files(pattern = ".h5")
    MassesFile <- GetPeakParametersFromH5(FileNames[i])
    PeakTable <- unlist(MassesFile$label, use.names = FALSE)
    MatrixInit <- matrix(data = 0, nrow = length(FileNames),
                         ncol = length(PeakTable))
    MatrixInit <- `colnames<-`(MatrixInit, PeakTable)
    assign(paste0("MeanIntensities_", DirNames[i]), MatrixInit)
  }
  Summary_MeanIntensities <- mget(ls(pattern = "MeanIntensities_"))
  rm(list = names(Summary_MeanIntensities))

  for (i in 1:length(DirNames)) {
    setwd(paste0(WorkDir, DirNames[i]))
    FileNames <- list.files(pattern = ".h5")
    for (j in 1:length(FileNames)) {
      InfoFile <- GetH5Descriptor(FileNames[j])
      NbrWaveforms <- InfoFile$nbrWaveforms
      NbrSegments <- InfoFile$nbrSegments
      NbrBufs <- InfoFile$nbrBufs
      NbrWrites <- InfoFile$nbrWrites
      CurData <- h5read(FileNames[j], "/PeakData/PeakData")
      for (k in 1:length(PeakTable)) {
        MeanIntensities <- mean(CurData[k, (1:NbrSegments), (1:NbrBufs),
                                    (1:NbrWrites)])*NbrWaveforms
        Summary_MeanIntensities[[i]][j,k] <- MeanIntensities
      }
    }
  }
  return(Summary_MeanIntensities)
}




# this function reads data series for a user-defined selection of
# isotopes from the output of the "readH5File" function and returns them
# in a data.frame object

extractDataSeries <- function(Data, IsotopeList){
  PeakTable <- colnames(Data[[1]])
  IndexIsotope <- match(IsotopeList, PeakTable)
  Amplitudes <- names(Data)

  Colnames <- vector(mode = "character", (length(IsotopeList)+2))
  for (i in 1:length(IsotopeList)) {
  Colnames[(i+2)] <- IsotopeList[i]
  }
  Colnames[1:2] <- c("MassAxis", "NotchAmplitude")

  LowerBound <- -(length(Data[[1]][,IndexIsotope[1]])-1)/2
  UpperBound <- (length(Data[[1]][,IndexIsotope[1]])-1)/2
  MassAxis <- rep(c(LowerBound:UpperBound), length(Amplitudes))

  AmplitudeList <- vector(mode = "character",
                          length = (length(Data[[1]][,1]))*
                          length(Amplitudes))
  for (i in 1:length(Amplitudes)) {
    for (j in 1:length(Data[[i]][,IndexIsotope[i]])) {
      AmplitudeList[(j+(length(Data[[i]][,IndexIsotope[i]])*(i-1)))] <-
      substring(Amplitudes[i], 17)
    }
  }
  DataMatrix <- matrix(data = 0, ncol = length(IsotopeList),
                       nrow = length(Data[[1]][,IndexIsotope[1]])*
                       length(Amplitudes))
  for (i in 1:length(Amplitudes)) {
    for (k in 1:length(IndexIsotope)) {
      for (j in 1:length(Data[[i]][,IndexIsotope[k]])) {
        DataMatrix[(j+(length(Data[[i]][,IndexIsotope[k]])*(i-1))),k] <-
        Data[[i]][j,IndexIsotope[k]]
      }
    }
  }

  DataToPlot <- cbind.data.frame(MassAxis, AmplitudeList, DataMatrix)
  colnames(DataToPlot) <- Colnames
  return(DataToPlot)
}






# this function reads data from the output of the "extractDataSeries"
# function and provides a graphical representation using functionalities
# of the "ggplot2" library

plotDataSeries <- function(Data, IsotopeOfInterest){
library("ggplot2")

ggplot(data = DataToPlot)+
  geom_point(mapping = aes(x = MassAxis,
                           y = DataToPlot[,IsotopeOfInterest],
                           color = NotchAmplitude), lwd = 1.2)+
  scale_y_continuous(trans = "log10", labels = scales::comma_format())+
  labs(x = "distance of notch position from target mass [m/z]",
       y = "signal intensity [cps]",
       title = NULL,
       subtitle = NULL) +
  theme(panel.background = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.95, 0.05),
        legend.justification = c("right", "bottom"),
        axis.line = element_line(),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 10,
                                   angle = 0, hjust = 0.5),
        axis.title = element_text(color = "black", size = 12))
}
