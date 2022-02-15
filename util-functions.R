## General utility functions for the analysis.

#' @title Plot of PCA results
#'
#' `plot_pca` is a simple utility function to plot the results from a PCA
#' analysis.
#'
#' @param pc the result from a principal component analysis (i.e. the result
#'     returned by `prcomp`.
#'
#' @param pch the point character. See [plot()] or [par()] for more information.
#'
#' @param col the color to be used for each data point/sample.
#'
#' @param pc_x `integer(1)` defining which principal component should be drawn
#'     on the x-axis.
#'
#' @param pc_y `integer(1)` defining the principal component to be drawn on the
#'     y-axis.
#'
#' @param main `character(1)` with the optional title of the plot.
#'
#' @param labels `character` with length equal to the number of samples. If
#'     provided, these will be displayed instead of data points.
#'
#' @param ... additional arguments to be passed to the [points()] or [text()]
#'     calls (if `labels = NULL` or not).
#'
#' @author Johannes Rainer
#'
#' @md
plot_pca <- function(pc, pch = 16, col = "#000000", pc_x = 1, pc_y = 2, 
                     main = "", labels = NULL, ...) {
    pcSummary <- summary(pc)
    plot(pc$x[, pc_x], pc$x[, pc_y], pch = NA, main = main,
         xlab = paste0("PC", pc_x, ": ",
                       format(pcSummary$importance[2, pc_x] * 100, 
                              digits = 3), " % variance"),
         ylab = paste0("PC", pc_y, ": ",
                       format(pcSummary$importance[2, pc_y] * 100, 
                              digits = 3), " % variance"))
    grid()
    if (!is.null(labels)) 
        text(pc$x[, pc_x], pc$x[, pc_y], labels = labels, col = col, 
             ...)
    else points(pc$x[, pc_x], pc$x[, pc_y], pch = pch, col = col, 
                ...)
}

#' @title Calculate relative standard deviations
#'
#' `rsd` and `rowRsd` are convenience functions to calculate the relative
#' standard deviation of a numerical vector or for rows of a numerical matrix,
#' respectively.
#'
#' @param x for `rsd` a `numeric`, for `rowRsd` a numerical `matrix`.
#'
#' @param na.rm `logical(1)` whether missing values (`NA`) should be removed
#' prior to the calculations.
#'
#' @author Johannes Rainer
#'
#' @md
rsd <- function(x, na.rm = TRUE) {
    sd(x, na.rm = na.rm) / abs(mean(x, na.rm = na.rm))
}
rowRsd <- function(x, na.rm = TRUE)
    apply(x, MARGIN = 1, rsd, na.rm = na.rm)

#' @title Determine the proportion of missing values
#'
#' `naProp` and `rowNaProp` determine the proportion of missing values in a
#' numeric vector or in rows of a numeric matrix
#'
#' @param x `numeric` or numeric `matrix`.
#'
#' @author Johannes Rainer
#'
#' @md
naProp <- function(x) {
    sum(is.na(x)) / length(x)
}
rowNaProp <- function(x) {
    apply(x, MARGIN = 1, naProp)
}

#'@title Extract run start time stamp
#'
#' @description
#'
#' Extract the start time stamps from all mzML files of an `MSnExp` object.
#'
#' @param x `MSnExp` object
#'
#' @param format `character(1)` defining the date/time format of the time
#'     stamp. If `NULL` the time stamp will be returned as a `character`.
#' 
#' @return `character` with the start time stamps.
extractTimeStamps <- function(x, format = "%Y-%m-%dT%H:%M:%S") {
    stopifnot(inherits(x, "MSnExp"))
    ts <- sapply(fileNames(x), function(z) {
        fl <- mzR::openMSfile(z)
        run_info <- mzR::runInfo(fl)
        mzR::close(fl)
        run_info$startTimeStamp
    })
    if (length(format))
        strptime(ts, format = format)
    else ts
}

ConfusionPlot <- function(Output) {
  
  #https://ragrawal.wordpress.com/2011/05/16/visualizing-confusion-matrix-in-r/
  
  #generate random data 
  data = data.frame(predict(Output$Model, Output$ModelParam), Output$ModelTrueSource)
  names(data) = c("Actual", "Predicted") 
   
  #compute frequency of actual categories
  actual = as.data.frame(table(data$Actual))
  names(actual) = c("Actual","ActualFreq")
   
  #build confusion matrix
  confusion = as.data.frame(table(data$Actual, data$Predicted))
  names(confusion) = c("Actual","Predicted","Freq")
   
  #calculate percentage of test cases based on actual frequency
  confusion = merge(confusion, actual, by=c("Actual"))
  confusion$Percent = confusion$Freq/confusion$ActualFreq*100
   
  #render plot
  # we use three different layers
  # first we draw tiles and fill color based on percentage of test cases
  tile <- ggplot() +
  geom_tile(aes(x=Actual, y=Predicted,fill=Percent),data=confusion, color="black",size=0.1) +
  labs(x="Actual",y="Predicted")
  tile = tile + 
  geom_text(aes(x=Actual,y=Predicted, label=sprintf("%.1f", Percent)),data=confusion, size=3, colour="black") +
  #scale_fill_gradient(low = "yellow", high = "blue")
    scale_fill_gradient(low = "#ffffbf", high = "#91bfdb")
   
  # lastly we draw diagonal tiles. We use alpha = 0 so as not to hide previous layers but use size=0.3 to highlight border
  tile = tile + 
  geom_tile(aes(x=Actual,y=Predicted),data=subset(confusion, as.character(Actual)==as.character(Predicted)), color="black",size=0.3, fill="black", alpha=0) 
   
  #render
  tile
}



    
    ConfusionPlotNew <- function(Actual,Predicted) {
  
  #https://ragrawal.wordpress.com/2011/05/16/visualizing-confusion-matrix-in-r/
  data <- data.frame(Actual,Predicted)
  #generate random data 
  names(data) = c("Actual", "Predicted") 
  
  #compute frequency of actual categories
  actual = as.data.frame(table(data$Actual))
  names(actual) = c("Actual","ActualFreq")
  
  #build confusion matrix
  confusion = as.data.frame(table(data$Actual, data$Predicted))
  names(confusion) = c("Actual","Predicted","Freq")
  
  #calculate percentage of test cases based on actual frequency
  confusion = merge(confusion, actual, by=c("Actual"))
  confusion$Percent = confusion$Freq/confusion$ActualFreq*100
  
  #render plot
  # we use three different layers
  # first we draw tiles and fill color based on percentage of test cases
  tile <- ggplot() +
    geom_tile(aes(x=Actual, y=Predicted,fill=Percent),data=confusion, color="black",size=0.1) +
    labs(x="Actual",y="Predicted")
  tile = tile + 
    geom_text(aes(x=Actual,y=Predicted, label=sprintf("%.1f", Percent)),data=confusion, size=3, colour="black") +
    #scale_fill_gradient(low = "yellow", high = "blue")
    scale_fill_gradient(low = "#ffffbf", high = "#91bfdb")
  
  # lastly we draw diagonal tiles. We use alpha = 0 so as not to hide previous layers but use size=0.3 to highlight border
  tile = tile + 
    geom_tile(aes(x=Actual,y=Predicted),data=subset(confusion, as.character(Actual)==as.character(Predicted)), color="black",size=0.3, fill="black", alpha=0) 
  
  #render
  tile
}



PLS_Training_Boot <- function(Data,Sources,indTest, ncomp = 10) {
#' A randomnly chosen subset of samples are used for training the samples
    
  
  
  #' The data for the  training and test set respectively is chosen along with the
  #' corresponding predicted variables
  training2 <- Data[,-indTest]
  test2 <- Data[,indTest]
  trainSource <- Sources[-indTest]
  testSource <- Sources[indTest]
  #' THe data is preprossesed (centered and scaled)
  preProcValues2 <- preProcess(t(training2))
  #' This preprossesing is applied to the data, note the confusing using of predict
  trainDescr2 <- predict(preProcValues2, t(training2))
  testDescr2 <- predict(preProcValues2, t(test2))
  #' A model is now trained using the centered data, and the predicted variables
  useSoftmax2 <- plsda(trainDescr2, trainSource, ncomp )
  Output <- list (Model = useSoftmax2,
                  ModelParam = testDescr2,
                  ModelTrueSource = testSource
                  )
}


PLS_Training <- function(Data,Sources,numCuts) {
#' A randomnly chosen subset of samples are used for training the samples
    indTrain <- sample(seq(along = Sources), numCuts)
  
  
  #' The data for the  training and test set respectively is chosen along with the
  #' corresponding predicted variables
  training2 <- Data[,indTrain]
  test2 <- assay(sexpFlt,"count")[,-indTrain]
  trainSource <- Sources[indTrain]
  testSource <- Sources[-indTrain]
  #' THe data is preprossesed (centered and scaled)
  preProcValues2 <- preProcess(t(training2))
  #' This preprossesing is applied to the data, note the confusing using of predict
  trainDescr2 <- predict(preProcValues2, t(training2))
  testDescr2 <- predict(preProcValues2, t(test2))
  #' A model is now trained using the centered data, and the predicted variables
  useSoftmax2 <- plsda(trainDescr2, trainSource, ncomp = 10)
  Output <- list (Model = useSoftmax2,
                  ModelParam = testDescr2,
                  ModelTrueSource = testSource,
                  indTrain = indTrain)
}

GGSaver <- function(PlotName, Samplename = SampleName){ 
  ggsave(filename = paste0("./Output/Figures/", Samplename,"/",Samplename,PlotName,".pdf")) }

NormFuncPrimi <- function(dta) {
  Dtacolmean <- apply(dta,1,mean)  
  DtacolSum <- dta/Dtacolmean
  Dtarowmean <- apply(DtacolSum,1,mean)
  DtaDeseqFinSum <- (DtacolSum)/Dtarowmean
}


FastaFormater <- function (X, ForwardConstant = "", ReverseConstant=""){
  Sequences <- X[,1]
  Names <- rownames(X)
Xfasta <- matrix(ncol=1, nrow = 2 * length(Sequences))
Xfasta[c(TRUE, FALSE)] <- paste0(">", Names)
Xfasta[c(FALSE, TRUE)] <- paste0 (ForwardConstant,Sequences,ReverseConstant)
Xfasta
}

ConfusionPlotFreq <- function(Actual,Predicted) {
  
  #https://ragrawal.wordpress.com/2011/05/16/visualizing-confusion-matrix-in-r/
  data <- data.frame(Actual,Predicted)
  #generate random data 
  names(data) = c("Actual", "Predicted") 
  
  #compute frequency of actual categories
  actual = as.data.frame(table(data$Actual))
  names(actual) = c("Actual","ActualFreq")
  
  #build confusion matrix
  confusion = as.data.frame(table(data$Actual, data$Predicted))
  names(confusion) = c("Actual","Predicted","Freq")
  
  #calculate percentage of test cases based on actual frequency
  confusion = merge(confusion, actual, by=c("Actual"))
  confusion$Percent = confusion$Freq/confusion$ActualFreq*100
  confusion$Label = paste0(confusion$Freq, "/", confusion$ActualFreq)
  #render plot
  # we use three different layers
  # first we draw tiles and fill color based on percentage of test cases
  tile <- ggplot() +
    geom_tile(aes(x=Actual, y=Predicted, fill=Percent),data=confusion, color="black",size=0.1) +
    labs(x="Actual",y="Predicted")
  tile = tile + 
    geom_text(aes(x=Actual,y=Predicted, label=Label), data=confusion, size=3, colour="black") +
    #geom_text(aes(x=Actual,y=Predicted, label=sprintf("%.1f", Percent)),data=confusion, size=3, colour="black") +
    #scale_fill_gradient(low = "yellow", high = "blue")
    scale_fill_gradient(low = "#ffffbf", high = "#31a354")
  
  # lastly we draw diagonal tiles. We use alpha = 0 so as not to hide previous layers but use size=0.3 to highlight border
  tile = tile + 
    geom_tile(aes(x=Actual,y=Predicted),
              data=subset(confusion, as.character(Actual)==
                            as.character(Predicted)), 
              color="black",size=0.3, fill="black", alpha=0) 
  
  #render
  tile
}

VariableWiseMean <- function(Data, Func = mean) {
  (apply(Data,2,"/",apply(Data,1,Func)))
}