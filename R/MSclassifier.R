#' Median-Supplement methods
#' 
#' MSclassifier implements a median-supplement approach to machine learning, supporting 
#' complete compliance efforts by never missing sensitive sub-datasets or allowing some 
#' sub-datasets to escape the classification process when balancing overall dataset as 
#' required in traditional classification models, for an automated and effective binary 
#' classification for optimal decision making. The median-supplement approaches to 
#' machine learning were first used to decipher hormone and HER2 receptor status phenotypes
#' in breast cancer.
#'  
#' @param X a data frame of values of attributes (e.g. gene expression levels) and 
#' classes (e.g. receptor status phenotypes in breast cancer).  Samples are in rows while attributes are in columns.  
#' The last column of X should have the classes for all instances in X (e.g. receptor status phenotypes of samples). 
#' This form the training set.
#
#' @param testset   the set of new instances to be classified. The Default is NULL. When set to NULL, the function 
#' returns only the model. To classify new instances, specify the data frame of the new instances as the test set. 
#' It should have the form (and attributes) of X.
#'   
#' @param method specifies whether to determine a median-supplement Random
#' Forest or median-supplement Naive Bayes. 
#' "MSRandomForest" infers median-supplement Random Forest. 
#' "MSNaiveBayes" applies the median-supplement Naive Bayes.
#' The default is median-supplement Random Forest. 

#' @note These methods are applied to binary classification problems such as 
#' finding receptor status phenotypes in breast cancer. 
#' They require unbalanced training data set
#' involving only two labels.
#' 
#' They enhance capacities of Naive Bayes and Random Forest
#' as predictive classifiers to determine the classes or labels of
#' test sets or new instances from a training set with unequal numbers 
#' of members of classes. For example, it can be used for identifying receptor status phenotypes
#' of breast cancer patients from gene expression
#' data with unequal numbers of members of receptor status
#' phenotypes. In this case, the method determines the presence or absence of
#' receptor status phenotype. 
#' 
#' The methods proceed by introducing a
#' predetermined number of supplementary instances based on
#' the median of each attribute (feature) of the 
#' training set. They then apply a Naive Bayes or a Random Forest
#' method to infer models for predicting classes of test sets or new instances of interest.
#' For example, the method could predict the class, receptor 
#' status phenotypes, of a breast cancer patient given a large gene expression
#' training data.

#' @author Adabor ES (emmanuelsadabor@gimpa.edu.gh), Acquaah-Mensah GK (George.Acquaah-Mensah@mcphs.edu), 
#' Mazandu GK (kuzamunu@aims.ac.za)
#' 
#' @references Adabor ES, Acquaah-Mensah GK. 
#' Machine learning approaches to decipher hormone
#' and HER2 receptor status phenotypes in breast cancer.
#' Brief Bioinf., 2017.doi: 10.1093/bib/bbx138

#' @examples 
#' data("her2")
#' data("testset")
#' Model <- MSclassifier(X=her2, testset=NULL, method="MSRandomForest")
#' predictions <- predict(Model, newdata=testset)
#' head(predictions)
#' table(predictions, testset$her2_status)
#' Predictions <- MSclassifier(her2,testset=tetset, method="MSRandomForest")
#' # MSNaiveBayes
#' Model <- MSclassifier(her2, testset=NULL, method="MSNaiveBayes")
#' predictions <- predict(Model, newdata=testset)
#' Predictions <- MSclassifier(her2,testset=tetset, method="MSNaiveBayes")
#' @export

#####################################################################################
MSclassifier <- function(X, testset=NULL, method="MSRandomForest"){
  # Testing the inputs
  
  if ( method == "MSRandomForest" | method == "MSNaiveBayes"){
    p=1
  }
  else { 
    stop("Please specify  MSRandomForest or MSNaiveBayes as method")}
  
  statuses = unique(as.vector(X[,ncol(X)]))
  if (length(statuses) !=2){
    
    stop("Data set should have two groups in order to 
         apply median-supplement method")
  }
  ############################################################################
  
  n = ncol(X) - 1  
  pos = which(as.vector(X[,ncol(X)]) == statuses[1])
  neg = which(as.vector(X[,ncol(X)]) == statuses[2])
  m = abs(length(pos) - length(neg))
  if ( m==0 ){
    stop("Data is a balanced set")}
  
  #######################################################################################
  ## Generating the random matrix
  set.seed(1); mylist=lapply(1:m, function(i) runif(n,0,1)) # set seed makes results reproducible
  mytab = do.call(rbind,mylist)
  #######################################################################################
  # Setting the labels of random matrix
  
  if (length(pos)<length(neg)){
    Label = statuses[1]
    newX = X[pos,]
  } else
  {Label = statuses[2]
  newX = X[neg,]}
  
  ########################################################################################
  # Median supplement data
  rdata=mytab
  for (i in 1:n){
    a = median(newX[,i])
    rdata[,i] = a*rdata[,i]
  }
  rdata = cbind(rdata,rep(Label,m))
  colnames(rdata) = colnames(X)
  rtX = data.frame(rbind(as.matrix(X),as.matrix(rdata)))  # median-supplement data
  colnames(rtX)[ncol(rtX)] = "Status"
  #######################################################################################
  # Inferring models
  #Alternative-- works well. Matrix equality have beeen proven
  rm1 = rtX[,-ncol(rtX)]
  crm1 = apply(rm1,2,as.numeric)
  rm2 =  rtX[,ncol(rtX)]
  training_set = data.frame(crm1,Status=rm2)
  
  # Median-supplement Naive Bayes
  if (method == "MSNaiveBayes"){
    # Median-supplement Naive Bayes
    
    Model = e1071::naiveBayes(Status~., data = training_set) # Data is a dataframe
  
  }
  else
  {  # Median-supplement Random forest
    
    Model = randomForest::randomForest(Status~., data = training_set)
    
  }
  
  results = Model
  
  if (!is.null(testset)){
    Group = predict(Model, newdata=testset[,-ncol(testset)])
    names(Group) = paste("Sample",1:length(Group),sep="")
    #results = table(py, testset[,ncol(testset)])
    results = Group
  }
  #return(Model)
  return(results)
  
  
  }


