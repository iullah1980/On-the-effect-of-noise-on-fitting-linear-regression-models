
# Load necessary libraries
library(SNPRelate) # For SNP-related data manipulation
library(plyr) # For data manipulation
library(tidyverse) # For data manipulation and visualization
library(snpStats) # For SNP statistics

# create download directory and set it
.exdir = 'tmp'
dir.create(.exdir)
.file = file.path(.exdir, 'tar.gz')


# download file
url = 'http://www.ricediversity.org/data/sets_hydra/HDRA-G6-4-RDP1-RDP2-NIAS-2.tar.gz'
options(timeout=500)
download.file(url, .file)

# untar it
untar(.file, exdir = path.expand(.exdir))

snpgdsBED2GDS("tmp/HDRA-G6-4-RDP1-RDP2-NIAS/HDRA-G6-4-RDP1-RDP2-NIAS.bed",
              "tmp/HDRA-G6-4-RDP1-RDP2-NIAS/HDRA-G6-4-RDP1-RDP2-NIAS.fam", 
              "tmp/HDRA-G6-4-RDP1-RDP2-NIAS/HDRA-G6-4-RDP1-RDP2-NIAS.bim", out="obGDS")

# Open the GDS file
genofile <- snpgdsOpen("obGDS")

# Extract genotype data from the GDS file
gds <- snpgdsGetGeno(genofile, with.id = TRUE)

# Close the GDS file to free up resources
closefn.gds(genofile)

# Prepare the genotype data for analysis
genotype_data <- gds$genotype
colnames(genotype_data) <- gds$snp.id
rownames(genotype_data) <- gds$sample.id

# Read PLINK files for additional analysis
path <- "tmp/HDRA-G6-4-RDP1-RDP2-NIAS"
ob.plink <- read.plink(file.path(path, "HDRA-G6-4-RDP1-RDP2-NIAS"))
names(ob.plink)

# Extract genotype data from the PLINK object
ob.geno <- ob.plink$genotypes
ob.geno

# Extract SNP annotations from the PLINK object
annotation <- ob.plink$map
head(annotation)

# download file
url = 'http://www.ricediversity.org/data/sets_hydra/phenoAvLen_G6_4_RDP12_ALL.tar.gz'
download.file(url, .file)

# untar it
untar(.file, exdir = path.expand(.exdir))

# Read phenotype data
path1 <- "tmp/phenoAvLen_G6_4_RDP12_ALL"
ob.pheno <- read.table(file.path(path1, "phenoAvLen_G6_4_RDP12_ALL.txt"),header = T)
head(ob.pheno)
colnames(ob.pheno)=c("ID","IID","Length")

# download file
.file = file.path(.exdir, 'FINAL6_Suppl_Table_1_Germplasm_20160725.xls')
url = 'http://www.ricediversity.org/data/sets_hydra/FINAL6_Suppl_Table_1_Germplasm_20160725.xls'
download.file(url, .file)

# Read germplasm data
# Before proceeding, ensure that the downloaded .xls file is saved as a .csv file.
path2 <- "tmp/"
ob.Germplasm <- read.csv(file.path(path2, "FINAL6_Suppl_Table_1_Germplasm_20160725.csv"), header = TRUE)
head(ob.Germplasm)

# Processing and merging phenotype and germplasm data based on ID matches
ID <- gsub('[^0-9]', '', ob.Germplasm$IRGC.ID)
ob.mod <- data.frame(ob.Germplasm,ID)
total <- merge(ob.mod,ob.pheno,by="ID")


# Filtering data based on genotype assay ID availability and subpopulation calls
filteredDataWithAssayID <- total[!is.na(total$HDRA.genotype.assay.ID),]
rownames(filteredDataWithAssayID) = filteredDataWithAssayID$HDRA.genotype.assay.ID

filteredDataNSFTV <- subset(ob.Germplasm, ob.Germplasm$Other.accession.ID %in% intersect(ob.pheno$ID, ob.Germplasm$Other.accession.ID))
NSFTV_MergedWithPheno <- data.frame(filteredDataNSFTV, ID=filteredDataNSFTV$Other.accession.ID)
totalMergedData <- merge(NSFTV_MergedWithPheno, ob.pheno, by="ID")

filteredDataWithAssayID2 <- totalMergedData[!is.na(totalMergedData$HDRA.genotype.assay.ID),]
rownames(filteredDataWithAssayID2) = filteredDataWithAssayID2$HDRA.genotype.assay.ID

# Combining filtered datasets and further filtering by subpopulation call
combinedFilteredData <- rbind(filteredDataWithAssayID, filteredDataWithAssayID2)
indicaSubpopData <- combinedFilteredData[combinedFilteredData$fastStructure.subpopulation.call %in% "indica",]

# Matching IDs between phenotype/genotype data and subsetting genotype data accordingly
ids <- intersect(rownames(indicaSubpopData), rownames(genotype_data))
finalGenotypeData <- genotype_data[ids, ]
indicaSubpopData <- indicaSubpopData[ids,]

# eno subset based on matching IDs
finalGenoSubset <- ob.geno[ids,]

# Performing summary statistics on the SNP data to determine call rate and MAF
info.snps <- col.summary(finalGenoSubset)

# Defining criteria for SNP quality control (QC): Call rate > 0.95 and MAF > 0.05
use <- info.snps$Call.rate > 0.95 & info.snps$MAF > 0.05
mask.snps <- use & !is.na(use)

# Applying QC filters to the genotype data
geno.qc.snps <- finalGenoSubset[ , mask.snps]

# Updating the annotation data to include only SNPs that passed QC
annotation.qc <- annotation[mask.snps, ]

# Calculating and displaying the number of SNPs removed due to QC criteria
# This includes those removed for bad call rate and those removed for low MAF
numSNPsRemovedCallRate <- sum(info.snps$Call.rate < 0.95)
numSNPsRemovedMAF <- sum(info.snps$MAF < 0.05, na.rm=TRUE)
totalSNPsNotPassQC <- sum(!mask.snps)


# Summarizing individual data based on the QC-passed genotype data
info.indv <- row.summary(geno.qc.snps)


# Plotting the distribution of grain length from the phenotype data
par(mar=c(4, 5, 2, 1))
hist(indicaSubpopData$Length, xlab="", ylab="", cex.axis=1.5, las=1, main="", breaks=30)
title(xlab="Average grain length (mm)", cex.lab=1.5, line = 2.8)
title(ylab="Frequency", cex.lab=1.5, line=3.8)

# Preparing the genotype data for regression analysis by applying QC filters
genoForAnalysis <- finalGenotypeData[, mask.snps]

# Scaling genotype data after centering SNP means and missing values imputation
CHR <- annotation.qc$chromosome
genoChromosome3 <- genoForAnalysis[, annotation.qc$chromosome == 3]
meanAdjustedGeno <- apply(genoChromosome3, 2, mean, na.rm=TRUE)
centeredGeno <- sweep(genoChromosome3, 2, meanAdjustedGeno)
centeredGeno[is.na(centeredGeno)] <- 0
scaledGeno <- scale(centeredGeno)

# Preparing phenotype data for regression analysis
scaledPhenotype <- scale(indicaSubpopData$Length, center = TRUE, scale = FALSE)

# Setting up the training data size and number of shuffles for the regression analysis
n <- 300 # size of training data
n_rpt <- 100 # number of shuffles

# Select every 5th marker to reduce high multicollinearity 
select_index <- seq(1, dim(scaledGeno)[2], 5)
selectedData <- data.frame(scaledGeno[,select_index])


# Calculate OLS and minimum norm OLS estimates, and evaluate error
ridgeless_reg_error <- function(xTrain, xTest, yTrain, yTest){
  # Calculate regression coefficients using the Moore-Penrose pseudoinverse
  # b <- ginv(crossprod(xTrain))%*%crossprod(xTrain, yTrain)
  b <- solve(crossprod(xTrain)+.00001*diag(ncol(xTrain)))%*%crossprod(xTrain, yTrain)
  # The constant .00001 is important as setting it to zero or a smaller value 
  # may lead to an unexpected jumps in the test error. This adjustment helps stabilize 
  # the calculation by ensuring the matrix inversion process is numerically stable.
  predTrain <- xTrain%*%b
  predTest <- xTest%*%b
  
  lossTrain <- sqrt(mean((yTrain - predTrain[,1])^2))
  lossTest <- sqrt(mean((yTest - predTest[,1])^2))
  list(lossTrain=lossTrain, lossTest=lossTest)
}


train_test_rmses <- function(n, x, y){
  
  # myData is the dataset prepared for analysis, currently set to 'x'. This placeholder allows for potential 
  # preprocessing or subsetting of the data before running the regression analysis.
  myData <- x
  
  # Define the number of columns in 'myData' to determine the range of features to evaluate in the regression model.
  n_cols_myData <- ncol(myData)
  
  # 'd_evaluate' defines the set of indices to evaluate in the regression model, starting from simple models 
  # with fewer predictors and gradually including more to see the impact on model performance.
  d_evaluate <- c(c(1,2,3,10,15,seq(20, 500, 10)), seq(525, dim(myData)[2], 25))
  
  # Creating a model matrix from 'myData', treating it as predictors in a linear model. This step prepares 
  # the data for regression analysis by encoding any categorical variables and adding an intercept term.
  X <- model.matrix(~., data=myData)
  
  # Selecting a subset of rows from 'X' randomly to serve as the training set, with 'n' specifying the 
  # size of this subset. This approach allows for validation of the model on unseen data.
  train_index <- sample(1:nrow(myData), n, replace = FALSE)
  xTrain = X[train_index, ]
  xTest = X[-train_index, ]
  
  # Subsetting 'y', the response variable, according to the indices chosen for the training and test sets.
  yTrain = y[train_index]
  yTest = y[-train_index]
  
  # Applying the 'ridgeless_reg' function over a range of model complexities ('d_evaluate') and 
  # computing the root mean square errors (RMSEs) for both the training and test sets. This process 
  # evaluates how well the models with different numbers of predictors perform.
  rmses <- map_df(d_evaluate,
                  ~ridgeless_reg_error(xTrain[, 1:.x, drop=FALSE],
                                 xTest[, 1:.x, drop=FALSE],
                                 yTrain,
                                 yTest)
  )
  
  # Returning the results as a tibble, with 'd_evaluate' indicating the number of predictors used in 
  # each model and 'rmses' storing the corresponding RMSEs for training and test sets.
  tibble(d_evaluate, rmses)
}



# Running the regression analysis and summarizing the results
results <- map_df(1:n_rpt, 
                  ~ train_test_rmses(n=n, x=selectedData, y=scaledPhenotype),
                  .id = "runs")

test_error_summary = results %>% 
  group_by(d_evaluate) %>% 
  dplyr::summarize(
    rmse_train = mean(lossTrain, na.rm = TRUE),
    rmse_test = mean(lossTest, na.rm = TRUE)
  ) %>% 
  pivot_longer(-d_evaluate, names_to = 'lossType', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(lossType, 'loss')) 

# Plotting the results
g <- results %>%
  pivot_longer(-c(runs, d_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "loss")) %>% 
  ggplot(aes(x=d_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, linewidth=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n, color="gray",lwd=1) +
  geom_hline(yintercept=sqrt(mean(scaledPhenotype^2)), color="black",lwd=1) +
  coord_cartesian(ylim = c(0, 5)) +#max(test_error_summary$RMSE) + .5)) +
  labs(x = 'Number of explanatory variables (d)', y = 'Error', title = "") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))+
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = test_error_summary %>% filter(lossType == 'rmse_test') %>% 
      filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), 
        plot.caption = element_text(hjust = 0),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14))
g


png(file="applicationGWASplot.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()


