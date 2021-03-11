#!/usr/bin/env Rscript

### use these attributes from VCF file to make a model
# attributes = c("QD","ReadPosRankSum", "FS", "VPR", "VAF", "VMMF",
#                "SPLICEADJ", "RPT", "Homopolymer", "Entropy", "RNAEDIT", "INDEL")

attributes = c("DJ","ReadPosRankSum","QD","FS","ED","PctExtPos","RS")

args = commandArgs(TRUE)
if(length(args) < 2 || length(args) > 3) {
    stop(sprintf("\n\nIncorrect number of arguments. \nUSAGE: RVboost.R input_matrix output_file [attrs=\"%s\"] \n\n", paste(attributes, collapse=",")))
}

###arguments
input_matrix = args[1]
output = args[2]


if (length(args) == 3) {
    attributes = strsplit(args[3], ",")[[1]]
}

if (! 'RS' %in% attributes) {
    message("Adding RS to attributes list, as it is essential")
    attributes = c(attributes, 'RS')  # required
}


message("Using attribute list: ", paste(attributes, collapse=","))

library(gbm)

data = read.table(input_matrix, header=T, row.names=1, stringsAsFactors = FALSE)

# Remove GT_1/2 if it is not present in the features
if (!("GT_1/2" %in% colnames(data))){
    attributes = attributes[attributes != "GT_1/2"]
}

if (! all(attributes %in% colnames(data))) {
    missing_atts = attributes[ ! attributes %in% colnames(data) ]
    stop(paste("Error, missing attributes in data matrix:", missing_atts, sep=","))
}


data = data[, attributes, drop=F] # restrict to what we want to analyze here and reorder columns

RS_col = which(colnames(data) %in% "RS")
RS = data[, RS_col, drop=T]
RS = ifelse(is.na(RS), 0, 1)

data = data[,-RS_col, drop=F]

## reset attributes sans RS
attributes = colnames(data)


###########################
## adjust data where needed.
if ("ReadPosRankSum" %in% attributes) {
    ReadPosRankSum_col = which(attributes %in% "ReadPosRankSum")
    data[,ReadPosRankSum_col] = abs(data[,ReadPosRankSum_col])
    is_NA_booleans = (is.na(data[,ReadPosRankSum_col]) | is.null(data[,ReadPosRankSum_col]) )
    data[is_NA_booleans, ReadPosRankSum_col] = median(data[ ! is_NA_booleans, ReadPosRankSum_col])
}

if ("SPLICEADJ" %in% attributes) {
    SPLICEADJ_col = which(attributes %in% "SPLICEADJ")
    data[,SPLICEADJ_col][is.na(data[,SPLICEADJ_col])] = -1
}

if ("RPT" %in% attributes) {
    RPT_col = which(attributes %in% "RPT")
    is_NA_rpt = is.na(data[,RPT_col])

    data[is_NA_rpt, RPT_col] = 0
    data[! is_NA_rpt, RPT_col] = 1
}

if ("Homopolymer" %in% attributes) {
    Homopolymer_col = which(attributes %in% "Homopolymer")
    is_NA_homopolymer = is.na(data[,Homopolymer_col])
    data[is_NA_homopolymer, Homopolymer_col] = 0
}

if ("RNAEDIT" %in% attributes) {
    RNAEDIT_col = which(attributes %in% "RNAEDIT")
    is_NA_rnaedit = is.na(data[, RNAEDIT_col])
    data[,RNAEDIT_col] = 0
    data[! is_NA_rnaedit,RNAEDIT_col] = 1

    if (sum(data[,RNAEDIT_col]) == 0) {
        message("warning, no RNAEDIT events assigned. Removing column")
        data = data[,-RNAEDIT_col]
    }
}

# Set values with NA to the median value 
for(j in 1:ncol(data)){
    is_na <- which(is.na(data[ ,j]))
    
    if(length(is_na) > 0){
        # get the median of the non NA values 
        not_na <- which(!is.na(data[ ,j]))
        median_value <- median(data[not_na,j])
        # set the Na's to the median 
        data[is_na,j] <- median_value
    }
}
############################
## Run adaboost

message("Running adaboost - rvboost-style")

NUMTREES = 2e4  #rvboost defaults

data = data.matrix(data)

res = gbm.fit(x=data,
              y=RS,
              n.trees=NUMTREES,
              interaction.depth=2,
              distribution="adaboost",
              verbose=FALSE)


## as per rvboost:
## convert it to the 0-1 scale since the adaboost method gives the predictions on logit scale.
## http://stats.stackexchange.com/questions/37497/how-to-use-r-gbm-with-distribution-adaboost

fitted_values <-  plogis(res$fit)

ecdf_func <- ecdf(fitted_values[which(RS==1)])

# apply all scores to ecdf_func
fitted_value_scores = ecdf_func(fitted_values)

result_table = data.frame(variants=rownames(data), RVBfitval=fitted_values, RVBscore=fitted_value_scores)

write.table(result_table, file=output, quote=F, row.names=F, sep="\t")
