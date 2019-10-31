#!/usr/bin/env Rscript

### use these attributes from VCF file to make a model
attributes = c("QD","ReadPosRankSum", "QD", "FS", "ED", "VPR", "VAF", "VMMF", "SPLICEADJ", "RPT", "Homopolymer", "Entropy", "RNAEDIT")

stdin = commandArgs(TRUE)
if(length(stdin) < 2 || length(stdin) > 3) {
    stop(sprintf("\n\nIncorrect number of arguments. \nUSAGE: RVboost.R input_matrix output_file [attrs=\"%s\"] \n\n", paste(attributes, collapse=",")))
}

###arguments
input_matrix = stdin[1]
output = stdin[2]




if (length(stdin) == 3) {
    attributes = strsplit(stdin[3], ",")[[1]]
}
message("Using attribute list: ", paste(attributes, collapse=","))

library(gbm)

data = read.table(input_matrix, header=T, row.names=1, stringsAsFactors = FALSE)

if (! all(attributes %in% colnames(data))) {
    missing_atts = attributes[ ! attributes %in% colnames(data) ]
    stop(paste("Error, missing attributes in data matrix:", missing_atts, sep=","))
}


RS_col = which(colnames(data) %in% "RS")
RS = data[, RS_col, drop=T]
RS = ifelse(RS=="Y", 1, 0)

data = data[,-RS_col]

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
    data[is_NA_rnaedit,RNAEDIT_col] = 0

    if (sum(data[,RNAEDIT_col]) == 0) {
        message("warning, no RNAEDIT events assigned. Removing column")
        data = data[,-RNAEDIT_col]
    }
}


############################
## Run adaboost

message("Running adaboost - rvboost-style")

NUMTREES = 5e3  #rvboost defaults

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

write.table(result_table, file=output, quote=F, row.names=F)
