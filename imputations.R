#### Load Data #####
authorships <- read.csv("../data/authorships_with_genderizer.csv")
fieldTable <- read.csv("../data/fieldChar_v2_post60.csv", header = T)

# Get list of terminal fields
terminalFieldsTotal <- fieldTable$cluster[which(fieldTable$termField == 1)]




# each authorship to their terminal field
splits <- strsplit(as.character(authorships$cluster), ":")
getSubs <- function(x){sapply(1:length(x), function(z){paste(x[1:z], collapse = ":")} )}
termIndex <- sapply(sapply(sapply(splits, getSubs), match, terminalFieldsTotal), max, na.rm = T)
terminalFieldsAuthorships <- terminalFieldsTotal[termIndex]
indicesToImpute <- which(authorships$gender == 10 & !is.na(terminalFieldsAuthorships))
indicesInFieldTable <- match(terminalFieldsAuthorships[indicesToImpute], fieldTable$cluster)


### Low Homophily Scenario ###
propForAuthorshipsLow <- fieldTable$pctFem[indicesInFieldTable] / (1 - fieldTable$pctUnest[indicesInFieldTable])

set.seed(10101)
for(i in 1:10){
  imputedGenders <- rbinom(length(propForAuthorshipsLow), size = 1, prob = propForAuthorshipsLow)
  genderAll <- authorships$gender
  genderAll[indicesToImpute] <- imputedGenders
  write.csv(genderAll, paste("lowHomophily_", i, ".csv", sep = ""))
}



### High Homophily Scenario ###
# all paperIDs with authorships to imput
pIdMissing <- unique(authorships$pID[indicesToImpute])

# If paper has at least 1 estimated author, get the m/f split
probs <- aggregate(authorships$gender ~ authorships$pID, FUN = mean, subset = authorships$pID %in% pIdMissing & authorships$gender!= 10)
# indices corresponding to at least 1 estimated author
regularInd <- indicesToImpute[which(authorships$pID[indicesToImpute] %in% probs[,1 ])]
regularProbs <- probs[match(authorships$pID[regularInd], probs[, 1]), 2]

# which indices correspond to sampling with same prob as TF and setting all to same
irregularInd <- indicesToImpute[which(!(authorships$pID[indicesToImpute] %in% probs[, 1]))]
# which papers correspond to sampling with same prob as TF and setting all to same
irregularPID <- unique(authorships$pID[irregularInd])
irregularIndicesInFieldTable <- match(terminalFieldsAuthorships[match(irregularPID, authorships$pID)], fieldTable$cluster)

irregularProbs <- fieldTable$pctFem[irregularIndicesInFieldTable] / (1 - fieldTable$pctUnest[irregularIndicesInFieldTable])


set.seed(10101)
for(i in 1:10){
  regularImpute <- rbinom(length(regularProbs), size = 1, prob = regularProbs)
  irregularImpute <- rbinom(length(irregularProbs), size = 1, prob = irregularProbs)
  genderAll <- authorships$gender
  genderAll[regularInd] <- regularImpute
  genderAll[irregularInd] <- irregularImpute[match(authorships$pID[irregularInd], irregularPID)]
  
  write.csv(genderAll, paste("highHomophily_", i, ".csv", sep = ""))
}
