authorships_raw <- read.csv("../data/authorships_with_genderizer.csv", header = T)
authorships_raw <- authorships_raw[authorships_raw$year >= 1960,]
lined_data <- read.csv("../data/citation_between_smallest_cluster.txt", header = F)
sampled.alpha.mix <- read.csv("../data/toplevel_022019.csv", header = T)

## Update ##
p.vals <- apply(sampled.alpha.mix[c(1, 5001:35000),], MAR = 2, FUN = function(x){mean(x[-1] >= x[1], na.rm = T)})
p.val.names <- gsub("X", "", gsub(".", ":", names(p.vals), fixed = T), fixed = T)

field.char <- read.csv("../data/field_characteristics.csv")


fields <- read.table("../data/fields.txt", header = T, sep = "\t")




fieldTable <- data.frame(cluster = field.char$cluster, 
                         pctUnest = rep(0, dim(field.char)[1]),
                         pctFem = rep(0, dim(field.char)[1]),
                         numAuth = rep(0, dim(field.char)[1]),
                         papers = rep(0, dim(field.char)[1]),
                         pctSingleFem = rep(0, dim(field.char)[1]),
                         pctMultiFem = rep(0, dim(field.char)[1]),
                         pVal = p.vals[match(field.char$cluster, p.val.names)],
                        pValBY = p.adjust(p.vals[match(field.char$cluster, p.val.names)], method = "BY"),
                        pValBH = p.adjust(p.vals[match(field.char$cluster, p.val.names)], method = "BH"),
                         expAlpha = colMeans(sampled.alpha.mix[-c(1:1000), ], na.rm = T)[match(field.char$cluster, p.val.names)],
                         obsAlpha = rep(0, length(field.char$cluster)), 
                         termField = field.char$term.field)

levelClust <- matrix(0,nrow =dim(authorships_raw)[1], ncol = 6)

for(level in 1:6){
  levelClust[, level] <- sapply(strsplit(as.character(authorships_raw$cluster), ":"), function(x){paste(x[1:min(level, length(x))], collapse = ":")})
}



getInfo <- function(clusterID){
  level <- length(strsplit(clusterID, ":")[[1]])
  tab <- authorships_raw[which(levelClust[, level] == clusterID), ]
  size <- dim(tab)[1]
  unest <- mean(tab$gender == 10)
  female <- mean(tab$gender == 1)
  papers <- length(unique(tab$pID))
  pctSingleFem <- mean(tab$gender[which(tab$authorCount ==1)] == 1)
  pctMultiFem <- mean(tab$gender[which(tab$authorCount >1)] == 1)
  return(c(unest, female, size, papers, pctSingleFem, pctMultiFem))
}

getInfoJSTOR <- function(){
  tab <- authorships_raw
  size <- dim(tab)[1]
  unest <- mean(tab$gender == 10)
  female <- mean(tab$gender == 1)
  papers <- length(unique(tab$pID))
  pctSingleFem <- mean(tab$gender[which(tab$authorCount ==1)] == 1)
  pctMultiFem <- mean(tab$gender[which(tab$authorCount >1)] == 1)
  return(c(unest, female, size, papers, pctSingleFem, pctMultiFem))
}


output <- sapply(as.character(field.char$cluster), getInfo)
fieldTable[, 2:7] <- t(output)
fieldTable[1765, 2:7] <- getInfoJSTOR()
outputObsAlpha <- sapply(as.character(field.char$cluster), FUN = function(x){
  if(x %in% p.val.names){return(sampled.alpha.mix[1, match(x, p.val.names)])}
  else{return(NA)}})
fieldTable$obsAlpha <- outputObsAlpha
head(fieldTable)



write.csv(fieldTable, "fieldChar_v2_post60.csv", row.names = F)

fieldTable <- read.csv("fieldChar_v2_post60.csv", header = T)




sum(fieldTable$pValBY[fieldTable$termField == 1] < .05, na.rm = T)
sum(fieldTable$pValBY[fieldTable$termField == 0 & !(fieldTable$cluster %in% c(1:25))] < .05, na.rm = T)
sum(fieldTable$pValBY[fieldTable$termField == 0 & (fieldTable$cluster %in% c(1:25))] < .05, na.rm = T)

# Testing at terminal fields
fieldTableTerm <- fieldTable[which(fieldTable$termField == 1 ), ]
ratio <- fieldTableTerm$pctSingleFem / fieldTableTerm$pctMultiFem 
majorityFemale <- fieldTableTerm$pctFem / (1- fieldTableTerm$pctUnest) > .5
pctFemale <- c(fieldTableTerm$pctFem / (1- fieldTableTerm$pctUnest))

ids <- sapply(strsplit(as.character(fieldTableTerm$cluster), ":"), function(x){x[1]})
gee_mod <- gee::gee(c(fieldTableTerm$pValBY < .05) ~ log(fieldTableTerm$numAuth)
    + pctFemale * majorityFemale + ratio, id = ids, family = binomial, silent = T)
summary(gee_mod)$co

pnorm(abs(c(summary(gee_mod)$coef[,5])), lower.tail = F) * 2
xtable::xtable(cbind(summary(gee_mod)$co[, -c(2:3)], round(pnorm(abs(c(summary(gee_mod)$coef[,5])), lower.tail = F) * 2,3)))
 

glm_mod <- glm(c(fieldTableTerm$pValBY < .05) ~ log(fieldTableTerm$numAuth)
                    + pctFemale * majorityFemale + ratio, family = binomial)
summary(glm_mod)
