#### Checking Stability #####
# field characteristics
field.char <- read.csv("../data/field_characteristics.csv")
# Field cluster and labels
fields <- read.table("../data/fields.txt", header = T, sep = "\t")
# samples from null distribution. First row corresponds to observed alpha
sampled.alpha.mix <- read.csv("../data/toplevel_022019.csv", header = T)
# samples from distribution without compositional homophily
level_1 <- read.csv("../data/newAlphas_level_1.csv", header = T)

fdr.rate <- .05

# p-values for total and subsets of 10K. First ?? are discarded as burn-in
## Update ##
p.vals <- apply(sampled.alpha.mix[c(1, 5001:35000),], MAR = 2, FUN = function(x){mean(x[-1] >= x[1], na.rm = T)})


# cluster labels for each p-value
p.val.fields <- gsub(".", ":",gsub("X", "", names(p.vals), fixed = T), fixed = T)
# 1 indicates a terminal field; 0 indicates non-terminal field
p.val.term <- field.char$term.field[match(p.val.fields, field.char$cluster)]

# p-values adjusted by the BY method
# adjustment is for all tests together, but we then subset into top level fields
# terminal fields and composite fields which do not include top level fields 
adj.p.by <- p.adjust(p.vals, method = "BY")
adj.p.by.top <- adj.p.by[which(gsub(".", ":", gsub("X", "", names(adj.p.by)), fixed = T) %in% c(1:25)[-24])]
adj.p.by.term <- adj.p.by[which(gsub(".", ":", gsub("X", "", names(adj.p.by)), fixed = T) %in% p.val.fields[p.val.term == 1])]
adj.p.by.comp <- adj.p.by[which(!(gsub(".", ":", gsub("X", "", names(adj.p.by)), fixed = T) %in% c(p.val.fields[p.val.term == 1], c(1:25)[-24], "JSTOR")))]



# Total number of significant fields for total and each subset of samples 
## Values for 2nd paragraph 2B ##
sum(adj.p.by.top < fdr.rate)    
sum(adj.p.by.term < fdr.rate, na.rm = T)
sum(!is.na(adj.p.by.term))
sum(adj.p.by.comp < fdr.rate, na.rm = T)
sum(!is.na(adj.p.by.comp))


##### Table 2 and 3 #####
# top level field corresponding to each terminal and composite fields 
top.level.term <- sapply(strsplit(gsub("X", "", names(adj.p.by.term)),".", fixed = T), function(x){x[1]})
top.level.comp <- sapply(strsplit(gsub("X", "", names(adj.p.by.comp)),".", fixed = T), function(x){x[1]})

# number of terminal and composite fields with well defined observed alpha for each top level field
num.term <- sapply(c(1:25)[-24], function(x){sum(top.level.term == x & !is.na(adj.p.by.term))})
num.comp <- sapply(c(1:25)[-24], function(x){sum(top.level.comp == x & !is.na(adj.p.by.comp))})

# number of significant terminal and composite fields for each top level field
num.sig.term.by <- sapply(c(1:25)[-24], function(x){sum(adj.p.by.term[which(top.level.term == x)] < fdr.rate, na.rm = T)})
num.sig.comp.by <- sapply(c(1:25)[-24], function(x){sum(adj.p.by.comp[which(top.level.comp == x)] < fdr.rate, na.rm = T)})


rec <- matrix(0, ncol = 6, nrow = 25)
for(t in c(1:25)[-24]){
  # Observed
  rec[t, 1] <- sampled.alpha.mix[1, paste("X", t, sep = "")]
  # mean of null distribution
  ## Update ##
  rec[t, 2] <- mean(sampled.alpha.mix[c(1, 5001:35000), paste("X", t, sep = "")])
  # mean of distribution without preserving compoisitional. First 1K is burn-in
  rec[t, 3] <- mean(level_1[-c(1:1000), paste("X", t, sep = "")])
}

colnames(rec) <- c("obs", "comp", "struct", "obs-mf", "comp-mf", "struct-mf")
rownames(rec) <- fields$label[match(1:25, fields$cluster)]

# Function to map alpha back to mf papers
getPropPapers <- function(alpha, males, females){
  mm = (alpha * males * females + males^2) / (2 * (females + males))
  mf = (1 - alpha) * males * females / (females + males)
  ff = (males + females) / 2 - mm - mf
  return(list(mm = mm,
              mf = mf,
              ff = ff))
}

for(i in c(1:25)[-24]){
  rec[i, 4:6] <- c(getPropPapers(rec[i, 1], (1 - field.char$pct.coauth.women[which(field.char$cluster == i)]) * 200, field.char$pct.coauth.women[which(field.char$cluster == i)] * 200)$mf,
                   getPropPapers(rec[i, 2], (1 - field.char$pct.coauth.women[which(field.char$cluster == i)]) * 200, field.char$pct.coauth.women[which(field.char$cluster == i)] * 200)$mf,
                   getPropPapers(rec[i, 3], (1 - field.char$pct.coauth.women[which(field.char$cluster == i)]) * 200, field.char$pct.coauth.women[which(field.char$cluster == i)] * 200)$mf)
}

rec <- rec[-24,]

tab2 <- data.frame(label = rep("", 24),
                   alpha = apply(weights::rd(rec[, 1:3], digits = 2, add = F), MAR = 1, FUN = paste, collapse = "/"),
                   alphaObs = weights::rd(rec[, 1], digits = 2, add = F),
                   alphaNull = weights::rd(rec[, 2], digits = 2, add = F),
                   alphaNullStruct = weights::rd(rec[, 3], digits = 2, add = F),
                   FM1 = round(rec[, 4], digits = 1),
                   FM2 = round(rec[, 5], digits = 1),
                   FM3 = round(rec[, 6], digits = 1),
                   p.val = rep(0, 24),
                   numSigTermBY = paste(num.sig.term.by, "/", num.term, sep = ""),
                   numSigAggBY = paste(num.sig.comp.by, "/", num.comp, sep = ""),
                   size = field.char$num.authors[match(c(1:25)[-24], field.char$cluster)],
                   pctWomen = field.char$pct.women[match(c(1:25)[-24], field.char$cluster)])

tab2$label <-fields$label[match(c(1:25)[-24], fields$cluster)]
tab2$p.val <- weights::rd(adj.p.by.top[paste("X", c(1:25)[-24], sep = "")],digits = 2, add = F)


colnames(tab2) <- c("Field", "Alpha","Obs", "Null", "Struct", "Papers (Obs)", "Papers (Null)","Papers (Struct)",  "P-value", "Term", "Comp", "size", "pctWomen")
tab2.order <- tab2[order(tab2$Field), ]

short_names <- c("Anthro",
                 "Classics",
                 "Cog Sci",
                 "Demography",
                 "Eco/Evol",
                 "Economics",
                 "Education",
                 "History",
                 "Law",
                 "Math",
                 "Mol/Cell Bio",
                 "Mycology",
                 "Opr Res",
                 "Org/mkt",
                 "Philosophy",
                 "Phys Anthro",
                 "Plant Phys",
                 "Intl Poli Sci",
                 "US Poli Sci",
                 "Occ Health",
                 "Prob/Stat",
                 "Radiation",
                 "Sociology",
                 "Vet Med")


tab2.order[,1] <- short_names
## Update ##
print(xtable::xtable(rbind(c("JSTOR", weights::rd(sampled.alpha.mix[1, "JSTOR"],2), weights::rd(mean(sampled.alpha.mix[c(5001:35000), "JSTOR"]), 2), ".00", ".00",
                             paste(sum(adj.p.by.term < fdr.rate, na.rm = T), "/",sum(!is.na(adj.p.by.term)), sep = ""),
                             paste(sum(adj.p.by.comp < fdr.rate, na.rm = T), "/", sum(!is.na(adj.p.by.comp)), sep = "")), 
                           as.matrix(tab2.order[order(tab2.order$size, decreasing = T), -c(2, 6:8, 12, 13)])), include.rownames = F))

tab3 <- cbind(tab2.order[, 1], weights::rd(tab2.order[, 6], digits = 1), weights::rd(tab2.order[, 7], digits = 1), weights::rd(tab2.order[, 8], digits = 1),
              weights::rd((tab2.order[, 7] - tab2.order[, 6]) , digits = 1))
print(xtable::xtable(tab3[order(tab3[, 5], decreasing = T), ], digits = 1), include.rownames = F)

## Update ##
c(getPropPapers(sampled.alpha.mix[1, "JSTOR"], (1 - field.char$pct.coauth.women[which(field.char$cluster == "JSTOR")]) * 200, field.char$pct.coauth.women[which(field.char$cluster == "JSTOR")] * 200)$mf,
getPropPapers(mean(sampled.alpha.mix[c(1, 5001:35000), "JSTOR"]), (1 - field.char$pct.coauth.women[which(field.char$cluster == "JSTOR")]) * 200, field.char$pct.coauth.women[which(field.char$cluster == "JSTOR")] * 200)$mf,
getPropPapers(rec[i, 3], (1 - field.char$pct.coauth.women[which(field.char$cluster == "JSTOR")]) * 200, field.char$pct.coauth.women[which(field.char$cluster == "JSTOR")] * 200)$mf)


##### Table 3 #####
fieldTable <- read.csv("../data/fieldChar_v2_post60.csv", header = T)



getSigs <- function(x){
  fields <- gsub(".", ":", gsub("X", "", names(x)), fixed = T)
  term <- mean(x[which(fieldTable$termField[match(fields, fieldTable$cluster)] == 1)] < .05, na.rm = T)
  comp <- mean(x[which(fieldTable$termField[match(fields, fieldTable$cluster)] == 0 & !(fields %in% c(1:25)))] < .05, na.rm = T)
  top <- mean(x[which(fields %in% c(1:25))] < .05, na.rm = T)
  c(term, comp, top)
}

pValSen <- list()
for(impNum in 1:10){
  imputedAlpha <- read.csv(paste("results/lowHomophily_", impNum,".csv", sep = ""), header = T)
  pValSen[[impNum]] <- p.adjust(apply(imputedAlpha[-c(2:1000),], MAR = 2, FUN = function(x){mean(x[-1] >= x[1], na.rm = T)}), method = "BY")
}

low <- sapply(pValSen, getSigs)

pValSen <- list()
for(impNum in 1:10){
  imputedAlpha <- read.csv(paste("results/highHomophily_", impNum,".csv", sep = ""), header = T)
  pValSen[[impNum]] <- p.adjust(apply(imputedAlpha[-c(2:1000),], MAR = 2, FUN = function(x){mean(x[-1] >= x[1], na.rm = T)}), method = "BY")
}

high <- sapply(pValSen, getSigs)

xtable(rbind(rowMeans(low),
             rowMeans(high)))

