### Table 2 ####
field.char <- read.csv("../data/fieldChar_v2.csv")
fields <- read.table("../data/fields.txt", header = T, sep = "\t")
sampled.alpha.mix <- read.csv("../data/newAlphasParl.csv", header = T)
fdr.rate <- .05
level_1 <- read.csv("../data/newAlphas_level_1.csv", header = T)

top.level.term <- sapply(strsplit(as.character(field.char$cluster[which(field.char$termField == 1)]), ":"), function(x){x[1]})
top.level.comp <- sapply(strsplit(as.character(field.char$cluster[which(field.char$termField == 0)]), ":"), function(x){x[1]})

num.term <- sapply(c(1:25)[-24], function(x){sum(top.level.term == x)})
num.comp <- sapply(c(1:25)[-24], function(x){sum(top.level.comp == x)})

num.sig.term.by <- sapply(c(1:25)[-24], function(x){sum(field.char$pValBY[which(field.char$termField == 1)][which(top.level.term == x)] < fdr.rate, na.rm = T)})
num.sig.comp.by <- sapply(c(1:25)[-24], function(x){sum(field.char$pValBY[which(field.char$termField == 0)][which(top.level.comp == x)] < fdr.rate, na.rm = T)}) 

tops <- c(1:25)
rec <- matrix(0, ncol = 6, nrow = 25)
for(t in c(1:25)[-24]){
  rec[t, 1] <- sampled.alpha.mix[1, paste("X", t, sep = "")]
  rec[t, 2] <- mean(sampled.alpha.mix[-c(2:1000), paste("X", t, sep = "")])
  rec[t, 3] <- mean(level_1[-c(2:1000), paste("X", t, sep = "")])
}

colnames(rec) <- c("obs", "comp", "struct", "obs-mf", "comp-mf", "struct-mf")
rownames(rec) <- fields$label[match(1:25, fields$cluster)]

getPropPapers <- function(alpha, males, females){
  mm = (alpha * males * females + males^2) / (2 * (females + males))
  mf = (1 - alpha) * males * females / (females + males)
  ff = (males + females) / 2 - mm - mf
  return(list(mm = mm,
              mf = mf,
              ff = ff))
}

for(i in c(1:25)[-24]){
  rec[i, 4:6] <- c(getPropPapers(rec[i, 1], (1 - field.char$pctMultiFem[which(field.char$cluster == i)]) * 200, field.char$pctMultiFem[which(field.char$cluster == i)] * 200)$mf,
                   getPropPapers(rec[i, 2], (1 - field.char$pctMultiFem[which(field.char$cluster == i)]) * 200, field.char$pctMultiFem[which(field.char$cluster == i)] * 200)$mf,
                   getPropPapers(rec[i, 3], (1 - field.char$pctMultiFem[which(field.char$cluster == i)]) * 200, field.char$pctMultiFem[which(field.char$cluster == i)] * 200)$mf)
}

rec <- rec[-24,]

tab2 <- data.frame(label = rep("", 24),
                   alpha = apply(weights::rd(rec[, 1:3], digits = 2, add = F), MAR = 1, FUN = paste, collapse = "/"),
                   alphaObs = weights::rd(rec[, 1], digits = 2, add = F),
                   alphaNull = weights::rd(rec[, 2], digits = 2, add = F),
                   alphaNullStruct = weights::rd(rec[, 3], digits = 3, add = F),
                   FM1 = round(rec[, 4], digits = 1),
                   FM2 = round(rec[, 5], digits = 1),
                   FM3 = round(rec[, 6], digits = 1),
                   p.val = rep(0, 24),
                   numSigTermBY = paste(num.sig.term.by, "/", num.term, sep = ""),
                   numSigAggBY = paste(num.sig.comp.by, "/", num.comp, sep = ""),
                   size = field.char$numAuth[match(c(1:25)[-24], field.char$cluster)],
                   pctWomen = field.char$pctFem[match(c(1:25)[-24], field.char$cluster)])

tab2$label <-fields$label[match(c(1:25)[-24], fields$cluster)]
tab2$p.val <- weights::rd(field.char$pValBY[match(c(1:25)[-24], field.char$cluster)],digits = 2, add = F)


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
print(xtable::xtable(tab2.order[order(tab2.order$size, decreasing = T), -c(2, 6:8, 12, 13)]), include.rownames = F)


tab3 <- cbind(tab2.order[, 1], weights::rd(tab2.order[, 6], digits = 1), weights::rd(tab2.order[, 7], digits = 1), weights::rd(tab2.order[, 8], digits = 1),
              weights::rd((tab2.order[, 7] - tab2.order[, 6]) , digits = 1))
print(xtable::xtable(tab3[order(tab3[, 5], decreasing = T), ], digits = 1), include.rownames = F)