library(doParallel)
library(igraph)
library(cppSampler)

imputationNum <- 1
##### Read raw data with genders from genderizer included #####
# authorships_raw <- read.csv("~/authorships/data/authorships_with_genderizer.csv", header = T)
# imputedGenders <- read.csv(paste("~/authorships/lowHomophily_",imputationNum,".csv", sep =""), header = T)[,2]
authorships_raw <- read.csv("../data/authorships_with_genderizer.csv", header = T)
imputedGenders <- read.csv(paste("lowHomophily_",imputationNum,".csv", sep =""), header = T)[,2]


authorships_raw$gender <- imputedGenders
##### Read citation flows #####
# lined_data <- read.csv("~/authorships/data/citation_between_smallest_cluster.txt", header = F)
lined_data <- read.csv("../data/citation_between_smallest_cluster.txt", header = F)
cross_citation_cutoff <- .05

###### Prune out authorships which do no qualify ######
## Remove authorships for which we are still unable to determine gender
## Remove any papers before 1960
authorships_noAmbig <- authorships_raw[-which(authorships_raw$gender == 10 | authorships_raw$year < 1960), ]
dim(authorships_noAmbig)


## Remove single author papers
## Note by the ordering of the removals, if there is a 2 author paper for which
## we cannot determine the gender; we end up removing that paper entirely
author_count <- table(authorships_noAmbig$pID)
single_author_papers <- as.numeric(names(author_count))[author_count == 1]
authorships_cleaned <- authorships_noAmbig[-which(authorships_noAmbig$pID %in% single_author_papers),]
rm(authorships_raw, authorships_noAmbig)


##### authorships_cleaned now holds data of multi-author papers and only authors with an identified gender ###
##### Process at citation graph #####
unique_mini_fields <- unique(lined_data$V1)
cluster.as.char <- as.character(unique_mini_fields)

# given a raw cluster id, it returns the clean  mini.field to which the author row belongs
# These are the mini-fields that appear in our citation network graph
# So if 10:10 is in the mini-field citation network, a paper with clustering 10:10:1:2, will
# be mapped into 10:10. Note, we match to the most detailed mini-field possible
get.minifield.id <- function(cluster.id){
  split.up <- unlist(strsplit(cluster.id, split = ":"))
  for(i in c(length(split.up):1)){
    mini_field_index <- which(cluster.as.char == paste(split.up[1:i], collapse = ":"))
    if(length(mini_field_index) > 0){
      return(cluster.as.char[mini_field_index])
    }
  }
  return("0:0")
}

# Authorships not in a mini-field in the citation network are denoted 0:0
mini_fields_raw <- sapply(as.character(authorships_cleaned$cluster), get.minifield.id)

original_mini_field <- mini_fields_raw[which(mini_fields_raw != "0:0")]
# Only relevant observational units are stored in the following vectors
field.key <- data.frame(orig_mf = original_mini_field,
                   numer_mf = as.numeric(as.factor(original_mini_field)))


gender <- authorships_cleaned$gender[which(mini_fields_raw != "0:0")]
paper_id <- authorships_cleaned$pID[which(mini_fields_raw != "0:0")]

# The mini-fields left in our data
mf_left <- unique(field.key)
mf_left <- mf_left[order(mf_left$numer_mf, decreasing = F),]

#### Get cross citation graphs ####

## Get adjacency matrix ##
## Note that this currently contains all mini-fields, but some have been dropped out of our
## data
adj_matrix <- matrix(lined_data[,3], nrow = length(unique_mini_fields), ncol = length(unique_mini_fields))
colnames(adj_matrix) <- rownames(adj_matrix) <- as.character(unique_mini_fields)

# now only include mini-fields left
adj_matrix <- adj_matrix[as.character(mf_left$orig_mf), as.character(mf_left$orig_mf)]

## original matrix has counts, make into proportions
adj_matrix <- adj_matrix / rowSums(adj_matrix)
## Take out anything under the cutoff, then renormalize
adj_matrix <- ifelse(adj_matrix > cross_citation_cutoff, adj_matrix, 0)
adj_matrix <- adj_matrix / rowSums(adj_matrix)
## ensure symmetric support
adj_matrix <- (adj_matrix + t(adj_matrix)) / 2
adj_matrix <- adj_matrix / rowSums(adj_matrix)



#### Get connected components of the graph ###
edge_set <- c()
for(j in 1:dim(adj_matrix)[1]){
  tuples <- paste(j, which(adj_matrix[j,] > 0) )
  pairs <- as.numeric(unlist(strsplit(tuples, " ")))
  edge_set <- c(edge_set, pairs)
}

g <- make_undirected_graph(edge_set)
clust <- clusters(g)
cluster.assignment <- clust$membership[field.key$numer_mf]

## list of mini-fields in each connected component
connected_components <- list()
for(i in 1:max(clust$membership)){
  connected_components[[i]] <- mf_left$numer_mf[which(clust$membership == i)]
}


##############################################################################
authorshipsTableList <- vector("list", length(connected_components))
count <- 1
for(i in 1:length(connected_components)){
    ind <- which(cluster.assignment == i)
    mat <- as.matrix(cbind(orig_mf = field.key$numer_mf[ind],
                      current_mf = field.key$numer_mf[ind],
                      gender = gender[ind],
                      paper_id = paper_id[ind]))
    authorshipsTableList[[i]] <- mat
}







split.clust <- strsplit(as.character(mf_left$orig_mf), ":")
## Gets all possible splits of the cluster
all_subsets <- function(x){
  depth <- length(x)
  ret <- rep("", depth)
  for(i in 1:depth){
    ret[i] <- paste(x[1:i], collapse = ":")
  }
  return(ret)
}

all_relevant_units <- unique(unlist(sapply(split.clust, all_subsets)))



create.relField.object <- function(relField.id){
  depth <- length(unlist(strsplit(relField.id, ":")))
  names.at.depth <- sapply(split.clust,
                           function(x){paste(x[1:min(depth, length(x))], collapse = ":")})
  indices <- which(names.at.depth == relField.id)
  field_names <- mf_left$numer_mf[indices]
  return(field_names)
}

relField_list <- lapply(all_relevant_units, create.relField.object)
names(relField_list) <- all_relevant_units

splitter <- function(x){
  return(split(as.data.frame(x), x[, 2]))
}


###########################################################################################
sim.size <- 10000

# Nslots <- as.numeric(Sys.getenv("NSLOTS"))
Nslots <- 8
print(sprintf("%d slots were allocated", Nslots))
cl <- makeCluster(6)

registerDoParallel(cl)
clusterEvalQ(cl, library(cppSampler))

singleTab <- unlist(parLapply(cl = cl, authorshipsTableList, splitter), recursive = F)
# singleTab <- unlist(lapply(authorshipsTableList, splitter), recursive = F)
singleTab <- singleTab[order(as.numeric(names(singleTab)))]
singleTab <- sapply(singleTab, data.matrix)

# Calculate Alpha for observed data
alpha <- matrix(parLapply(cl = cl, relField_list, calcAlpha, singleTab), nrow =1)
colnames(alpha) <- names(relField_list)

write.table(alpha, file = paste("~/authorships/results/lowHomophily_",imputationNum,".csv", sep = ""), append = F, row.names = F,
            sep = ",", col.names = T)

out <- authorshipsTableList



for(s in 1:sim.size){
  # Permute Everything
  # Out will be the permuted components
  out <- parLapply(cl = cl, out, permuteField, transMat = adj_matrix, cycleProb = .5, max_cycle = 100, thinning = 1e5)


  # singleTab is a list of data matrices
  # that is split by current mini-field
  singleTab <- unlist(parLapply(cl = cl, out, splitter), recursive = F)
  singleTab <- singleTab[order(as.numeric(names(singleTab)))]
  singleTab <- sapply(singleTab, data.matrix)

  # Calculate Alpha for everything
  alpha <- matrix(unlist(parLapply(cl = cl, relField_list, calcAlpha, singleTab)), nrow = 1)
  colnames(alpha) <- names(relField_list)

  write.table(alpha, file = paste("~/authorships/results/lowHomophily_",imputationNum,".csv", sep = ""), append = T, row.names = F,
              sep = ",", col.names = F)

}

stopCluster(cl)
