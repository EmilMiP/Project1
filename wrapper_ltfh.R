'Liability Threshold model conditional on Family History (LT-FH) 
Usage:
wrapper_ltfh.R (-h | --help)
wrapper_ltfh.R --popPrev <popPrev> --sampPrev <sampPrev> --obsHeri <obsHeri> --phenoFile <phenoFile> --outDist <outDist> [--debug] [--fixedHeri <fixedHeri>] [--identifier <identifier>] [--childCol <childCol>] [--dadCol <dadCol>] [--momCol <momCol>] [--sibCol <sibCol>] 
Options:
-h --help                         This print
--popPrev <popPrev>               Population prevalence estimate
--sampPrev <sampPrev>             Sample prevalence
--obsHeri <obsHeri>               Observed scale heritability
--phenoFile <phenoFile>           Phenotype file 
--outDist <outDist>               Root directory for the output, file name will be constructed based off of provided phenotype file name
--debug                           Prints all arguements, debugging purpose 
--fixedHeri <fixedHeri>           Fixed value of liability scale heritability to be used instead of calculated
--identifier <identifier>         String snipet to be attacched to the output name, to easier identify, e.g. fixed heritability from calculated heritability output
--childCol <childCol>             [default: CHILD_STATUS]
--dadCol <dadCol>                 [default: P1_STATUS]
--momCol <momCol>                 [default: P2_STATUS]
--sibCol <sibCol>                 [default: SIB_STATUS]
' -> doc
library(docopt)
args = docopt(doc, version = "alpha")

if (args$'--debug') {
  print(args)
}




cat("            ___________________________________ \n")
cat("           |                                   | \n")
cat("           |                                   | \n")
cat("           |                                   | \n")
cat("           |                                   | \n")
cat("           |       LT-FH Command Line          | \n")
cat("           |                                   | \n")
cat("           |                                   | \n")
cat("           |                                   | \n")
cat("           |___________________________________| \n")
cat("\n")
cat("\n")



setwd("/home/emp/faststorage/project1/ltfh/software v2")
source("assign_ltfh.R")

library(data.table)

lia_h2 = function(obs_h2, P, K) {
  #P: fraction of cases in data
  #K: prevalens in population
  z = qnorm(K, lower.tail = F)
  return((K*(1 - K))^2 / (dnorm(z)^2 * (P * (1 - P))) * obs_h2)
}
#lia_h2(0.11, 0.05, 0.05): 0.4912118

if (args$'--outDist' == F) {
  cat("Please Supply an output destination. \n")
  break
}


cat(" ---- LTFH w/ these parameters ---- \n")
cat(paste("Phenotype File:", args$'phenoFile', "\n"))
samp_prev = as.numeric(args$'sampPrev')
child_sib_T = qnorm(samp_prev, lower.tail = F)
pop_prev = as.numeric(args$'popPrev')
parents_T = qnorm(pop_prev, lower.tail = F)

h2 = lia_h2(as.numeric(args$'ohsHeri'), P = samp_prev, K = pop_prev)

cat(paste("Sample prevalens: ", samp_prev, "\n"))
cat(paste("Population prevalens: ", pop_prev, "\n"))
cat(paste("Liability scale heritability: ", h2, "\n"))
cat(paste("Observed scale heritability: ", as.numeric(args$'ohsHeri'), "\n"))

if (args$'childCol' != "CHILD_STATUS") {
  cat("Offspring column name used:", args$'childCol', "\n")
}

if (args$'dadCol' != "P1_STATUS") {
  cat("Father column name used:", args$'dadCol', "\n")
}

if (args$'momCol' != "P2_STATUS") {
  cat("Mother column name used:", args$'momCol', "\n")
}

if (args$'sibCol' != "SIB_STATUS") {
  cat("Sibling column name used:", args$'sibCol', "\n")
}



if (!is.null(args$'--fixedHeri')) {
  h2 = as.numeric(args$'fixedHeri')
  cat(" ----------> FIXED HERITABILITY IS SET TO ", h2, "<---------- \n" )
}

distName = unlist(strsplit(args$'phenoFile', "/"))
distName = distName[length(distName)]
outRoot = args$'outDist' #unlist(strsplit(args$'outDist', "/"))

if (!is.null(args$'--identifier')) {
  outDist = paste(outRoot, "LTFH", args$'identifier', distName, sep = "")
} else {
  outDist = paste(outRoot, "LTFH", distName, sep = "")
}

cat("Printing output to:", outDist, "\n")
pheno = as.data.frame(fread(args$'phenoFile'))

childCol = args$'childCol'
dadCol = args$'dadCol'
momCol = args$'momCol'
sibCol = args$'sibCol'

new_pheno = create_pheno(data = pheno,
                         trait_h2 = h2,
                         T_val_child = child_sib_T,
                         T_val_parent = parents_T,
                         relevant_trait_child = childCol,
                         relevant_trait_dad = dadCol,
                         relevant_trait_mom = momCol,
                         number_siblings_col = "NUM_SIBS",
                         relevant_trait_sib = sibCol,
                         maximum_siblings_to_compute = 4)
new_pheno = cbind(new_pheno[,1], new_pheno)
colnames(new_pheno)[1] = "FID"

fwrite(x = as.data.frame(new_pheno), file = outDist, sep = " ", quote = F, na = "NA")

