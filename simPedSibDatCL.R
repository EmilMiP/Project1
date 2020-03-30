'Generate Genotype Data in .ped/.map format
Usage:
simPedSibDat.R (-h | --help)
simPedSibDat.R --noChildren <noChildren> --M <M> --out <out> --v <v> [--C <C>] [--lia-h2 <lia-h2>] [--K <K>] [--nsib <nsib>] [--nthreads <nthreads>] [--fixedEffect <fixedEffect>] [--baseName <baseName>] [--sizeName <sizeName>] [--identifier <identifier>] [--debug]
Options:
  -h --help                       This print
  --noChildren <noChildren>       Number of children, i.e. individuals to simulate
  --M <M>                         Number of SNPs to simulate per individual
  --out <out>                     path to output destination, file name constructed from settings
  --v <v>                         Number denoting the current version, e.g. v between 1 and 10 for 10 different sets of data.
  --C <C>                         Number of causal SNPs, e.g. M * 0.02 
  --lia-h2 <lia-h2>               heritability on liability scale [default: 0.5]
  --K <K>                         prevalens in the population [default: 0.05]
  --nthreads <nthreads>           Number of threads to use for simulations [default: 10]
  --nsib <nsib>                   Number of siblings per individual [default: 2]
  --fixedEffect <fixedEffect>     Fixed effect sizes(betas) ? [default: FALSE]
  --baseName <baseName>           Base of output file name [default: sibs2]
  --sizeName <sizeName>           Size of genotype matrix [default: 10kx10k]
  --identifier <identifier>       Identifying name, used when all other parameters are the same
' -> doc
library(docopt)
args = docopt(doc, version = "alpha")

if (args$'--debug') {
  print(args)
}
library(MASS)
library(stringr)
library(doSNOW)
library(foreach)
library(parallel)
library(progress)
library(data.table)


#function to generate binary phenotype for parents and offspring
#and unscaled genotype for the offspring.
generateOffspring = function(parents.halved, nsib, M, MAF, environ, lia.T, lia.beta) {
  offspring.geno.unscaled = replicate(1 + nsib, rowSums(round(parents.halved + rnorm(2*M, sd = 0.02))))
  offspring.geno.scaled = sweep(offspring.geno.unscaled, 1, 2*MAF, FUN = "-" ) 
  offspring.geno.scaled = sweep(offspring.geno.scaled, 1, sqrt(2*MAF*(1 - MAF)), FUN = "/")
  offspring.geno.lia.first = lia.beta %*% offspring.geno.scaled
  offspring.lia.first = offspring.geno.lia.first  + environ[-(2:3)]
  return(list(
    "offspring.cc" = as.numeric( offspring.lia.first > lia.T ) ,
    "offspring.geno.unscaled" = offspring.geno.unscaled[,1], # I do not need the scaled offsprings genotype after this point.
    "offspring.geno.lia" = offspring.geno.lia.first,
    "offspring.lia" = offspring.lia.first
  ))
}

generateChildren = function(
  NoChildren = 10000,
  M = 10000,
  MAF,
  C = M*0.01,
  lia.h2 = 0.5,
  K = 0.05,
  nsib = 2,
  nthreads = 6,
  allele.mat = matrix(c("A","A","A","G", "G","G"), nrow = 3, ncol = 2, byrow = T),
  fixed.effect = T,
  out = "C:\\Users\\FIUN7293\\CommandStationBeta\\EphemeralFolder\\Results\\ph"
) {
  #liability threshold:
  lia.T = qnorm(1 - K) 
  #liability scale effect sizes:
  if ( fixed.effect == F) {
    nonzero.betas = rnorm(n = C, mean = 0, sd = sqrt(lia.h2/C))
  } else {
    cat("OBS: currently using FIXED effect sizes! \n")
    nonzero.betas = sqrt(lia.h2/C)
  }
  
  #vector of liability scale betas:
  lia.beta = rep(0, length.out = M)
  #positions of nonzero betas:
  beta.pos = sample(1:M, size = C, replace = FALSE)
  
  #assigning nonzero betas to their positions:
  lia.beta[beta.pos] = nonzero.betas
  
  #covariance matrix for environmen:  
  covarMat = diag(1 - lia.h2, 3 + nsib) #this can be extended by eq (4) from the article to include environmental correlation.
  
  #mean value vector for the environment variable:
  muVec = rep(0,3 + nsib)
  
  ##### doing first indiv outside parallel loop to start the output file:  
  
  #generates a matrix containing genotype [0/1/2] for each potential parent:
  parents =  replicate(n = 2, rbinom(n = M, size = 2, prob = MAF))
  #normalising genotypes:
  parents.scaled = sweep(parents, 1, 2*MAF, FUN = "-") ###### THIS!!
  parents.scaled = sweep(parents.scaled, 1, sqrt(2*MAF*(1 - MAF)), FUN = "/")###### THIS!!
  #environment for parents and offspring:
  environ = mvrnorm(n = 1, mu = muVec, Sigma = covarMat) 
  
  #assign parental liabilities:
  parents.geno.lia.first = lia.beta %*%  parents.scaled
  parents.lia = parents.geno.lia.first + environ[2:3] 
  #assign parental phenotypes:
  parents.cc = as.numeric(parents.lia > lia.T)
  ##generating child geno and phenotype:
  #Other Method:
  parents.halved = parents/2 #same parents for all siblings, no need to recalculate this
  
  offspring = generateOffspring(parents.halved, nsib, M, MAF, environ, lia.T, lia.beta)
  
  
 # lock = tempfile()
  
  fwrite(as.data.table(matrix(c(1,1,0,0,0, #FID, IID, Father, Mother, Sex
                                "offspring_cc" = offspring$offspring.cc[1],
                                "offspring_geno" = t(allele.mat[offspring$offspring.geno.unscaled + 1,]) ), nrow = 1)),
         file = paste(out,".ped", sep = ""), sep = " ",
         row.names = F, col.names = F, quote = F)
  #here "offspring.geno" is technically a 2xM matrix, however when it is being written to the file it is collapsed
  #into a vetor where each column comes one after the other.
  
  
 # cat("starting parallelization backend with", nthreads, "threads for generation of children:\n")
#  cl = makeCluster(nthreads, type = "SOCK")
#  registerDoSNOW(cl)
 # iterations = NoChildren
  
#  pb = progress_bar$new(
#    format = "[:bar] :percent",
#    total = iterations,
#    width = 100)
  
  #  progress_num = 1:iterations
  #  progress = function(n){
  #    pb$tick(tokens = list(letter = progress_num[n]))
  #}
  
  #opts = list(progress = progress)
  
  ph = foreach(i = 2:NoChildren, .packages = c("MASS", "stringr", "flock", "data.table"), .export = c("generateOffspring")) %do% { #.options.snow = opts,
    #generates a matrix containing genotype [0/1/2] for each potential parent:
    parents =  replicate(n = 2, rbinom(n = M, size = 2, prob = MAF))
    #normalising genotypes:
    parents.scaled = sweep(parents, 1, 2*MAF, FUN = "-") ###### THIS!!
    parents.scaled = sweep(parents.scaled, 1, sqrt(2*MAF*(1 - MAF)), FUN = "/")###### THIS!!
    #environment for parents and offspring:
    environ = mvrnorm(n = 1, mu = muVec, Sigma = covarMat) 
    
    #assign parental liabilities:
    parents.geno.lia = lia.beta %*%  parents.scaled
    parents.lia = parents.geno.lia + environ[2:3] 
    #assign parental phenotypes:
    parents.cc = as.numeric(parents.lia > lia.T)
    ##generating child geno and phenotype:
    #Other Method:
    parents.halved = parents/2 #same parents for all siblings, no need to recalculate this
    
    offspring = generateOffspring(parents.halved, nsib, M, MAF, environ, lia.T, lia.beta)
    mat.ph = as.data.table(matrix(c(i,i,0,0,0, #FID, IID, Father, Mother, Sex
                                    "offspring_cc" = offspring$offspring.cc[1],
                                    "offspring_geno" = t(allele.mat[offspring$offspring.geno.unscaled + 1,]) ), nrow = 1))
    
    # locked = flock::lock(lock)
    fwrite(mat.ph,
           file = paste(out,".ped", sep = ""), sep = " ",
           append = T)
    #   flock::unlock(locked)
    
    c("offspring_cc" = offspring$offspring.cc[1], "parents_cc" = parents.cc, "offspring_lia" = offspring$offspring.lia[1], "parents_lia" = parents.lia, "offspring_geno_lia" =  offspring$offspring.geno.lia[1], "parents_geno_lia" = parents.geno.lia, "siblings_cc" = offspring$offspring.cc[-1], "siblings_lia" = offspring$offspring.lia[-1], "siblings_geno_lia" = offspring$offspring.geno.lia[-1] )
  }
  # close(pb)
  # stopCluster(cl)
  ph = matrix(unlist(ph), ncol = 9 + 3*nsib, byrow = T)
  return(list("offspring_cc" = c(offspring$offspring.cc[1], ph[,1]), 
              "parents_cc" = rbind(parents.cc,ph[,2:3]), 
              "offspring_lia" = c(offspring$offspring.lia[1], ph[,4]),
              "parents_lia" = rbind(parents.lia, ph[,5:6]),
              "offspring_geno_lia" = c(offspring$offspring.geno.lia[1], ph[,7]),
              "parents_geno_lia" = rbind(parents.geno.lia.first, ph[,8:9]),
              "siblings_cc" = rbind(offspring$offspring.cc[-1], ph[,9 + 1:nsib]),
              "siblings_lia" = rbind(offspring$offspring.lia[-1], ph[,9 + nsib + 1:nsib]),
              "siblings_geno_lia" = rbind(offspring$offspring.geno.lia[-1], ph[, 9 + 2*nsib + 1:nsib]),
              "lia_betas" = lia.beta,
              "lia_threshold" = lia.T))
}





#No. snps:
M = as.numeric(args$'M')
#No. Children:
NoChildren = as.numeric(args$'noChildren')
#simulated maf values:
MAF = runif(n = M, min = 0.01, max = 0.5)
#No. causal snps:
C = as.numeric(args$'C')
#liability scale heritability:
lia.h2 = as.numeric(args$'lia-h2')
#prevalence in parents:
K = as.numeric(args$'K')
#out = "C:\\Users\\FIUN7293\\CommandStationBeta\\EphemeralFolder\\Results\\sim100kx100k_v10"
nsib = as.numeric(args$'nsib')

#for (out in paste("C:\\Users\\FIUN7293\\CommandStationBeta\\EphemeralFolder\\Results\\sibs2_10k2_F_C",C,"_v", 1:10,"_NDS", sep = "")) {

if (as.logical(args$'fixedEffect')) {
  NFF = "F"
} else {
  NFF = "NF"
}
if (is.null(args$'C')) {
  C = M*0.02
}


base = args$'baseName'
size = args$'sizeName'
identifier = args$'identifier'
v = args$'v'
nthreads = args$'nthreads' 

out = paste(args$'out', "/", base, "_", size,"_", NFF , "_C", C, "_v",v,identifier, sep = "")
cat("\n:-================================================================================-:\n")
cat("\nworking on:\n", out, "\n")

children = generateChildren(MAF = MAF,
                            lia.h2 = lia.h2, 
                            M = M, 
                            C = C, 
                            nthreads = nthreads, 
                            NoChildren = NoChildren,
                            K = K,
                            out = out,
                            fixed.effect = as.logical(args$'fixedEffect'),
                            nsib = nsib)


#save input file for LTFH, i.e., case/ctrl status for offspring and parents. sibs not included here:
fwrite(list("FID" = 1:NoChildren,
            "IID" = 1:NoChildren,
            "CHILD_STATUS" = children$offspring_cc,
            "P1_STATUS" = children$parents_cc[,1],
            "P2_STATUS" = children$parents_cc[,2],
            "NUM_SIBS" = rep(nsib, NoChildren),
            "SIB_STATUS" = ifelse(rowSums(children$siblings_cc) > 0, 1, 0)),
       file = paste(out,".phen", sep = ""), sep = " ")




#save df.map:
fwrite(list("CHR" = rep(1,M),
            "SNP" = paste("rs", 1:M, sep = ""),
            "cM" = rep(0, M),
            "BP" = 1:M), 
       paste(out,".map", sep = ""), sep = " ", col.names = F)


#df.true = data.frame("lia.betas" = children$lia.betas,
#                     "lia.threshold" = children$lia.threshold,
#                     "prevalens" = K,
#                     "No. Causal SNPs" = C)
fwrite(data.frame("FID" = 1:NoChildren,
                  "IID" = 1:NoChildren,
                  "lia_betas" = children$lia_betas,
                  "offspring_lia" = children$offspring_lia,
                  "offspringgeno_lia" = children$offspring_geno_lia,
                  "parents_lia" = children$parents_lia,
                  "parents_geno_lia" = children$parents_geno_lia,
                  "lia_threshold" = children$lia_threshold,
                  "prevalens" = K,
                  "No_Causal_SNPs" = C,
                  "maf" = MAF,
                  "sex" = as.numeric(runif(NoChildren) >= .5),
                  "sibling_cc" = children$siblings_cc,
                  "sibling_lia" = children$siblings_lia, 
                  "sibling_geno_lia" = children$siblings_geno_lia,
                  "sibling_sex" = replicate(nsib, as.numeric(runif(NoChildren) >= .5))),
       paste(out,".true", sep = ""), sep = " ", col.names = T)




