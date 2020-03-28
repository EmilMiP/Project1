library(data.table)
library(stringr)

#get list of all the summary stats to use, and their true values.
fileRoot = "/home/emp/faststorage/project1/simulatedData"
#fileRoot = "D:\\Work\\Project1\\SimlatedData"
base = "sibs2"
sword = "10kx10k"
cword = "C200"
bword = "_NF"
identifier = ""

true.files = list.files(path = fileRoot, pattern = "\\.true", full.names = T)
true.files = str_subset(true.files, pattern = base)
true.files = str_subset(true.files, pattern = sword)
true.files = str_subset(true.files, pattern = cword)
true.files = str_subset(true.files, pattern = bword)
true.files = str_subset(true.files, pattern = identifier)


keeplist.files = gsub("true", "keeplist", true.files)


phen.files = list.files(path = fileRoot, pattern = "\\.phen", full.names = T)
phen.files = str_subset(phen.files, pattern = base)
phen.files = str_subset(phen.files, pattern = sword)
phen.files = str_subset(phen.files, pattern = cword)
phen.files = str_subset(phen.files, pattern = bword)
phen.files = str_subset(phen.files, pattern = identifier)

#phen.files = str_subset(phen.files, pattern = "/sibs")


ExcludeListGen = function(keeplist.files, K_child) {
  for (i in seq_along(keeplist.files)) {
    phen = fread(phen.files[i])
    FID = 1:dim(phen)[1]
    phen.vec = str_subset(colnames(phen), "CHILD_STATUS")
    for (j in seq_along(phen.vec)) {
      ctrls = sort(sample(FID[phen[[phen.vec[j]]] == 0], size = dim(phen)[1] - sum(phen[[phen.vec[j]]] == 1) - 5000))
      removers = data.frame("FID" = ctrls, "IID" = ctrls)
      fwrite(removers, paste(gsub("keeplist", "excludeList", keeplist.files[i]), sep = ""), sep = " ", quote = F)
    }
  }
}