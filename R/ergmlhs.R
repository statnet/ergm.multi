.combine_ergmlhs <- function(nwl, ignore.settings=c()){
  ergml <- nwl %>% map(`%n%`, "ergm") %>% map(NVL, list())
  Reduce(function(l1,l2){
    settings <- union(names(l1),names(l2))
    for(setting in settings){
      if(!identical(l1[[setting]],l2[[setting]])) stop("Two subnetworks have inconsistent ERGM settings.")
      l1
    }
  }, ergml)
}
