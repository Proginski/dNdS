# home <- "/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski"
# setwd(paste(home,"/Eukaryotes",sep=""))

# In order to pass arguments with bash
args <- commandArgs(trailingOnly = TRUE)

tree_path= args[1]
out_path = args[2]

# tree_path <- file("PHYLO/SCER_NCBI_NEXT/test.nwk","r")
tree_text <- readLines(tree_path,n=1)
print(paste("Orignal tree :",tree_text))


strip.nodelabels<-function(text){
  obj<-strsplit(text,"")[[1]]
  cp<-grep(")",obj)
  csc<-c(grep(":",obj),length(obj))
  exc<-cbind(cp,sapply(cp,function(x,y) y[which(y>x)[1]],y=csc))
  exc<-exc[(exc[,2]-exc[,1])>1,]
  inc<-rep(TRUE,length(obj))
  if(nrow(exc)>0) for(i in 1:nrow(exc)) 
    inc[(exc[i,1]+1):(exc[i,2]-1)]<-FALSE
  paste(obj[inc],collapse="")
}

tree_text <- strip.nodelabels(tree_text)
print(paste("New tree :",tree_text))

out <-file(out_path)
writeLines(tree_text, out)
close(out)

