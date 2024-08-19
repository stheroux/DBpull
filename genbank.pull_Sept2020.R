library(rentrez)
library(reutils)
options(error=utils::recover)
api_key="b2ccdcd6619a29d3bf31de74e7cde9a1c209"  

fetch.gb<-function(a,b,c) {
  
  for(i in a)
  {
    dir.create("fasta")
    query<-paste('"',i,'"', "AND", b, "AND", c)
    query.num<-entrez_search(db="nuccore", query)
    count<-query.num$count
    print(paste(i, "count is", count))
    
    #  
    # if count < 1   
    #
    if (count < 1) {
      print("count equals zero") 
      no.ref.seq.taxa<-i
      write(no.ref.seq.taxa, file=paste(b,".no.ref.seq.taxa.txt"), append=T) }
    
    else
      
      #  
      # if count =< 50  
      # 
      if ((count > 0) & (count <= 50)) {
        print("count is less than or equal to 50")
        
        search.out <- entrez_search(db="nuccore", query, use_history=TRUE, retmax=count)
        fetch.out <- entrez_fetch(db="nuccore", web_history=search.out$web_history,
                                  rettype="fasta", retmax=count)
        
        write(fetch.out, file=paste("fasta/",i,"fetch.out.fasta"), append=TRUE) }
    
    #
    # if count > 50 and divisable by 50   
    # 
    if ((count > 50) & (count%%50 == 0) & (count <= 2000)) {
      print("count is greater than 50 and divis by 50")
      
      for(seq_start in seq(1,count,50)) {
        search.out <- entrez_search(db="nuccore", query, use_history=TRUE, retmax=count)
        fetch.out <- entrez_fetch(db="nuccore", web_history=search.out$web_history,
                                  rettype="fasta", retmax=50, retstart=seq_start)
        write(fetch.out, file=paste("fasta/",i,"fetch.out.fasta"), append=TRUE) } }
    
    else
      
      #
      # if count > 50 and NOT divisable by 50   
      # 
      if ((count > 50) & (count%%50 != 0) & (count < 2000))  {
        
        print("count is greater than 50 and NOT divis by 50") 
        
        remainder<-(count%%50)
        
        search.out <- entrez_search(db="nuccore", query, use_history=TRUE, retmax=count)
        fetch.out <- entrez_fetch(db="nuccore", web_history=search.out$web_history,
                                  rettype="fasta", retmax=remainder, retstart=1)
        write(fetch.out, file=(paste("fasta/",i,"fetch.out.fasta")))
        
        restart<-(remainder+1)
        
        for(seq_start in seq(restart,count,50)) {
          search.out <- entrez_search(db="nuccore", query, use_history=TRUE, retmax=count)
          fetch.out <- entrez_fetch(db="nuccore", web_history=search.out$web_history,
                                    rettype="fasta", retmax=50, retstart=seq_start)
          write(fetch.out, file=(paste("fasta/",i,"fetch.out.fasta")), append=T) }
        
      }         
    
    else 
      
      if (count > 1000)  {
        
        print("count is greater than 1000, skip")
        omit.taxa<-i
        write(omit.taxa, file=paste(b,".omit.taxa.txt"), append=T)
      }
    
  }
  
  # optional quick summary of output
  #zed<-read.table("no.ref.seq.taxa.txt", header=F, sep='\t')
  #zed<-t(zed)
  #perc.no.ref.seq<-(length(zed))/(length(TAXA))*100
  #print(paste("percent no ref seq is", perc.no.ref.seq))
  #write(perc.no.ref.seq, file="perc.no.ref.seq.txt")
  
  
  
  }
