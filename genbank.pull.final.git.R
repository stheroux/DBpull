#####################################################################
#                                                                   #
# Purpose:       Download genbank sequences in fasta format         #
#                                                                   #
# Author:        Susie Theorux                                      #
# Contact:       susannat@sccwrp.org                                #
# Client:        SCCWRP                                             #
#                                                                   #
# Code created:  2016-6-01                                          #
# Last updated:  2016-06-05                                         #
# Source:        github.com/stheroux/genbank.pull                   #
#                                                                   #
# Requires:       set working directory (setwd)                     #
#                 "taxa.txt" tab-delim taxa name list in directory  #
#                                                                   #
# Usage:  Edit setwd                                                #
#         Source file                                               #
#         fetch.gb(TAXA, "18S[ALL]")                                #
#         fetch.gb(TAXA, "rbcL[ALL]")                               #
#                                                                   #
# ###################################################################


setwd("~/")

#install.packages("rentrez")
#install.packages("reutils")

library(rentrez)
library(reutils)
options(error=utils::recover)

TAXA<-read.table("taxa.txt", sep='\t')
TAXA<-t(TAXA)  

#######################
###  genbank.pull #####
#######################


#where a is the taxa list 
#where b is the gene "18S[ALL]" or "rbcL[ALL]"

fetch.gb<-function(a,b) {

#TAXA<-read.table(a, sep='\t')
#TAXA<-t(TAXA)  
  
for(i in a)
  {
  query<-paste(i, "AND", b, "AND 200:1800[SLEN]")
  query.num<-entrez_search(db="nuccore", query)
  count<-query.num$count
  print(paste(i, "count is", count))

#  
# if count < 1   
#
    if (count < 1) {
    print("count equals zero") 
    no.ref.seq.taxa<-i
    write(no.ref.seq.taxa, file="no.ref.seq.taxa.txt", append=T) }

else

#  
# if count =< 50  
# 
   if ((count > 0) & (count <= 50)) {
    print("count is less than or equal to 50")
    
    search.out <- entrez_search(db="nuccore", query, use_history=TRUE, retmax=count)
    fetch.out <- entrez_fetch(db="nuccore", web_history=search.out$web_history,
                              rettype="fasta", retmax=count)
  
    write(fetch.out, file=paste(query,"fetch.out.fasta"), append=TRUE) }

#
# if count > 50 and divisable by 50   
# 
   if ((count > 50) & (count%%50 == 0) & (count <= 2000)) {
    print("count is greater than 50 and divis by 50")
   
    for(seq_start in seq(1,count,50)) {
    search.out <- entrez_search(db="nuccore", query, use_history=TRUE, retmax=count)
    fetch.out <- entrez_fetch(db="nuccore", web_history=search.out$web_history,
                                  rettype="fasta", retmax=50, retstart=seq_start)
    write(fetch.out, file=paste(query,"fetch.out.fasta"), append=TRUE) } }

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
    write(fetch.out, file=(paste(query,"fetch.out.fasta")))
    
    restart<-(remainder+1)
  
    for(seq_start in seq(restart,count,50)) {
      search.out <- entrez_search(db="nuccore", query, use_history=TRUE, retmax=count)
      fetch.out <- entrez_fetch(db="nuccore", web_history=search.out$web_history,
                                rettype="fasta", retmax=50, retstart=seq_start)
      write(fetch.out, file=(paste(query,"fetch.out.fasta")), append=T) }
    
  }         

else 
  
  if (count > 2000)  {
    
    print("count is greater than 2000, skip")
    omit.taxa<-i
    write(omit.taxa, file="omit.taxa.txt", append=T)
    }
  
}}
    
      
  
