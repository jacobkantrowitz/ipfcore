# basic functions for general use

writeDelim = function(dat,filename,row.names=F) {
 write.table(dat,filename,sep="\t",quote=F,na="",row.names=row.names)
}

readDelim = function(filename,stringsAsFactors=F) {
 read.delim(filename,stringsAsFactors=stringsAsFactors)
}

