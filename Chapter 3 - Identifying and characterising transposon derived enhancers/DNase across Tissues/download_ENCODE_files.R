

try(system('xargs -n 1 curl -O -L < files.txt'))

bed = read.delim('bed_metadata.tsv',as.is=T)
fname = paste(bed$File.accession,'bed.gz',sep='.')
new.name = paste(gsub(' ','_',bed$Biosample.term.name),fname,sep='_')

file.rename(fname,new.name)
