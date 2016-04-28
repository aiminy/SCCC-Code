#to use
#R --slave --args input.txt pair_name output.txt < ~code/code_R/GBS_QC.R

myarg <- commandArgs()
cat(myarg,"\n");
m=length(myarg)
cat(m,"\n");

f_in<-myarg[4:4]
p_name<-myarg[5:5]
f_out<-myarg[6:6]

data.gbs<-read.csv(f_in,sep="\t",header=F)

print(dim(data.gbs)[1])


m=dim(data.gbs)[1];
#n=dim(data.gbs)[2];
#nn=n-1;


cat("Pairs","\t","Chr","\t","Pos","\t","REF","\t","ALT","\t","normal_REF","\t","normal_ALT","\t","tumor_REF","\t","tumor_ALT","\n",file=f_out,append=FALSE)


for(i in 1:m){

#REF
if(data.gbs[i,3]=="A"){
tumor=unlist(strsplit(as.character(data.gbs[i,6]),split=","))[2]
normal=unlist(strsplit(as.character(data.gbs[i,5]),split=","))[2]
tumor_REF=unlist(strsplit(as.character(tumor),split=":"))[1]
normal_REF=unlist(strsplit(as.character(normal),split=":"))[1]
#cat(data.gbs[i,1],"\t",data.gbs[i,2],"\t",tumor1,"\t",normal1,"\n")
}

if(data.gbs[i,3]=="C"){
tumor=unlist(strsplit(as.character(data.gbs[i,6]),split=","))[3]
normal=unlist(strsplit(as.character(data.gbs[i,5]),split=","))[3]
tumor_REF=unlist(strsplit(as.character(tumor),split=":"))[1]
normal_REF=unlist(strsplit(as.character(normal),split=":"))[1]
#cat(data.gbs[i,1],"\t",data.gbs[i,2],"\t",tumor1,"\t",normal1,"\n")
#cat(tumor1,"\t",normal1,"\n")
}

if(data.gbs[i,3]=="G"){
tumor=unlist(strsplit(as.character(data.gbs[i,6]),split=","))[4]
normal=unlist(strsplit(as.character(data.gbs[i,5]),split=","))[4]
tumor_REF=unlist(strsplit(as.character(tumor),split=":"))[1]
normal_REF=unlist(strsplit(as.character(normal),split=":"))[1]
#cat(data.gbs[i,1],"\t",data.gbs[i,2],"\t",tumor1,"\t",normal1,"\n")
#cat(tumor1,"\t",normal1,"\n")
}

if(data.gbs[i,3]=="T"){
tumor=unlist(strsplit(as.character(data.gbs[i,6]),split=","))[5]
normal=unlist(strsplit(as.character(data.gbs[i,5]),split=","))[5]
tumor_REF=unlist(strsplit(as.character(tumor),split=":"))[1]
normal_REF=unlist(strsplit(as.character(normal),split=":"))[1]
#cat(data.gbs[i,1],"\t",data.gbs[i,2],"\t",tumor1,"\t",normal1,"\n")
#cat(tumor1,"\t",normal1,"\n")
}

#ALT
if(data.gbs[i,4]=="A"){
tumor=unlist(strsplit(as.character(data.gbs[i,6]),split=","))[2]
normal=unlist(strsplit(as.character(data.gbs[i,5]),split=","))[2]
tumor_ALT=unlist(strsplit(as.character(tumor),split=":"))[1]
normal_ALT=unlist(strsplit(as.character(normal),split=":"))[1]
#cat(data.gbs[i,1],"\t",data.gbs[i,2],"\t",tumor1,"\t",normal1,"\n")
}


if(data.gbs[i,4]=="C"){
tumor=unlist(strsplit(as.character(data.gbs[i,6]),split=","))[3]
normal=unlist(strsplit(as.character(data.gbs[i,5]),split=","))[3]
tumor_ALT=unlist(strsplit(as.character(tumor),split=":"))[1]
normal_ALT=unlist(strsplit(as.character(normal),split=":"))[1]
}

if(data.gbs[i,4]=="G"){
tumor=unlist(strsplit(as.character(data.gbs[i,6]),split=","))[4]
normal=unlist(strsplit(as.character(data.gbs[i,5]),split=","))[4]
tumor_ALT=unlist(strsplit(as.character(tumor),split=":"))[1]
normal_ALT=unlist(strsplit(as.character(normal),split=":"))[1]

}

if(data.gbs[i,4]=="T"){
tumor=unlist(strsplit(as.character(data.gbs[i,6]),split=","))[5]
normal=unlist(strsplit(as.character(data.gbs[i,5]),split=","))[5]
tumor_ALT=unlist(strsplit(as.character(tumor),split=":"))[1]
normal_ALT=unlist(strsplit(as.character(normal),split=":"))[1]
}


cat(p_name,"\t",as.character(data.gbs[i,1]),"\t",data.gbs[i,2],"\t",as.character(data.gbs[i,3]),"\t",as.character(data.gbs[i,4]),"\t",normal_REF,"\t",normal_ALT,"\t",tumor_REF,"\t",tumor_ALT,"\n",file=f_out,append=TRUE)

}

quit("yes")