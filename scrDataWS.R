

##### This is a function to create an SCR object of class scrobj. ##### 
##### This is necessary for correct formatting of data            #####
##### In order to aid users, there are some error checks          #####
##### incorporated here. Users are encouraged add more such error #####
##### trapping routines                                           #####
############################### MODIFIED 8 MAY 2019 ############################
scrData <-
function(traps,captures,statespace,alive=NULL,
         Xsex=NULL,
         Xd=NULL,Xeff=NULL, Xeff1=NULL, Xeff2=NULL, Xeff3=NULL, Ytel=NULL, Xtel=NULL){
			 
################################################################################

if(ncol(captures)==3){
session=rep(1,nrow(captures))
individual <- captures[,"ANIMAL_ID"]
occasion<-captures[,"SO"]
trapid<- captures[,"LOC_ID"] 
}else{
 session<- captures[,"session"]
individual<-captures[,"ANIMAL_ID"]
occasion<-captures[,"SO"]
trapid<- captures[,"LOC_ID"]}

captures<-cbind(session=session,individual=individual,occasion=occasion,trapid=trapid)

traplocs<-traps[,2:3]
MASK<-as.matrix(traps[,4:ncol(traps)])
nind<-max(individual)
K<-dim(MASK)[2]
ntraps<-nrow(traplocs)

alive=matrix(1,nrow=length(unique(captures[,"individual"])),ncol=ncol(MASK))
Xd<- rep(1,nrow(statespace))



checkdata.fn<-function(Y){
encfreqs<- table(table(Y[,2]))
totrecaps<- sum(encfreqs*(as.numeric(names(encfreqs)) -1))
trapfreqs<- table(apply(table(Y[,2],Y[,1])>0,1,sum))
multitraps<- sum(trapfreqs*(as.numeric(names(trapfreqs))-1))
c(totrecaps,multitraps)
}
			 
## Some minimum data requirements are used here. This is completely arbitrary but
## it is meant to FORCE the user to proceed AT THEIR OWN RISK knowing
## that they have limited data and thinking about the potential consequences
			 
if(1==2){
smydata<-checkdata.fn(Y)
cat("You had ",smydata[1]," recaptures in your study",fill=TRUE)
cat("You had ",smydata[2]," multi-trap recaptures in your study",fill=TRUE)
if( (smydata[1]< 10) | (smydata[2] < 5) ){
 cat("You should go into the field and obtain more data.",fill=TRUE)
 cat("Try again later.",fill=TRUE)
 return(NULL)
}
}


check<-cbind(captures[,"trapid"],captures[,"occasion"])
check2<- MASK[check]
if(any(check2 != 1)){
    cat("Some trap operation data != 1. Check the following record numbers on capture history file:",fill=TRUE)
    cat(which(check2 != 1))
    return(NULL)
}



if(!is.null(Xsex)){

if(!is.numeric(Xsex)){
cat("Error, sex needs to be binary",fill=TRUE)
return(NULL)
}


}


############################## MODIFIED 8 MAY 2019 ###########################
obj<-list(
    traps=traps,
    captures=captures,
    statespace=statespace,
    alive=alive,
    Xd=Xd,
    Xsex=Xsex,
    Xeff=Xeffort, Xeff1=Xeffort1, Xeff2=Xeffort2, Xeff3=Xeffort3)
##############################################################################

class(obj)<-c("scrdata","list")

return(obj)

}
