library("RGeostats")
library("maptools")
newfullpoly=readRDS("newfullpoly.rds")

areas=read.csv("gridareas.csv")

fishdata2014=read.csv("fishdata2014.csv")[,2:4]
fishdata2015=read.csv("fishdata2015.csv")[,4:6]
fishdata2016=read.csv("fishdata2016.csv")[,2:4]

colnames(fishdata2014)=c("lon","lat","NASC")
colnames(fishdata2015)=c("lon","lat","NASC")
colnames(fishdata2016)=c("lon","lat","NASC")


############## sig_bs proportion simulation ############

sigpropdata2016=read.table("2016spratsigprop.csv",header=T,sep=",")
sigpropdata2015=read.table("2015spratsigprop.csv",header=T,sep=",")
sigpropdata2014=read.table("2014spratsigprop.csv",header=T,sep=",")


#2016
projec.toggle(0)
sigpropdata2016.db=db.create(sigpropdata2016,flag.grid=F,ndim=2,autoname=F)
db.data1=db.polygon(sigpropdata2016.db,newfullpoly)

projec.define(projection="mean", db=db.data1,flag.update=T) 
model.anam1=anam.fit(db.data1,type="gaus",draw=T,nbpoly=8,title="Anamophosis",xlim=c(-3,3))
db.data.trans1=anam.z2y(db.data1,anam=model.anam1)
vario.data1=vario.calc(db.data.trans1,lag=8,nlag=5)
model.vario1 <- model.auto(vario.data1,struc=melem.name(c(1,2)),draw=T)
grid.simu1=db.create(flag.grid=T,x0=c(-5.9,54),dx=c(0.015839,0.008971),nx=c(100,600))
grid.simu1 <- db.polygon(grid.simu1,newfullpoly);beep(2)
neigh.simu <- neigh.create(ndim=2,type=0)
grid.simu1 <- simtub(dbin=db.data.trans1, dbout=grid.simu1, model=model.vario1, 
                     neigh=neigh.simu, uc = "", mean = 0, nbsimu = 100, 
                     nbtuba = 1000, radix = "Simu",modify.target = TRUE); beep(2)
grid.simu1 <- anam.y2z(grid.simu1,name="Simu*",anam=model.anam1)
names=vector(,100)
for(i in 1:100)
{names[i]=paste("Raw.Simu.Gaussian.sigprop.S",i,sep="")}
x2016=matrix(,3710,100)
for(i in 1:100){
  x2016[,i]=db.extract(grid.simu1,names=names[i],flag.compress=T)
}
hist(x2016)
x2016[x2016>1]=1
x2016[x2016<0]=0
saveRDS(grid.simu1,"2016_sigprop_anam_sims_db.rds")
saveRDS(x2016,"2016sigprop_anam_sims.rds")


#2015

projec.toggle(0)
sigpropdata2015.db=db.create(sigpropdata2015,flag.grid=F,ndim=2,autoname=F)
db.data1=db.polygon(sigpropdata2015.db,newfullpoly)
projec.define(projection="mean", db=db.data1,flag.update=T) 
model.anam1=anam.fit(db.data1,type="gaus",draw=T,nbpoly=14,title="Anamophosis",xlim=c(-5,5))
db.data.trans1=anam.z2y(db.data1,anam=model.anam1)
vario.data1=vario.calc(db.data.trans1,lag=8,nlag=5)
model.vario1 <- model.auto(vario.data1,struc=melem.name(c(1,2)),draw=T)
grid.simu1=db.create(flag.grid=T,x0=c(-5.9,54),dx=c(0.015839,0.008971),nx=c(100,600))
grid.simu1 <- db.polygon(grid.simu1,newfullpoly)
neigh.simu <- neigh.create(ndim=2,type=0)
grid.simu1 <- simtub(dbin=db.data.trans1, dbout=grid.simu1, model=model.vario1, 
                     neigh=neigh.simu, uc = "", mean = 0, nbsimu = 100, 
                     nbtuba = 1000, radix = "Simu",modify.target = TRUE); beep(2)
grid.simu1 <- anam.y2z(grid.simu1,name="Simu*",anam=model.anam1)
names=vector(,100)
for(i in 1:100)
{names[i]=paste("Raw.Simu.Gaussian.sigprop.S",i,sep="")}
x2015=matrix(,3710,100)
for(i in 1:100){
  x2015[,i]=db.extract(grid.simu1,names=names[i],flag.compress=T)
}
hist(x2015)
x2015[x2015>1]=1
x2015[x2015<0]=0

saveRDS(grid.simu1,"2015_sigprop_anam_simugrid.rds")
saveRDS(x2014,"2015sigprop_anam_sims.rds")


#2014

projec.toggle(0)
sigpropdata2014.db=db.create(sigpropdata2014,flag.grid=F,ndim=2,autoname=F)
db.data1=db.polygon(sigpropdata2014.db,newfullpoly)
projec.define(projection="mean", db=db.data1,flag.update=T) 
model.anam1=anam.fit(db.data1,type="gaus",draw=T,nbpoly=14,title="Anamophosis",xlim=c(-5,5))
db.data.trans1=anam.z2y(db.data1,anam=model.anam1)
vario.data1=vario.calc(db.data.trans1,lag=8,nlag=5)
model.vario1 <- model.auto(vario.data1,struc=melem.name(c(1)),draw=T)
grid.simu1=db.create(flag.grid=T,x0=c(-5.9,54),dx=c(0.015839,0.008971),nx=c(100,600))
grid.simu1 <- db.polygon(grid.simu1,newfullpoly);beep(2)
neigh.simu <- neigh.create(ndim=2,type=0)
for(i in 1:100){
grid.simu1 <- simtub(dbin=db.data.trans1, dbout=grid.simu1, model=model.vario1, 
                     neigh=neigh.simu, uc = "", mean = 0, nbsimu = 1,seed=sample(1:100000,1,replace=F), 
                     nbtuba = 1000, radix = paste("Simu_",i,sep=""),modify.target = TRUE)
}
grid.simu1 <- anam.y2z(grid.simu1,name="Simu*",anam=model.anam1)
names=vector(,100)
for(i in 1:100)
{names[i]=paste("Raw.Simu_",i,".Gaussian.sigprop.S1",sep="")}
x2014=matrix(,3710,100)
for(i in 1:100){
  x2014[,i]=db.extract(grid.simu1,names=names[i],flag.compress=T)
}
hist(x2014)
x2014[x2014>1]=1
x2014[x2014<0]=0

saveRDS(grid.simu1,"2014_sigprop_anam_simugrid.rds")
saveRDS(x2014,"2014sigprop_anam_sims.rds")




############## 2016 fish NASC sim #####################

data=fishdata2016
data.db1=db.create(data,flag.grid=F,ndim=2,autoname=F)
projec.toggle(0)
db.data1=db.polygon(data.db1,newfullpoly)
projec.define(projection="mean", db=db.data1,flag.update=T) 
model.anam1=anam.fit(db.data1,type="emp",ndisc=db.data1$nech,draw=T,title="Anamophosis")
db.data.trans1=anam.z2y(db.data1,anam=model.anam1)
db.data.trans1 <- db.rename(db.data.trans1,name="Gaussian.NASC",newname="Yp")

ycut1 <- qnorm(sum(db.extract(db.data.trans1,"NASC") == 0) / length(db.extract(db.data.trans1,name="NASC",flag.compress=T)))
Y1<- db.extract(db.data.trans1,"Yp",flag.compress=T)
NASC1 <- db.extract(db.data.trans1,"NASC")
Y1[NASC1 == 0] <- ycut1
db.data.trans1 <- db.replace(db.data.trans1,"Yp",Y1)
n.H <- 50
vario.Yp1 <- vario.calc(db.data.trans1,lag=1,nlag=40)
vario.Y1  <- vario.trans.cut(vario.Yp1,ycut1,n.H)
model.vario.Y1 <- model.auto(vario.Y1,struc=melem.name(c(1,2)),draw=T)
Ymax1 <- db.extract(db.data.trans1,name="Yp",flag.compress=F)
Ymin1 <- db.extract(db.data.trans1,name="Yp",flag.compress=F)
Ymin1[Ymin1 <= ycut1] <- -10
db.data.trans1<-db.add(db.data.trans1,Ymax1)
db.data.trans1<-db.locate(db.data.trans1,db.data.trans1$natt,"upper")
db.data.trans1<-db.add(db.data.trans1,Ymin1)
db.data.trans1<-db.locate(db.data.trans1,db.data.trans1$natt,"lower")

db.data.trans1 <-gibbs(db = db.data.trans1, model = model.vario.Y1, seed = 2322, 
                       nboot = 10, niter = 500, flag.norm=FALSE, toleps = 1,
                       radix = "Gibbs", modify.target = TRUE)

db.data.trans1<-db.rename(db.data.trans1,"Gibbs.G1","Y2")

grid.simu2=db.create(flag.grid=T,x0=c(-5.9,54),dx=c(0.015839,0.008971),nx=c(100,600))
grid.simu2 <- db.polygon(grid.simu2,newfullpoly)
neigh.simu <- neigh.create(ndim=2,type=2,nmini=5,nmaxi=100)
grid.simu2 <- simtub(dbin=db.data.trans1, dbout=grid.simu2, model=model.vario.Y1, 
                     neigh=neigh.simu, uc = "", mean = 0, nbsimu = 100, 
                     nbtuba = 1000, radix = "Simu",modify.target = TRUE)

grid.simu1 <- anam.y2z(grid.simu1,name="Simu*",anam=model.anam1)

area1data=matrix(,100,length(db.extract(grid.simu.mean1,"Simu.Y2.S2")))
area1datanames=vector(,100)
for(i in 1:100)
{area1datanames[i]=paste("Raw.Simu.Y2.S",i,sep="")}

for(i in 1:100)
{area1data[i,]=db.extract(grid.simu.mean1,names=area1datanames[i])}
area1data[area1data<0]=0

saveRDS(area1data,"2016fishNASCsims.rds")
saveRDS(grid.simu.mean1,"2016fishNASCdb.rds")



############## 2015 fish NASC sim #####################



data=fishdata2015
data.db1=db.create(data,flag.grid=F,ndim=2,autoname=F)
projec.toggle(0)
db.data1=db.polygon(data.db1,newfullpoly)
projec.define(projection="mean", db=db.data1,flag.update=T) 
model.anam1=anam.fit(db.data1,type="emp",ndisc=db.data1$nech,draw=T,title="Anamophosis")
db.data.trans1=anam.z2y(db.data1,anam=model.anam1)
db.data.trans1 <- db.rename(db.data.trans1,name="Gaussian.NASC",newname="Yp")
ycut1 <- qnorm(sum(db.extract(db.data.trans1,"NASC") == 0) / length(db.extract(db.data.trans1,name="NASC",flag.compress=T)))
Y1<- db.extract(db.data.trans1,"Yp",flag.compress=T)
NASC1 <- db.extract(db.data.trans1,"NASC")
Y1[NASC1 == 0] <- ycut1
db.data.trans1 <- db.replace(db.data.trans1,"Yp",Y1)
n.H <- 50
vario.Yp1 <- vario.calc(db.data.trans1,lag=1,nlag=40)
vario.Y1  <- vario.trans.cut(vario.Yp1,ycut1,n.H)
model.vario.Y1 <- model.auto(vario.Y1,struc=melem.name(c(1,2)),draw=T)
Ymax1 <- db.extract(db.data.trans1,name="Yp",flag.compress=F)
Ymin1 <- db.extract(db.data.trans1,name="Yp",flag.compress=F)
Ymin1[Ymin1 <= ycut1] <- -10
db.data.trans1<-db.add(db.data.trans1,Ymax1)
db.data.trans1<-db.locate(db.data.trans1,db.data.trans1$natt,"upper")
db.data.trans1<-db.add(db.data.trans1,Ymin1)
db.data.trans1<-db.locate(db.data.trans1,db.data.trans1$natt,"lower")

db.data.trans1 <-gibbs(db = db.data.trans1, model = model.vario.Y1, seed = 2322, 
                       nboot = 10, niter = 500, flag.norm=FALSE, toleps = 1,
                       radix = "Gibbs", modify.target = TRUE)

db.data.trans1<-db.rename(db.data.trans1,"Gibbs.G1","Y2")

grid.simu1=db.create(flag.grid=T,x0=c(-5.9,54),dx=c(0.015839,0.008971),nx=c(100,600))
grid.simu1 <- db.polygon(grid.simu1,newfullpoly)
neigh.simu <- neigh.create(ndim=2,type=2,nmini=5,nmaxi=100)
grid.simu1 <- simtub(dbin=db.data.trans1, dbout=grid.simu1, model=model.vario.Y1, 
                     neigh=neigh.simu, uc = "", seed=10, mean = 0, nbsimu = 100, 
                     nbtuba = 1000, radix = "Simu" ,modify.target = TRUE)

grid.simu1 <- anam.y2z(grid.simu.mean1,name="Simu*",anam=model.anam1)

area1data=matrix(,100,length(db.extract(grid.simu1,"Simu.Y2.S2")))
area1datanames=vector(,100)
for(i in 1:100)
{area1datanames[i]=paste("Raw.Simu.Y2.S",i,sep="")}

for(i in 1:100)
{area1data[i,]=db.extract(grid.simu.mean1,names=area1datanames[i])}
area1data[area1data<0]=0

saveRDS(area1data,"2015fishNASCsims.rds")
saveRDS(grid.simu.mean1,"2015fishNASCdb.rds")


############## 2014 fish NASC sim #####################


data=fishdata2014
data.db1=db.create(data,flag.grid=F,ndim=2,autoname=F)
projec.toggle(0)
db.data1=db.polygon(data.db1,newfullpoly)
projec.define(projection="mean", db=db.data1,flag.update=T) 
model.anam1=anam.fit(db.data1,type="emp",ndisc=db.data1$nech,draw=T,title="Anamophosis")
db.data.trans1=anam.z2y(db.data1,anam=model.anam1)
db.data.trans1 <- db.rename(db.data.trans1,name="Gaussian.NASC",newname="Yp")

ycut1 <- qnorm(sum(db.extract(db.data.trans1,"NASC") == 0) / length(db.extract(db.data.trans1,name="NASC",flag.compress=T)))
Y1<- db.extract(db.data.trans1,"Yp",flag.compress=T)
NASC1 <- db.extract(db.data.trans1,"NASC")
Y1[NASC1 == 0] <- ycut1
db.data.trans1 <- db.replace(db.data.trans1,"Yp",Y1)
n.H <- 50
vario.Yp1 <- vario.calc(db.data.trans1,lag=1,nlag=40)
vario.Y1  <- vario.trans.cut(vario.Yp1,ycut1,n.H)
model.vario.Y1 <- model.auto(vario.Y1,struc=melem.name(c(1,2)),draw=T)
Ymax1 <- db.extract(db.data.trans1,name="Yp",flag.compress=F)
Ymin1 <- db.extract(db.data.trans1,name="Yp",flag.compress=F)
Ymin1[Ymin1 <= ycut1] <- -10

db.data.trans1<-db.add(db.data.trans1,Ymax1)
db.data.trans1<-db.locate(db.data.trans1,db.data.trans1$natt,"upper")
db.data.trans1<-db.add(db.data.trans1,Ymin1)
db.data.trans1<-db.locate(db.data.trans1,db.data.trans1$natt,"lower")

db.data.trans1 <-gibbs(db = db.data.trans1, model = model.vario.Y1, seed = 23722, 
                       nboot = 10, niter = 500, flag.norm=FALSE, toleps = 1,
                       radix = "Gibbs", modify.target = TRUE)
db.data.trans1<-db.rename(db.data.trans1,"Gibbs.G1","Y2")
grid.simu1=db.create(flag.grid=T,x0=c(-5.9,54),dx=c(0.015839,0.008971),nx=c(100,600))
grid.simu1 <- db.polygon(grid.simu1,newfullpoly)
neigh.simu <- neigh.create(ndim=2,type=2,nmini=5,nmaxi=100)

grid.simu1 <- simtub(dbin=db.data.trans1, dbout=grid.simu1, model=model.vario.Y1, 
                     neigh=neigh.simu, uc = "", mean = 0, seed = 27852, nbsimu = 100, 
                     nbtuba = 1000, radix = "Simu",modify.target = TRUE)

grid.simu1 <- anam.y2z(grid.simu1,name="Simu*",anam=model.anam1)

area1data=matrix(,100,length(db.extract(grid.simu1,"Simu.Y2.S2")))
area1datanames=vector(,100)
for(i in 1:100)
{area1datanames[i]=paste("Raw.Simu.Y2.S",i,sep="")}
for(i in 1:100)
{area1data[i,]=db.extract(grid.simu.mean1,names=area1datanames[i])}
area1data[area1data<0]=0

saveRDS(area1data,"2014fishNASCsims.rds")
saveRDS(grid.simu.mean1,"2014fishNASCdb.rds")


############## 2016 length/sigbs sim ################


lengthdata=read.table("1916Aclupeidmeanlengths.csv",header=T,sep=",")
lengthdata_sprat=lengthdata[which(lengthdata$n.sprat>19),c(1,2,5)]
lengthdata_herring=lengthdata[which(lengthdata$n.herring>19),c(1,2,7)]
projec.toggle(0)
lengthdata.db=db.create(lengthdata_sprat,flag.grid=F,ndim=2,autoname=F)

projec.toggle(0)
db.data1=db.polygon(lengthdata.db,newfullpoly)
projec.define(projection="mean", db=db.data1,flag.update=T) 
model.anam1=anam.fit(db.data1,type="gaus",nbpoly=12,draw=T,title="Anamophosis",xlim=c(-2,2))
db.data.trans1=anam.z2y(db.data1,anam=model.anam1)
vario.data1=vario.calc(db.data.trans1,lag=6,nlag=8)
model.vario1 <- model.auto(vario.data1,struc=melem.name(c(1,3)),draw=T)

grid.simu1=db.create(flag.grid=T,x0=c(-5.9,54),dx=c(0.015839,0.008971),nx=c(100,600))
grid.simu1 <- db.polygon(grid.simu1,newfullpoly)
neigh.simu <- neigh.create(ndim=2,type=0)
grid.simu1 <- simtub(dbin=db.data.trans1, dbout=grid.simu1, model=model.vario1, 
                     neigh=neigh.simu, uc = "", mean = 0, nbsimu = 100, 
                     nbtuba = 1000, radix = "Simu",modify.target = TRUE)

grid.simu1 <- anam.y2z(grid.simu1,name="Simu*",anam=model.anam1)

lengthnames=vector(,100)
for(i in 1:100)
{lengthnames[i]=paste("Raw.Simu.Gaussian.herring.length.S",i,sep="")}

x2014=matrix(,3710,100) # 500 columns
for(i in 1:100){
  x2014[,i]=db.extract(grid.simu1,names=lengthnames[i],flag.compress=T)
}

saveRDS(x2014,"2016herring_length_anam_sims.rds")
saveRDS(grid.simu1,"2016sprat_length_anam_sims_db.rds")



############## 2015 length/sigbs sim ################



lengthdata=read.table("1615Aclupeidmeanlengths.csv",header=T,sep=",")
lengthdata_sprat=lengthdata[which(lengthdata$n.sprat>19),c(1,2,5)]
lengthdata_herring=lengthdata[which(lengthdata$n.herring>19),c(1,2,7)]
projec.toggle(0)
plot(newfullpoly)
lengthdata.db=db.create(lengthdata_herring,flag.grid=F,ndim=2,autoname=F)
db.data1=db.polygon(lengthdata.db,newfullpoly)
projec.define(projection="mean", db=db.data1,flag.update=T) 
model.anam1=anam.fit(db.data1,type="gaus",nbpoly=12,draw=T,title="Anamophosis")
db.data.trans1=anam.z2y(db.data1,anam=model.anam1)
vario.data1=vario.calc(db.data.trans1,lag=8,nlag=5)
model.vario1 <- model.auto(vario.data1,struc=melem.name(c(2)),draw=T)
grid.simu1=db.create(flag.grid=T,x0=c(-5.9,54),dx=c(0.015839,0.008971),nx=c(100,600))
grid.simu1 <- db.polygon(grid.simu1,newfullpoly)
neigh.simu <- neigh.create(ndim=2,type=0)
grid.simu1 <- simtub(dbin=db.data.trans1, dbout=grid.simu1, model=model.vario1, 
                     neigh=neigh.simu, uc = "", mean = 0, nbsimu = 100, 
                     nbtuba = 1000, radix = "Simu",modify.target = TRUE)
grid.simu1 <- anam.y2z(grid.simu1,name="Simu*",anam=model.anam1)
lengthnames=vector(,100)
for(i in 1:100)
{lengthnames[i]=paste("Raw.Simu.Gaussian.herring.length.S",i,sep="")}
x2014=matrix(,3710,100)
for(i in 1:100){
  x2014[,i]=db.extract(grid.simu1,names=lengthnames[i],flag.compress=T)
}

saveRDS(x2014,"2015herring_length_anam_sims.rds")
saveRDS(grid.simu.length,"2015length_anam_sims_db.rds")


############## 2014 length/sigbs sim ################ 

lengthdata=read.table("2014Aclupeidmeanlengths.csv",header=T,sep=",")
lengthdata_sprat=lengthdata[which(lengthdata$n.sprat>19),c(1,2,5)]
lengthdata_herring=lengthdata[which(lengthdata$n.herring>19),c(1,2,7)]
projec.toggle(0)
lengthdata.db=db.create(lengthdata_herring,flag.grid=F,ndim=2,autoname=F)
db.data1=db.polygon(lengthdata.db,newfullpoly)
projec.define(projection="mean", db=db.data1,flag.update=T) 
model.anam1=anam.fit(db.data1,type="gaus",nbpoly=10,draw=T,title="Anamophosis")
db.data.trans1=anam.z2y(db.data1,anam=model.anam1)
vario.data1=vario.calc(db.data.trans1,lag=4,nlag=10)
model.vario1 <- model.auto(vario.data1,struc=melem.name(c(1)),draw=T)
grid.simu1=db.create(flag.grid=T,x0=c(-5.9,54),dx=c(0.015839,0.008971),nx=c(100,600))
grid.simu1 <- db.polygon(grid.simu1,newfullpoly)
neigh.simu <- neigh.create(ndim=2,type=0)
for(i in 1:100){
grid.simu1 <- simtub(dbin=db.data.trans1, dbout=grid.simu1, model=model.vario1, 
                     neigh=neigh.simu, uc = "", mean = 0, nbsimu = 1,seed=sample(1:100000,1), 
                     nbtuba = 1000, radix = paste("Simu_",i,sep=""),modify.target = TRUE)
}

grid.simu1 <- anam.y2z(grid.simu1,name="Simu*",anam=model.anam1)
lengthnames=vector(,100)
for(i in 1:100)
{lengthnames[i]=paste("Raw.Simu_",i,".Gaussian.herring.length.S1",sep="")}
x2014=matrix(,3710,100) # 500 columns
for(i in 1:100){
  x2014[,i]=db.extract(grid.simu1,names=lengthnames[i],flag.compress=T)
}

saveRDS(x2014,"2014herring_length_anam_sims.rds")
saveRDS(grid.simu1,"2014length_anam_sims_db.rds")



############## 2016 sprat numbers ############## ####from here####

y=2016
nascdata=readRDS(paste(y,"fishNASCsims.rds",sep=""))
sigprop=readRDS(paste(y,"sigprop_anam_sims.rds",sep=""))
lengths=readRDS(paste(y,"sprat_length_anam_sims.rds",sep=""))
sigbs=(10^-7.1666)*(lengths^2)

spratnasc=matrix(,3710,10000) # nasc for sprat only
for(i in 1:100)
{
  for(j in 1:100)
  {
    spratnasc[,(i+(100*(j-1)))]=nascdata[i,]*sigprop[,j]
  }
}
rm(nascdata,sigprop)

spratdens1=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens1[,((100*(i-1))+(j))]=spratnasc[,i]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens1,"2016spratdensities1.rds")
rm(spratdens1)

spratdens2=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens2[,((100*(i-1))+(j))]=spratnasc[,(i+2500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens2,"2016spratdensities2.rds")
rm(spratdens2)

spratdens3=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens3[,((100*(i-1))+(j))]=spratnasc[,(i+5000)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens3,"2016spratdensities3.rds")
rm(spratdens3)

spratdens4=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens4[,((100*(i-1))+(j))]=spratnasc[,(i+7500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens4,"2016spratdensities4.rds")
rm(spratdens4)

############## 2015 sprat numbers ##############

y=2015
nascdata=readRDS(paste(y,"fishNASCsims.rds",sep=""))
sigprop=readRDS(paste(y,"sigprop_anam_sims.rds",sep=""))
lengths=readRDS(paste(y,"sprat_length_anam_sims.rds",sep=""))
sigbs=(10^-7.1666)*(lengths^2)

spratnasc=matrix(,3710,10000) # nasc for sprat only
for(i in 1:100)
{
  for(j in 1:100)
  {
    spratnasc[,(i+(100*(j-1)))]=nascdata[i,]*sigprop[,j]
  }
}
rm(nascdata,sigprop)

spratdens1=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens1[,((100*(i-1))+(j))]=spratnasc[,i]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens1,"2015spratdensities1.rds")
rm(spratdens1)

spratdens2=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens2[,((100*(i-1))+(j))]=spratnasc[,(i+2500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens2,"2015spratdensities2.rds")
rm(spratdens2)

spratdens3=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens3[,((100*(i-1))+(j))]=spratnasc[,(i+5000)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens3,"2015spratdensities3.rds")
rm(spratdens3)

spratdens4=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens4[,((100*(i-1))+(j))]=spratnasc[,(i+7500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens4,"2015spratdensities4.rds")
rm(spratdens4)

############## 2014 sprat numbers ##############

y=2014
nascdata=readRDS(paste(y,"fishNASCsims.rds",sep=""))
sigprop=readRDS(paste(y,"sigprop_anam_sims.rds",sep=""))
lengths=readRDS(paste(y,"sprat_length_anam_sims.rds",sep=""))
sigbs=(10^-7.1666)*(lengths^2)

spratnasc=matrix(,3710,10000) # nasc for sprat only
for(i in 1:100)
{
  for(j in 1:100)
  {
    spratnasc[,(i+(100*(j-1)))]=nascdata[i,]*sigprop[,j]
  }
}
rm(nascdata,sigprop)

spratdens1=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens1[,((100*(i-1))+(j))]=spratnasc[,i]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens1,"2014spratdensities1.rds")
rm(spratdens1)

spratdens2=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens2[,((100*(i-1))+(j))]=spratnasc[,(i+2500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens2,"2014spratdensities2.rds")
rm(spratdens2)

spratdens3=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens3[,((100*(i-1))+(j))]=spratnasc[,(i+5000)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens3,"2014spratdensities3.rds")
rm(spratdens3)

spratdens4=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    spratdens4[,((100*(i-1))+(j))]=spratnasc[,(i+7500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(spratdens4,"2014spratdensities4.rds")
rm(spratdens4); beep(2)


############## 2016 herring numbers ##############

y=2016
nascdata=readRDS(paste(y,"fishNASCsims.rds",sep=""))
sigprop=readRDS(paste(y,"sigprop_anam_sims.rds",sep=""))
hsigprop=1-sigprop
lengths=readRDS(paste(y,"herring_length_anam_sims.rds",sep=""))
sigbs=(10^-7.1666)*(lengths^2)

herringnasc=matrix(,3710,10000) # nasc for herring only
for(i in 1:100)
{
  for(j in 1:100)
  {
    herringnasc[,(i+(100*(j-1)))]=nascdata[i,]*(hsigprop[,j])
  }
}
rm(nascdata,sigprop)

herringdens1=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens1[,((100*(i-1))+(j))]=herringnasc[,i]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens1,"2016herringdensities1.rds")
rm(herringdens1)

herringdens2=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens2[,((100*(i-1))+(j))]=herringnasc[,(i+2500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens2,"2016herringdensities2.rds")
rm(herringdens2)

herringdens3=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens3[,((100*(i-1))+(j))]=herringnasc[,(i+5000)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens3,"2016herringdensities3.rds")
rm(herringdens3)

herringdens4=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens4[,((100*(i-1))+(j))]=herringnasc[,(i+7500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens4,"2016herringdensities4.rds")
rm(herringdens4)

############## 2015 herring numbers ##############
y=2015

nascdata=readRDS(paste(y,"fishNASCsims.rds",sep=""))
sigprop=readRDS(paste(y,"sigprop_anam_sims.rds",sep=""))
hsigprop=1-sigprop
lengths=readRDS(paste(y,"herring_length_anam_sims.rds",sep=""))
sigbs=(10^-7.1666)*(lengths^2)

herringnasc=matrix(,3710,10000) # nasc for herring only
for(i in 1:100)
{
  for(j in 1:100)
  {
    herringnasc[,(i+(100*(j-1)))]=nascdata[i,]*(hsigprop[,j])
  }
}
rm(nascdata,sigprop)

herringdens1=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens1[,((100*(i-1))+(j))]=herringnasc[,i]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens1,"2015herringdensities1.rds")
rm(herringdens1)

herringdens2=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens2[,((100*(i-1))+(j))]=herringnasc[,(i+2500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens2,"2015herringdensities2.rds")
rm(herringdens2)

herringdens3=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens3[,((100*(i-1))+(j))]=herringnasc[,(i+5000)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens3,"2015herringdensities3.rds")
rm(herringdens3)

herringdens4=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens4[,((100*(i-1))+(j))]=herringnasc[,(i+7500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens4,"2015herringdensities4.rds")
rm(herringdens4)



############## 2014 herring numbers ##############

y=2014
nascdata=readRDS(paste(y,"fishNASCsims.rds",sep=""))
sigprop=readRDS(paste(y,"sigprop_anam_sims.rds",sep=""))
hsigprop=1-sigprop
lengths=readRDS(paste(y,"herring_length_anam_sims.rds",sep=""))
sigbs=(10^-7.1666)*(lengths^2)

herringnasc=matrix(,3710,10000) # nasc for herring only
for(i in 1:100)
{
  for(j in 1:100)
  {
    herringnasc[,(i+(100*(j-1)))]=nascdata[i,]*hsigprop[,j]
  }
}
rm(nascdata,sigprop)

herringdens1=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens1[,((100*(i-1))+(j))]=herringnasc[,i]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens1,"2014herringdensities1.rds")
rm(herringdens1)

herringdens2=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens2[,((100*(i-1))+(j))]=herringnasc[,(i+2500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens2,"2014herringdensities2.rds")
rm(herringdens2)


herringdens3=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens3[,((100*(i-1))+(j))]=herringnasc[,(i+5000)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens3,"2014herringdensities3.rds")
rm(herringdens3)

herringdens4=matrix(,3710,250000) # numbers per square kmMETRE
for(i in 1:2500)
{
  for(j in 1:100)
  {
    herringdens4[,((100*(i-1))+(j))]=herringnasc[,(i+7500)]/(sigbs[,j]*(4*pi*1.852^2))
  }
}
saveRDS(herringdens4,"2014herringdensities4.rds")
rm(herringdens4);beep(2)


############## 2016-14 sprat + herring densities to abundances ##############

y="2016sprat"
for(i in 1:4){
z=i
spratdens=readRDS(paste(y,"densities",z,".rds",sep=""))

for(j in 1:ncol(spratdens))
{spratdens[,j]=spratdens[,j]*areas[,7]}
saveRDS(spratdens,paste(y,"abundances",z,".rds",sep=""))
}

y="2015sprat"
for(i in 1:4){
  z=i
  spratdens=readRDS(paste(y,"densities",z,".rds",sep=""))

  for(j in 1:ncol(spratdens))
  {spratdens[,j]=spratdens[,j]*areas[,7]}
  saveRDS(spratdens,paste(y,"abundances",z,".rds",sep=""))
}

y="2014sprat"
for(i in 1:4){
  z=i
  spratdens=readRDS(paste(y,"densities",z,".rds",sep=""))

  for(j in 1:ncol(spratdens))
  {spratdens[,j]=spratdens[,j]*areas[,7]}
  saveRDS(spratdens,paste(y,"abundances",z,".rds",sep=""))
}


y="2016herring"
for(i in 1:4){
  z=i
  herringdens=readRDS(paste(y,"densities",z,".rds",sep=""))
  
  for(j in 1:ncol(herringdens))
  {herringdens[,j]=herringdens[,j]*areas[,7]}
  saveRDS(herringdens,paste(y,"abundances",z,".rds",sep=""))
}

y="2015herring"
for(i in 1:4){
  z=i
  herringdens=readRDS(paste(y,"densities",z,".rds",sep=""))

  for(j in 1:ncol(herringdens))
  {herringdens[,j]=herringdens[,j]*areas[,7]}
  saveRDS(herringdens,paste(y,"abundances",z,".rds",sep=""))
}

y="2014herring"
for(i in 1:4){
  z=i
  herringdens=readRDS(paste(y,"densities",z,".rds",sep=""))

  for(j in 1:ncol(herringdens))
  {herringdens[,j]=herringdens[,j]*areas[,7]}
  saveRDS(herringdens,paste(y,"abundances",z,".rds",sep=""))
}


############## collating sprat abundances #############


y="2016sprat"
for(z in 1:4){
  spratabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  temp=colSums(spratabund)
  saveRDS(temp,paste(y,"abundance_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"abundance_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_abundance_estimates_all",".rds",sep=""))

y="2015sprat"
for(z in 1:4){
  spratabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  temp=colSums(spratabund)
  saveRDS(temp,paste(y,"abundance_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"abundance_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_abundance_estimates_all",".rds",sep=""))

y="2014sprat"
for(z in 1:4){
  spratabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  temp=colSums(spratabund)
  saveRDS(temp,paste(y,"abundance_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"abundance_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_abundance_estimates_all",".rds",sep=""))

y="2016herring"
for(z in 1:4){
  herringabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  temp=colSums(herringabund)
  saveRDS(temp,paste(y,"abundance_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"abundance_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_abundance_estimates_all",".rds",sep=""))

y="2015herring"
for(z in 1:4){
  herringabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  temp=colSums(herringabund)
  saveRDS(temp,paste(y,"abundance_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"abundance_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_abundance_estimates_all",".rds",sep=""))

y="2014herring"
for(z in 1:4){
  herringabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  temp=colSums(herringabund)
  saveRDS(temp,paste(y,"abundance_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"abundance_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_abundance_estimates_all",".rds",sep=""))


############## 2016-4 sprat abundances to biomasses #############

y="2016sprat"
x="2016"
fishlengths=readRDS(paste(x,"sprat_length_anam_sims.rds",sep=""))
lengthcode=sort(rep(1:100,2500))
spratweights=0.0023*(fishlengths^3.469)

for(i in 1:4){
  z=i
  spratabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  for(j in 1:250000){
  spratabund[,j]=spratabund[,j]*spratweights[,lengthcode[j]]
  }
  saveRDS(spratabund,paste(y,"biomasses",z,".rds",sep=""))
}

y="2015sprat"
x="2015"
fishlengths=readRDS(paste(x,"sprat_length_anam_sims.rds",sep=""))
lengthcode=sort(rep(1:100,2500))
spratweights=0.0024*(fishlengths^3.4341)

for(i in 1:4){
  z=i
  spratabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  for(j in 1:250000){
    spratabund[,j]=spratabund[,j]*spratweights[,lengthcode[j]]
  }
  saveRDS(spratabund,paste(y,"biomasses",z,".rds",sep=""))
}


y="2014sprat"
x="2014"
fishlengths=readRDS(paste(x,"sprat_length_anam_sims.rds",sep=""))
lengthcode=sort(rep(1:100,2500))
spratweights=0.0025*(fishlengths^3.3849)

for(i in 1:4){
  z=i
  spratabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  for(j in 1:250000){
    spratabund[,j]=spratabund[,j]*spratweights[,lengthcode[j]]
  }
  saveRDS(spratabund,paste(y,"biomasses",z,".rds",sep=""))
}


############## 2016-4 herring abundances to biomasses #############


y="2016herring"
x="2016"
fishlengths=readRDS(paste(x,"herring_length_anam_sims.rds",sep=""))
lengthcode=sort(rep(1:100,2500))
herringweights=0.004*(fishlengths^3.1709)

for(i in 1:4){
  z=i
  herringabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  for(j in 1:250000){
    herringabund[,j]=herringabund[,j]*herringweights[,lengthcode[j]]
  }
  saveRDS(herringabund,paste(y,"biomasses",z,".rds",sep=""))
}

y="2015herring"
x="2015"
fishlengths=readRDS(paste(x,"herring_length_anam_sims.rds",sep=""))
lengthcode=sort(rep(1:100,2500))
herringweights=0.0023*(fishlengths^3.3722)

for(i in 1:4){
  z=i
  herringabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  for(j in 1:250000){
    herringabund[,j]=herringabund[,j]*herringweights[,lengthcode[j]]
  }
  saveRDS(herringabund,paste(y,"biomasses",z,".rds",sep=""))
}

y="2014herring"
x="2014"
fishlengths=readRDS(paste(x,"herring_length_anam_sims.rds",sep=""))
lengthcode=sort(rep(1:100,2500))
herringweights=0.0034*(fishlengths^3.2443)
for(i in 1:4){
  z=i
  herringabund=readRDS(paste(y,"abundances",z,".rds",sep=""))
  for(j in 1:250000){
    herringabund[,j]=herringabund[,j]*herringweights[,lengthcode[j]]
  }
  saveRDS(herringabund,paste(y,"biomasses",z,".rds",sep=""))
}

############## collating sprat biomass estimates ###########

y="2016sprat"
for(z in 1:4){
spratabund=readRDS(paste(y,"biomasses",z,".rds",sep=""))
temp=colSums(spratabund)
saveRDS(temp,paste(y,"biomass_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"biomass_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_biomass_estimates_all",".rds",sep=""))

y="2015sprat"
for(z in 1:4){
  spratabund=readRDS(paste(y,"biomasses",z,".rds",sep=""))
  temp=colSums(spratabund)
  saveRDS(temp,paste(y,"biomass_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"biomass_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_biomass_estimates_all",".rds",sep=""))

y="2014sprat"
for(z in 1:4){
  spratabund=readRDS(paste(y,"biomasses",z,".rds",sep=""))
  temp=colSums(spratabund)
  saveRDS(temp,paste(y,"biomass_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"biomass_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_biomass_estimates_all",".rds",sep=""))


############## collating herring biomass estimates ###########

y="2016herring"
for(z in 1:4){
  herringabund=readRDS(paste(y,"biomasses",z,".rds",sep=""))
  temp=colSums(herringabund)
  saveRDS(temp,paste(y,"biomass_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"biomass_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_biomass_estimates_all",".rds",sep=""))

y="2015herring"
for(z in 1:4){
  herringabund=readRDS(paste(y,"biomasses",z,".rds",sep=""))
  temp=colSums(herringabund)
  saveRDS(temp,paste(y,"biomass_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"biomass_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_biomass_estimates_all",".rds",sep=""))

y="2014herring"
for(z in 1:4){
  herringabund=readRDS(paste(y,"biomasses",z,".rds",sep=""))
  temp=colSums(herringabund)
  saveRDS(temp,paste(y,"biomass_estimates",z,".rds",sep=""))
}
for(z in 1:4){
  assign(paste("temp_",z,sep=""),readRDS(paste(y,"biomass_estimates",z,".rds",sep="")))
}  
saveRDS(c(temp_1,temp_2,temp_3,temp_4),paste(y,"_biomass_estimates_all",".rds",sep=""))


############ all biomass estimates ###############

a="herring"
for(z in 2014:2016){
  assign(paste(a,"_biomass_",z,sep=""),readRDS(paste("Final estimates/",z,a,"_biomass_estimates_all",".rds",sep="")))
}  

a="sprat"
for(z in 2014:2016){
  assign(paste(a,"_biomass_",z,sep=""),readRDS(paste("Final estimates/",z,a,"_biomass_estimates_all",".rds",sep="")))
}  

a="herring"
for(z in 2014:2016){
  assign(paste(a,"_abundance_",z,sep=""),readRDS(paste("Final estimates/",z,a,"_abundance_estimates_all",".rds",sep="")))
}  

a="sprat"
for(z in 2014:2016){
  assign(paste(a,"_abundance_",z,sep=""),readRDS(paste("Final estimates/",z,a,"_abundance_estimates_all",".rds",sep="")))
}  


bio_lists=c("herring_biomass_2014","sprat_biomass_2014","herring_biomass_2015","sprat_biomass_2015","herring_biomass_2016","sprat_biomass_2016")
results=matrix(,3,6)
colnames(results)=bio_lists
rownames(results)=c("mean","lowCI","highCI")

for(i in 1:6){
  bio_list=get(bio_lists[i])
  results[1,i] = mean(bio_list)/1000000000
  results[2,i] = (sort(bio_list)[round(length(bio_list)*.025)])/1000000000
  results[3,i] = (sort(bio_list)[round(length(bio_list)*.975)])/1000000000
}

results


abun_lists=c("herring_abundance_2014","sprat_abundance_2014","herring_abundance_2015","sprat_abundance_2015","herring_abundance_2016","sprat_abundance_2016")
results_abund=matrix(,3,6)
colnames(results_abund)=abun_lists
rownames(results_abund)=c("mean","lowCI","highCI")

for(i in 1:6){
  abun_list=get(abun_lists[i])
  results_abund[1,i] = mean(abun_list)/1000000000
  results_abund[2,i] = (sort(abun_list)[round(length(abun_list)*.025)])/1000000000
  results_abund[3,i] = (sort(abun_list)[round(length(abun_list)*.975)])/1000000000
}


results_abund



