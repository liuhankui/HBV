
Args<-commandArgs()

library(grid)
library(ggplot2)
all<-read.table(Args[6],sep="\t",head=F,colC=c("character","numeric","numeric"))
names(all)<-c("chr","site","p")
all$p[all$p==0]<-NA
chrname<-data.frame(name=unique(all$chr))
qq<-data.frame(obs=-log10(sort(all$p)))
qq$exp<- -log10(seq(length(qq$obs))/(length(qq$obs)+1))

add<-0
interval<-max(all$site)/20
xbreak<-vector()
stat<-data.frame()
for(i in chrname$name){
    temp<-subset(all,chr==i,select=c(chr,site,p))
    temp$xlab<-temp$site+add+interval
    add<-max(temp$xlab)
    xbreak<-c(xbreak,sum(range(temp$xlab))/2)
    stat<-rbind(stat,temp)
}
chrname$xbreak<-xbreak
stat$chr<-factor(stat$chr,levels=chrname$name,order=T)
pdf(Args[7],height=5,width=10)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,8)))
p<-ggplot(stat,aes(xlab,-log10(p)))+
geom_point(aes(colour=chr))+
scale_x_continuous("Chromosome",breaks=chrname$xbreak,labels=chrname$name,expand=c(0.005,0.005))+
scale_y_continuous(expression(paste(-log[10],"P")),expand=c(0.005,0.01),limits=c(0,7))+
scale_colour_manual(values=rep(c("olivedrab","orange1"),len=length(chrname$name)))+
labs(title="Manhattan plot")+geom_hline(yintercept=-log10(0.05/nrow(stat)),colour="red")+
theme(panel.grid=element_blank(),panel.background=element_blank(),axis.line=element_line(),axis.text=element_text(colour="black"),axis.text.x=element_text(angle=-45,hjust=0),legend.position="none")
q<-ggplot(qq,aes(exp,obs))+
geom_point()+
geom_abline(intercept=0,slope=1,colour="lightblue")+
scale_x_continuous(expression(paste("Expected ",-log[10],"P")),expand=c(0,0.1))+
scale_y_continuous(expression(paste("Observed ",-log[10],"P")),expand=c(0,0.1))+
labs(title="QQ plot")+
theme(panel.grid=element_blank(),panel.background=element_blank(),axis.line=element_line(),axis.text=element_text(colour="black"),legend.position="none")
print(p,vp=viewport(layout.pos.row=1:2,layout.pos.col=1:6))
print(q,vp=viewport(layout.pos.row=2,layout.pos.col=7:8))

all<-read.table(Args[8],head=F,sep="\t")
all<-data.frame(phenotype=factor(all[,6],order=T))
pq<-ggplot(all,aes(x=phenotype))+
geom_bar(width=0.5)+
labs(title="Histogram")+
scale_x_discrete(breaks=c(2,1),labels=c("case","control"))+
theme(panel.grid=element_blank(),panel.background=element_blank(),axis.line=element_line(),axis.text=element_text(colour="black"))
print(pq,vp=viewport(layout.pos.row=1,layout.pos.col=7:8))
dev.off()
print(warnings())
q(save="no")
