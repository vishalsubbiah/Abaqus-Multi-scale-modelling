library('R.matlab')
FinResults=readMat('FinResults1.mat')
FinResults=FinResults$FinResults
FinResults[,1]=FinResults[,1]*180/pi
MisorientFin=FinResults[,1]
EngFin=FinResults[,2]
tFin=cor(MisorientFin,EngFin)
yFin=ccf(MisorientFin,EngFin,type="correlation",ylab="Correlation",xlab="lag",cex.lab=1.3,main=" ")
#dev.copy(jpeg,filename="PA_CCF10.jpeg");
#dev.off ();
boxplot(EngFin,ylab="Energy (J)", cex.lab=1.3)
dev.copy(jpeg,filename="PA_BoxPlotEng.jpeg")
dev.off();

Results=readMat('FinalResults1.mat')
Results=Results$Results
Results[,1]=Results[,1]*180/pi
Misorient=Re(Results[,1])
Eng=Re(Results[,2])
t=cor(Misorient,Eng)
y=ccf(Misorient,Eng,type="correlation",ylab="Correlation",xlab="lag",cex.lab=1.3,main=" ")
#dev.copy(jpeg,filename="BA_CCF10.jpeg");
#dev.off ();
boxplot(Eng,ylab="Energy (J)", cex.lab=1.3)
dev.copy(jpeg,filename="BA_BoxPlotEng.jpeg")
dev.off();