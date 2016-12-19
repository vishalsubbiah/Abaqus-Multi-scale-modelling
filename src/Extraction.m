function [FinalData]=Extraction()
clear;
clc;

ElementNode=importdata('ElementNode');
Data=importdata('Data');

Elements=unique(Data(:,1));
Nodes=ElementNode.data;

[NoElements,~]=size(Elements);
NoIntPoints=8;
EngAvg=zeros(NoElements,1);
Theta=zeros(NoElements,3);

for i=1:NoElements
    VolAvg=0;
    for j=1:NoIntPoints
        EngAvg(i)=EngAvg(i)+energy(Data(NoIntPoints*(i-1)+j,3:8),Data(NoIntPoints*(i-1)+j,10:15),Data(NoIntPoints*(i-1)+j,9));
        VolAvg=VolAvg+Data(NoIntPoints*(i-1)+j,9);
    end
    EngAvg(i)=EngAvg(i)/VolAvg;
end

for i=1:NoElements
    Theta(i,1:3)=Data(NoIntPoints*(i-1)+1,16:18);
end

FinalData=[Elements,Nodes,EngAvg,Theta];
save('EngInput1','FinalData');

end

function [Eng] = energy(E,S,Vol)
    Eng=0.5*(E(1)*S(1)+E(2)*S(2)+E(3)*S(3)+2*(E(4)*S(4)+E(5)*S(5)+E(6)*S(6)))*Vol;
end
        

