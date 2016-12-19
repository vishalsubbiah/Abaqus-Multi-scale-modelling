function EngGrad()
tic;
load('EngInput1');
[N,~]=size(FinalData);
counter=1;

for i=2:N/10
    for j=1:i-1
        [~,n]=size(intersect(FinalData(i,2:9),FinalData(j,2:9)));
        if(n==4)
            EngGrad(counter,1)=abs(FinalData(i,10)-FinalData(j,10));
            MisAngle(counter,1)=Misorient(FinalData(i,11:13),FinalData(j,11:13));
            counter=counter+1;
        end
        
    end
end
Results=[MisAngle,EngGrad];
save('FinalResults1-10','Results');
time=toc;
save('Time10','time');
end

function [MisAngle] = Misorient(E1,E2)
M1=eulrotm(E1);
M2=eulrotm(E2);
R=M1*pinv(M2);
MisAngle=acos((R(1,1)+R(2,2)+R(3,3)-1)/2);
%A=vrrotmat2vec(R);
%MisAngle=A(4);
end

function [M]= eulrotm(E)
c1=cosd(E(1));
c=cosd(E(2));
c2=cosd(E(3));
s1=sind(E(1));
s=sind(E(2));
s2=sind(E(3));
M=[c1*c2-s1*s2*c,s1*c2+c1*s2*c,s2*s;-c1*s2-s1*c2*c, -s1*s2+c1*c2*c, c2*s; s1*s, -c1*s,c];
end
