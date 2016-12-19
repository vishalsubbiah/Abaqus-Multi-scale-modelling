clear;
clc;
load('FinalResults1');
[n,~]=size(Results);
Results=sort(Results);
counter=1;
t=1;
for i=2:n
    if(Results(i-1,1)~=Results(i,1))
        FinResults(counter,:)=[Results(t,1);mean(Results(t:i-1,2))];
        counter=counter+1;
        t=i;
    end
end
save('FinResults1','FinResults');    
        
        
    