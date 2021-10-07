function [bestx,fval]=ECA(fun,runfornum,Gperson,Pperson,iter_max,sigma,R_range,LB,UB)   %竞选算法
NC=size(LB,1);
if length(Pperson)==1
    Pperson=Pperson.*ones(1,runfornum);
end
if length(Pperson)<runfornum
     n=runfornum-length(Pperson);
     Pperson=[Pperson Pperson(end).*ones(1,n)];
elseif  length(Pperson)>runfornum
    Pperson=Pperson(1:runfornum);
end
if length(R_range)==1
    R_range=[R_range R_range];
end
if length(sigma)==1
    sigma=[sigma sigma];
end
for i=1:runfornum
   runforLocation(i,:)=LB'+(UB-LB)'.*rand(1,NC);   %竞选人
   runforRestige(i,1)=fun(runforLocation(i,:));
end
[a,b]=sort(runforRestige,'ascend');
pmax=a(1);pmin=a(end);
bestx=runforLocation(b(1),:);
fval=a(1);
voterNum=Gperson+sum(Pperson);
for iter=1:iter_max
    for i=1:runfornum
        runforRange(i,1)=(runforRestige(i,1)-pmin)*(R_range(2)-R_range(1))/(pmax-pmin)+R_range(1);
        sigma(i,1)=sigma(2)-(runforRestige(i,1)-pmin)/(pmax-pmin)*(sigma(2)-sigma(1));
    end
    for i=1:Gperson  %全局选民位置
        voterLocation(i,:)=LB'+(UB-LB)'.*rand(1,NC);
        voterFitness(i,1)=fun(voterLocation(i,:));
    end
    for i=1:runfornum     %局部选民
       for j=1:Pperson(i)
           L(j,:)=runforLocation(i,:)+sigma(i,1)*sin(2*pi*rand)*sqrt(-2*log(rand));
           L(j,:)=boundtest(L(j,:),LB,UB);
           F1(j,1)=fun(L(j,:)); 
           if F1(j,1)<runforRestige(i,1)   %局部选民转换成竞选者
               temp1=runforLocation(i,:);
               runforLocation(i,:)=L(j,:);
               L(j,:)=temp1;
               temp2=runforRestige(i,1);
               runforRestige(i,1)=F1(j,1);
               F1(j,1)=temp2;
           end
       end
       voterLocation=[voterLocation;L];
       voterFitness=[voterFitness;F1];
    end   
    voterDist=zeros(voterNum,runfornum);
    for i=1:voterNum   %竞选者与选民的距离
        for j=1:runfornum
            voterDist(i,j)=norm(voterLocation(i,:)-runforLocation(j,:));
        end
    end
    voterF=zeros(voterNum,runfornum);
    total_voterF=zeros(voterNum,1);
    for i=1:voterNum   %每个选民受到的影响力
        for j=1:runfornum 
            if runforRange(j,1)>=voterDist(i,j)
                voterF(i,j)=(runforRange(j,1)-voterDist(i,j))*runforRestige(j,1)/runforRange(j,1);
            end
        end
        total_voterF(i,1)=sum(voterF(i,:));
    end
    voterS=zeros(voterNum,runfornum);
    for i=1:voterNum   %选民的支持力
        for j=1:runfornum
            if total_voterF(i,1)~=0
              voterS(i,j)=voterF(i,j)*voterFitness(i,1)/total_voterF(i,1);
            end
        end 
    end
    runforS=sum(voterS);%竞选者获得的支持
    voterQ=zeros(voterNum,runfornum);
    for i=1:voterNum   %选民的贡献
        for j=1:runfornum
            if runforS(j)~=0
               voterQ(i,j)=voterS(i,j)/runforS(j);
            end
        end
    end
    for i=1:runfornum  %竞选者的重心
        runforLocation(i,:)=zeros(1,NC);
        for j=1:voterNum
           runforLocation(i,:)=runforLocation(i,:)+voterLocation(j,:)*voterQ(j,i);
        end
        runforLocation(i,:)=boundtest(runforLocation(i,:),LB,UB);
        runforRestige(i,1)=fun(runforLocation(i,:));
    end
    for i=1:Gperson  %全局选民转换为竞选者
        for j=1:runfornum
            if voterFitness(i,1)<runforRestige(j,1)
               temp=runforLocation(j,:);
               runforLocation(j,:)=voterLocation(i,:);
               voterLocation(i,:)=temp;
               temp1=runforRestige(j,1);
               runforRestige(j,1)=voterFitness(i,1);
               voterFitness(i,1)=temp1;
            end
        end
    end
    [a,b]=sort(runforRestige,'ascend');
    favg=mean(runforRestige);
    best=runforLocation(b(1),:);
    if max(abs(runforRestige-favg))>=1
        f=1;
    else
        f=max(abs(runforRestige-favg));
    end
    sigma1=sum(((runforRestige-favg)./f).^2)/runfornum;
   %px=(pmax-pmin)*sigma/iter_max+(pmin-pmax)*2*sigma/iter_max+pmax;
    px=0.7*(1-sigma1/runfornum);
    if rand<px
        for i=1:runfornum
            runforLocation(i,:)=runforLocation(i,:)+(best-runforLocation(i,:))*rand;
            runforLocation(i,:)=boundtest(runforLocation(i,:),LB,UB);
            runforRestige(i,1)=fun(runforLocation(i,:));
        end
    end
    [a,b]=sort(runforRestige,'ascend');
    pmax=a(1);pmin=a(end);
    if a(1)<fval
        fval=a(1);
        bestx=runforLocation(b(1),:);
    end 
end
    
               
               
        
            
        
    
        
    
    
    
    
            
            
        
    
    
                 


