function [bestx,fval]=SSO(fun,popsize,iter_max,pm,LB,UB)  %蜘蛛算法
NC=size(LB,1);
Nf=floor((0.9-rand*0.25)*popsize);  %雌性个体
Nm=popsize-Nf;                    %雄性个体
pf=exp(-(0.1:(3-0.1)/(iter_max-1):3));
for i=1:popsize
    spider(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    r1(i,1)=(max(spider(i,:))-min(spider(i,:)))/(2*popsize);
    fitness(i,1)=fun(spider(i,:));
end
d=squareform(pdist(spider));   %距离矩阵
[a,b]=sort(fitness,'ascend');
bestx=spider(b(1),:);
fval=a(1);
best=a(1);worst=a(end);
for iter=1:iter_max
    %w=wmax-iter/(iter_max-iter)*(wmax-wmin);
    y1=0;y2=0;
    for i=1:Nf+Nm
        w(i,1)=0.001+(fitness(i,1)-worst)/(best-worst);
        if i>Nf
           y1=y1+spider(i,:)*w(i,1);
           y2=y2+w(i,1);
        end
    end
    MalespiderCenter=y1./y2;
    [a1,b1]=sort(w,'descend');
    sb=spider(b1(1),:);   %全局最优蜘蛛
    [a2,b2]=sort(w(Nf+1:end,1),'descend');
    MalespiderW=w(Nf+b2,1);  %雄性蜘蛛权重
    spider(Nf+1:end,:)=spider(Nf+b2,:);
    MalespiderW_middle=MalespiderW(ceil(Nm/2),1);   %中间权重
    for i=1:Nf+Nm 
        vibb=a1(1)*exp(-d(i,b1(1))^2);
        if i<=Nf    %雌性个体
           temp=d(i,1:Nf);
           temp(i)=inf;
           [a3,b3]=sort(temp,'ascend');  
           for j=1:Nf
               if w(b3(j),1)>w(i,1)
                  sc=spider(b3(j),:);
                  vibc=w(b3(j),1)*exp(-d(i,b3(j))^2);  %最近的雌性蜘蛛
                  break
               end
           end
           if exist('sc','var')==0
                sc=spider(b3(1),:);
                vibc=w(b3(1),1)*exp(-d(i,b3(1))^2);
           end
           if rand<pf(iter)
               new(i,:)=spider(i,:)+rand*vibc*(sc-spider(i,:))+rand*vibb*(sb-spider(i,:))+rand*(rand-0.5);
           else
               new(i,:)=spider(i,:)-rand*vibc*(sc-spider(i,:))-rand*vibb*(sb-spider(i,:))+rand*(rand-0.5);
           end
        elseif i>Nf  %雄性个体
            temp1=d(i,1:Nf);
            [a4,b4]=min(temp1);
            sf=spider(b4,:);
            vibf=w(b4,1)*exp(-d(i,b4)^2);
            if w(i,1)>MalespiderW_middle
                new(i,:)=spider(i,:)+rand*vibf*(sf-spider(i,:))+rand*(rand-0.5);
            else
                new(i,:)=spider(i,:)+rand*(MalespiderCenter-spider(i,:));
            end
        end
        new(i,:)=boundtest(new(i,:),LB,UB);
        new_r(i,1)=(max(new(i,:))-min(new(i,:)))/(2*popsize);
        new_F(i,1)=fun(new(i,:));
    end
    pop1=[spider(1:Nf,:);new(1:Nf,:)];
    pop2=[spider(Nf+1:end,:);new(Nf+1:end,:)];
    f1=[fitness(1:Nf,1);new_F(1:Nf,1)];
    f2=[fitness(Nf+1:end,1);new_F(Nf+1:end,1)];
    R1=[r1(1:Nf,1);new_r(1:Nf,1)];
    R2=[r1(Nf+1:end,1);new_r(Nf+1:end,1)];
    [a4,b4]=sort(f1,'ascend');
    Femalespider_minF=a4(Nf);
    spider(1:Nf,:)=pop1(b4(1:Nf),:);
    fitness(1:Nf,1)=f1(b4(1:Nf),1);
    r1(1:Nf,1)=R1(b4(1:Nf),1);
    [a5,b5]=sort(f2,'ascend');
    Malespider_minF=a5(Nm);
    spider(Nf+1:end,:)=pop2(b5(1:Nm),:);
    fitness(Nf+1:end,1)=f2(b5(1:Nm),1);
    r1(Nf+1:end,1)=R2(b5(1:Nm),1);
    r=sum(r1);
    for i=Nf+1:Nf+ceil(Nm/2)
        temp3=d(i,1:Nf);
        Tg=find(temp3<=r);
        if ~isempty(Tg)
            ps=w(Tg,1)./sum(w(Tg,1));
            for j=1:NC
                flag=0;
                for k=1:length(Tg)
                    PS=sum(ps(1:k));
                    if rand<PS
                        new1(j)=spider(Tg(k),j);
                        flag=1;
                        break;
                    end 
                end        
                if flag==0
                    new1(j)=LB(j,1)+(UB(j,1)-LB(j,1))*rand;
                end
            end
            new1=boundtest(new1,LB,UB);
            new1_F=fun(new1);
            [a6,b6]=max([Femalespider_minF Malespider_minF]);
            if new1_F<a6
                if b6==1
                    spider(Nf,:)=new1;
                    fitness(Nf,1)=new1_F;
                    r1(Nf,1)=(max(new1)-min(new1))/(2*popsize);
                else
                    spider(end,:)=new1;
                    fitness(end,1)=new1_F;
                    r1(end,1)=(max(new1)-min(new1))/(2*popsize);
                end
            end
        end
    end
    [a7,b7]=min(fitness);
    if a7<fval
        fval=a7;
        bestx=spider(b7,:);
    end
    favg=mean(fitness);
    if max(abs(fitness-favg))>=1
        f=1;
    else
        f=max(abs(fitness-favg));
    end
    sigma=sum(((fitness-favg)./f).^2)/popsize;
   % px=(pmax-pmin)*sigma/iter_max+(pmin-pmax)*2*sigma/iter_max+pmax;
    pm=pm*(1-sigma/popsize); 
    for i=1:popsize   %差分变异
        if rand<pm
            new2=spider(i,:);
            pos=sort(ceil((Nf+Nm)*rand(1,2)));
            while pos(1)==pos(2)||pos(1)==i||pos(2)==i
               pos=sort(ceil((Nf+Nm)*rand(1,2)));
            end
            new2=new2+rand*(bestx-new2)+rand*(spider(pos(1),:)-spider(pos(2),:));
            new2=boundtest(new2,LB,UB);
            new2_F=fun(new2);
            if new2_F<fitness(i,1)
                spider(i,:)=new2;
                fitness(i,1)=new2_F;
                r1(i,1)=(max(new2)-min(new2))/(2*popsize);
            end
        end
    end
    d=squareform(pdist(spider));   %距离矩阵
    [a,b]=sort(fitness,'ascend');
    best=a(1);worst=a(end);
    if a(1)<fval
       bestx=spider(b(1),:);
       fval=a(1);
    end 
end
    
    
                
                
                
    
    
        
            
                
        
        
        
   
        
        
            
               
               
           
           
    
    
    
    

