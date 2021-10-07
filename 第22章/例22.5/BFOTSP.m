function [best_route,best_f]=BFOTSP(city,bacterialnum,nc,nre,ned,pmax)%细菌觅食算法,求TSP
[citynum,d]=city2d(city);
for i=1:bacterialnum
    bacterial(i,:)=randperm(citynum);
    fitness(i,1)=value(bacterial(i,:),d);
end
[a,b]=sort(fitness,'ascend');
best_route=bacterial(b(1),:);
best_f=a(1);
for l=1:ned      %迁徙
     for n=1:nre  %繁殖
        for m=1:nc   %趋向
            k=randperm(bacterialnum);
            g1=bacterial(k(1:ceil(bacterialnum/2)),:);  %分成两组
            fitness1=fitness(k(1:ceil(bacterialnum/2)),1);
            [a,b]=sort(fitness1,'ascend');
            g1=g1(b,:);
            g2=bacterial(k(ceil(bacterialnum/2)+1:end),:);
            fitness2=fitness(k(ceil(bacterialnum/2)+1:end),:);
            [a,b]=sort(fitness2,'descend');
            g2=g2(b,:);
            new=[];f=[];
            for i=1:bacterialnum/2
                [new1,new2,f1,f2]=greedycross(g1(i,:),g2(i,:),d);
                new=[new;new1;new2];
                f=[f;f1;f2];
            end
            bacterial=new;
            fitness=f;
        end
        [a,b]=sort(fitness,'ascend');    %以适应度衡量
        bacterial=bacterial(b,:);
        sr=bacterialnum*uint16((11-n)*(101-l)/(ned*100))/2;
        for i=1:sr
            bacterial(bacterialnum-i+1,:)=bacterial(i,:);
        end  
    end
    for i=1:bacterialnum
        fitness(i,1)=value(bacterial(i,:),d);
    end
    [a,b]=sort(fitness,'ascend');
    favg=mean(fitness);
    best_route=bacterial(b(1),:);
    fmax=a(1);
    for i=1:bacterialnum
        if fitness(i,1)>favg
            ped=pmax*exp(((fitness(i,1)-favg)*(l-ned))/((ned-1)*(fmax-favg+1))); 
        else
            ped=pmax;
        end
        if ped>rand
           point=ceil(rand*citynum);
           n_retail=uint16(citynum*(51-l)/100);
           if citynum-point>n_retail
               new=best_route(point+1:point+n_retail);
               new1=[best_route(1:point) best_route(point+n_retail+1:end)];
           else
               new=best_route(point+1:end);
               a1=n_retail-length(new);
               new=[new best_route(1:a1)];
               new1=best_route(a1+1:point);
           end
           a1=length(new1);
           k=randperm(a1);
           bacterial(i,:)=[new1(k) new];
        end
    end
    fitness(i,1)=value(bacterial(i,:),d);
end
[a,b]=min(fitness);
if a<best_f
   best_f=a;
   best_route=bacterial(b,:);
end
TSPplot(city,best_route);










