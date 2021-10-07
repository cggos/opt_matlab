function cbest=memeticTSP(city,popsize,searchnum,iter_max)%memetic算法求解TSP
%myval=[80 3 0.6 0.05 100 8];
%popsize=myval(1);searchnum=myval(2);pc=myval(3);pm=myval(4);iter_max=myval(5);num=myval(6);
%LB=[-1;-1];UB=[2;2];eps1=1e-3;
[m,dis]=city2d(city);  
for i=1:popsize
    pop(i).route=randperm(m); 
    pop(i).fitness=value(pop(i).route,dis);
    f(i)=pop(i).fitness;
end
[a,b]=sort(f);
cbest=pop(b(1));
for i=1:iter_max
    new=crossover_memetic(pop,f,1);
    pop=local_search_memetic(new,dis,searchnum,4);
    for j=1:popsize
       pop(j).fitness=value(pop(j).route,dis);
       f(j)=pop(j).fitness;
    end
    new=crossover_memetic(pop,f,2);
    pop=local_search_memetic(new,dis,searchnum,4);
    pop=mutation_memetic(pop,dis,f,3);
    pop=local_search_memetic(pop,dis,searchnum,1);
    for j=1:popsize
       pop(j).fitness=value(pop(j).route,dis);
       f(j)=pop(j).fitness;
    end
    [a,b]=sort(f);
    pop=selection_MTC(pop);
    pop(b(end))=pop(b(1));
    cbest1=pop(b(1));
    if cbest1.fitness<cbest.fitness     %求极小值
         cbest=cbest1;
    end
end
Route=cbest.route;        %搜索到的最优解的环游路径
TSPplot(city,Route);      %画图
gname;

function new=crossover_memetic(pop,f,type)  %交叉算子，顺序交叉
popsize=size(pop,2);
NC=length(pop(1).route);
[a,b]=sort(f);
if type==1
    M=popsize/2;
elseif type==2
    M=popsize;
end
for i=1:M
    if type==1
       a=(i*2)-1;
       xpo=a+1;
    elseif type==2
       a=i;
       xpo=b(1);
    end
    x_max=max(pop(a).fitness,pop(xpo).fitness);
    if x_max>=mean(f)
        pc=0.8*(f(b(end))-x_max)/(f(b(end))-mean(f));
    else
        pc=0.8;
    end
    if rand<pc
        idex=[a xpo];
        pop3=[pop(a) pop(xpo)];
        for j=1:2
           point=sort(ceil(rand(1,2)*NC));   
           if point(1)==point(2)
              point=sort(ceil(rand(1,2)*NC));
           end
           new(idex(j)).route=zeros(1,NC);
           new(idex(j)).route(point(1):point(2))=pop3(j).route(point(1):point(2));
           temp1=pop3(j).route(point(1):point(2));
           if j==1
               temp=pop3(j+1).route;
           else
               temp=pop3(j-1).route;
           end
           [y1,y2]=mycompare1(temp,temp1);
           temp=redu(temp,y2,'c');
           if point(1)==1
               new(idex(j)).route(point(2)+1:end)=temp;
           elseif point(2)==NC
               new(idex(j)).route(1:point(1)-1)=temp;
           else
               new(idex(j)).route(1:point(1)-1)=temp(1:point(1)-1);
               new(idex(j)).route(point(2)+1:end)=temp(point(1):end);
           end
           new(idex(j)).fitness=0;
        end
     else
        new(a)=pop(a);
        new(xpo)=pop(xpo);
    end  
end

function pop=local_search_memetic(pop,dis,searchnum,type)   %memetic算法的局部搜索算子
NC=length(pop(1).route);
popsize=size(pop,2);
y=exhaus({'a';'b';'c'});
y=redu(y,findpos(y,[1 2 3],'r'),'r');
[r,c]=size(y);
for i=1:popsize
   for t=1:searchnum       %爬山法   
      pop_search=pop(i);
      if type==1
         k=ceil(rand(1,2)*NC);
         while k(1)==k(2)
              k=ceil(rand(1,2)*NC);
          end
         temp=pop_search.route(k(1));
         pop_search.route(k(1))=pop_search.route(k(2));
         pop_search.route(k(2))=temp;  
      elseif type==2    %贪婪倒位变异算子
          pop_search=pop(i);
          point=ceil(rand*NC);
          while point==NC
              point=ceil(rand*NC);
          end
          temp=pop_search.route;
          if point==1
              temp=redu(temp,point:point+1,'c');
          else
              temp=redu(temp,point-1:point+1,'c');
          end
          minDis=dis(pop_search.route(point),temp(1));
          point1=temp(1);
          for k=2:length(temp)
             if dis(pop_search.route(point),temp(k))<minDis
                point1=temp(k);
             end
          end
          point1=find(pop_search.route==point1);
          if point1<point
             pop_search.route(point1:point-1)=fliplr(pop(i).route(point1:point-1));
          else
             pop_search.route(point+1:point1)=fliplr(pop(i).route(point+1:point1));
          end
      elseif type==3    %递归插孤算子
          num=0;
          while num<5
             point=sort(ceil(rand(1,2)*NC));
             while 1
                if point(2)-point(1)<2||NC-(point(2)-point(1))<=2
                   point=sort(ceil(rand(1,2)*NC));
                else
                   break
                end
             end
             d1=dis(pop_search.route(point(1)),pop_search.route(point(1)+1));   %a
             d2=dis(pop_search.route(point(2)),pop_search.route(point(2)-1));   %b
             d3=dis(pop_search.route(point(1)),pop_search.route(point(1)));
             delta1=d1+d2-d3;
             if point(2)<NC
                 temp=[pop_search.route(point(2)+1:end) pop_search.route(1:point(1)-1)];
             else
                 temp=[pop_search.route(1:point(1)-1) pop_search.route(point(2)+1:end)];
             end 
             flag=0;
             for j=1:length(temp)
                 for k=1:length(temp)
                     if k~=j
                        d1=dis(temp(j),pop_search.route(point(1)+1));   %m
                        d2=dis(temp(k),pop_search.route(point(2)-1));   %n   
                        d3=dis(pop_search.route(temp(j)),pop_search.route(temp(k)));
                        delta2=d1+d2-d3;
                        if delta2<delta1
                            m=[find(temp(j)==pop_search.route) find(temp(k)==pop_search.route)];
                            flag=1;
                            break;
                        end
                     end
                 end
                 if flag==1
                     break
                 end
             end
             if exist('m','var')
                 temp1=[];
                 if m(1)>point(2)
                    temp1=[temp1 pop_search.route(point(1)) pop_search.route(point(2):m(1))];
                 else
                    temp1=[temp1 pop_search.route(point(1)) pop_search.route(point(2):end) pop_search.route(1:m(1))];
                 end
                 temp1=[temp1 pop_search.route(point(1)+1:point(2)-1)];
                 temp1=[temp1 pop_search.route(m(2):point(1)-1)];
                 temp1=[temp1 pop_search.route(m(1)+1:m(2)-1)];
                 new(i).route=temp1;
                 break
              else
                 num=num+1;
              end
          end
      elseif type==4   %三个体变异
          point=sort(ceil(rand(1,3)*NC));
          while 1
            if point(2)==point(1)||point(3)==point(1)||point(2)==point(3)
               point=sort(ceil(rand(1,3)*NC));
            else
               break
            end
          end
          for j=1:r
              new1(j)=pop_search;
              for k=1:c
                new1(j).route(point(k))=pop_search.route(point(y(j,k)));
              end
              newF(j)=value(new1(j).route,dis);
          end
          [a,b]=min(newF);
          pop_search=new1(b);
      end
      pop_search.fitness=value(pop_search.route,dis);
      if pop_search.fitness<pop(i).fitness
         pop(i)=pop_search;  
      end
   end
end

function new=mutation_memetic(pop,dis,f,type) %变异算子
NC=length(pop(1).route);
popsize=size(pop,2);
for i=1:popsize
   new(i)=pop(i);
   if type==1
      point=sort(ceil(rand(1,3)*NC));
      while 1
         if point(2)==point(1)||point(3)==point(1)||point(2)==point(3)
             point=sort(ceil(rand(1,3)*NC));
         else
             break
         end
      end
      y=exhaus({'a';'b';'c'});
      y=redu(y,findpos(y,[1 2 3],'r'),'r');
      [r,c]=size(y);
      for j=1:r
         new1(j,:)=new(i).route;
         for k=1:c
             new1(j,point(k))=new(i).route(point(y(j,k)));
         end
         newF(j)=value(new1(j,:),dis);
      end
      [a,b]=min(newF);
      new(i).route=new1(b,:);
   else  
      [a,b]=sort(f);
      if pop(i).fitness>=mean(f)
         pm=0.1*(f(b(end))-pop(i).fitness)/(f(b(end))-mean(f));
      else
         pm=0.1;
      end
      test=rand;
      if test<pm
         if type==2    %一点交换
            point=sort(ceil(rand(1,2)*NC));
            while point(1)==point(2)
                point=sort(ceil(rand(1,2)*NC));
            end
            temp=new(i).route(point(1));
            new(i).route(point(1))=new(i).route(point(2));
            new(i).route(point(2))=temp;
         elseif type==3     %二点交换      
            point=sort(ceil(rand(1,2)*NC));
            while point(1)==point(2)||point(2)-point(1)<=2
                point=sort(ceil(rand(1,2)*NC));
            end
            temp=new(i).route(point(1):point(1)+1);
            new(i).route(point(1):point(1)+1)=new(i).route(point(2)-1:point(2));
            new(i).route(point(2)-1:point(2))=temp;
         end
      end
   end 
end




