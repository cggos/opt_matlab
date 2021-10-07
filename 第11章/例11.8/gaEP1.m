function [maxx,maxf]=gaEP1(fun,popsize,iterm_max)          %进化规划算法求背包问题
a=[198 30 167 130 35 20 105 196 94 126];
b=546;
nvar=length(a);
chome=zeros(popsize,nvar);
q=round(0.9*popsize);
for i=1:popsize
   chome(i,:)=rand(1,nvar)<0.5;       %初始化变量
   chomeV(i,1)=fun(chome(i,:));
end
[a1,b1]=sort(chomeV,'descend');
chome=chome(b1,:);%已经排序
chomeV=chomeV(b1);
maxx=chome(1,:);
maxf=chomeV(1);
for i=1:iterm_max
       for j=1:popsize
           newchome(j,:)=chome(j,:);
           b1=ceil(nvar*rand);
           if b>=a(b1)&chome(j,b1)==0
              newchome(j,b1)=1;
           else
               for k=1:nvar
                  s=0;
                  if b-chome(j,k)*a(k)<a(k)
                     s=s+1;
                  end
               end
               if s==nvar
                  b3=ceil(nvar*rand);
                  newchome(j,b3)=0;
               end
           end
           newchomeV(j,1)=fun(newchome(j,:));
       end
       chome=EP_select2(chome,newchome,chomeV,newchomeV,q);%选择子代
       for j=1:popsize
           chomeV(j,1)=fun(chome(j,:)); 
       end
       [a,b]=sort(chomeV,'descend');
       chome=chome(b,:);%已经排序
       chomeV=chomeV(b);
       if maxf<chomeV(1)
          maxx=chome(1,:);
          maxf=chomeV(1);
       end
end

function chome=EP_select2(old,new,oldF,newF,q)
num=size(old,1);
total_chome=[old;new];
total_F=[oldF;newF];
competitor=randperm(2*num);
for i=1:2*num
    score=0;
    for j=1:q
        if total_F(i)>total_F(competitor(j))
           score=score+1;
        end
    end
    temp1(i)=score;
end
[a,b]=sort(temp1,'descend');
total_chome=total_chome(b,:);
chome=total_chome(1:num,:);

