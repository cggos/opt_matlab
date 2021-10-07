function val_f=gaDE1(fun,popsize,iterm_max,LB,UB,type)    %差分进化算法求积分
%popsize为种群大小,NC为维数,pc为交叉概率,F为放大因子
if nargin==5
    type=[];
end
NC=size(LB,1);
chome=myinitialize(popsize,NC,LB,UB);
for i=1:popsize
   [chomeV1(i,1),chomeV2(i,1)]=fun(chome(i,:));
end
[a,index]=sort(chomeV1);   %求极小值
val_f=chomeV2(index(1),1);
chome=chome(index,:);
chomeV1=chomeV1(index,1);
for i=1:iterm_max
    alpha=(iterm_max-i)/iterm_max;
    beta=1-alpha;
    pc=0.1+0.8*(i/iterm_max)^2;
    F=0.5*(exp(iterm_max/(iterm_max+i))-1);
    if isempty(type)
        m=ceil(rand*7);
        new=mutation_DE(chome,F,index(1),alpha,beta,LB,UB,m);
    else
        new=mutation_DE(chome,F,index(1),alpha,beta,LB,UB,type);
    end
    new=crossover_DE(new,chome,index(1),pc);
    for j=1:popsize
        newV(j,1)=fun(chome(j,:));
    end
    chome=selection_DE(new,chome,chomeV1,newV);
    for j=1:popsize
       [chomeV1(j,1),chomeV2(j,1)]=fun(chome(j,:));
    end
    [a,index]=sort(chomeV1);
    chome(index(end),:)=chome(index(1),:);
    if val_f>chomeV1(index(1),1)
       val_f=chomeV2(index(1),1);
    end
   chome=chome(index,:);
   chomeV1=chomeV1(index,1);
end


function old=myinitialize(popsize,NC,LB,UB)    %初始化
old=zeros(popsize,NC);
for i=1:NC
  old(:,i)=LB(i)+(UB(i)-LB(i)).*rand(popsize,1);
end
old=sort(old,2);



function new=mutation_DE(chome,F,index,alpha,beta,LB,UB,type)   %差分
[r,c]=size(chome);
if type==7
    M=round(r/2);
else
    M=r;
end
for i=1:M
  k=randperm(r);
  k=redu(k,find(i==k),'c');     %扣除当前个体
  switch type
      case 1   %随机向量差分法
         new(i,:)=chome(i,:)+F.*(chome(k(1),:)-chome(k(2),:));
      case 2   %最优解
         new(i,:)=chome(index,:)+F.*(chome(k(1),:)-chome(k(2),:));
      case 3  %最优解加多个个体
         new(i,:)=chome(index,:)+F.*(chome(k(1),:)+chome(k(2),:)-chome(k(3),:)-chome(k(4),:));
      case 4
         new(i,:)=chome(i,:)+F.*(chome(index,:)-chome(k(1),:)+chome(k(2),:)-chome(k(3),:));
      case 5
         new(i,:)=alpha.*chome(k(1),:)+beta.*chome(index,:)+F.*(chome(k(2),:)-chome(k(3),:));
      case 6
          if rand<0.99
             if rand>0.5
                new(i,:)=chome(index,:)+rand.*(chome(k(1),:)-chome(k(2),:)); 
             else
                new(i,:)=chome(k(1),:)+rand.*(chome(k(2),:)-chome(k(3),:)); 
             end
          else
             new(i,:)=LB'+(UB'-LB').*rand(1,c);
          end
      case 7   
         if i+round(r/2)>r
             berak
         end
         new(i,:)=chome(index,:)+F.*(chome(k(1),:)-chome(k(2),:));     
         new(i+round(r/2),:)=(chome(k(3),:)+chome(k(4),:))/2;
  end
  if type==7
     for j=1:c
        new(i,j)=boundtest(new(i,j),LB(j),UB(j));
        new(i+round(r/2),j)=boundtest(new(i+round(r/2),j),LB(j),UB(j));
     end
  else
     for j=1:c
        new(i,j)=boundtest(new(i,j),LB(j),UB(j));
     end
  end
end
new=sort(new,2);

function new=crossover_DE(new,old,index,pc)
[r,c]=size(old);
for i=1:r
   for j=1:c
      if rand<pc 
         new(i,j)=old(i,j);
      else
         new(i,j)=old(i,j)+rand*(old(index,j)-old(i,j));
      end
   end
end
new=sort(new,2);


function  old=selection_DE(new,old,fitness,newV) %已排序
r=size(new,1);
for i=1:r
    if fitness(i,1)>newV(i,1)   
        old(i,:)=new(i,:);
    end
end
old=sort(old,2);