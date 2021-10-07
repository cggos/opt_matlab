function [val_x,val_f]=gaES(fun,popsize,iterm_max,LB,UB)   %进化策略算法
nvar=size(LB,1);   %变量的个数
chome=zeros(popsize,nvar);
chomesigma=zeros(popsize,nvar);
lamda=7*popsize;
for i=1:popsize
    for j=1:nvar
       chome(:,j)=unifrnd(LB(j),UB(j),popsize,1);       %初始化变量
       chomesigma(:,j)=3.0.*ones(popsize,1);
    end
    chomeV(i,1)=fun(chome(i,:));
end
[a,b]=sort(chomeV);
val_x=chome(b(1),:);
val_f=chomeV(b(1),1);
for i=1:iterm_max
    [newchome,newsigma]=recombination1(chome,chomesigma,lamda,1);
    [newchome,newsigma]=ES_mutation(newchome,newsigma,LB,UB);
    for j=1:popsize
        newchomeV(j,1)=fun(newchome(j,:));
    end
    [chome,chomesigma]=ES_select1(chome,newchome,chomeV,newchomeV,chomesigma,newsigma,1);
    for j=1:popsize
        chomeV(j,1)=fun(chome(j,:));
    end
    [a,b]=sort(chomeV);
    if val_f>chomeV(b(1),1)
       val_x=chome(b(1),:);
       val_f=chomeV(b(1),1);
    end
end

function [new,newsigma]=recombination1(oldchome,oldsigma,lamda,type)
[num,nvar]=size(oldchome);
temp=ceil(num*rand);   %选择一个父代
for i=1:lamda
    temp1=ceil(num.*rand(1,2));
    switch type
       case 1   %离散重组          
           mask=round(rand(1,nvar));
           if any(mask-1)||any(mask)
               mask=round(rand(1,nvar));
           end
            for j=1:nvar  %变量数目
               if mask(j)==0
                   new(i,j)=oldchome(temp1(1),j);
                   newsigma(i,j)=oldsigma(temp1(1),j);
               else
                   new(i,j)=oldchome(temp1(2),j);
                   newsigma(i,j)=oldsigma(temp1(2),j);
               end
            end 
    case 2   %中间重组
        for j=1:nvar
            new(i,j)=(oldchome(temp1(1),j)+oldchome(temp1(2),j))/2;
            newsigma(i,j)=(oldsigma(temp1(1),j)+oldsigma(temp1(2),j))/2;
        end    
    case 3   %混杂重组
       temp2=ceil(num*rand);
       while temp2==temp
           temp2=ceil(num*rand);
       end
       for j=1:nvar
           new(i,j)=(oldchome(temp,j)+oldchome(temp2,j))/2;
           newsigma(i,j)=(oldsigma(temp,j)+oldsigma(temp2,j))/2;
       end     
    end      
end

function [new,newsigma]=ES_mutation(oldchome,oldsigma,LB,UB)
[popsize,nvar]=size(oldchome);
for j=1:popsize
     a=randn;
     for k=1:nvar
        newsigma(j,k)=oldsigma(j,k)*exp(a+randn);
        new(j,k)=oldchome(j,k)+normrnd(0,newsigma(j,k),1,1);
        new(j,k)=boundtest(new(j,k),LB(k),UB(k));
     end 
end

function [chome,newsigma]=ES_select1(old,new,oldF,newF,oldsigma,newsigma,type)
num=size(old,1);
total_chome=[old;new];
total_F=[oldF;newF];
total_sigma=[oldsigma;newsigma];
if type==1   %(N+λ）方式
   [a,b]=sort(total_F);
   total_chome=total_chome(b,:);
   total_sigma=total_sigma(b,:);
   chome=total_chome(1:num,:);
   newsigma=total_sigma(1:num,:);
elseif type==2   %(N,λ）方式
   if size(new,1)<num
       error('新个体的数目太少');
   end
   [a,b]=sort(newF);
   new=new(b,:);
   newsigma=newsigma(b,:);
   chome=new(1:num,:);
   newsigma=newsigma(1:num,:);
end
    
