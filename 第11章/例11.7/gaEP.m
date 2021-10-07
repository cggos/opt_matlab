function [minx,minf]=gaEP(fun,popsize,iterm_max,LB,UB,type)          %进化规划算法
nvar=size(LB,1);   %变量的个数
tao1=sqrt(2*nvar)^-1;
tao2=sqrt(2*sqrt(nvar))^-1;
q=round(0.9*popsize);
chome=zeros(popsize,nvar);
for i=1:popsize
    for j=1:nvar
       chome(:,j)=unifrnd(LB(j),UB(j),popsize,1);       %初始化变量
       chomeLamda(j)=(UB(j)-LB(j))/2;
       chomeVar(j)=var(chome(:,j));
    end
    chomeV(i,1)=fun(chome(i,:));
    chomeSigma(i,1)=sqrt(chomeV(i,1));
end
[a,b]=sort(chomeV);
chome=chome(b,:);%已经排序
chomeV=chomeV(b);
chomeSigma=chomeSigma(b);
minx=chome(1,:);
minf=chomeV(1);
for i=1:iterm_max
        if type==1    %标准 
           for j=1:popsize
              for k=1:nvar
                 newchome(j,k)=chome(j,k)+normrnd(0,chomeSigma(k),1,1);
                 newchome(j,k)=boundtest(newchome(j,k),LB(k),UB(k));
              end
           end
        elseif type==2     %自适应标准
             for j=1:popsize
                a=randn;
                for k=1:nvar
                   newchome(j,k)=chome(j,k)+randn*chomeVar(k);
                   newchome(j,k)=boundtest(newchome(j,k),LB(k),UB(k));
                   chomeVar(k)=chomeVar(k)*exp(tao1*a+tao2*randn);
                end
             end
        elseif type==3     %单点变异
             for j=1:popsize
                 b1=ceil(nvar*rand);
                 if chomeLamda(b1)<1e-4;
                    chomeLamda(b1)=(UB(b1)-LB(b1))/2;
                 end
                 newchome(j,b1)=chome(j,b1)+chomeLamda(b1)*randn;
                 newchome(j,b1)=boundtest(newchome(j,b1),LB(b1),UB(b1));
                 chomeLamda(b1)=chomeLamda(b1)*exp(-1.01);
             end
        end
        for j=1:popsize
           newchomeV(j,1)=fun(newchome(j,:));
        end
        chome=EP_select1(chome,newchome,chomeV,newchomeV,q);%选择子代
        for j=1:popsize
           chomeV(j,1)=fun(chome(j,:)); 
           chomeSigma(j,1)=sqrt(chomeV(j,1));
        end
        [a,b]=sort(chomeV);
        chome=chome(b,:);%已经排序
        chomeV=chomeV(b);
        chomeSigma=chomeSigma(b);
        if minf>chomeV(1)
           minx=chome(1,:);
           minf=chomeV(1);
        end
end

function chome=EP_select1(old,new,oldF,newF,q)
num=size(old,1);
total_chome=[old;new];
total_F=[oldF;newF];
competitor=randperm(2*num);
for i=1:2*num
    score=0;
    for j=1:q
        if total_F(i)<total_F(competitor(j))
           score=score+1;
        end
    end
    temp1(i)=score;
end
[a,b]=sort(temp1,'descend');
total_chome=total_chome(b,:);
chome=total_chome(1:num,:);

function y=boundtest(x,LB,UB)    %边界检测
if x>=UB
    y=LB+(x-UB);
elseif x<=LB
    y=UB-(LB-x);
else
    y=x;
end


