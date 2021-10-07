function [bestx,bestf]=COApso(fun,type)    %混沌粒子群算法求极值
prompt={'粒子数';'最大迭代数';'混沌迭代次数';'变量下界';'变量上界'};
name='输入算法各参数';
defaultanswer={'30','100','20','-inf','inf'};
answer=inputdlg(prompt,name,1,defaultanswer);
popsize=str2num(answer{1});
max_iterm=str2num(answer{2});
L=str2num(answer{3});
LB=str2num(answer{4});
UB=str2num(answer{5});
n=size(LB,1);    %确定维数
vmax=1;vmin=-1;
for i=1:2*popsize    %混沌初始化
    z=createCOA(L,n,'logistic');
    x1(i,:)=LB'+(UB-LB)'.*z;
    y1(i)=fun(x1(i,:));
end
[a,b]=sort(y1);
pop=x1(b(1:popsize),:);
fitness=y1(b(1:popsize));
v=rand(popsize,n);
[best_y best_index]=min(fitness);
z_best=pop(best_index,:);    %全局极值位置
g_best=pop;                  %个体极值位置
y_g_best=fitness;            %个体极值
y_z_best=best_y;             %全局极值
w_max=0.9;w_min=0.4;
for i=1:max_iterm
    w=w_max-i*(w_max-w_min)/max_iterm;   
    for j=1:popsize
        v(j,:)=w.*v(j,:)+2.05*rand.*(g_best(j,:)-pop(j,:))+2.05*rand.*(z_best-pop(j,:));
        v(j,find(v(j,:)>vmax))=vmax;
        v(j,find(v(j,:)<vmin))=vmin;
        pop(j,:)=pop(j,:)+v(j,:);
        if rand>0.9
            k=ceil(n*rand);
            pop(j,k)=LB(k)+(UB(k)-LB(k))*rand;
        end
        pop(j,:)=boundtest(pop(j,:),LB,UB);
        y(j)=fun(pop(j,:));
        if y(j)<y_g_best(j)
           g_best(j,:)=pop(j,:);
           y_g_best(j)=y(j);
        end
        if y(j)<y_z_best
           z_best=pop(j,:);
           y_z_best=y(j);
        end
    end
    switch type
        case 1       %扰动
            for j=1:popsize
               z=createCOA(L,n,'logistic');
               deltax=-1+2.*z;
               pop1=pop(j,:)+deltax;
               pop1=boundtest(pop1,LB,UB);
               y2=fun(pop1);
               if y2<y(j)
                  pop(j,:)=pop1;
                  y(j)=y2;
               end     
            end
        case 2   %替换粒子
            num=round(popsize*0.2);
            m=randperm(popsize);
           for j=1:num
              for k=1:n
                 z=(pop(m(j),k)-LB(k))/(UB(k)-LB(k));
                 for s=1:L
                    z=4.*z.*(1-z);
                 end
                 pop1(k)=LB(k)+(UB(k)-LB(k))*z;
              end
              y3=fun(pop1);
              if y3<y(m(j))
                 pop(m(j),:)=pop1;
                 y(m(j))=y3;
              end
           end
        case 3     %扰动局部最优
           for j=1:popsize
              for k=1:n
                 z=(g_best(j,k)-LB(k))/(UB(k)-LB(k));
                 for s=1:L
                    z=4.*z.*(1-z);
                 end
                 zg(j,k)=LB(k)+(UB(k)-LB(k))*z;
              end
              fg(j)=fun(zg);
           end
           [a,b]=min(fg);
           aa=ceil(popsize*rand);
           pop(aa,:)=zg(b,:);
           y(aa)=a;
        case 4     %扰动全局最优
          z=createCOA(L,n,'logistic');
          deltax=-1+2.*z;
          g=z_best+deltax;
          g=boundtest(g,LB,UB);
          y1=fun(g);
         if y1<y_z_best
            z_best=g;
            y_z_best=y1;
         end
    end
end
bestx=z_best;bestf=y_z_best;
