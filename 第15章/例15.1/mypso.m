function [bestx,bestf]=mypso(fun,type)    %粒子群算法求极值
if type==1   %连续问题解
   prompt={'粒子数';'最大迭代数';'变量下界';'变量上界'};
   name='输入算法各参数';
   defaultanswer={'30','200','-inf','inf'};
   answer=inputdlg(prompt,name,1,defaultanswer);
   popsize=str2num(answer{1});
   max_iterm=str2num(answer{2});
   LB=str2num(answer{3});
   UB=str2num(answer{4});
  [n,m]=size(LB);    %确定维数及组数(m)
  for i=1:popsize
     pop(i).x=LB'+(UB-LB)'.*rand(m,n);
     v(i).x=rand(m,n);
     fitness(i)=fun(pop(i).x);
  end
  vmax=1;vmin=-1;
  [best_y best_index]=min(fitness);
  z_best=pop(best_index).x;    %全局极值位置
  g_best=pop;                  %个体极值位置
  y_g_best=fitness;            %个体极值
  y_z_best=best_y;             %全局极值
  w_max=0.9;w_min=0.4;
  for i=1:max_iterm
    w=w_max-i*(w_max-w_min)/max_iterm;
    for j=1:popsize
        v(j).x=w.*v(j).x+2.05*rand.*(g_best(j).x-pop(j).x)+2.05*rand.*(z_best-pop(j).x);
        v(j).x(find(v(j).x>vmax))=vmax;
        v(j).x(find(v(j).x<vmin))=vmin;
        pop(j).x=pop(j).x+v(j).x;
        pop(j).x=boundtest(pop(j).x,LB,UB);
        if rand>0.9
            k=ceil(n*rand);
            for k1=1:m
               pop(j).x(k1,k)=LB(k,k1)+(UB(k,k1)-LB(k,k1))*rand;
            end
        end
        y(j)=fun(pop(j).x);
    end 
    if y(j)<y_g_best(j)
        g_best(j)=pop(j);
        y_g_best(j)=y(j);
    end
    if y(j)<y_z_best
        z_best=pop(j).x;
        y_z_best=y(j);
    end
  end
  bestx=z_best;bestf=y_z_best;
elseif type==2      %离散问题
   prompt={'粒子数';'最大迭代数';'数据集维数'};
   name='输入算法各参数';
   defaultanswer={'30','200','[]','[]'};
   answer=inputdlg(prompt,name,1,defaultanswer);
   popsize=str2num(answer{1});
   max_iterm=str2num(answer{2});
   dim=str2num(answer{3});
   for i=1:popsize
      pop(i,:)=rand(1,dim)<0.5;
      v(i,:)=rand(1,dim);
      fitness(i)=fun(pop(i,:));
   end
   [best_y best_index]=min(fitness);
   g_best=pop(best_index,:);    %全局极值位置
   p_best=pop;                  %个体极值位置
   pbest=fitness;            %个体极值
   gbest=best_y;             %全局极值
   alpha=0.3;beta=0.7;
   c1=0.4;c2=0.3;c3=0.3;
   for i=1:max_iterm
       for j=1:popsize
          Pbest=alpha*pbest(j)+beta*(1-pbest(j));
          Gbest=alpha*gbest+beta*(1-gbest);
          v1=c1.*v(j,:)+c2*Pbest+c3*Gbest;
          v1=1./(1+exp(-v1));
          for k=1:dim
              if v1(k)<rand
                  pop(j,k)=1;
              else
                  pop(j,k)=0;
              end
          end
          y(j)=fun(pop(j,:));
          if y(j)<pbest(j)
              p_best(j,:)=pop(j,:);
              pbest(j)=y(j);
          end
          if y(j)<gbest
             g_best=pop(j,:);
             gbest=y(j);
          end
       end
   end
   bestx=g_best;bestf=gbest;   
end