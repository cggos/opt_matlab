function [x_best,y_best]=antmin(fun,type)  
if type==1
    prompt={'蚂蚁数m';'最大迭代数';'小数点的位数';'信息素重要程度';'信息素挥发系数'};
    name='输入算法各参数';
    defaultanswer={'20','1000','7','0.9','0.9'};
    answer=inputdlg(prompt,name,1,defaultanswer);
    m=str2num(answer{1});
    iterm=str2num(answer{2});
    d=str2num(answer{3});
    alpha=str2num(answer{4});
    rho=str2num(answer{5});
    q0=0.9;
    tau0=0.01;
    tau1(1:10)=tau0;
    NC=1;
    for i=1:10
       for j=1:10
          tau(i,j)=tau0;
       end
    end
    tabu=zeros(m,d);    %初始信息素和记录表
    while NC<=iterm
       for i=1:m   %对每只蚂蚁
           for j=1:d  %对每层城市
               if j==1
                   for k=1:10
                       p(k)=tau1(k)/sum(tau1);
                   end    %起始城市与下个城市的概率
               else
                   for k=1:10    %中间层的城市概率
                       if tabu(i,j-1)==0
                          p(k)=tau(10,k)/sum(tau(10,:));
                       else
                          p(k)=tau(tabu(i,j-1),k)/sum(tau(tabu(i,j-1),:));
                       end
                   end
               end
               q=rand;
               if q<=q0   %路径选择规则之一
                  if j==1   %选路径，并更新
                     [aa,idex]=max(tau1);
                     tau1(idex)=(1-rho)*tau1(idex)+rho*tau0; %信息素
                  else
                     if tabu(i,j-1)==0
                        [aa,idex]=max(tau(10,:));
                        tau(10,idex)=(1-rho)*tau(10,idex)+rho*tau0;
                     else
                        [aa,idex]=max(tau(tabu(i,j-1),:));
                        tau(tabu(i,j-1),idex)=(1-rho)*tau(tabu(i,j-1),idex)+rho*tau0;
                     end
                  end
                  if idex==10
                      tabu(i,j)=0;
                  else
                      tabu(i,j)=idex;
                  end
               else   %路径选择规则之二
                    Pcum=cumsum(p); 
                    select=find(Pcum>=rand);
                    if isempty(select)
                         tabu(i,j)=round(1+9*rand);
                    else
                        tabu(i,j)=select(1);
                    end
                    if j==1
                       tau1(tabu(i,j))=(1-rho)*tau1(tabu(i,j))+rho*tau0; 
                    else
                        if tabu(i,j-1)==0
                            tau(10,tabu(i,j))=(1-rho)*tau(10,tabu(i,j))+rho*tau0;
                        else
                            tau(tabu(i,j-1),tabu(i,j))=(1-rho)*tau(tabu(i,j-1),tabu(i,j))+rho*tau0;
                        end
                    end
                    if tabu(i,j)==10
                         tabu(i,j)=0;
                    end
               end
           end
       end
       format long   %　计算函数值和最小值
       for i=1:m; 
           x(i)=0;
           for k=1:d
              x(i)=x(i)+tabu(i,k)*10^(-k);
           end
           y(i)=fun(x(i));
       end
       [y_min(NC),y_min_idex]=min(y);
        y_min(NC)=y_min(NC);
       x_min(NC)=0;
       for k=1:d
           x_min(NC)=x_min(NC)+tabu(y_min_idex,k)*10^(-k);
       end   %路径解码
       for k=1:d-1    %信息素全局更新
            if tabu(y_min_idex,k)==0
                tabu(y_min_idex,k)=10;
            end
            if tabu(y_min_idex,k+1)==0
               tabu(y_min_idex,k+1)=10;
            end
            tau(tabu(y_min_idex,k),tabu(y_min_idex,k+1))=(1-alpha)*tau(tabu(y_min_idex,k),tabu(y_min_idex,k+1))+alpha*y_min(NC)^(-1);
       end
       NC=NC+1;
    end
    [y_best,y_bestidex]=min(y_min);
    x_best=x_min(y_bestidex);  %输出结果
    x_best=vpa(x_best,d);
    y_best=vpa(y_best,d);
elseif type==2
    format short
    prompt={'蚂蚁数m';'最大迭代数';'信息素挥发系数';'自变量下界';'自变量上界'};
    name='输入算法各参数';
    answer=inputdlg(prompt,name,1);
    if isempty(answer{1})
        m=20;
    else
        m=str2num(answer{1});
    end
    if isempty(answer{2})
        iterm=200;
    else
        iterm=str2num(answer{2});
    end
    if isempty(answer{3})
        Rou=0.2;
    else
        Rou=str2num(answer{3});
    end
    LB=str2num(answer{4});
    UB=str2num(answer{5});
    n=size(LB,1);     %维数
    for i=1:m
       X(i,:)=LB'+(UB-LB)'.*rand(1,n); 
       Tau(i)=fun(X(i,:));
    end
    for k=1:iterm
        lamda=1/k;
        [Tau_Best,BestIndex]=max(Tau);
        for i=1:m
           P=(Tau(BestIndex)-Tau(i))/Tau(BestIndex);  
           if P<rand  
               temp=X(i,:)+(2*rand(1,n)-1)*lamda;  
           else  
               temp=X(i,:)+(UB-LB)'.*(rand-0.5);
           end
           temp=boundtest(temp,LB,UB);    
           if fun(temp)>fun(X(i,:))  
               X(i,:)=temp;
           end
           Tau(i)=(1-Rou)*Tau(i)+fun(X(i,:));  
        end
    end
    [max_value,max_index]=max(Tau);
    x_best=X(max_index,:);
    y_best=fun(X(max_index,:));
end


