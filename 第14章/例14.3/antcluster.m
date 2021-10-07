function [pattern_best,d_best]=antcluster(data,type)   %蚂蚁聚类
x=guiyi(data);   %归一化
num=size(x,1);
if type==1     %已知类别
    prompt={'蚂蚁数m';'最大迭代数';'信息素重要程度';'信息素挥发系数';'类别数'};
    name='输入算法各参数';
    answer=inputdlg(prompt,name,1);
    if isempty(answer{1})
        m=30;
    else
        m=str2num(answer{1});
    end
    if isempty(answer{2})
        iterm=2000;
    else
        iterm=str2num(answer{2});
    end
    if isempty(answer{3})
        alpha=0.9;
    else
        alpha=str2num(answer{3});
    end
    if isempty(answer{4})
        rho=0.9;
    else
        rho=str2num(answer{4});
    end
    c=str2num(answer{5});    %类别数
    q0=0.8;tau0=0.01;
    tau1(1:c)=tau0;
    tau2=ones(c,c)*tau0;%路径信息素
    tau3=ones(num,c)*tau0;%样品归属信息素
    tabu=zeros(m,num);%存储并记录路径的生成
    NC=1;%迭代计数器
    while NC<=iterm%停止条件之一：达到最大迭代次数
       for i=1:m
          for j=1:num
             if j==1
                for k=1:c
                    p(k)=(tau1(k)*tau3(j,k))/sum(tau1.*tau3(j,:));
                end
             else
                for k=1:c
                    p(k)=(tau2(tabu(i,j-1),k)*tau3(j,k))/sum((tau2(tabu(i,j-1),:).*tau3(j,:)));
                end
             end
             q=rand;
             if q<=q0
                if j==1
                   [aa,idex]=max(tau1.*tau3(j,:));
                   tau1(idex)=(1-rho)*tau1(idex)+rho*tau0;
                   %tau3(1,idex)=(1-alpha)*tau3(1,idex)+rho*tau0; 
                else
                   [aa,idex]=max(tau2(tabu(i,j-1),:).*tau3(j,:));
                   tau2(tabu(i,j-1),idex)=(1-rho)*tau2(tabu(i,j-1),idex)+rho*tau0;
                   %tau3(j,idex)=(1-alpha)*tau3(j,idex)+rho*tau0; 
                end
                tabu(i,j)=idex;
             else
                Pcum=cumsum(p);
                select=find(Pcum>=rand);
                if isempty(select)
                   tabu(i,j)=ceil(c*rand); 
                else
                   tabu(i,j)=select(1);   
                end 
                if j==1
                   tau1(tabu(i,j))=(1-rho)*tau1(tabu(i,j))+rho*tau0; 
                   %tau3(1,tabu(i,j))=(1-alpha)*tau3(1,tabu(i,j))+rho*tau0; 
                else
                   tau2(tabu(i,j-1),tabu(i,j))=(1-rho)*tau2(tabu(i,j-1),tabu(i,j))+rho*tau0;
                   %tau3(j,tabu(i,j))=(1-alpha)*tau3(j,tabu(i,j))+rho*tau0; 
                end
             end
          end
       end
       pattern=zeros(m,num);
       for i=1:m
           b1=0;k1=0;      %未知类别数时，找出每只蚂蚁走过的路径，并求出类别数，即按类别排列的模式
           for k=1:num
               b=find(tabu(i,:)==k);
               if isempty(b)
                  b1=b1+1;
                  continue
               else
                  k1=k1+1;
                  for kk=1:length(b)
                     pattern(i,b(kk))=k-b1;
                  end
                end
           end
           if length(unique(pattern(i,:)))==c
              y=cluster_center(x,pattern(i,:));           %求聚类中心
              d(i)=cluster_dis(x,tabu(i,:),y);              %求每只蚂蚁的函数值
           else
               d(i)=inf;
           end
       end
       [d_min(NC),d_min_idex]=min(d);
           %要记录此时的路径
       tabu_min(NC,:)=tabu(d_min_idex,:);
       for k=1:num-1
           tau2(tabu(d_min_idex,k),tabu(d_min_idex,k+1))=(1-alpha)*tau2(tabu(d_min_idex,k),tabu(d_min_idex,k+1))+alpha*d_min(NC)^(-1);     
       end
       for k=1:num
           tau3(k,tabu(d_min_idex,k))=(1-alpha)*tau3(k,tabu(d_min_idex,k))+alpha*d_min(NC)^(-1); 
       end
       NC=NC+1;
    end
    [d_best,d_bestidex]=min(d_min);
    pattern_best=tabu_min(d_bestidex,:);
elseif type==2
    prompt={'最大迭代数';'信息素重要程度';'能见度';'信息素挥发系数'};
    name='输入算法各参数';
    defaultanswer={'1000','0.8','1','0.2'};
    answer=inputdlg(prompt,name,1,defaultanswer);
    iterm=str2num(answer{1});
    alpha=str2num(answer{2});
    rho=str2num(answer{3});
    beta=str2num(answer{4});
    q0=0.8;
    d=squareform(pdist(x));
    for i=1:num-1
       temp1(i)=min(d(i,i+1:end));
       temp2(i)=max(d(i,i+1:end));
    end
    m_min=min(temp1);m_max=max(temp2);
    r=m_min+(m_max-m_min)/2;
    tau=zeros(num,num);
    h=ones(num,num);
    p=zeros(num,num);
    centernum=num;
    temp=zeros(1,centernum);
    m_pattern=1:num;
    m_center=cluster_center(x,m_pattern);
    NC=1;
    while NC<=iterm||centernum==1%停止条件之一：达到最大迭代次数
       for i=1:centernum-1
          temp(i)=0;
          for j=i+1:centernum     
             y(i,j)=center_dis(m_center(i:j,:));
             if y(i,j)<r
                tau(i,j)=1;
             else
               tau(i,j)=0;
             end
             h(i,j)=1;
             temp(i)=temp(i)+(tau(i,j)^alpha)*(h(i,j)^beta);
          end
       end
       flag=1;
       for i=1:centernum-1
          for j=i+1:centernum
              p(i,j)=(tau(i,j)^alpha)*h(i,j)^beta/temp(i);
              if p(i,j)>q0
                 for k=1:num
                    if m_pattern(k)==j
                       m_pattern(k)=i;
                   elseif m_pattern(k)>j
                       m_pattern(k)=m_pattern(k)-1;
                    end
                 end
                 b1=0;k1=0;    %同上
                 for k=1:num
                    b=find(m_pattern==k);
                    if isempty(b)
                       b1=b1+1;
                       continue
                    else
                       k1=k1+1;
                       for kk=1:length(b)
                         m_pattern(b(kk))=k-b1;
                       end
                    end
                 end
                 flag=0;
                 centernum=centernum-1;
                 break
              end
          end
          if flag==0
             break
          end
      end
      if flag==1
          break
      end 
      for i=1:centernum
         m_center=cluster_center(x,m_pattern);
      end  
    end
    pattern_best=m_pattern;
    y=cluster_center(x,pattern_best);
    d_best=cluster_dis(x,pattern_best,y); 
end
    
    