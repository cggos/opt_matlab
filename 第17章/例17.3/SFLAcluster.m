function [class,fval]=SFLAcluster(fun,frogNum,m,max_iterm,L,LB,UB,data) %��������㷨,classΪ���
[n,c]=size(LB);     %nά��m��
for i=1:frogNum
    frog(i).center=LB'+(UB'-LB').*rand(c,n);
    y(i)=fun(frog(i).center);
end
[a,b]=sort(y,'ascend');  %��Ӧ�Ⱥ���������
frog=frog(b);
y=y(b);
fval=y(1);best_x=frog(1).center;  %ȫ�����Ž⣬����Сֵ
for iter=2:max_iterm
    [new,fitness]=grouping1(frog,y,m);   %����
    [new,fitness,best_x,fval]=frog_alter_C(fun,new,fitness,m,L,best_x,fval,LB,UB);%����
    old=[];
    y=[];
    for i=1:m
        y=[y fitness(i).y];
        for j=1:length(new(i).center)
           old=[old;new(i).center{j}];   %����
        end 
    end
    for i=1:frogNum
        frog(i).center=old(c*(i-1)+1:c*i,:);
    end
    [a,b]=sort(y,'ascend');  %��Ӧ�Ⱥ���������
    frog=frog(b);
    y=y(b);
end
class=myclass(data,best_x);

function [new,fitness]=grouping1(old,y,m)   %����,���ھ���
N=length(old);
for i=1:m    %����
    j=1;
    while 1
        if m*(j-1)+i>N
            break
        else
           new(i).center{j}=old(m*(j-1)+i).center;   
           fitness(i).y(j)=y(m*(j-1)+i);
           j=j+1;
        end
    end
end


function [new,fitness,best_x,fval]=frog_alter_C(fun,new,fitness,m,L,best_x,fval,LB,UB)
[NC,c]=size(LB);
for i=1:m
  for j=1:L
    [y_best,b1]=min(fitness(i).y);
     x_best=new(i).center{b1}; %��������
    [y_worst,b2]=max(fitness(i).y);
     x_worst=new(i).center{b2}; %�������
     x_temp=x_worst+rand(c,NC).*(x_best-x_worst);
     y_temp=fun(x_temp);
    if y_temp<y_worst
       new(i).center{b2}=x_temp;
       fitness(i).y(b2)=y_temp;
       if y_temp<fval
          fval=y_temp;
          best_x=x_temp;
       end
    else
       x_best=best_x;
       x_temp=x_worst+rand(c,NC).*(x_best-x_worst);
       y_temp=fun(x_temp);
       if y_temp<y_worst
          new(i).center{b2}=x_temp;
          fitness(i).y(b2)=y_temp;
          if y_temp<fval
             fval=y_temp;
             best_x=x_temp;
          end        
       else
           temp1=LB'+(UB'-LB').*rand(c,NC);
           new(i).center{b2}=temp1;
           fitness(i).y(b2)=fun(new(i).center{b2});
           if fitness(i).y(b2)<fval
              fval=fitness(i).y(b2);
              best_x=new(i).center{b2};
           end
       end
    end
  end
end
















 
 
 