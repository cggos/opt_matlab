function [bestx,fval]=bitSFLA(fun,frogNum,m,max_iterm,L,dimension) 
%��������㷨����0-1�滮,dimensionΪά�����������Ŀ
for i=1:frogNum
    frog(i,:)=-5+10.*rand(1,dimension);
    c(i,:)=realbit(frog(i,:));   %������
    y(i)=fun(c(i,:));
end
[a,b]=sort(y,'ascend');  %��Ӧ�Ⱥ���������
frog=frog(b,:);y=y(b);c=c(b,:);
fval=y(1);best_x=frog(1,:);  %ȫ�����Ž⣬����Сֵ
bestx=c(1,:);
for iter=1:max_iterm
  [new,fitness]=grouping(frog,y,m);   %����
  [new,fitness,best_x,fval,bestx]=frog_alter2(fun,new,fitness,m,L,best_x,fval,bestx);%����
  old=[];
  y=[];
  for i=1:m
      old=[old;new(i).x];   %����
      y=[y fitness(i).y];
  end
  [a,b]=sort(y,'ascend');  %��Ӧ�Ⱥ���������
  frog=old(b,:);
  y=y(b);
end

function [new,fitness]=grouping(old,y,m)   %����
N=size(old,1);
for i=1:m    %����
    j=1;
    while 1
        if m*(j-1)+i>N
            break
        else
           new(i).x(j,:)=old(m*(j-1)+i,:);   
           fitness(i).y(j)=y(m*(j-1)+i);
           j=j+1;
        end
    end
end

function [old,fitness,best_x,fval,bestx]=frog_alter2(fun,old,fitness,m,L,best_x,fval,bestx)   %����
NC=length(best_x);
for i=1:m     %ÿ��ģ���������
  for j=1:L
     [y_best,b1]=min(fitness(i).y);
     x_best=old(i).x(b1,:); %��������
     [y_worst,b2]=max(fitness(i).y);
     x_worst=old(i).x(b2,:); %�������
     for k=1:NC
        d_temp(k)=rand*(x_best(k)-x_worst(k));
     end
    temp=x_worst+d_temp;%���º��λ��
    temp1=realbit(temp);
    y_temp=fun(temp1);
    if y_temp<y_worst
       old(i).x(b2,:)=temp;
       fitness(i).y(b2)=y_temp;
       if y_temp<fval
           best_x=temp;
           fval=y_temp;
           bestx=temp1;
       end
    else
       x_best=best_x;
       for k=1:NC
          d_temp(k)=0.729*rand*2.05*rand*(x_best(k)-x_worst(k));
       end
       temp=x_worst+d_temp;
       temp1=realbit(temp);
       y_temp=fun(temp1);
       if y_temp<y_worst
           old(i).x(b2,:)=temp;
           fitness(i).y(b2)=y_temp;
           if y_temp<fval
              best_x=temp;
              fval=y_temp;
              bestx=temp1;
           end
       else
           old(i).x(b2,:)=-5+10.*rand(1,NC);
           fitness(i).y(b2)=fun(old(i).x(b2,:));
           if fitness(i).y(b2)<fval
               best_x=old(i).x(b2,:);
               fval=fitness(i).y(b2);
               bestx=realbit(best_x);
           end 
        end
     end
  end
end

function b=realbit(x)   %ʵ��ת���ɶ�����
n=length(x);
b=zeros(1,n);
for i=1:n
    if 1/(1+exp(-x(i)))<=0.5
        b(i)=0;
    else
        b(i)=1;
    end
end



