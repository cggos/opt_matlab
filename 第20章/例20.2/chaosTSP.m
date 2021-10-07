function [best_route,fval]=chaosTSP(city,N1,pc,pr,gama)   %�����㷨���TSP
[NC,d]=city2d(city);
M=zeros(NC,NC);
for i=1:NC
    for j=1:NC
        if i==j
          M(i,j)=1;
        end
    end
end
fval=inf;
route=[];
for i=1:NC
   route=[route find(M(:,i)==1)];
end
for j=1:2
    z1(j)=rand;
    while z1(j)==0.25||z1(j)==0.5||z1(j)==0.75
       z1(j)=rand;
    end
    if j==2
        while z1(2)==z1(1)
           z1(2)=rand;
        end
    end
    x1(j)=4*z1(j)*(1-z1(j));
end
for i=1:N1
   mx=ceil(x1.*NC);%�������У��У���
   M1=M;
   if rand<pc     %  �н���
      temp=M(:,mx(1));   %����
      M(:,mx(1))=M(:,mx(2));
      M(:,mx(2))=temp;
   end
   if rand<pr     %�н���
      temp=M(mx(1),:);   %����
      M(mx(1),:)=M(mx(2),:);
      M(mx(2),:)=temp;
   end
   route=[];
   for j=1:NC
      route=[route find(M(:,j)==1)];
   end 
   y=value(route,d);
   if y<fval
      best_route=route;
      fval=y;
   elseif y>gama*fval
      M=M1;
   end
   x1=4*x1.*(1-x1);  %����n1ά�Ļ������
end
TSPplot(city,best_route);
