function [route,f]=hopfieldTSP(city,iter_max)
global A D
[N,distance]=city2d(city);
A=200;
D=100;
U0=0.01;
step=0.0001;
delta=2*rand(N,N)-1;
U=U0*log(N-1)+delta;
V=(1+tansig(U/U0))/2;
E=zeros(1,iter_max);
for k=1:iter_max
    dU=diff_u(V,distance);
    U=U+dU*step;
    V=(1+tanh(U/U0))/2;
    e=energy(V,distance);
    E(k)=e;
end
[r,c]=size(V);
V1=zeros(r,c);
[V_max,V_ind]=max(V);
for j=1:c
    V1(V_ind(j),j)=1;
end
C=sum(V1);
R=sum(V1');
flag=isequal(C,ones(1,N))&isequal(R,ones(1,N));
if flag==1
    [V1_max,V1_ind]=max(V1);
    city_end=city(V1_ind,:);
    Length_end=citydist(city_end(1,:),city_end(end,:)');
    for i=2:size(city_end,1)
        Length_end=Length_end+citydist(city_end(i-1,:),city_end(i,:)');
    end
    for i=1:N
        route(i)=find(V1(:,i)==1);
    end
    figure(1)
    TSPplot(city,route);
    f=value(route,distance);
    figure(2)
    plot(1:iter_max,E)
    title(['能量函数变化曲线(最优能量):' num2str(E(end)) ')']);
    xlabel('迭代次数');
    ylabel('能量函数');
else
    disp('寻优路径无效');
end



function du=diff_u(V,d)
global A D
n=size(V,1);
sum_x=repmat(sum(V,2)-1,1,n);
sum_i=repmat(sum(V,1)-1,n,1);
V_temp=V(:,2:n);
V_temp=[V_temp V(:,1)];
sum_d=d*V_temp;
du=-A*sum_x-A*sum_i-D*sum_d;

function E=energy(V,d)
global A D
n=size(V,1);
sum_x=sumsqr(sum(V,2)-1);
sum_i=sumsqr(sum(V,1)-1);
V_temp=V(:,2:n);
V_temp=[V_temp V(:,1)];
sum_d=d*V_temp;
sum_d=sum(sum(V.*sum_d));
E=0.5*(A*sum_x+A*sum_i+D*sum_d);

function d=citydist(city1,city2)
N=size(city1,1);
d=zeros(N,N);
for i=1:N
    for j=1:N
        if i~=j
           d(i,j)=norm(city1(i,:)'-city2(:,j));
        end
    end
end


