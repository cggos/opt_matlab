function [road,f0]=dyprogTSP(d,city)   %用动态规划求解TSP问题,指定开始地
if nargin==1
    city=1;
end
n=size(d,1);
city1=redu(1:n,city,'c');
for i=1:n
   f0(i,1)=d(i,city);
   f(i,:)=i;
end
f0=redu(f0,city,'r');
f=redu(f,city,'r'); 
m=nchoosek(city1,1);
r=size(m,1);
f1=[];
f2=[];
road1=[];
d1=[];
for i=1:r
   y1=mycompare(city1,m(i,:));
   num=length(y1);
   num2=length(m(i,:));
   for k1=1:num
     f2=[f2;y1(k1) m(i,:)];  %状态
     for j=1:num2
        y3=findpos1(f,m(i,j),'r',1);
        d1(k1,j)=d(y1(k1),m(i,j))+f0(y3);
     end
     [fi,id]=min(d1(k1,:));
     road1=[road1;y1(k1) m(i,:)];
     f1=[f1;fi];   %路径长
  end
end
f0=f1;
f=f2;
for k=2:n-1 
    m=nchoosek(city1,k);
    r=size(m,1);
    f1=[];
    f2=[];
    road2=[];
    y3=[];
    d1=[];
    for i=1:r
        if k==n-1
            y1=city;
        else
            y1=mycompare(city1,m(i,:));
        end
        num=length(y1);
        num2=length(m(i,:));
        m1=nchoosek(m(i,:),k-1);
        for k1=1:num
           f2=[f2;y1(k1) m(i,:)];  %状态
           for j=1:num2
                 y2=mycompare(m(i,:),m1(j,:));
                 y3(j,:)=[y2 m1(j,:)];
                 y4=findpos1(f,y3(j,:),'r',1);
                 d1(k1,j)=d(y1(k1),y2)+f0(y4);
           end
           [fi,id]=min(d1(k1,:));
           y5=findpos1(f,y3(id,:),'r',1);
           road2=[road2;y1(k1) road1(y5,:)];    %路径
           f1=[f1;fi];   %路径长
        end   
    end
    f0=f1;
    f=f2;
    road1=road2;
end
road=[road1 city];



        
    
    
