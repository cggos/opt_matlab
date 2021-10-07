function [new1,new2,d1,d2]=greedycross(old1,old2,d)    %贪心交叉求路径
citynum=size(d,1);
d1=value(old1,d);
d2=value(old2,d);
point=ceil(citynum*rand);
new1=zeros(1,citynum);
new2=zeros(1,citynum);
new1(1)=old1(point);    %向左
new2(1)=old1(point);    %向右
visit=old1(point);
point1=point;
for i=2:citynum
   while 1
        if point+1>citynum&&find(old2==old1(point))==citynum
            temp1=old1(1);
            temp2=old2(1);
        elseif point+1>citynum&&find(old2==old1(point))~=citynum
            temp1=old1(1);
            temp2=old2(find(old2==old1(point))+1);
        elseif find(old2==old1(point))==citynum
            temp1=old1(point+1);
            temp2=old2(1);
        else
            temp1=old1(point+1);
            temp2=old2(find(old2==old1(point))+1);
        end
        if isempty(mycompare1(visit,temp1))&&isempty(mycompare1(visit,temp2))
            if d(new1(i-1),temp1)<d(new1(i-1),temp2)
                new1(i)=temp1;
                visit=[visit new1(i)];
            else
                new1(i)=temp2;
                visit=[visit new1(i)];
            end
        elseif ~isempty(mycompare1(visit,temp1))||~isempty(mycompare1(visit,temp2))
              temp=d(new1(i-1),:);
              temp(visit)=inf;
              [a,new1(i)]=min(temp);    %与一个城市最近的城市
        end
        break
   end
   point=find(old1==new1(i));
end
point=point1;
visit=old1(point);
for i=2:citynum
    while 1
         if point-1<1&&find(old2==old1(point))==1
             temp1=old1(end);
             temp2=old2(end);
         elseif point-1<1&&find(old2==old1(point))~=1
             temp1=old1(end);
             temp2=old2(find(old2==old1(point))-1);
         elseif find(old2==old1(point))==1
             temp1=old1(point-1);
             temp2=old2(end);
         else
             temp1=old1(point-1);
             temp2=old2(find(old2==old1(point))-1);
         end
         if isempty(mycompare1(visit,temp1))&&isempty(mycompare1(visit,temp2))
             if d(new2(i-1),temp1)<d(new2(i-1),temp2)
                new2(i)=temp1;
                visit=[visit new2(i)];
             else
                new2(i)=temp2;
                visit=[visit new2(i)];
             end
         elseif ~isempty(mycompare1(visit,temp1))||~isempty(mycompare1(visit,temp2))
             temp=d(new2(i-1),:);
             temp(visit)=inf;
             [a,new2(i)]=min(temp);    %与一个城市最近的城市
         end
         break
    end
    point=find(old1==new2(i));
end
d3=value(new1,d);
d4=value(new2,d);
dd=[d1 d2 d3 d4];
new=[old1;old2;new1;new2];
[a,b]=sort(dd,'ascend');
new1=new(b(1),:);
new1=isin_TSP(new1,citynum);
new2=new(b(2),:);
new2=isin_TSP(new2,citynum);
d1=dd(b(1));
d2=dd(b(2));
