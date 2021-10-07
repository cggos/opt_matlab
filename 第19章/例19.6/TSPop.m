function pop=TSPop(varargin)     %TSP各种算子
type=length(varargin{end});
if findstr(type,'c')
    pop1=varargin{1};pop2=varargin{2};
    pop=[pop1;pop2];
elseif strcmp(type,'tm')
    pop=varargin{1};
    pm=varargin{2};
elseif strcmp(type,'imm')
    pop=varargin{1};
    d=varargin{2};
else
    pop1=varargin{1};
    pop=pop1;
end
switch type
    case 'cn'    %Non-ABEL法
         pop(1,:)=pop1(pop2);
         pop(2,:)=pop2(pop1);
    case 'cpmx'   %部分匹配
         point=sort(ceil(rand(1,2)*n));
         while point(1)==point(2)||point(2)-point(1)==1
              point=sort(ceil(rand(1,2)*n));
         end
         pop(1,point(1):point(2))=fliplr(pop2(point(1):point(2)));
         pop(2,point(1):point(2))=fliplr(pop1(point(1):point(2)));      
         pop(1,:)=isin_TSP(pop(1,:));
         pop(2,:)=isin_TSP(pop(2,:));
    case 'cox1'   %OX算子
         point=sort(ceil(rand(1,2)*n));
         while point(1)==point(2)||point(1)==1||point(2)==n
             point=sort(ceil(rand(1,2)*n));
         end
         temp1=pop1(point(1):point(2));
         pop(1,point(1):point(2))=temp1;
         temp2=pop2(point(1):point(2));
         pop(2,point(1):point(2))=temp2;
         temp3=[pop1(point(2)+1:end) pop1(1:point(1)-1) temp1];
         [y1,y2]=mycompare1(temp3,temp2);
         temp3=redu(temp3,y2,'c');
         pop(2,1:point(1)-1)=temp3(1:point(1)-1);
         pop(2,point(2)+1:end)=temp3(point(1):end);
         temp3=[pop2(point(2)+1:end) pop2(1:point(1)-1) temp2];
         [y1,y2]=mycompare1(temp3,temp1);
         temp3=redu(temp3,y2,'c');
         pop(1,1:point(1)-1)=temp3(1:point(1)-1);
         pop(1,point(2)+1:end)=temp3(point(1):end);
    case 'cox2'   %类OX算子
         point=sort(ceil(rand(1,2)*n));
         while point(1)==point(2)
             point=sort(ceil(rand(1,2)*n));
         end
         temp1=pop1(point(1):point(2));
         pop(1,point(1):point(2))=temp1;
         temp2=pop2(point(1):point(2));
         pop(2,point(1):point(2))=temp2;
         temp3=[temp2 pop1];
         pop(1,:)=isin_TSP(temp3,n);
         temp3=[temp1 pop2];
         pop(2,:)=isin_TSP(temp3,n);
    case 'cospx'    %单点顺序 
         point=ceil(rand*n);
         while point==n
             point=ceil(rand*n);
         end
         temp1=[pop2(point+1:end) pop1];
         temp2=[pop1(point+1:end) pop2];
         pop(1,:)=isin_TSP(temp1,n);
         pop(2,:)=isin_TSP(temp2,n);
    case 'remm'    %对换变异
         point=sort(ceil(rand(1,2)*n));
         while point(1)==point(2)
             point=sort(ceil(rand(1,2)*n));
         end
         pop(point(1))=pop1(point(2));
         pop(point(2))=pop1(point(1));
    case 'im'  %插入变异
        point1=ceil(rand*n);
        point2=ceil(rand*n);
        while point2==point1||abs(point2-point1)==1
           point2=ceil(rand*n);
        end
        pop1=redu(pop,point1,'c');
        pop1=matrixinsert(pop1,pop(point1),point2-1);
        pop=pop1;
    case 'tm'
        m1=[];m2=[];
        for i=1:n
            if rand<pm
              m1=[m1 pop(i)];
              m2=[m2 i];
            end
        end
        m1=fliplr(m1);
        pop(m2)=m1;
        for i=1:length(m1)-1
            for j=i+1:length(m1)
                pop(m2(i)+1:m2(j)-1)=fliplr(pop(m2(i)+1:m2(j)-1));
            end
        end
    case 'em'   %进化逆转
        point=sort(ceil(rand(1,2)*n));
        while point(1)==point(2)
             point=sort(ceil(rand(1,2)*n));
        end
        pop(point(1)+1:point(2)-1)=fliplr(pop(point(1)+1:point(2)-1));
    case 'imm'   %免疫算子
        point=ceil(rand*n);
        pop=greed(route,d,point);        
end
         