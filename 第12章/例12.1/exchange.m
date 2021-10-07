function [df,route_after]=exchange(route,d,flag,type)
n=length(route);
switch type
    case 1   %交换1,2个城市
        loc=sort(ceil(rand(1,2)*n));
        if flag==1    %固定起点
            while loc(1)==loc(2)||loc(1)==1
               loc=sort(ceil(rand(1,2)*n));
            end
        elseif flag==0
            while loc(1)==loc(2)
              loc=sort(ceil(rand(1,2)*n));
            end
        end
        d1=d(route(loc(1)-1),route(loc(1)))+d(route(loc(1)+1),route(loc(1)));
        if loc(2)~=n
           d1=d1+d(route(loc(2)-1),route(loc(2)))+d(route(loc(2)+1),route(loc(2)));
        else
           d1=d1+d(route(loc(2)-1),route(loc(2)));
        end
        temp=route(loc(1));
        route(loc(1))=route(loc(2));
        route(loc(2))=temp;
        route_after=route;
        d2=d(route_after(loc(1)-1),route_after(loc(1)))+d(route_after(loc(1)+1),route_after(loc(1)));
        if loc(2)~=n
           d2=d2+d(route_after(loc(2)-1),route_after(loc(2)))+d(route_after(loc(2)+1),route_after(loc(2)));
        else
           d2=d2+d(route_after(loc(2)-1),route_after(loc(2)));
        end
        df=d2-d1;
    case 2   %将一段路程插入另一城市后
        loc=sort(ceil(rand(1,3)*n));
        if flag==1    %固定起点
            while loc(1)==loc(2)||loc(1)==1||loc(2)==loc(3)
               loc=sort(ceil(rand(1,3)*n));
            end
        elseif flag==0
            while loc(1)==loc(2)||loc(2)==loc(3)
              loc=sort(ceil(rand(1,3)*n));
            end
        end       
        temp=matrixinsert(route,{route(loc(1):loc(2))},loc(3));
        route_after=redu(temp,loc(1):loc(2),'c');
        d1=d(route(loc(1)-1),route(loc(1)))+d(route(loc(2)+1),route(loc(2)));
        if loc(3)~=n
          d1=d1+d(route(loc(3)+1),route(loc(3)));
        end
        d2=d(route(loc(1)-1),route(loc(2)+1))+d(route(loc(3)),route(loc(1)));
        if loc(3)~=n
           d2=d2+d(route(loc(2)),route(loc(3)+1));           
        end 
        df=d2-d1;     
    case 3
       route_after=route;
       loc=ceil(rand(1,2)*n);
       if loc(1)>=loc(2)    %中间倒置
           loc=sort(loc);
           while loc(1)==loc(2)||loc(2)-loc(1)<=2
              loc=ceil(rand(1,2)*n);
           end
           route_after(loc(1)+1:loc(2)-1)=fliplr(route_after(loc(1)+1:loc(2)-1));
           d1=d(route(loc(1)+1),route(loc(1)))+d(route(loc(2)-1),route(loc(2)));
           d2=d(route_after(loc(1)+1),route_after(loc(1)))+d(route_after(loc(2)-1),route_after(loc(2)));
           df=d2-d1;      
       else               %两端倒置
           if flag==1
              loc=sort(loc);
              while loc(1)<3||loc(2)>n-2
                loc=sort(ceil(rand(1,2)*n));
              end
              route_after(2:loc(1)-1)=fliplr(route_after(2:loc(1)-1));
              route_after(loc(2)+1:end)=fliplr(route_after(loc(2)+1:end));
           elseif flag==0
               loc=sort(loc);
              while loc(1)<2||loc(2)>n-2
                loc=sort(ceil(rand(1,2)*n));
              end
              route_after(1:loc(1)-1)=fliplr(route_after(loc(1)-1));
              route_after(loc(2)+1:end)=fliplr(route_after(loc(2)+1:end));
           end
           d1=d(route(loc(1)-1),route(loc(1)))+d(route(loc(2)+1),route(loc(2)));
           d2=d(route_after(loc(1)-1),route_after(loc(1)))+d(route_after(loc(2)+1),route_after(loc(2)));
           df=d2-d1;          
       end
 end
