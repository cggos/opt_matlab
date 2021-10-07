function varargout=redu(x,m,type)%删除指定的行或列
temp=0;y=x;m=sort(m);
switch type
    case 'c'
       y1=[];
       for i=1:length(m)   %删除列
         y1=[y1 y(:,m(i)-temp)];
         y(:,m(i)-temp)=[];
         temp=temp+1;
       end
    case 'r'
       y1=[];
       for i=1:length(m)   %删除行
          y1=[y1;y(m(i)-temp,:)];
          y(m(i)-temp,:)=[];
          temp=temp+1;
       end
    case 'rc'     %删除矩阵
        y1={};
        for i=1:length(m)
           y1=[y1;y(m(i)-temp,:)];
           y(m(i)-temp,:)=[];
           y1=[y1;y(:,m(i)-temp)];
           y(:,m(i)-temp)=[];
           temp=temp+1;
        end        
end
if nargout==2
    varargout={y,y1};%y为剩余的矩阵，y1为删除的矩阵
else
    varargout={y};
end