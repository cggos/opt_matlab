function varargout=redu(x,m,type)%ɾ��ָ�����л���
temp=0;y=x;m=sort(m);
switch type
    case 'c'
       y1=[];
       for i=1:length(m)   %ɾ����
         y1=[y1 y(:,m(i)-temp)];
         y(:,m(i)-temp)=[];
         temp=temp+1;
       end
    case 'r'
       y1=[];
       for i=1:length(m)   %ɾ����
          y1=[y1;y(m(i)-temp,:)];
          y(m(i)-temp,:)=[];
          temp=temp+1;
       end
    case 'rc'     %ɾ������
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
    varargout={y,y1};%yΪʣ��ľ���y1Ϊɾ���ľ���
else
    varargout={y};
end