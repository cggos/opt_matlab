function y=findpos1(x,x1,dim,type)     %��x1��x�е�λ��
[r,c]=size(x);
y=[];
if strcmp(dim,'r') %���� 
    for j=1:r
        if type==1  %��ȷ
            if isequal(x(j,:),x1)
                y=[y j];
            end
        elseif type==2  %ֻҪ���ַ���
            if ismember(x1,x(j,:))
                y=[y j];
            end
        end
    end
elseif strcmp(dim,'c')   %����
   for j=1:c
        if type==1  %��ȷ
            if isequal(x(:,j),x1)
                y=[y j];
            end
        elseif type==2
            if ismember(x1,x(:,j));
                y=[y j];
            end
        end
    end
end