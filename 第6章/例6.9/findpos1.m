function y=findpos1(x,x1,dim,type)     %找x1在x中的位置
[r,c]=size(x);
y=[];
if strcmp(dim,'r') %按行 
    for j=1:r
        if type==1  %精确
            if isequal(x(j,:),x1)
                y=[y j];
            end
        elseif type==2  %只要数字符合
            if ismember(x1,x(j,:))
                y=[y j];
            end
        end
    end
elseif strcmp(dim,'c')   %按列
   for j=1:c
        if type==1  %精确
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