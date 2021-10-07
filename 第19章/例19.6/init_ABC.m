function y=init_ABC(beenum,NC,type)   %离散ABC算法的初始化函数
for i=1:beenum
    if strcmp(type,'tsp')
       y(i,:)=randperm(NC);
    elseif strcmp(type,'bit')
       y(i,:)=rand(1,NC)<0.5;
    end
end
