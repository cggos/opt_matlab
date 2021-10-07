function y=crossvalindtest(data,k,num,target)   %k-折交叉验证,num为重复次数
%target为样品的分类等级
M=size(data,1);
y=0;
for j=1:num
  indices=crossvalind('Kfold',M,k);  %进行随机分包
  y2=0;
  for i=1:k               %交叉验证,k为折数
    test=(indices==i);   %获得test集元素在数据集中对应的单元编号
    train=~test;         %train集元素的编号为非test元素的编号
    train_data=data(train,:);   %从数据集中划分出train样本的数据
    train_target=target(train);   %获得样本集的测试目标
    test_data=data(test,:);         %test样本集
    test_target=target(test);
    knnmd1=fitcknn(train_data,train_target);
    y1=predict(knnmd1,test_data);
    %y1=knnclassify(test_data,train_data,train_target);
    y2=y2+correct(y1,test_target);
  end
  y2=y2/k;
  y=y+y2;
end
y=y/num;      %正确率