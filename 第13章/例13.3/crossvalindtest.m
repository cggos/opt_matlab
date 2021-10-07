function y=crossvalindtest(data,k,num,target)   %k-�۽�����֤,numΪ�ظ�����
%targetΪ��Ʒ�ķ���ȼ�
M=size(data,1);
y=0;
for j=1:num
  indices=crossvalind('Kfold',M,k);  %��������ְ�
  y2=0;
  for i=1:k               %������֤,kΪ����
    test=(indices==i);   %���test��Ԫ�������ݼ��ж�Ӧ�ĵ�Ԫ���
    train=~test;         %train��Ԫ�صı��Ϊ��testԪ�صı��
    train_data=data(train,:);   %�����ݼ��л��ֳ�train����������
    train_target=target(train);   %����������Ĳ���Ŀ��
    test_data=data(test,:);         %test������
    test_target=target(test);
    knnmd1=fitcknn(train_data,train_target);
    y1=predict(knnmd1,test_data);
    %y1=knnclassify(test_data,train_data,train_target);
    y2=y2+correct(y1,test_target);
  end
  y2=y2/k;
  y=y+y2;
end
y=y/num;      %��ȷ��