function [MinValue,MinFounction]=gaexhause(fun,LB,UB,eps)        %�Ŵ��㷨�ı����㷨
%epsΪ����
n=ceil(log2((UB-LB)/eps+1));
Individual={};
for i=1:n
    Individual=[Individual;[0 1]];
end
IndividualCode=exhaus(Individual);
%Ⱦɫ������Ӧ��10������
IndividualCodeValue=zeros(1,2^n);
for i=1:2^n
   for j=1:n
       IndividualCodeValue(i)=IndividualCodeValue(i)+IndividualCode(i,j)*(2^(n-j));
   end
end
%Ⱦɫ������Ӧ�����������ʵ��
CodeValue=zeros(1,2^n);
for k=1:2^n
   CodeValue(k)=CodeValue(k)+ LB+IndividualCodeValue(k)*(UB-LB)/(2^n-1);
end
MinFounction=fun(CodeValue(1)); %����Ŀ�꺯������ʼ�����һ��ֵΪ��Сֵ
MinIndividual=IndividualCode(1);
MinValue=CodeValue(1);
for i=2:2^n
   FounctionValue(i)=fun(CodeValue(i));
   if FounctionValue(i)<MinFounction
       MinFounction=FounctionValue(i);
       MinIndividual=IndividualCode(i);
       MinValue=CodeValue(i);
   end
end
