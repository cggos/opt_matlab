function [JbValue,A_Matrix,Center] = FCMCluster(X,CN)   %ģ������
%====�������====
   %X            �������������
   %CN           ���������Ŀ
%====�������====
   %JbValue      ��Ŀ�꺯��ֵ
   %A_Matrix     �����������
   %Center       �������������
figure(1)
plot(X(:,1),X(:,2),'o')
xlabel('������X');ylabel('������Y');title('��������')
options=[3,20,1e-6,0];
[Center,U,obj_fcn] = fcm(X,CN,options);
JbValue = obj_fcn(end);
maxU = max(U);
Aindex = cell(CN,1);
for i = 1:CN
    Aindex{i} = find(U(i,:) == maxU);
end
Newindex = Aindex;
for i = 1:CN
    figure(i+1)
    subplot(1,2,2)
    line(X(Newindex{i},1), X(Newindex{i},2), 'linestyle', 'none','marker', 'x', 'color', 'r');
    hold on
    plot(Center(i,1),Center(i,2),'diamond','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
    text(Center(i,1),Center(i,2),'����λ��');
    title(['FCM����' num2str(i) '�ķֲ�ͼ']);
    xlabel('������');ylabel('������');    
end
A_Matrix0 = cell(CN,1);
for i = 1:CN
    A_Matrix0{i} = [X(Newindex{i},1), X(Newindex{i}, 2)];
end
A_Matrix = A_Matrix0;