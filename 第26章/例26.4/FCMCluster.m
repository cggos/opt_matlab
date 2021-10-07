function [JbValue,A_Matrix,Center] = FCMCluster(X,CN)   %模糊聚类
%====输入参数====
   %X            ：待聚类的数据
   %CN           ：聚类的数目
%====输出参数====
   %JbValue      ：目标函数值
   %A_Matrix     ：各聚类矩阵
   %Center       ：各聚类的中心
figure(1)
plot(X(:,1),X(:,2),'o')
xlabel('横坐标X');ylabel('纵坐标Y');title('样本数据')
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
    text(Center(i,1),Center(i,2),'中心位置');
    title(['FCM聚类' num2str(i) '的分布图']);
    xlabel('横坐标');ylabel('纵坐标');    
end
A_Matrix0 = cell(CN,1);
for i = 1:CN
    A_Matrix0{i} = [X(Newindex{i},1), X(Newindex{i}, 2)];
end
A_Matrix = A_Matrix0;