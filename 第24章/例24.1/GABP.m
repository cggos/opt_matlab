function [F1,F2,D1,D2] = GABP(datax,datay,iterm_max,sizepop,pc,pm,net_opt)
%====输入参数====
   %input        ：输入数据
   %output       ：输出数据
   %maxgen       ：最大进化代数
   %sizepop      ：种群规模
   %pc       ：交叉概率选择，在0和1之间取值
   %pm    ：变异概率选择，在0和1之间取值
%====输出参数====
   %F1           ：未经优化的BP神经网络的预测误差
   %F2           ：经优化的BP神经网络的预测误差
   %D1           ：未经优化的BP神经网络的均方差
   %D2           ：经优化的BP神经网络的均方差   
inputnum = size(datax,2);   %输入层节点数
hiddennum = 5;  %隐层节点数
outputnum = 1;  %输出层节点数
if nargin==6
    epochs=100;
    goal=1e-6;
    lr=0.1;
else
    epochs=net_opt(1);
    lr=net_opt(2);
    goal=net_opt(3); 
end

%（2）设置训练数据和预测数据
Num1=ceil(0.85*size(datax,1));
input_train = datax(1:Num1,:)';
input_test =datax(Num1+1:end,:)';
output_train = datay(1:Num1,:)';
output_test = datay(Num1+1:end,:)';
LB=min(min(datax,[],1));
UB=max(max(datax,[],1));


%（3）选连样本输入输出数据归一化
[inputn,inputps] = mapminmax(input_train);      %对训练输入数据进行归一化
[outputn,outputps] = mapminmax(output_train);   %对训练输出数据进行归一化

%（4）生成BP神经网络
net = newff(inputn,outputn,hiddennum);

%（5）用未经优化的BP网络进行拟合并输出误差，以与经遗传算法优化后的BP网络的预测结果进行对比
   %网络进化参数
   net.trainParam.epochs = epochs;    %最大训练次数
   net.trainParam.lr =lr;        %自学习的学习效率
   net.trainParam.goal=goal;    %网络训练目标
   
   %网络训练
   [net,per2] = train(net,inputn,outputn);
   
   % BP网络预测
   %数据归一化与反归一化处理
   inputn_test = mapminmax('apply',input_test,inputps);  %测试数据归一化
   an = sim(net,inputn_test);                            %网络利用归一化后的测试数据进行预测
   test_simu1= mapminmax('reverse',an,outputps);         %对预测数据进行反归一化
   error = test_simu1- output_test;                      %预测数据与期望数据之间的差
   F1 = error;  %未经优化的BP神经网络的预测误差
   D1 = var(F1);  
   
%权值与阈值总数
numsum = inputnum*hiddennum + hiddennum + hiddennum*outputnum + outputnum;

%% 2）遗传算法寻优
lenchrom = ones(1,numsum);                     %染色体个体长度      
bound = [LB.*ones(numsum,1) UB.*ones(numsum,1)];  %染色体个体范围

%（1）初始化种群
individuals = struct('fitness',zeros(1,sizepop), 'chrom',[]);  %将种群信息定义为一个结构体
for i = 1:sizepop
    %个体初始化：随机产生一个种群，并赋予个体
    individuals.chrom(i,:) = Code(lenchrom,bound);    %染色体编码
    x = individuals.chrom(i,:);
    %调用子函数fun计算个体适应度
    individuals.fitness(i) = fun(x,inputnum,hiddennum,outputnum,inputn,outputn);   %染色体的适应度
end

%（2）找最好的染色体
[bestfitness,bestindex] = min(individuals.fitness);
bestchrom = individuals.chrom(bestindex,:);     %最好的染色体
avgfitness = sum(individuals.fitness)/sizepop;  %染色体的平均适应度

%（3）记录每一代进化中最好的适应度和平均适应度
%avgfitness：每一代种群的平均适应度
%bestfitness：每一代种群的最佳适应度
trace = [avgfitness,bestfitness]; 
 
%（4）迭代求解最佳初始阀值和权值
%进化开始
for i = 1:iterm_max
    %选择：调用子函数Select进行选择
    individuals = Select(individuals,sizepop); 
    %交叉：调用子函数Cross进行交叉选择
    individuals.chrom = Cross(pc,lenchrom,individuals.chrom,sizepop,bound);
    %变异：调用子函数Mutation进行变异选择
    individuals.chrom = Mutation(pm,lenchrom,individuals.chrom,sizepop,i,iterm_max,bound);
    %计算适应度 
    for j = 1:sizepop
        x = individuals.chrom(j,:);    %染色体解码
        individuals.fitness(j) = fun(x,inputnum,hiddennum,outputnum,inputn,outputn);   
    end
    %找到最小和最大适应度的染色体及它们在种群中的位置
    [newbestfitness,newbestindex] = min(individuals.fitness);   %最小的适应度及其位置
    [worestfitness,worestindex] = max(individuals.fitness);     %最大的适应度及其位置
    %用新的最好的染色体替换上一次进化中最好的染色体
    if bestfitness > newbestfitness
        bestfitness = newbestfitness;
        bestchrom = individuals.chrom(newbestindex,:);          %适应度最好的染色体
    end
    individuals.chrom(worestindex,:) = bestchrom;
    individuals.fitness(worestindex) = bestfitness;    
    avgfitness = sum(individuals.fitness)/sizepop;    
    trace = [trace;avgfitness bestfitness];          %记录每一代进化中最好的适应度和平均适应度
end
%（5）遗传算法结果分析 
 figure(1)    %设定绘制各代平均适应度的图形窗口
[r,c] = size(trace);
plot([1:r]',trace(:,1),'b--diamond');      %绘制各代的平均适应度
hold on
plot([1:r]',trace(:,2),'r-*');             %绘制各代的最好适应度
title(['适应度曲线（终止代数=' num2str(iterm_max) '）']);
xlabel('进化代数');ylabel('适应度');
legend('平均适应度','最佳适应度');

%% 3)把最优初始阀值权值赋予BP神经网络预测
% （1）用遗传算法优化的BP网络进行值预测
x = bestchrom;
w1 = x(1:inputnum*hiddennum);                                %输入层与隐层权值
B1 = x(inputnum*hiddennum+1:inputnum*hiddennum+hiddennum);   %隐层节点阈值
w2 = x(inputnum*hiddennum+hiddennum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum);  %隐层与输出层权值
B2 = x(inputnum*hiddennum+hiddennum+hiddennum*outputnum+1:numsum); %输出层节点阈值

net.iw{1,1} = reshape(w1,hiddennum,inputnum);
net.lw{2,1} = reshape(w2,outputnum,hiddennum);
net.b{1} = reshape(B1,hiddennum,1);
net.b{2} = B2;

%（2）BP网络训练
%BP神经网络进化参数
net.trainParam.epochs = 100;
net.trainParam.lr = 0.1;
net.trainParam.goal = 0.00001;

%（3）BP神经网络训练
[net,per2] = train(net,inputn,outputn);

%（4）BP网络预测
%数据归一化与反归一化
inputn_test = mapminmax('apply',input_test,inputps);   %对经过优化后的测试数据归一化
an = sim(net,inputn_test);                             %对归一化的数据进行网络预测
test_simu2= mapminmax('reverse',an,outputps);          %对网络预测数据进行反归一化
error = test_simu2-output_test;                       %网络预测值与期望值之间的差
F2 = error;   %经优化后的BP神经网络的预测误差
D2 = var(F2);
figure(2)
plot(input_test,test_simu1,'P');hold on
plot(input_test,test_simu2,'*');
plot(input_test,output_test,'or');
legend('BP','GABP','origin');

%% 1)适应度子函数
function error = fun(x,inputnum,hiddennum,outputnum,inputn,outputn)
%该函数用来计算适应度值
%输入参数
    %x          个体
    %inputnum   输入层节点数
    %outputnum  隐含层节点数
    %net        网络
    %inputn     训练输入数据
    %outputn    训练输出数据
%输出参数   
    %error      个体适应度值
%提取权值与阈值
numsum = inputnum*hiddennum + hiddennum + hiddennum*outputnum + outputnum;
w1 = x(1:inputnum*hiddennum);
B1 = x(inputnum*hiddennum+1:inputnum*hiddennum+hiddennum);
w2 = x(inputnum*hiddennum+hiddennum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum);
B2 = x(inputnum*hiddennum+hiddennum+hiddennum*outputnum+1:numsum);

%生成BP神经网络
net = newff(inputn,outputn,hiddennum);

%网络进化参数
net.trainParam.epochs = 20;
net.trainParam.lr = 0.1;
net.trainParam.goal = 0.00001;
net.trainParam.show = 100;
net.trainParam.showWindow = 0;
 
%网络权值赋值
net.iw{1,1} = reshape(w1,hiddennum,inputnum);
net.lw{2,1} = reshape(w2,outputnum,hiddennum);
net.b{1} = reshape(B1,hiddennum,1);
net.b{2} = B2;

%网络训练
net = train(net,inputn,outputn);
an = sim(net,inputn);
error = sum(abs(an-outputn));   %预测误差和作为个体适应度值

%% 2)选择操作子函数
function ret = Select(individuals,sizepop)
% 该函数用于进行选择操作
%输入参数
   % individuals 种群信息
   % sizepop     种群规模
%输出参数
   % ret         选择后的新种群

%求适应度值倒数
fitness1 = 10./individuals.fitness; %individuals.fitness为个体适应度值

%个体选择概率
sumfitness = sum(fitness1);
sumf = fitness1./sumfitness;

%采用轮盘赌法选择新个体
index = []; 
for i = 1:sizepop   %sizepop为种群数
    pick = rand;
    while pick == 0    
        pick = rand;        
    end
    for j = 1:sizepop    
        pick = pick-sumf(j);        
        if pick < 0        
            index = [index j];            
            break;  
        end
    end
end

%新种群
individuals.chrom = individuals.chrom(index,:);   %individuals.chrom为种群中个体
individuals.fitness = individuals.fitness(index);
ret = individuals;   %选择的染色体

%% 3)交叉操作子函数
function ret = Cross(pc,lenchrom,chrom,sizepop,bound)
%本函数完成交叉操作
%输入参数
    %pcorss    交叉概率
    %lenchrom  染色体的长度
    %chrom     染色体群
    %sizepop   种群规模
    %bound     个体上界和下界
%输出参数
    % ret       交叉后的染色体
    
 for i = 1:sizepop  %每一轮for循环中，可能会进行一次交叉操作，染色体是随机选择的，交叉位置也是随机选择的，
     %但该轮for循环中是否进行交叉操作则由交叉概率决定（continue控制）     
     pick=rand(1,2);    %随机选择两个染色体进行交叉
     while prod(pick)==0
         pick=rand(1,2);
     end
     index = ceil(pick.*sizepop);
     pick = rand;       %交叉概率决定是否进行交叉
     while pick == 0
         pick = rand;
     end
     if pick > pc
         continue;
     end
     flag = 0;
     while flag == 0
         % 随机选择交叉位置
         pick = rand;
         while pick == 0
             pick = rand;
         end
         pos = ceil(pick.*sum(lenchrom)); %随机选择进行交叉的位置，即选择第几个变量进行交叉，注意：两个染色体交叉的位置相同
         pick = rand;   %开始个体交叉
         v1 = chrom(index(1),pos);
         v2 = chrom(index(2),pos);
         chrom(index(1),pos) = pick*v2+(1-pick)*v1;
         chrom(index(2),pos) = pick*v1+(1-pick)*v2;       %交叉结束
         flag1 = test(lenchrom,bound,chrom(index(1),:));  %检验染色体1的可行性
         flag2 = test(lenchrom,bound,chrom(index(2),:));  %检验染色体2的可行性
         if   flag1*flag2 == 0
             flag = 0;
         else flag = 1;
         end    %如果两个染色体不是都可行，则重新交叉
     end
 end
ret = chrom;    %交叉操作后的染色体

%% 4)变异操作子函数
function ret = Mutation(pm,lenchrom,chrom,sizepop,num,iterm_max,bound)
% 本函数完成变异操作
%输入参数
   %pcorss    变异概率
   %lenchrom  染色体长度
   %chrom     染色体群
   %sizepop   种群规模
   %opts      变异方法的选择
   %pop       当前种群的进化代数和最大的进化代数信息
   %bound     每个个体的上届和下届
   %maxgen    最大迭代次数
   %num       当前迭代次数
%输出参数
   %ret       变异后的染色体

for i=1:sizepop   %每一轮for循环中，可能会进行一次变异操作，染色体是随机选择的，变异位置也是随机选择的，
    %但该轮for循环中是否进行变异操作则由变异概率决定（continue控制）
    pick = rand;  %随机选择一个染色体进行变异
    while pick == 0
        pick = rand;
    end
    index = ceil(pick*sizepop);
    pick = rand;  %变异概率决定该轮循环是否进行变异
    if pick > pm
        continue;
    end
    flag = 0;
    while flag == 0
        %随机选择变异位置
        pick = rand;
        while pick == 0      
            pick = rand;
        end
        pos = ceil(pick*sum(lenchrom));  %随机选择了染色体变异的位置，即选择了第pos个变量进行变异
        %开始个体变异
        pick = rand;      
        fg = (rand*(1-num/iterm_max))^2;
        if pick > 0.5
            chrom(index,pos) = chrom(index,pos)+(bound(pos,2)-chrom(index,pos))*fg;
        else
            chrom(index,pos) = chrom(index,pos)-(chrom(index,pos)-bound(pos,1))*fg;
        end
        %变异结束
        flag = test(lenchrom,bound,chrom(i,:));     %检验染色体的可行性
    end
end
ret = chrom;   %变异操作后的染色体

%% 5)编码子函数
function ret = Code(lenchrom,bound)
%本函数将变量编码成染色体，用于随机初始化一个种群
%输入参数
   % lenchrom  染色体长度
   % bound     变量的取值范围
%输出参数
   % ret       染色体的编码值
flag = 0;
while flag == 0
    pick = rand(1,length(lenchrom));
    ret = bound(:,1)'+(bound(:,2)-bound(:,1))'.*pick;  %线性插值，编码结果以实数向量存入ret中
    flag = test(lenchrom,bound,ret);                   %检验染色体的可行性
end

%% 6)测试子函数
function flag = test(lenchrom,bound,code)
%输入参数
   % lenchrom   染色体长度
   % bound      变量的取值范围
%输出参数
   % code       染色体的编码值
x = code; %先解码
flag = 1;