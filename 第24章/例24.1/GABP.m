function [F1,F2,D1,D2] = GABP(datax,datay,iterm_max,sizepop,pc,pm,net_opt)
%====�������====
   %input        ����������
   %output       ���������
   %maxgen       ������������
   %sizepop      ����Ⱥ��ģ
   %pc       ���������ѡ����0��1֮��ȡֵ
   %pm    ���������ѡ����0��1֮��ȡֵ
%====�������====
   %F1           ��δ���Ż���BP�������Ԥ�����
   %F2           �����Ż���BP�������Ԥ�����
   %D1           ��δ���Ż���BP������ľ�����
   %D2           �����Ż���BP������ľ�����   
inputnum = size(datax,2);   %�����ڵ���
hiddennum = 5;  %����ڵ���
outputnum = 1;  %�����ڵ���
if nargin==6
    epochs=100;
    goal=1e-6;
    lr=0.1;
else
    epochs=net_opt(1);
    lr=net_opt(2);
    goal=net_opt(3); 
end

%��2������ѵ�����ݺ�Ԥ������
Num1=ceil(0.85*size(datax,1));
input_train = datax(1:Num1,:)';
input_test =datax(Num1+1:end,:)';
output_train = datay(1:Num1,:)';
output_test = datay(Num1+1:end,:)';
LB=min(min(datax,[],1));
UB=max(max(datax,[],1));


%��3��ѡ����������������ݹ�һ��
[inputn,inputps] = mapminmax(input_train);      %��ѵ���������ݽ��й�һ��
[outputn,outputps] = mapminmax(output_train);   %��ѵ��������ݽ��й�һ��

%��4������BP������
net = newff(inputn,outputn,hiddennum);

%��5����δ���Ż���BP���������ϲ���������뾭�Ŵ��㷨�Ż����BP�����Ԥ�������жԱ�
   %�����������
   net.trainParam.epochs = epochs;    %���ѵ������
   net.trainParam.lr =lr;        %��ѧϰ��ѧϰЧ��
   net.trainParam.goal=goal;    %����ѵ��Ŀ��
   
   %����ѵ��
   [net,per2] = train(net,inputn,outputn);
   
   % BP����Ԥ��
   %���ݹ�һ���뷴��һ������
   inputn_test = mapminmax('apply',input_test,inputps);  %�������ݹ�һ��
   an = sim(net,inputn_test);                            %�������ù�һ����Ĳ������ݽ���Ԥ��
   test_simu1= mapminmax('reverse',an,outputps);         %��Ԥ�����ݽ��з���һ��
   error = test_simu1- output_test;                      %Ԥ����������������֮��Ĳ�
   F1 = error;  %δ���Ż���BP�������Ԥ�����
   D1 = var(F1);  
   
%Ȩֵ����ֵ����
numsum = inputnum*hiddennum + hiddennum + hiddennum*outputnum + outputnum;

%% 2���Ŵ��㷨Ѱ��
lenchrom = ones(1,numsum);                     %Ⱦɫ����峤��      
bound = [LB.*ones(numsum,1) UB.*ones(numsum,1)];  %Ⱦɫ����巶Χ

%��1����ʼ����Ⱥ
individuals = struct('fitness',zeros(1,sizepop), 'chrom',[]);  %����Ⱥ��Ϣ����Ϊһ���ṹ��
for i = 1:sizepop
    %�����ʼ�����������һ����Ⱥ�����������
    individuals.chrom(i,:) = Code(lenchrom,bound);    %Ⱦɫ�����
    x = individuals.chrom(i,:);
    %�����Ӻ���fun���������Ӧ��
    individuals.fitness(i) = fun(x,inputnum,hiddennum,outputnum,inputn,outputn);   %Ⱦɫ�����Ӧ��
end

%��2������õ�Ⱦɫ��
[bestfitness,bestindex] = min(individuals.fitness);
bestchrom = individuals.chrom(bestindex,:);     %��õ�Ⱦɫ��
avgfitness = sum(individuals.fitness)/sizepop;  %Ⱦɫ���ƽ����Ӧ��

%��3����¼ÿһ����������õ���Ӧ�Ⱥ�ƽ����Ӧ��
%avgfitness��ÿһ����Ⱥ��ƽ����Ӧ��
%bestfitness��ÿһ����Ⱥ�������Ӧ��
trace = [avgfitness,bestfitness]; 
 
%��4�����������ѳ�ʼ��ֵ��Ȩֵ
%������ʼ
for i = 1:iterm_max
    %ѡ�񣺵����Ӻ���Select����ѡ��
    individuals = Select(individuals,sizepop); 
    %���棺�����Ӻ���Cross���н���ѡ��
    individuals.chrom = Cross(pc,lenchrom,individuals.chrom,sizepop,bound);
    %���죺�����Ӻ���Mutation���б���ѡ��
    individuals.chrom = Mutation(pm,lenchrom,individuals.chrom,sizepop,i,iterm_max,bound);
    %������Ӧ�� 
    for j = 1:sizepop
        x = individuals.chrom(j,:);    %Ⱦɫ�����
        individuals.fitness(j) = fun(x,inputnum,hiddennum,outputnum,inputn,outputn);   
    end
    %�ҵ���С�������Ӧ�ȵ�Ⱦɫ�弰��������Ⱥ�е�λ��
    [newbestfitness,newbestindex] = min(individuals.fitness);   %��С����Ӧ�ȼ���λ��
    [worestfitness,worestindex] = max(individuals.fitness);     %������Ӧ�ȼ���λ��
    %���µ���õ�Ⱦɫ���滻��һ�ν�������õ�Ⱦɫ��
    if bestfitness > newbestfitness
        bestfitness = newbestfitness;
        bestchrom = individuals.chrom(newbestindex,:);          %��Ӧ����õ�Ⱦɫ��
    end
    individuals.chrom(worestindex,:) = bestchrom;
    individuals.fitness(worestindex) = bestfitness;    
    avgfitness = sum(individuals.fitness)/sizepop;    
    trace = [trace;avgfitness bestfitness];          %��¼ÿһ����������õ���Ӧ�Ⱥ�ƽ����Ӧ��
end
%��5���Ŵ��㷨������� 
 figure(1)    %�趨���Ƹ���ƽ����Ӧ�ȵ�ͼ�δ���
[r,c] = size(trace);
plot([1:r]',trace(:,1),'b--diamond');      %���Ƹ�����ƽ����Ӧ��
hold on
plot([1:r]',trace(:,2),'r-*');             %���Ƹ����������Ӧ��
title(['��Ӧ�����ߣ���ֹ����=' num2str(iterm_max) '��']);
xlabel('��������');ylabel('��Ӧ��');
legend('ƽ����Ӧ��','�����Ӧ��');

%% 3)�����ų�ʼ��ֵȨֵ����BP������Ԥ��
% ��1�����Ŵ��㷨�Ż���BP�������ֵԤ��
x = bestchrom;
w1 = x(1:inputnum*hiddennum);                                %�����������Ȩֵ
B1 = x(inputnum*hiddennum+1:inputnum*hiddennum+hiddennum);   %����ڵ���ֵ
w2 = x(inputnum*hiddennum+hiddennum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum);  %�����������Ȩֵ
B2 = x(inputnum*hiddennum+hiddennum+hiddennum*outputnum+1:numsum); %�����ڵ���ֵ

net.iw{1,1} = reshape(w1,hiddennum,inputnum);
net.lw{2,1} = reshape(w2,outputnum,hiddennum);
net.b{1} = reshape(B1,hiddennum,1);
net.b{2} = B2;

%��2��BP����ѵ��
%BP�������������
net.trainParam.epochs = 100;
net.trainParam.lr = 0.1;
net.trainParam.goal = 0.00001;

%��3��BP������ѵ��
[net,per2] = train(net,inputn,outputn);

%��4��BP����Ԥ��
%���ݹ�һ���뷴��һ��
inputn_test = mapminmax('apply',input_test,inputps);   %�Ծ����Ż���Ĳ������ݹ�һ��
an = sim(net,inputn_test);                             %�Թ�һ�������ݽ�������Ԥ��
test_simu2= mapminmax('reverse',an,outputps);          %������Ԥ�����ݽ��з���һ��
error = test_simu2-output_test;                       %����Ԥ��ֵ������ֵ֮��Ĳ�
F2 = error;   %���Ż����BP�������Ԥ�����
D2 = var(F2);
figure(2)
plot(input_test,test_simu1,'P');hold on
plot(input_test,test_simu2,'*');
plot(input_test,output_test,'or');
legend('BP','GABP','origin');

%% 1)��Ӧ���Ӻ���
function error = fun(x,inputnum,hiddennum,outputnum,inputn,outputn)
%�ú�������������Ӧ��ֵ
%�������
    %x          ����
    %inputnum   �����ڵ���
    %outputnum  ������ڵ���
    %net        ����
    %inputn     ѵ����������
    %outputn    ѵ���������
%�������   
    %error      ������Ӧ��ֵ
%��ȡȨֵ����ֵ
numsum = inputnum*hiddennum + hiddennum + hiddennum*outputnum + outputnum;
w1 = x(1:inputnum*hiddennum);
B1 = x(inputnum*hiddennum+1:inputnum*hiddennum+hiddennum);
w2 = x(inputnum*hiddennum+hiddennum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum);
B2 = x(inputnum*hiddennum+hiddennum+hiddennum*outputnum+1:numsum);

%����BP������
net = newff(inputn,outputn,hiddennum);

%�����������
net.trainParam.epochs = 20;
net.trainParam.lr = 0.1;
net.trainParam.goal = 0.00001;
net.trainParam.show = 100;
net.trainParam.showWindow = 0;
 
%����Ȩֵ��ֵ
net.iw{1,1} = reshape(w1,hiddennum,inputnum);
net.lw{2,1} = reshape(w2,outputnum,hiddennum);
net.b{1} = reshape(B1,hiddennum,1);
net.b{2} = B2;

%����ѵ��
net = train(net,inputn,outputn);
an = sim(net,inputn);
error = sum(abs(an-outputn));   %Ԥ��������Ϊ������Ӧ��ֵ

%% 2)ѡ������Ӻ���
function ret = Select(individuals,sizepop)
% �ú������ڽ���ѡ�����
%�������
   % individuals ��Ⱥ��Ϣ
   % sizepop     ��Ⱥ��ģ
%�������
   % ret         ѡ��������Ⱥ

%����Ӧ��ֵ����
fitness1 = 10./individuals.fitness; %individuals.fitnessΪ������Ӧ��ֵ

%����ѡ�����
sumfitness = sum(fitness1);
sumf = fitness1./sumfitness;

%�������̶ķ�ѡ���¸���
index = []; 
for i = 1:sizepop   %sizepopΪ��Ⱥ��
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

%����Ⱥ
individuals.chrom = individuals.chrom(index,:);   %individuals.chromΪ��Ⱥ�и���
individuals.fitness = individuals.fitness(index);
ret = individuals;   %ѡ���Ⱦɫ��

%% 3)��������Ӻ���
function ret = Cross(pc,lenchrom,chrom,sizepop,bound)
%��������ɽ������
%�������
    %pcorss    �������
    %lenchrom  Ⱦɫ��ĳ���
    %chrom     Ⱦɫ��Ⱥ
    %sizepop   ��Ⱥ��ģ
    %bound     �����Ͻ���½�
%�������
    % ret       ������Ⱦɫ��
    
 for i = 1:sizepop  %ÿһ��forѭ���У����ܻ����һ�ν��������Ⱦɫ�������ѡ��ģ�����λ��Ҳ�����ѡ��ģ�
     %������forѭ�����Ƿ���н���������ɽ�����ʾ�����continue���ƣ�     
     pick=rand(1,2);    %���ѡ������Ⱦɫ����н���
     while prod(pick)==0
         pick=rand(1,2);
     end
     index = ceil(pick.*sizepop);
     pick = rand;       %������ʾ����Ƿ���н���
     while pick == 0
         pick = rand;
     end
     if pick > pc
         continue;
     end
     flag = 0;
     while flag == 0
         % ���ѡ�񽻲�λ��
         pick = rand;
         while pick == 0
             pick = rand;
         end
         pos = ceil(pick.*sum(lenchrom)); %���ѡ����н����λ�ã���ѡ��ڼ����������н��棬ע�⣺����Ⱦɫ�彻���λ����ͬ
         pick = rand;   %��ʼ���彻��
         v1 = chrom(index(1),pos);
         v2 = chrom(index(2),pos);
         chrom(index(1),pos) = pick*v2+(1-pick)*v1;
         chrom(index(2),pos) = pick*v1+(1-pick)*v2;       %�������
         flag1 = test(lenchrom,bound,chrom(index(1),:));  %����Ⱦɫ��1�Ŀ�����
         flag2 = test(lenchrom,bound,chrom(index(2),:));  %����Ⱦɫ��2�Ŀ�����
         if   flag1*flag2 == 0
             flag = 0;
         else flag = 1;
         end    %�������Ⱦɫ�岻�Ƕ����У������½���
     end
 end
ret = chrom;    %����������Ⱦɫ��

%% 4)��������Ӻ���
function ret = Mutation(pm,lenchrom,chrom,sizepop,num,iterm_max,bound)
% ��������ɱ������
%�������
   %pcorss    �������
   %lenchrom  Ⱦɫ�峤��
   %chrom     Ⱦɫ��Ⱥ
   %sizepop   ��Ⱥ��ģ
   %opts      ���췽����ѡ��
   %pop       ��ǰ��Ⱥ�Ľ������������Ľ���������Ϣ
   %bound     ÿ��������Ͻ���½�
   %maxgen    ����������
   %num       ��ǰ��������
%�������
   %ret       ������Ⱦɫ��

for i=1:sizepop   %ÿһ��forѭ���У����ܻ����һ�α��������Ⱦɫ�������ѡ��ģ�����λ��Ҳ�����ѡ��ģ�
    %������forѭ�����Ƿ���б���������ɱ�����ʾ�����continue���ƣ�
    pick = rand;  %���ѡ��һ��Ⱦɫ����б���
    while pick == 0
        pick = rand;
    end
    index = ceil(pick*sizepop);
    pick = rand;  %������ʾ�������ѭ���Ƿ���б���
    if pick > pm
        continue;
    end
    flag = 0;
    while flag == 0
        %���ѡ�����λ��
        pick = rand;
        while pick == 0      
            pick = rand;
        end
        pos = ceil(pick*sum(lenchrom));  %���ѡ����Ⱦɫ������λ�ã���ѡ���˵�pos���������б���
        %��ʼ�������
        pick = rand;      
        fg = (rand*(1-num/iterm_max))^2;
        if pick > 0.5
            chrom(index,pos) = chrom(index,pos)+(bound(pos,2)-chrom(index,pos))*fg;
        else
            chrom(index,pos) = chrom(index,pos)-(chrom(index,pos)-bound(pos,1))*fg;
        end
        %�������
        flag = test(lenchrom,bound,chrom(i,:));     %����Ⱦɫ��Ŀ�����
    end
end
ret = chrom;   %����������Ⱦɫ��

%% 5)�����Ӻ���
function ret = Code(lenchrom,bound)
%�����������������Ⱦɫ�壬���������ʼ��һ����Ⱥ
%�������
   % lenchrom  Ⱦɫ�峤��
   % bound     ������ȡֵ��Χ
%�������
   % ret       Ⱦɫ��ı���ֵ
flag = 0;
while flag == 0
    pick = rand(1,length(lenchrom));
    ret = bound(:,1)'+(bound(:,2)-bound(:,1))'.*pick;  %���Բ�ֵ����������ʵ����������ret��
    flag = test(lenchrom,bound,ret);                   %����Ⱦɫ��Ŀ�����
end

%% 6)�����Ӻ���
function flag = test(lenchrom,bound,code)
%�������
   % lenchrom   Ⱦɫ�峤��
   % bound      ������ȡֵ��Χ
%�������
   % code       Ⱦɫ��ı���ֵ
x = code; %�Ƚ���
flag = 1;