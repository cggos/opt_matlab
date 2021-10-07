function [cbest_var,cbest_f]=qGA1(fun)   %С���������Ŵ��㷨
%numΪÿ������εĳ���
prompt={'������';'��������';'�����½�';'�����Ͻ�';'�������';'�Ա�����ɢ����'};
name='�����㷨������';
defaultanswer={'20','200','-inf','inf','0.02','1e-10'};
answer=inputdlg(prompt,name,1,defaultanswer);
popsize=str2num(answer{1});
iter_max=str2num(answer{2});
LB=str2num(answer{3});
UB=str2num(answer{4});
pm=str2num(answer{5});
eps=str2num(answer{6});
NC=size(LB,1);
%num=NC*max(ceil(log2((UB-LB)/eps + 1)));      %�����Ա�����ɢ���ȣ�ȷ�������Ʊ���λ���ĳ���
num=max(ceil(log2((UB-LB)/eps + 1))); 
for i=1:popsize
    pop(i).p=zeros(1,NC*num);
    pop(i).fai=pi/4*ones(1,NC*num);
    pop(i).q=zeros(2,NC*num);
    for j=1:NC*num
        pop(i).q(1,j)=cos(pop(i).fai(j));
        pop(i).q(2,j)=sin(pop(i).fai(j)); 
    end
    pop(i).fitness=0;
end
for i=1:popsize
    for j=1:NC*num
        p=rand;
        if p>0.5
            pop(i).p(1,j)=1;
        else
            pop(i).p(1,j)=0;
        end
    end
end
pop=Dec1(pop,LB,UB,NC*num,1);%����,1Ϊ����������  (pop,LB,UB,CodeLen,n) 
for i=1:popsize
    pop(i).fitness=fun(pop(i).var);
    f(i)=pop(i).fitness;
end
[a,b]=sort(f,'ascend');   %��С����
cbest1=pop(b(1));
table=[0 0 0 0 0;...
       0 0 0 0 0;...
       0 0 0 0 0;...
       0.05*pi -1 1 round(rand*2-1) 0;
       0.01*pi -1 1 round(rand*2-1) 0;
       0.025*pi 1 -1 0 round(rand*2-1);
       0.005*pi 1 -1 0 round(rand*2-1);
       0.025*pi 1 -1 0 round(rand*2-1)];
[new,new1,k1]=niching(pop,LB,UB,6);     %newΪС����,new1Ϊ����
for i=1:k1+1
    cbest(i)=cbest1;
end
for iter=2:iter_max
    for i=1:k1
        [new{i},cbest(i),favg(i),cworst(i)]=rotationq(fun,new{i},cbest(i),table,pm,LB,UB);
    end
    [new1,cbest(end),favg(k1+1),cworst(k1+1)]=rotationq(fun,new1,cbest(end),table,pm,LB,UB);
    f_avg=zeros(k1,k1);   %��С����ƽ������ֵ�Ĳ�ֵ
    for i=1:k1
        for j=1:k1
            if i~=j
              f_avg(i,j)=abs(favg(i)-favg(j));
            end
        end
    end
    [y,idx]=trimatrix(f_avg,'min');
    [a,b]=sort([cbest(idx(1)).fitness cbest(idx(2)).fitness],'ascend');
    new{idx(b(end))}(cworst(idx(b(end))))=cbest(idx(b(1)));   
end
for i=1:k1+1
    cbest_var(i,:)=cbest(i).var;
    cbest_f(i)=cbest(i).fitness;
end







            

