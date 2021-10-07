function [bestx,fval]=SAPSO(fun,popsize,c1,c2,lamda,iter_max,LB,UB)
D=size(LB,1);
for i = 1:popsize
    x(i,:) =LB'+(UB-LB)'.*rand(1,D);  %�����ʼ��λ��
    v(i,:) = randn(1,D);  %�����ʼ���ٶ�
    p(i)=fun(x(i,:));
    y(i,:)=x(i,:);
end
[a,b]=min(p);
pg=x(b,:);             %PgΪȫ������
fval=a;
T=fun(pg)/log(5);   %��3���������ʼ�¶�
for iter = 1:iter_max
    groupFit = fun(pg);    
    for i = 1:popsize          %��4������ǰ�¶�T�¸���pi����Ӧֵ
        Tfit(i)=exp(-(p(i)-groupFit)/T);       
    end    
    SumTfit=sum(Tfit);    
    Tfit=Tfit/SumTfit;    
    pBest=rand();    
    for i=1:popsize          %��5�����������̶Ĳ���ȷ��ȫ�����ŵ�ĳ�����ֵ        
        ComFit(i)=sum(Tfit(1:i));        
        if pBest<=ComFit(i)            
            pg_plus=x(i,:);            
            break;            
        end        
    end    
    C=c1+c2;    
    ksi=2/abs(2-C-sqrt(C^2-4*C));     %�ٶ�ѹ������   
    for i=1:popsize                               %���¸�΢�����ٶȺ�λ��
        v(i,:)=ksi*(v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg_plus-x(i,:)));
        x(i,:)=x(i,:)+v(i,:);
        if fun(x(i,:))<p(i)   %���¸�΢����piֵ��Ⱥ���pgֵ
            p(i)=fun(x(i,:));
            y(i,:)=x(i,:);
        end
        if p(i)<fval
            pg=y(i,:);
            fval=p(i);
        end
    end
    T=T*lamda;       %����   
end                 
bestx=pg;
