function [max_x,maxfval]=myga1(fun,LB,UB,popsize,iterm_max,px,pm,eps1)
  %fun       �����Ż���Ŀ�꺯��
  %LB            ���Ա����½�
  %UB            ���Ա����Ͻ�
  %popsize       ����Ⱥ��С
  %iterm_max      ����������
  %px            ���ӽ�����
  %pm            ���������
  %eps1           ���Ա�����ɢ����
if isempty(pm)
    pm=0.1;
end
if isempty(px)
    px=0.90;
end
if isempty(iterm_max)
    iterm_max=8000;
end
if isempty(popsize)
    popsize=50;
end
if isempty(LB)&&isempty(UB)
    nvar=input('�����������Ŀnvar�� ');
else
    nvar=size(LB,1);     %������
end
CodeLen=nvar*max(ceil(log2((UB-LB)/eps1 + 1)));      %�����Ա�����ɢ���ȣ�ȷ�������Ʊ���λ���ĳ���
x=zeros(popsize,CodeLen);                      %��Ⱥ����ĳ�ʼֵ
for i=1:popsize
    x(i,:)=Initial(CodeLen);              %�����ӳ���Initial��ʼ����Ⱥ
    FitValue(i)=fun(Dec1(LB,UB,x(i,:),CodeLen));    %������Ӧֵ
end
for i=1:iterm_max
    SumFitValue=sum(FitValue);     %���и�����Ӧֵ֮��
    Ave_x=FitValue/SumFitValue;    %���и�����Ӧֵ��ƽ��ֵ
    Prob_Ave_x=0;
    Prob_Ave_x(1)=Ave_x(1);
    for j=2:popsize                    %�������̶Ĳ��Եĸ����ۼ�
        Prob_Ave_x(j)=Prob_Ave_x(j-1)+Ave_x(j);
    end
    for k=1:popsize
        sita=rand();
        for n=1:popsize
            if sita<=Prob_Ave_x(n)
                FatherSelection=n;   %�������̶Ĳ���ѡ�񸸴�
                break;
            end
        end
        MotherSelection=floor(rand()*(popsize-1))+1;      %���ȷ��ĸ��
        PosCutPoint=floor(rand()*(CodeLen-2))+1;     %���ȷ�����㽻����е�λ��
        r1=rand();
        if r1<=px     %���н������
            nx(k,1:PosCutPoint)=x(FatherSelection,1:PosCutPoint);
            nx(k,(PosCutPoint+1):CodeLen)=x(MotherSelection,(PosCutPoint+1):CodeLen);
            r2=rand();
            if r2<=pm %���б������
                PosMutPoint=round(rand()*(CodeLen-1)+1);   %���ȷ������Ԫ�ص�λ��
                nx(k,PosMutPoint)=~nx(k,PosMutPoint);
            end
        else
            nx(k,:)=x(FatherSelection,:);
        end
    end
    x=nx;
    for m=1:popsize
        FitValue(m)=fun(Dec1(LB,UB,x(m,:),CodeLen));    %�Ӵ�������Ӧֵ�������Ӻ���Dec
    end
    [a,b]=max(FitValue);
    best=x(b,:);
    x(popsize,:)=best;
end
max_x=Dec1(LB,UB,best,CodeLen);
maxfval=fun(max_x);


function result=Initial(length)           %��ʼ������
for i=1:length
    r=rand();
    result(i)=round(r);
end


function y=Dec1(LB,UB,x,CodeLen)           %������ת��Ϊʮ���Ʊ���
nvar=size(LB,1);
sublen=CodeLen/nvar;
pow_two=2.^(0:sublen)';
maxintval=((2^sublen))-1;
range=UB'-LB';
start=1;
fin=sublen;
for j=1:nvar
    tvars(1:sublen)=x(start:fin);
    start=start+sublen;
    fin=fin+sublen;
    temp1=0;
    for k=1:sublen
         temp1=temp1+pow_two(k)*tvars(sublen-k+1);
    end
    y(j)=(range(j)*(temp1/maxintval))+LB(j);
end