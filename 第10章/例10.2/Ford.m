function [f,nf,wf,nwf,No]=Ford(C,Flow)   %ͼ���е�Ford�㷨��������
   %f         �������������
   %nf        �������������
   %wf         ���������
   %nwf        ���������
   %No                ����ź������ɵ���С��
   %C                 �����ϵ�����
   %Flow              �����ϵ�������������Ĭ�������Ϊ0

n=size(C,2);
if nargin == 1   
    f=zeros(n,n);
else
    f=Flow;
end

%��2�����б�Ź�������
No=zeros(1,n);
d=zeros(1,n);
while 1
    No(1)=n + 1;
    d(1)=inf;         %��Դ���
    while 1
        pd=1;         %��Ź���
        for i=1:n
            if No(i)    %ѡ��һ���ѱ�ŵĵ�
                for j=1:n
                    if No(j)==0 && f(i,j)<C(i,j)  %����δ��ŵĵ�j��������i,j���Ǳ���ʱ
                        No(j)=i;
                        d(j)=C(i,j)-f(i,j);
                        pd=0;
                        if d(j)>d(i)
                            d(j)=d(i);
                        end
                    elseif No(j)==0 && f(j,i)>0    %����δ��ŵĵ�j��������i,j��������ʱ
                        No(j)=-i;
                        d(j)=f(j,i);
                        pd=0;
                        if d(j)>d(i)
                            d(j)=d(i);
                        end
                    end
                end
            end
        end
        if (No(n)||pd)
            break;
        end
    end                 %����õ���Ż��޷����ʱ������ֹ��Ź���
    if (pd)
        break;
    end                 %����δ�õ���ţ�Flow�Ѿ�Ϊ����������㷨��ֹ
    %���������������̣�dvt��ʾ������
    dvt=d(n);
    t=n;
    while 1
        if No(t)> 0
            f(No(t),t)=f(No(t),t) + dvt;     %ǰ�򻡵���
        elseif No(t)<0
            f(-No(t),t)=f(-No(t),t) - dvt;   %���򻡵���
        end
        if No(t)==1
            for i=1:n
                No(i)=0;
                d(i)=0;
            end
            break
        end                   %��t�ı��ΪԴʱ������ֹ��������
        t=No(t);
    end
end                           %��������ǰһ�λ��ϵ���Flow
%��3���������������
if nargin==1
   nf=f;
else
    nf=f-Flow;  %�����������
end
%��4�������������
wf=0;
nwf=0;
for j=1:n
    wf=wf+f(1,j);
    nwf=nwf+nf(1,j);
end