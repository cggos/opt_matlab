function [f,wf,zwf]=Busacker(C,b)    %ͼ���е�Busacker_Gowan�����㷨
%====�������====
   %f     ����С�������������
   %wf              ���������
   %zwf              ����С����
%====�������====
   %C                    ��������
   %b                    �����ϵ�λ�����ķ���
n=size(C,2);                    %ͼ�Ľڵ����
wf=0;wf0=Inf;      %�������MaxFlow��Ԥ��������ֵMaxFlow0
f=zeros(n,n);    %ȡ��ʼ������Ϊ��������
while 1
    a=ones(n,n)*inf;
    for i=1:n
        a(i,i)=0;
    end
    %��1������������Ȩͼ
    for i=1:n
        for j=1:n
            if(C(i,j)>0 &&f(i,j)==0)
                a(i,j)=b(i,j);
            elseif(C(i,j)>0 && f(i,j)==C(i,j))
                a(j,i)=-b(i,j);
            elseif(C(i,j)>0)
                a(i,j)=b(i,j);
                a(j,i)=-b(i,j);
            end
        end
    end
    for i=2:n
        p(i)=inf;
        s(i)=i;
    end
    %��2���������·����������Ȩͼ����㵽�յ�����·
    for k=1:n
        pd=1;      
        for i=2:n
            for j=1:n
                if(p(i)>p(j)+a(j,i))
                    p(i)=p(j)+a(j,i);
                    s(i)=j;
                    pd=0;
                end
            end
        end
        if(pd)
            break;
        end;
    end
    
    if(p(n)==inf)
        break;
    end          %��������㵽�յ�����·, �㷨��ֹ��
    dvt=inf;
    t=n;
    %��3���������������,dvt��ʾ������
    while 1    %���������
        if(a(s(t),t)>=0)
            dvtt = C(s(t),t)-f(s(t),t); %ǰ�򻡵�����
        elseif(a(s(t),t)<0)
            dvtt=f(t,s(t));
        end                                            %���򻡵�����
        if(dvt>dvtt)
            dvt=dvtt;
        end
        if(s(t)==1)
            break;
        end      %��t�ı��Ϊ���ʱ, ��ֹ���������
        t=s(t);
    end
    %��������ǰһ�λ��ϵ���f
    pd=0;
    if(wf+dvt>=wf0)
        dvt=wf0-wf;
        pd=1;
    end
    t=n;
    while 1
        if(a(s(t),t)>=0)
            f(s(t),t)=f(s(t),t)+dvt;     %ǰ�򻡵���
        elseif(a(s(t),t) < 0)
            f(t,s(t))=f(t,s(t))-dvt;
        end                                                              %���򻡵���
        if(s(t)==1)
            break;
        end   %��t�ı��Ϊ���ʱ, ��ֹ��������
        t=s(t);
    end
    if(pd)
        break;
    end       %�����������ﵽԤ��������ֵ
    %�����������
    wf=0;
    for j=1:n
        wf=wf+f(1,j);
    end;
end
%������С����
zwf=0;
for i=1:n
    for j=1:n
        zwf=zwf+b(i,j)*f(i,j);
    end
end