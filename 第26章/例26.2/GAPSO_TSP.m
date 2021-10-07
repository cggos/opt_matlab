function [MinDistance,Path] = GAPSO_TSP(city,popsize,MaxGen)
[n,d]=city2d(city); %n��ʾ����Ĺ�ģ�����и�����     
individual=zeros(popsize,n);
for i=1:popsize
    individual(i,:)=randperm(n);    %�����������λ��
    indiFit(i,1)=value(individual(i,:),d);
end
[value1,index]=min(indiFit);
tourPbest=individual;                                    %��ǰ��������
tourGbest=individual(index,:) ;                          %��ǰȫ������
xnew1=individual;
trace=zeros(1,MaxGen);   %��Ÿ�������·���ĳ��ȼ�����·����ƽ������
for N = 1:MaxGen
    for i = 1:popsize
       % ��������Ž��н���
        c1 = unidrnd(n-1);    %��������λ
        c2 = unidrnd(n-1);    %��������λ
        while c1 == c2
            c1 = round(rand*(n-2))+1;
            c2 = round(rand*(n-2))+1;
        end
        chb1 = min(c1,c2);
        chb2 = max(c1,c2);
        cros = tourPbest(i,chb1:chb2);
        ncros = size(cros,2);      
        %ɾ���뽻��������ͬԪ��
        for j = 1:ncros
            for k = 1:n
                if xnew1(i,k) == cros(j)
                    xnew1(i,k) = 0;
                    for t = 1:n-k
                        temp = xnew1(i,k+t-1);
                        xnew1(i,k+t-1) = xnew1(i,k+t);
                        xnew1(i,k+t) = temp;
                    end
                end
            end
        end
        %���뽻������
        xnew1(i,n-ncros+1:n) = cros;
        %��·�����ȱ�������
        Dist = 0;
        for j=1:n-1
            Dist = Dist+d(xnew1(i,j),xnew1(i,j+1));
        end
        Dist = Dist+d(xnew1(i,1),xnew1(i,n));
        if indiFit(i) > Dist
            individual(i,:)=xnew1(i,:);
            indiFit(i,1)=Dist;
        end
        
        % ��ȫ�����Ž��н���
        c1 = round(rand*(n-2))+1;  %��������λ
        c2 = round(rand*(n-2))+1;  %��������λ
        while c1 == c2
            c1 = round(rand*(n-2))+1;
            c2 = round(rand*(n-2))+1;
        end
        chb1 = min(c1,c2);
        chb2 = max(c1,c2);
        cros = tourGbest(chb1:chb2); 
        ncros = size(cros,2);      
        %ɾ���뽻��������ͬԪ��
        for j = 1:ncros
            for k = 1:n
                if xnew1(i,k) == cros(j)
                    xnew1(i,k) = 0;
                    for t = 1:n-k
                        temp = xnew1(i,k+t-1);
                        xnew1(i,k+t-1) = xnew1(i,k+t);
                        xnew1(i,k+t) = temp;
                    end
                end
            end
        end
        %���뽻������
        xnew1(i,n-ncros+1:n) = cros;
        %��·�����ȱ�������
        Dist = 0;
        for j = 1:n-1
            Dist = Dist+d(xnew1(i,j),xnew1(i,j+1));
        end
        Dist = Dist+d(xnew1(i,1),xnew1(i,n));
        if indiFit(i) > Dist
            individual(i,:) = xnew1(i,:);
            indiFit(i,1)=Dist;
        end
        
       % �������
        c1 = round(rand*(n-1))+1;   %��������λ
        c2 = round(rand*(n-1))+1;   %��������λ
        while c1 == c2
            c1 = round(rand*(n-2))+1;
            c2 = round(rand*(n-2))+1;
        end
        temp = xnew1(i,c1);
        xnew1(i,c1) = xnew1(i,c2);
        xnew1(i,c2) = temp;
        
        %��·�����ȱ�������
        Dist = 0;
        for j = 1:n-1
            Dist = Dist + d(xnew1(i,j),xnew1(i,j+1));
        end
        Dist = Dist + d(xnew1(i,1),xnew1(i,n));      %���㻷��·���ĳ���
        if indiFit(i) > Dist
            individual(i,:) = xnew1(i,:);
            indiFit(i,1)=Dist;
        end
    end
    [value2,index] = min(indiFit);
    if value2<value1
       trace(2,N) = indiFit(index);    %����·���ĳ���
       tourGbest = individual(index,:); 
       value1=value2;
    end  
end
MinDistance = value1;
Path = tourGbest;
TSPplot(city,Path);