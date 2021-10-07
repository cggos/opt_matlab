function [k,C]=colorcodf(W,type)      %图的染色问题
%k为染色数，C为染色方案
if nargin==1
    type='p';
end
if strcmp(type,'p')    %点染色
    G=W;
    n=size(G,1);
    k=1;
    C=zeros(1,n);
    Z=1:n;
    while sum(find(C==0))
       tcol=find(C==0);
       m=sum(G(tcol,:),2);
       minm=min(m);
       k1=min(find(m==minm));
       c=G(tcol(k1),:);
       c(1,tcol(k1))=1;
       C(tcol(k1))=k;
       Sn=find(c~=0);
       flag=1;
       while flag
          tc=setdiff(Z,Sn);
          if isempty(tc)
             flag=0;
             k=k+1;
          else
             c=G(tc(1),:);
             c(1,tc(1))=1;
             C(tc(1))=k;
             Sn1=find(c~=0);
             Sn=union(Sn,Sn1);
          end
       end
       trow=find(C==k-1);
       G(:,trow)=1;
    end
    k=k-1;
elseif strcmp(type,'s')     %边染色
    [G,id]=graphpts(W);
    [k,C1]=colorcodf(G);
    C(1,:)=C1;
    C(2:3,:)=id';
elseif isnumeric(type)     %二分图的边染色，此时type为染色数
    [G,id]=graphpts(W);
    G1=G;
    n1=sum(G1,2);    %各点的权值
    num=size(G,1);       %边的数目
    k1=max(sum(W));           %最小的染色数
    if type<k1
        error('染色数太小，请重新输入染色数');
    end
    n=1:num;
    n2=floor(num/type);
    n4=rem(num,type);
    c=n2.*ones(type,1);
    c(1:n4)=c(1:n4)+1;
    c2=cell(1,type);
    for i=1:type
       [a,b]=sort(n1);
       G1=G1(b,:);
       n=n(b);
       c1=n(1);m=1;
       for j=num:-1:2
           if j~=n(1)
             if i~=type
                 if length(c1)==c(i)
                     break  
                 else
                    if G(n(1),n(j))==0&&all(G(c1(1:end),n(j))~=1)
                      c1=[c1 n(j)]; 
                      m=[m j];
                    end
                 end
             else
                if G(n(1),n(j))==0&&all(G(c1(1:end),n(j))~=1)
                    c1=[c1 n(j)]; 
                    m=[m j];
                end 
             end
           end
       end
       c2{i}=c1;
       c3(i)=length(c1);    %实际取值的长度
       G1=redu(G1,m,'r');
       G1=redu(G1,m,'c');
       n=redu(n,m,'c');
       num=size(G1,2);
       n1=sum(G1,2);
       flag=0;
    end
    while 1
        if isempty(n)
            break
        end
        if flag==0
           a=find((c3(1:end-1)-c(1:end-1)')~=0);
        elseif  flag==1
            a=find(c3-c'~=0);
        end
        if ~isempty(a)     %实际取边数少于理论值  
            for i=1:length(a)
                if all(G(c2{a(i)},n)==0)
                    c2{a(i)}=[c2{a(i)} n];
                 end
            end
         else
             for i=1:length(n)      %有多少个没有归位的边
                for j=1:type-1
                     if all(G(c2{j},n(i))==0)   %与所有的边都没有矛盾
                         b=randperm(length(c2{j}));
                         temp=c2{j}(b(1));
                         c2{j}(b(1))=n(i);
                         c2{type}=[c2{type} temp];
                     end
                end
             end
         end
         n=[];flag=1;
         for i=1:type
             for j=1:c(i)
                 if (G(c2{i}(1),c2{i}(j))~=0)
                     c3(i)=redu(c3(i),j,'c');
                     c2{i}=redu(c2{i},j,'c');
                     n=[n c2{i}(j)];
                 end
             end
         end
    end
    k=type;
    for i=1:type
        for j=1:length(c2{i})
           id(c2{i}(j),3)=i;
        end
    end
    C=id;
    C(:,2)=C(:,2)-C(end,1);
    C=C';
end
    

