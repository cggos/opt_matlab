function  [best_fval,best_x]=MainAneal(fun,LB,UB,t0,L)
if nargin==4
    L=30000;
elseif nargin==3
    L=30000;
    t0=1;    %��ʼ�¶�
end
rand('state',sum(clock));
f0=inf;
[n,c]=size(LB);    %nΪ������ά����cΪ������������ڷ���
for i=1:10
  x_temp=LB'+(UB'-LB').*rand(c,n);      %���ʼ��
  f=fun(x_temp);
  if f<f0
      f0=f;
      x0=x_temp;
  end
end
best_fval=f0;
best_x=x0;
t=t0;
for i=1:L
    for j=1:10
      f1=inf;
      for k=1:5
         x1=x0+unifrnd(-1,1,c,n);
         for k1=1:c
            x1(k1,:)=boundtest(x1(k1,:),LB(:,c),UB(:,c));   %����߽�
         end
         temp=fun(x1);
         if temp<f1
           f1=temp;
           x2=x1;
         end
      end
      df=f1-f0;
      if df<0
         x0=x2;
         f0=fun(x0);
      elseif exp(-df/t)>rand
         x0=x2;
         f0=fun(x0);
      end
      if f0<best_fval
         best_x=x0;
         best_fval=f0;
      end
    end
    if i<16
        t=t0/log(i);
    else
        t=t0*0.99^i;
    end
    if t<0.1^30  %����¶� 
         break
    end
end
