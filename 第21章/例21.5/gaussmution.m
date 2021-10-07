function chome=gaussmution(chome,chome_var,LB,UB,type)    %��˹����,�ṹ����ʽ
popsize=size(chome,2); 
nvar=size(LB,1);   %�����ĸ���
tao1=sqrt(2*nvar)^-1;
tao2=sqrt(2*sqrt(nvar))^-1;
chomeLamda=(UB-LB)'./2;
if type==1    %��׼ 
   for j=1:popsize
       for k=1:nvar
           chome(j).x(k)=chome(j).x(k)+normrnd(0,chome(j).sigma,1,1);
           chome(j).x(k)=boundtest(chome(j).x(k),LB(k),UB(k));
       end
   end
elseif type==2     %����Ӧ��׼
   for j=1:popsize
       a=randn;
       for k=1:nvar
           chome(j).x(k)=chome(j).x(k)+randn*chome_var(k);
           chome(j).x(k)=boundtest(chome(j).x(k),LB(k),UB(k));
           chome_var(k)=chome_var(k)*exp(tao1*a+tao2*randn);
       end
   end
elseif type==3     %�������
   for j=1:popsize
       b=ceil(nvar*rand);
       if chomeLamda(b)<1e-4;
          chomeLamda(b)=(UB(b)-LB(b))/2;
       end
       chome(j).x(b)=chome(j).x(b)+chomeLamda(b)*randn;
       chome(j).x(b)=boundtest(chome(j).x(b),LB(b),UB(b));
       chomeLamda(b)=chomeLamda(b)*exp(-1.01);
   end
end