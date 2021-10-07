function chome=gaussmution(chome,chome_var,LB,UB,type)    %高斯变异,结构体形式
popsize=size(chome,2); 
nvar=size(LB,1);   %变量的个数
tao1=sqrt(2*nvar)^-1;
tao2=sqrt(2*sqrt(nvar))^-1;
chomeLamda=(UB-LB)'./2;
if type==1    %标准 
   for j=1:popsize
       for k=1:nvar
           chome(j).x(k)=chome(j).x(k)+normrnd(0,chome(j).sigma,1,1);
           chome(j).x(k)=boundtest(chome(j).x(k),LB(k),UB(k));
       end
   end
elseif type==2     %自适应标准
   for j=1:popsize
       a=randn;
       for k=1:nvar
           chome(j).x(k)=chome(j).x(k)+randn*chome_var(k);
           chome(j).x(k)=boundtest(chome(j).x(k),LB(k),UB(k));
           chome_var(k)=chome_var(k)*exp(tao1*a+tao2*randn);
       end
   end
elseif type==3     %单点变异
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