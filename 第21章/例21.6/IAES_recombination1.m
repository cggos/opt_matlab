function [new,newsigma]=IAES_recombination1(oldchome,oldsigma,lamda,type)
num=size(oldchome,2);
nvar=length(oldchome(1).x);
temp=ceil(num*rand);   %选择一个父代
for i=1:lamda
    temp1=ceil(num.*rand(1,2));
    switch type
       case 1   %离散重组          
           mask=round(rand(1,nvar));
           if any(mask-1)||any(mask)
               mask=round(rand(1,nvar));
           end
            for j=1:nvar  %变量数目
               if mask(j)==0
                   new(i).x(j)=oldchome(temp1(1)).x(j);
                   newsigma(i,j)=oldsigma(temp1(1),j);
               else
                   new(i).x(j)=oldchome(temp1(2)).x(j);
                   newsigma(i,j)=oldsigma(temp1(2),j);
               end
            end 
    case 2   %中间重组
        for j=1:nvar
            new(i).x(j)=(oldchome(temp1(1)).x(j)+oldchome(temp1(2)).x(j))/2;
            newsigma(i,j)=(oldsigma(temp1(1),j)+oldsigma(temp1(2),j))/2;
        end    
    case 3   %混杂重组
       temp2=ceil(num*rand);
       while temp2==temp
           temp2=ceil(num*rand);
       end
       for j=1:nvar
           new(i).x(j)=(oldchome(temp).x(j)+oldchome(temp2).x(j))/2;
           newsigma(i,j)=(oldsigma(temp,j)+oldsigma(temp2,j))/2;
       end     
    end      
end