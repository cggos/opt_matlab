function [xmin,minf]=sumt1(fun,phi,Pphi,x0,x_syms,esp,type)   %罚函数法求解约束最优化
%%fun为目标函数，phi为包括目标函数、不等式和等式约束在内的罚函数，x0为初始点，
%phi为含有罚因子在内的符号变量表达式,Pphi为罚函数的表达式
if nargin==4
    type='Of'; 
    esp=1e-6;
end
syms m
M=10;C=10;
n=length(x_syms);
if strcmp(type,'Of') ||strcmp(type,'If')    % 函数法
   while 1
      phi1=subs(phi,m,M);
      xmin=DFP1([],phi1,x0,x_syms,esp);
      P=subs(Pphi,m,M);
      if eval(subs(P,x_syms,xmin))<=esp
          break
      else
          if strcmp(type,'Of')       %外点法
              M=C*M;
          elseif strcmp(type,'If')   %内点法
              M=M/C;
          end
      end
   end
elseif strcmp(type,'Oj')||strcmp(type,'Ij')
    z=myjacobian1(phi,x_syms);
    ans1=solve(z,x_syms);
    ans2=struct2cell(ans1);
    for i=1:n
       if strcmp(type,'Oj')
           xmin(i)=limit(ans2{i}(1),m,inf);
       elseif strcmp(type,'Ij')
           xmin(i)=limit(ans2{i}(1),m,0);
       end
    end
end
if isa(fun,'function_handle')
   minf=fun(xmin);
elseif isa(fun,'sym')
    minf=eval(subs(fun,x_syms,xmin));
end












    


