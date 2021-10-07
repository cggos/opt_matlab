function [xmin,minf]=sumt1(fun,phi,Pphi,x0,x_syms,esp,type)   %�����������Լ�����Ż�
%%funΪĿ�꺯����phiΪ����Ŀ�꺯��������ʽ�͵�ʽԼ�����ڵķ�������x0Ϊ��ʼ�㣬
%phiΪ���з��������ڵķ��ű������ʽ,PphiΪ�������ı��ʽ
if nargin==4
    type='Of'; 
    esp=1e-6;
end
syms m
M=10;C=10;
n=length(x_syms);
if strcmp(type,'Of') ||strcmp(type,'If')    % ������
   while 1
      phi1=subs(phi,m,M);
      xmin=DFP1([],phi1,x0,x_syms,esp);
      P=subs(Pphi,m,M);
      if eval(subs(P,x_syms,xmin))<=esp
          break
      else
          if strcmp(type,'Of')       %��㷨
              M=C*M;
          elseif strcmp(type,'If')   %�ڵ㷨
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












    


