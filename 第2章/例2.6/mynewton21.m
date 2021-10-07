function [xmin,minf]=mynewton21(fun,phi,x0,x_syms,type,esp)   %ţ�ٷ����Ԫ������ֵ
syms lamda
if nargin<7
    esp=1e-6;
end
z=myjacobian1(phi,x_syms);
H=hessian(phi,x_syms);
n=length(x0);
if strcmp(type,'xz')
    tau=input('�������= ');   
    if isempty(tau)
        tau=0.1;
    end
end
k=0;
while k<5000
    f=eval(subs(z,x_syms,x0));
    h=eval(subs(H,x_syms,x0));
    s=-inv(h)*f';
    if norm(f)<=esp
         xmin=x0;
         if ~isempty(fun)
             minf=fun(xmin);
         else
             minf=eval(subs(phi,x_syms,xmin));
         end
         break
    end
    if strcmp(type,'nt')   %ţ�ٷ�
        xk=x0+s';
    elseif strcmp(type,'zn')    %����ţ�ٷ�
        if ~isempty(fun)
            x3=mysearch1(fun,phi,x0,x_syms,s','d');  %һά�����󲽳�����20��~��31���������������� 
        else
            x3=mysearch1([],phi,x0,x_syms,s','d');
        end
        xk=x0+x3*s';
    elseif strcmp(type,'ng')      %ţ�٣��ݶȷ�
        if f*s<0
            if ~isempty(fun)
               x3=mysearch1(fun,phi,x0,x_syms,s','d');  %һά�����󲽳�����20��~��31���������������� 
            else
               x3=mysearch1([],phi,x0,x_syms,s','d');
            end
           % x3=mysearch1(fun,phi,x0,x_syms,xsyms,s','d');
            xk=x0+x3*s';
        else
            if ~isempty(fun)
               x3=mysearch1(fun,phi,x0,x_syms,-f,'d');  %һά�����󲽳�����20��~��31���������������� 
            else
               x3=mysearch1([],phi,x0,x_syms,-f,'d');
            end
            %x3=mysearch1(fun,phi,x0,x_syms,xsyms,-f,'d');
            xk=x0-x3*f;
        end
    elseif strcmp(type,'xz')     %����ţ�ٷ�
        miuk=(norm(f))^(1+tau);
        %ak=h+miuk*eye(n);
        %s=-ak\f';       %������������
        s=-(inv(h+miuk*eye(n)))'*f';
        if ~isempty(fun)
            x3=mysearch1(fun,phi,x0,x_syms,s','d');  %һά�����󲽳�����20��~��31���������������� 
        else
            x3=mysearch1(fun,phi,x0,x_syms,s','d');
        end
       % x3=mysearch1(fun,phi,x0,x_syms,xsyms,s','d');
        xk=x0+x3*s';  
    end
    x0=xk;
    k=k+1;
end