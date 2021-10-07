function [y,x]=mysearch1(fun,phi,x0,x_syms,d0,type)    %不精确搜索求步长,d0为方向,x_syms为变量
if isempty(d0)
    z=myjacobian1(phi,x_syms);
    g=subs(z,x_syms,x0);
    d0=-g;
end
switch type
    case 'd'   %直接法
        z=myjacobian1(phi,x_syms);
        a=0;b=inf;lamda=1;y=lamda;m=0;c1=0.1;c2=0.5;
        if ~isempty(fun)
            f=fun(x0);
        else
            f=eval(subs(phi,x_syms,x0));
        end
        g=subs(z,x_syms,x0);
        t=g*d0';
        while 1
           if m>20
               break
           end
           x1=x0+lamda*d0;
           if ~isempty(fun)
              f1=fun(x1);
           else
             f1=eval(subs(phi,x_syms,x1));
           end
           %f1=fun(x1);
           g1=subs(z,x_syms,x1);
           t1=g1*d0';
           a1=f-f1;
           a2=-c1*lamda*t;
           a3=c2*t;
           if (a1>=a2)&&(t1>=a3)
               y=lamda;
               x=x1;
               break;
           elseif (a1>=a2)&&t1<a3
               a=lamda;
               lamda=min(2*lamda,(lamda+b)/2);
           elseif a1<a2
               b=lamda;
               lamda=(lamda+a)/2;
           end
           m=m+1;
        end
    case 'a'   %Armijo规则
        z=myjacobian1(phi,x_syms);
        m=0;maxm=20;beta=0.5;sigma=0.2;
        y1=0;
        while m<=maxm
            if ~isempty(fun)
              x1=fun(x0+beta^m*d0);
            else
              x1=eval(subs(phi,x_syms,x0+beta^m*d0));
            end
          %x1=fun(x0+beta^m*d0);
            if ~isempty(fun)
               x2=fun(x0)+sigma*beta^m*sbs(z,x_syms,x0)*d0';
            else
               x2=eval(subs(phi,x_syms,x0))+sigma*beta^m*sbs(z,x_syms,x0)*d0';
            end
         %  x2=fun(x0)+sigma*beta^m*sbs(z,x_syms,x0)*d0';
           if x1<=x2
               y1=m;  
               break;
           end
           m=m+1;
        end
        y=beta^y1;
        x=x0+y*d0;
    case 'equ'     %解方程求步长
        syms al x 
        x1=x0+al*d0;      
        phi1=subs(phi,x_syms,x1);
        phi=(subs(phi1,al,'x'));
        ds=diff(phi);
        x2=solve(ds);   %解方程求步长
        for i=1:length(x2)
           if isreal(x2(i))
              y1=mymin(phi,x2(i),'min');
                if y1==1
                  y=eval(x2(i));
                  break
                end 
           end
        end
        if exist('y','var')==1
            x=x0+y*d0;
        else
            y=[];x=[];
        end
    case 's'    %一维搜索法
        syms al x
        x1=x0+al*d0;
        phi1=(subs(phi,x_syms,x1));
        phi=subs(phi1,al,'x');
        y=goldcut1(phi,1);
        x=x0+y*d0;  
end



       
        