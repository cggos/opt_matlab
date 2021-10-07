function [x0,minf]=interpolation1(fun,phi,varargin)    %二点三次插值法求极值
if nargin<5
    esp=1e-6;
end
if nargin==4
    alpha=varargin{2};   %初始步长
    esp=1e-6;
elseif nargin==5
    alpha=varargin{2};
    esp=varargin{3};
end
df=mydiff1(phi);
switch length(varargin{1})
    case 2
      a=varargin{1}(1);b=varargin{1}(2);
      if ~isempty(fun)
          fa=fun(a);
      else
          fa=eval(subs(phi,a));
      end
      dfa=eval(subs(df,a));
      if ~isempty(fun)
          fb=fun(b);
      else
          fb=eval(subs(phi,b));
      end
      dfb=eval(subs(df,b));
      n=1;
      while 1
        if n>10000
             error('此区间内没有极值')
        end
        if abs(a-b)<=esp
            x0=(a+b)/2;
            break
        else
            u=(fb-fa)/(b-a)-dfa;
            v=dfb-dfa;
            beta=(3*u-v)/(b-a);
            alpha=(v-2*u)/(b-a)^2;
            if beta^2-3*alpha*dfa<0
                error('此区间内没有极值')
            end
            xp=a-dfa/(beta+sqrt(beta^2-3*alpha*dfa));
            fp=eval(subs(df,xp));
            if fp>0
                b=xp;
                fb=fp;
            else
                a=xp;
                fa=fp;
            end
        end
        n=n+1;
     end
    case 1
        x1=varargin{1};
        if ~isempty(fun)
          f1=fun(x1);
        else
          f1=eval(subs(phi,x1));
        end
        df1=eval(subs(df,x1));
        if abs(df1)<=esp
            x0=x1;
        else
            n=1;dfu=1;
            while abs(dfu)>esp
                n=n+1;
                if n>10000
                    error('没有极值');
                end
                if df1>0
                   alpha=-abs(alpha);
                else
                   alpha=abs(alpha); 
                end
                x2=x1+abs(alpha);
                if ~isempty(fun)
                   f2=fun(x2);
                else
                   f2=eval(subs(phi,x2));
                end
                df2=eval(subs(df,x2));
                if abs(df2)<=esp
                    x0=x2;
                    break
                else
                    m=1;
                    while df1*df2>0
                      m=m+1;
                      if m>10000
                          error('没有极小');
                      end
                      alpha=2*alpha;
                      x1=x2;
                      f1=f2;
                      df1=df2;
                      x2=x1+alpha;
                      if ~isempty(fun)
                         f2=fun(x2);
                      else
                         f2=eval(subs(phi,x2));
                      end
                      df2=eval(subs(df,x2));
                      if abs(df2)<=esp
                        x0=x2;
                        break
                      end
                    end
                    k=exist('x0','var');
                    if k==0
                        w=3*(f2-f1)/(x2-x1);
                        v=w-df1-df2;
                        u=sign(x2-x1)*sqrt(v^2-df1*df2);
                        x=x1+(u-v-df1)*(x2-x1)/(df2-df1+2*u);
                        if ~isempty(fun)
                           f=fun(x);
                        else
                           f=eval(subs(phi,x));
                        end
                        dfu=eval(subs(df,x));
                        if abs(dfu)<=esp
                           x0=x;
                           break
                        else
                           alpha=alpha/10;
                           x1=x;
                           f1=f;
                           df1=dfu;
                        end
                    else
                        break
                    end
                end
            end
        end      
end
if ~isempty(fun)
   minf=fun(x0);
else
   minf=eval(subs(phi,x0));
end

            
        