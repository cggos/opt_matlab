function varargout=interval1(opfun,x0,lamda,type)  %进退法求区间
%当用导数时，函数式为原函数的导数
switch type
    case 'f'   %函数比较
        if isa(opfun,'function_handle')
            phai0=opfun(x0);
        elseif isa(opfun,'sym')
            phai0=eval(subs(opfun,x0));
        end   
        t2=x0+lamda;
        if isa(opfun,'function_handle')
            phai2=opfun(t2);
        elseif isa(opfun,'sym')
            phai2=eval(subs(opfun,t2));
        end   
        if phai2>phai0
            lamda=-lamda;
            t1=x0+lamda;
            if isa(opfun,'function_handle')
               phai1=opfun(t1);
            elseif isa(opfun,'sym')
               phai1=eval(subs(opfun,t1));
           end   
        else
            t1=x0+lamda;
            if isa(opfun,'function_handle')
               phai1=opfun(t1);
            elseif isa(opfun,'sym')
               phai1=eval(subs(opfun,t1));
            end   
        end       
        while phai1<=phai0
            lamda=2*lamda;
            t2=x0;
            x0=t1;
            phai0=phai1;
            t1=x0+lamda;
            if isa(opfun,'function_handle')
               phai1=opfun(t1);
            elseif isa(opfun,'sym')
               phai1=eval(subs(opfun,t1));
            end   
        end
        a=min(t1,t2);
        b=max(t1,t2);    
    case 'd'  %导数比较
        if isa(opfun,'sym')
            opfun1=diff(opfun);
        end
        if isa(opfun,'function_handle')
            phai=opfun(x0);
        elseif isa(opfun,'sym')
            phai=eval(subs(opfun1,x0));
        end   
        if phai==0
            value=x0;
        elseif phai<0
              x1=x0+lamda;
              if isa(opfun,'function_handle')
                  phai1=opfun(x1);
              elseif isa(opfun,'sym')
                  phai1=eval(subs(opfun1,x1));
              end   
              if phai1<0
                  while phai1<0
                      x2=x1+2*lamda;
                      if isa(opfun,'function_handle')
                         phai2=opfun(x2);
                      elseif isa(opfun,'sym')
                         phai2=eval(subs(opfun1,x2));
                      end   
                      if phai2>0
                         a=x1;
                         b=x2;
                         break
                      else
                          x1=x2;
                      end
                  end  
              elseif phai1==0
                  value=x1;
              else
                  a=x0;
                  b=x1;
              end
        else
              x1=x0-lamda;
              if isa(opfun,'function_handle')
                  phai1=opfun(x1);
              elseif isa(opfun,'sym')
                  phai1=eval(subs(opfun1,x1));
              end   
              if phai1>0
                  while phai1>0
                      x2=x1-2*lamda;
                      if isa(opfun,'function_handle')
                         phai2=opfun(x2);
                      elseif isa(opfun,'sym')
                         phai2=eval(subs(opfun1,x2));
                      end   
                      if phai2<0
                         a=x1;
                         b=x2;
                         break
                      else
                          x1=x2;
                      end
                  end  
              elseif phai1==0
                  value=x1;
              else
                  a=x1;
                  b=x0;
              end
        end        
end
if exist('value','var')
    varargout={value};
elseif nargout==1
    varargout={[a,b]};
else
    varargout={a,b};
end

