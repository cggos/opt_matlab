function out=mymultifun1(phi,x_range,x0,x_syms,xsyms,type)   %求多元函数的极值
syms lamda
switch type
    case 'u_L'  %根据特征值求极值
        num=length(xsyms);
        for i=1:num
          d{i}=diff(phi,xsyms{i});
          for j=1:num
            H(i,j)=diff(d{i},xsyms{j});
          end
        end
        z=det(H-lamda.*eye(num));
        lamda1=solve(z,lamda);
        fun=[];
        for i=1:num
            fun=[fun d{i}];
        end
        x0=solve(fun);   %驻点
        names=fieldnames(x0)';
        for i=1:num
           x_value(:,i)=getfield(x0,names{i});
        end
        n=size(x_value,1);
        for i=1:n
           for j=1:num
              lam(i,j)=subs(lamda1(j),names,x_value(i,:));
              if eval(lam(i,j))>0
                xx(i,j)=1;
              elseif eval(lam(i,j))<0
                xx(i,j)=-1;
              else
               xx(i,j)=0;
              end   
           end
       end
       k=0;
       for i=1:n
          if any(xx(i,:)-ones(1,num))==0
             k=k+1;
             out.x(k,:)=eval(x_value(i,:));
             out.pb{k}='min';
             out.value(k)=eval(subs(phi,names,out.x(k,:)));
          elseif any(xx(i,:)+ones(1,num))==0
             k=k+1;
             out.x(k,:)=eval(x_value(i,:));
             out.pb{k}='max';
             out.value(k)=eval(subs(phi,names,out.x(k,:)));
          end
       end
    case 'u_c'   %求条件极值
       x_syms1=x_syms{1};
       x_syms2=x_syms{2};
       num1=length(x_syms1);
       num2=length(x_syms2);
       num=num1+num2;
       for i=1:num1
          d{i}=diff(phi,x_syms1{i});
       end
       for i=1:num2
         d{i+num1}=diff(phi,x_syms2{i});
       end
       for i=1:num1
         for j=1:num1
           H(i,j)=diff(d{i},x_syms1{j});
         end
       end
       z=det(H-lamda.*eye(num1));
       lamda1=solve(z,lamda);
       fun=[];
       for i=1:num
           fun=[fun d{i}];
       end
       x0=solve(fun);   %驻点
       names=fieldnames(x0)';
       for i=1:num
          x_value(:,i)=getfield(x0,names{i});
       end
       n=size(x_value,1);
       for i=1:n
         for j=1:num1
           lam(i,j)=subs(lamda1(j),names,x_value(i,:));
           if eval(lam(i,j))>0
              xx(i,j)=1;
           elseif eval(lam(i,j))<0
              xx(i,j)=-1;
           else
              xx(i,j)=0;
           end   
         end
       end
       k=0;
       for i=1:n
         if any(xx(i,:)-ones(1,num1))==0
           k=k+1;
           out.x(k,:)=eval(x_value(i,:));
           out.pb{k}='min';
           out.value(k)=eval(subs(phi,names,out.x(k,:)));
         elseif any(xx(i,:)+ones(1,num1))==0
           k=k+1;
           out.x(k,:)=eval(x_value(i,:));
           out.pb{k}='max';
           out.value(k)=eval(subs(phi,names,out.x(k,:)));
         end
       end
    case 'u_u'
      df=myjacobian1(phi,x_syms);
      num=length(x_syms);
      for i=1:num
         s{i}=char(df(i));
         if findstr(s{i},'*')
           s{i}=strrep(s{i},'*','.*');
         end
         if findstr(s{i},'^')
           s{i}=strrep(s{i},'^','.^');
         end
         if findstr(s{i},'/')
           s{i}=strrep(s{i},'/','./');
         end
         if findstr(s{i},'exp')
           s{i}=strrep(s{i},'exp','l');
         end
         for j=1:num
           if findstr(s{i},xsyms{j})
             s{i}=strrep(s{i},xsyms{j},['x(',num2str(j),')']);
           end
         end
         if findstr(s{i},'l')
           s{i}=strrep(s{i},'l','exp');
         end    
      end     
      fun=['[',s{1},','];
      for i=2:num-1
          fun=strcat(fun,s{i},',');
      end
      fun=strcat(fun,s{num},']');
      options=optimset('fsolve');
      options=optimset('display','off');
      options.TolFun=1e-8;
      n=input('区间等分数=');
      dx(1:num)=(max(x_range)-min(x_range))/n;
      npoint=size(x0,1);
      for k=1:npoint
        xk=fsolve(fun,x0(k,:),options); % 求解方程组
        for j=1:num
          out.x(k,j)=xk(j);
          xk1(j)=xk(j)+dx(j);
        end
        if (subs(phi,x_syms,xk)>subs(phi,x_syms,xk1))
           out.pb{k}='max';
           out.value{k}=eval(subs(phi,x_syms,xk));
        else
          out.pb{k}='min';
          out.value{k}=eval(subs(phi,x_syms,xk));
        end
      end
end






      