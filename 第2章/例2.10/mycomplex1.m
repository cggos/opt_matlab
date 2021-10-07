function [minx,minf]=mycomplex1(fun,phi,x0,type,esp)   %���η�������ֵ,phiΪԼ������,x0Ϊ��ʼ�㣬��Ԫ����ʽ����
%%��һ�ַ��������Ϊ2*varnum���㣬�ڶ��ַ���Ϊn+1����,phiԼ������
if nargin==5
    esp=1e-6;
end
switch type
    case 1
         y=fun(x0);
         if isempty(phi)
             flag=0;
         else
             flag=1;
         end
         while 1
            alpha=1.3;
            [a,b]=max(y);
            xl=x0{b};  %���
            b1=redu(x0,b,'c');
            a1=mean(cell2mat(b1')); 
            if alpha*norm(a1-xl)<=1e-8  %ֹͣ����һ
                [minf,b]=min(y);
                minx=x0{b}; 
                break;
            end
            while 1
               x_alpha=a1+alpha*(a1-xl);     %���µĵ�
               if flag==1
                  y1=phi(x_alpha);   %Լ������Ҫ����ʵ�ʱ�д������ʽ
                  if y1==1
                     break
                  else
                     alpha=alpha/2;
                  end
               else
                  break;
               end
            end
            ynew=fun(x_alpha);
            if ynew<a
               x0{b}=x_alpha;
               y(b)=ynew;
            elseif alpha>0.01 
               while 1
                   while 1
                      alpha=alpha/2;
                      if alpha<=0.01
                         x0=redu(x0,b,'c');
                         break;
                      else
                         x_alpha=a1+alpha*(a1-xl);
                         if flag==1
                             y1=phi(x_alpha);   
                             if y1==1
                                break
                             end
                         else
                             break;
                         end
                      end
                   end
                   y2=fun(x_alpha);
                   if y2<a
                      x0{b}=x_alpha;
                      y(b)=y2;
                      break
                   end 
              end
            else
               x0=redu(x0,b,'c');   
            end
            total=0;
            n1=size(x0,2);
            for i=1:n1
               total=total+(fun(x0{i})-a)^2;
            end
            if total/n1<esp
              [minf,b]=min(y);
              minx=x0{b};  
              break;
            end
         end
    case 2
      n=size(x0,2);
      y=fun(x0);
      k=1;
      while k<2000 
         lamda=0.75;mu=2;
         [a1,b1]=max(y);
         xh=x0{b1};  %���
         fh=a1;
         b=redu(x0,b1,'c');
         xc=mean(cell2mat(b'));   %ʣ��������
         [a2,b2]=min(y);   
         xl=x0{b2};       %��õ�
         xr=2*xc-xh;      %�����
         fr=fun(xr);
         if fr>=fh
            xs=xh+lamda*(xr-xh);
            fs=fun(xs);
         else
            xe=xh+mu*(xr-xh);
            fe=fun(xe);
            if fe<=fr
               xs=xe;
               fs=fe;
            else
               xs=xr;
               fs=fr;
            end
         end
         if  fs<fh
            x0{b1}=xs;
            y(b1)=fs;
         else
            total=0;
            for i=1:n
                total=total+(fun(x0{i})-a1)^2;
            end
            if total/n<esp
               minx=x0{b2};
               minf=a2;
               break
            else
               for i=1:n
                  x0{i}=(x0{i}+xl)/2;
               end
               y=fun(x0);
            end
         end
         k=k+1;
      end
end



       
           
       
           
   
   



   
      










       




