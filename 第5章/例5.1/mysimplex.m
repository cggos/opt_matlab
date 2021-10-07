function [minx,minf,A2,b,s]=mysimplex(A,b,c)     %�����η������Թ滮��c,ΪĿ�꺯����
%AΪ��׼ʽ��Լ��������,bΪ������,nΪ�����������û�л����⣬����������⣬��ʱtype=2
n=find(c==0);
[r1,c1]=size(A);
m=findeye(A);
minx1=zeros(1,c1);
if isempty(m)||size(m,1)<r1
    c_L=length(c);
    c_temp=c;
    c=[zeros(1,c_L) ones(1,r1)];
    A=[A eye(r1)];
    n1=c1+1:c1+r1;
    c1=size(A,2);
    m=findeye(A);
    [a,b6]=sort(m(:,1));
    m=m(b6,:);
    m=m(end-r1+1:end,:);
    type=2;
else
    type=1;
end
cI=c(m(:,1));
for i=1:c1
    sigma(i)=cI*A(:,i)-c(i);
end
if isempty(find(sigma>0))
    [minx,b]=dualsimplex(A,b,c);   %�ö�ż�����
    minf=minx.f;
    A2=minx.A;
    s=[];
else
   f0=cI*b;
   flag=0;
   [A,b,sigma,f0]=f_simplex(A,b,sigma,f0,flag);
   if type==2
      A=redu(A,n1,'c');
      for i=1:r1
         a2=findzeros(A(i,:),0); 
         if strcmp(a2,'all')
            A=redu(A,i,'r');  %ɾ������
            b=redu(b,i,'r');
         end
      end
      m=findeye(A);
      if isempty(m)
          error('�����Ž�');
      else
          cI=c_temp(m(:,1));
          for i=1:size(A,2)
             sigma1(i)=cI*A(:,i)-c_temp(i);
          end
          f0=cI*b;
          flag=0;
          aa=find(sigma1>0);
          if ~isempty(aa)
            [A,b,sigma,f0]=f_simplex(A,b,sigma1,f0,flag);
          else
             sigma=sigma1;
          end
      end
   end
   y=findeye(A);
   minx1(y(:,1)')=b(y(:,2));
   a=findzeros(sigma,0);
   if length(a)>r1
      s={'���������'};
   else
      s=[];
   end
   if ~isempty(n)
      minx2=redu(minx1,n,'c');
   else
      minx2=[];
   end
   minf=f0;
   A2=[A b;sigma minf];
   minx=struct;
   minx.x=minx1;
   minx.x1=minx2; 
end

        
function [A,b,sigma,f0]=f_simplex(A,b,sigma,f0,flag)   %�������е���Ԫ
[r,c1]=size(A);
while 1
    A_temp=A;
    b_temp=b;
    sigma1=sigma;
    a=find(sigma>0);
    if isempty(a)          %�н�
          break
    else
       for i=1:length(a)  %�޽�
          a1=find(A_temp(:,a(i))<0);
          if ~isempty(a1)&&length(a1)==r
              flag=1;
              break;
          end
       end
       switch flag
           case 1
              error('�����Ž�');
           case 0
              [a,id]=sort(sigma1,'descend');%�����
              b1=find(A_temp(:,id(1))>0);
              [a4,ser]=redu(1:r,b1,'c');
              [a,a3]=sort(b(b1)./A_temp(b1,id(1)));
              b3=find(a==a(1));
              idex=ser(a3(b3(end)));
              zuyan=A_temp(idex,id(1));   %��Ԫ
              for i=1:r
                for j=1:c1 
                  if i==idex
                     A(i,j)=A_temp(idex,j)/zuyan;
                  else           
                     A(i,j)=A_temp(i,j)-A_temp(idex,j)*A_temp(i,id(1))/zuyan;
                  end
                end
              end
              for i=1:r
                if i==idex
                  b(i,1)=b_temp(idex,1)/zuyan;
                else           
                  b(i,1)=b_temp(i,1)-b_temp(idex,1)*A_temp(i,id(1))/zuyan;
                end 
              end
              for i=1:c1
                sigma(i)=sigma1(i)-A_temp(idex,i)*sigma1(id(1))/zuyan;
              end
              f0=f0-b_temp(idex,1)*sigma1(id(1))/zuyan;
        end
    end
end

        
        




