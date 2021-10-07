function [p_opt,fval]=dyprog2(x1,x2,decisfun,subfun,trafun,objfun)   %��ά��̬�滮
[k1,k]=size(x1);[k2,k]=size(x2);%kΪ�׶���
x1_isnan=~isnan(x1);
x2_isnan=~isnan(x2);
f_opt=nan*ones(k1,k2,k);       %f_opt(i,j,m)Ϊ��m�׶γ�״ֵ̬Ϊ(x1(i,m),x2(j,m))�µ�����ֵ,��ֵΪ����
d_opt1=f_opt;d_opt2=f_opt;     %(d opt1(i,j,m),d opt2(i,j,m))Ϊ��m�׶γ�״ֵ̬Ϊ(x1(i,m),x2(j,m))�µ����ž���ֵ,��ֵΪ����
tmp11=find(x1_isnan(:,k));
tmp21=find(x2_isnan(:,k));
%�ҳ���k�׶�״ֵ̬(���Ƿ���)���±�
tmp12=length(tmp11);
tmp22=length(tmp21);
for i=1:tmp12
  for t=1:tmp22
     [u1,u2]=decisfun(k,x1(tmp11(i),k),x2(tmp21(t),k));%�����Ӧ�������������,�����߱���Ϊһά,��ô�ڶ���DecisFun��
     tmp13=length(u1);tmp14=length(u2);t_vubm=inf;
     %��������for���������k�׶γ�״ֵ̬Ϊx1-(tmp11(i),k),x2(tmp21(t),k)ʱ�����ź���ֵ�����ž���ֵ
     for  j=1:tmp13
        for l=1:tmp14
           tmp=subfun(k,x1(tmp11(i),k),x2(tmp21(t),k),u1(j),u2(l));
           if  tmp<=t_vubm
               f_opt(tmp11(i),tmp21(t),k)=tmp;
               d_opt1(tmp11(i),tmp21(t),k)=u1(j);
               d_opt2(tmp11(i),tmp21(t),k)=u2(l);t_vubm=tmp;
           end
        end
     end
  end
end
for  ii=k-1:-1:1  %�Ӻ���ǰ��������f_opt�Լ�d_opt
    tmp011=find(x1_isnan(:,ii));tmp021=find(x2_isnan(:,ii)); 
    tmp012=length(tmp011);
    tmp022=length(tmp021);
    for  i=1:tmp012
       for t=1:tmp022
          [u1,u2]=decisfun(ii,x1(tmp011(i),ii),x2(tmp021(t),ii));
          tmp013=length(u1);
          tmp014=length(u2);
          t_vubm=inf;
          for  j=1:tmp013
              for l=1:tmp014
                   tmp000=subfun(ii,x1(tmp011(i),ii),x2(tmp021(t),ii),u1(j),u2(l));
                   tmp100=trafun(ii,x1(tmp011(i),ii),x2(tmp021(t),ii),u1(j),u2(l));
                   %�ɸý׶ε�״ֵ̬����Ӧ�ľ���ֵ�����һ�׶ε�״ֵ̬
                   tmp200=x1(:,ii+1)-tmp100(1);
                   tmp300=x2(:,ii+1)-tmp100(2);
                   tmp400=find(tmp200==0);%�ҳ���һ�׶ε�״ֵ̬��x1(:,ii+ 1)���±�
                   tmp500=find(tmp300==0);%�ҳ���һ�׶ε�״ֵ̬��x2(:,ii+ 1)���±�
                   if ~isempty(tmp400)&&~isempty(tmp500)
                       if  nargin<6
                          tmp000=tmp000+f_opt(tmp400(1),tmp500(1),ii+ 1);
                       else
                          tmp000=objfun(tmp000,f_opt(tmp400(1),tmp500(1),ii+ 1));
                       end
                       if  tmp000<t_vubm
                           f_opt(tmp011(i),tmp021(t),ii)=tmp000;
                           d_opt1(tmp011(i),tmp021(t),ii)=u1(j);
                           d_opt2(tmp011(i),tmp021(t),ii)=u2(l);
                           t_vubm=tmp000;
                       end
                   end
              end
          end
        end
    end
end
fval=f_opt(x1_isnan(:,1),x2_isnan(:,1),1);
p_opt=[];tmpx1=[];tmpx2=[];tmpd1=[];tmpf=[];tmpd2=[];
tmp11=find(x1_isnan(:,1));tmp01=length(tmp11);
tmp12=find(x2_isnan(:,1));tmp02=length(tmp12);
for  i=1:tmp01
        q=(i-1)*k*tmp02;
        for j=1:tmp02
            t=k*(j-1);
            t=q+t;
            tmpd1(i)=d_opt1(tmp11(i),tmp12(j),1);
            tmpd2(j)=d_opt2(tmp11(i),tmp12(j),1);%�����һ�׶εľ���ֵ
            tmpx1(i)=x1(tmp11(i),1);
            tmpx2(j)=x2(tmp12(j),1);%�����һ�׶ε�״ֵ̬
            tmpf(i,j)=subfun(1,tmpx1(i),tmpx2(j),tmpd1(i),tmpd2(j));
            %�����һ�׶ε�ָ�꺯��ֵ
            p_opt(t+1,[1 2 3 4 5 6])=[1,tmpx1(i),tmpx2(j),tmpd1(i),tmpd2(j),tmpf(i,j)];
            for  ii=2:k
               %��˳��������׶εľ���ֵ��״ֵ̬�Լ�ָ�꺯��ֵ
               u=trafun(ii-1,tmpx1(i),tmpx2(j),tmpd1(i),tmpd2(j));
               tmpx1(i)=u(1);
               tmpx2(j)=u(2);
               tmp1=x1(:,ii)-tmpx1(i);
               tmp2=x2(:,ii)-tmpx2(j);
               tmp3=find(tmp1==0);
               tmp4=find(tmp2==0);
               if ~isempty(tmp3)&&~isempty(tmp4)
                   tmpd1(i)=d_opt1(tmp3(1),tmp4(1),ii);
                   tmpd2(j)=d_opt2(tmp3(1),tmp4(1),ii);
               end
               tmpf(i,j)=subfun(ii,tmpx1(i),tmpx2(j),tmpd1(i),tmpd2(j));
               p_opt(t+ii,[ 1 2 3 4 5 6])=[ii,tmpx1(i),tmpx2(j),tmpd1(i),tmpd2(j),tmpf(i,j)];
            end
       end
end 



