function [old,newsigma]=IAES_select1(old,new,oldsigma,newsigma,type)  %ES��ѡ�����ӣ��ṹ�壩
num=size(old,2);
total_chome=[];
total_F=[];
for i=1:num
    total_F=[total_F;old(i).fitness;new(i).fitness];
    total_chome=[total_chome;old(i).x;new(i).x];
end
total_sigma=[oldsigma;newsigma];
if type==1   %(N+�ˣ���ʽ
   [a,b]=sort(total_F);
   total_chome=total_chome(b,:);
   total_sigma=total_sigma(b,:);
   for i=1:num
     old(i).x=total_chome(i,:);
   end
   newsigma=total_sigma(1:num,:);
elseif type==2   %(N,�ˣ���ʽ
   if size(new,2)<num
       error('�¸������Ŀ̫��');
   end
   newF=[];
   for i=1:num
      newF=[newF;new(i).fitness];
   end
   [a,b]=sort(newF);
   new=new(b);
   newsigma=newsigma(b,:);
   for i=1:num
      old(i).x=new(i).x;
   end
   newsigma=newsigma(1:num,:);
end
