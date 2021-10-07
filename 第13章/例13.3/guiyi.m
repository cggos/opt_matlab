function [y,s]=guiyi(x,y1)   %归一化
if iscell(x)
    x=cell2mat(x);
end
if nargin==1
    y1=[];
else 
   s_mean=y1.mean;
   s_std=y1.std;
end
r=size(x,1);
a1=mean(x);stdr=std(x);
if isempty(y1)  %归一化
    y=(x-a1(ones(r,1),:))./stdr(ones(r,1),:);
    s.mean=a1;s.std=stdr;
else
    y=s_std(ones(r,1),:).*x+s_mean(ones(r,1),:);   %反归一化
end
        
