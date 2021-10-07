function afish=fishswarm3_1(fun,afish0,afish1,parameter)   %Èº¾Û
fishnum=size(afish1,1);
NC=length(afish0);
n=0;
center=zeros(1,NC);
for i=1:fishnum
    if (fishdstc3_1(afish1(i,:),afish0)<parameter.visual)
        n=n+1;
        center=center+afish1(i,:);      
    end
end
if n~=0
   center=center/n<0.5;    %ÖÐÐÄ
   if (fun(center)/n>fun(afish0)*parameter.delta)
      afish=afish0;
      m=find((center-afish0)~=0);     
      n=randperm(length(m));
      if length(m)>parameter.step
         n=n(1:parameter.step);
      end
      afish(m(n))=center(m(n));
      return
   end
end
afish=fishprey3_1(fun,afish0,parameter);

