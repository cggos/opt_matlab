function [xmin,minf]=coord1(phi,fun,x0,esp)   %×ø±êÂÖ»»·¨
if nargin<4
    esp=1e-6;
end
k=0;a=mean(x0(2,:));
y1=phi(mean(x0,2));
while k<2000
   phi1=subs(fun,'y',a);
   x1=goldcut1(phi1,x0(1,1),x0(1,2),esp);
   phi1=subs(fun,'x',x1);
   if find(phi1=='y')
       phi1=subs(phi1,'y','x');
   end
   x2=goldcut1(phi1,x0(2,1),x0(2,2),esp);
   y2=phi([x1 x2]);
   if abs(y2-y1)<esp
      xmin=[x1 x2];
      minf=y2;
      break;
   else
      a=x2; 
      y1=y2;
   end
   k=k+1;
end
        
        
        
