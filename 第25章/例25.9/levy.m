function y=levy(varargin)  %levy·ÉÐÐ  
beta=1.5;
a1=mfun('gamma',beta+1);
a2=mfun('gamma',(1+beta)/2);
a3=beta*2^((beta-1)/2);
sigmau=(a1*sin(pi*beta/2)/(a2*a3))^(1/beta);
if nargin==0
     v=abs(normrnd(0,1));
    % u=normrnd(0,sigmau^2);
     u=normrnd(0,1)*sigmau;
     y=u/(v)^(1/beta);
else
    m=varargin{1};n=varargin{2};
    y=zeros(m,n);
    v=abs(normrnd(0,1,m,n));
  % u=normrnd(0,sigmau^2,m,n);
    u=normrnd(0,1,m,n).*sigmau;
    y=u./(v).^(1/beta);
end

