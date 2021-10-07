function y=optifun68(x,data)
pattern=myclass(data,x);    %x为类中心
y=cluster_dis(data,pattern,x);

