function y=optifun76(x,data)
m=cluster_center(data,x);
y=dis2(data,m,x);
y=-y;