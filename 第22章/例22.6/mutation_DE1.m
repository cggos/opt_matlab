function new=mutation_DE1(old,F,index,LB,UB,type)   %多类聚类中心的差分
r=size(old,2);
for i=1:r
  k=randperm(r);
  k=my_delrow(k,i);%扣除当前个体
  switch type
      case 1   %随机向量差分法
         new(i).center=old(i).center+F.*(old(k(1)).center-old(k(2)).center);
      case 2   %最优解
         new(i).center=old(index).center+F.*(old(k(1)).center-old(k(2)).center);
      case 3  %最优解加多个个体
         new(i).center=old(index).center+F.*(old(k(1)).center+old(k(2)).center-old(k(3)).center-old(k(4)).center);
      case 4
         new(i).center=old(i).center+F.*(old(index).center-old(k(1)).center+old(k(2)).center-old(k(3)).center);
  end
  new(i).center=boundtest(new(i).center,LB,UB);
end
