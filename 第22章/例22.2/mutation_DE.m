function new=mutation_DE(old,F,index,val_bound,type)
[r,c]=size(old);
for i=1:r
  k=randperm(r);
  k=my_delrow(k,i);%�۳���ǰ����
  switch type
      case 1   %���������ַ�
         new(i,:)=old(i,:)+F.*(old(k(1),:)-old(k(2),:));
      case 2   %���Ž�
         new(i,:)=old(index,:)+F.*(old(k(1),:)-old(k(2),:));
      case 3  %���Ž�Ӷ������
         new(i,:)=old(index,:)+F.*(old(k(1),:)+old(k(2),:)-old(k(3),:)-old(k(4),:));
      case 4
         new(i,:)=old(i,:)+F.*(old(index,:)-old(k(1),:)+old(k(2),:)-old(k(3),:));
  end
  for j=1:c
    if new(i,j)>val_bound(j,2)||new(i,j)<val_bound(j,1)
        new(i,j)=val_bound(j,1)+(val_bound(j,2)-val_bound(j,1))*rand;
    end
  end
end
