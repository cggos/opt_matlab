function [x,val,C]=Hungarian(C)      %匈牙利算法解指派问题
[r,c]=size(C);
C1=C;
while 1
  y1=find(C==0);    %找零
  if isempty(y1)
     y1=min(C,[],2);  %每行的最小值
     C=C-repmat(y1,1,c);
     y1=min(C,[],1);
     C=C-repmat(y1,c,1);
  end
  C2=C;
  C=lattice(C);
  [Q1,Q2]=find(C==-inf);
  num=length(Q1);
  if num==r
      x=[Q1 Q2];   %行及列序号
      val=0;
      for i=1:num
          val=val+C1(x(i,1),x(i,2));
      end
      break;
  else
     [y1,y7]=find(C==-inf);
     y1=redu(1:r,y1,'c');   %找-inf所在的行
     y2=find(C(y1,:)==inf); %找inf所在的列
     m1=[];
     for i=1:length(y2)
        a=find(C(:,y2(i))==-inf);   %找-inf所在的行
        if ~isempty(a)
            m1=[m1 a];
        end
     end
     m2=[];
     for i=1:length(m1)
         a=find(C(m1(i),:)==inf);
         if ~isempty(a)
             m2=[m2 a];
         end
     end
     n3=[y1 m1];
     n1=redu(1:r,n3,'c');
     C=redu(C,n1,'r');
     n2=[y2 m2];
     n4=redu(1:c,n2,'c');
     C=redu(C,n2,'c');
     sita=min(min(C));
     [r2,c2]=size(C);
     if c2>=r2   %列数多
         C2(n3,:)=C2(n3,:)-sita.*ones(length(n3),c);   %每行减sita
         [y1,y2]=find(C2<0);
         C2(:,unique(y2))=C2(:,unique(y2))+sita.*ones(r,length(unique(y2)));
     else
         C2(:,n4)=C2(:,n4)-sita.*ones(r,length(n4));   %每列减sita
         [y1,y2]=find(C2<0);
         C2(unique(y1),:)=C2(unique(y1),:)+sita.*ones(length(unique(y1)),c);
     end    
  end
  C=C2;
end




  
  

   
