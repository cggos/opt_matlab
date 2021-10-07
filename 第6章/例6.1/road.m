function [y,rd2]=road(x,d,str)    %Çó×î¶ÌÂ·¾¶
if nargin==2
    str=[];
end
yr=exhaus(x);
r1=size(d,1);
r2=size(yr,1);
y=inf;
a=zeros(1,r2);
for i=1:r2
    for j=1:r1
        [y1,y2]=mycompare1(yr(i,j),x{j});
        [y3,y4]=mycompare1(yr(i,j+1),x{j+1});      
        a(i)=a(i)+d{j}(y2,y4);
    end
    if a(i)<=y
       y=a(i);
    end
end
rd1=yr(find(a<=y),:);
if ~isempty(str)
    r3=size(str,1);
    num=size(rd1,1);
    rd2=cell(1,num);
    for i=1:num
        for j=1:r3
           [syms,n]=findletter(str{j},2);
           [y1,y2]=mycompare1(rd1(i,j),x{j});
           if j==r3
              if n==1
                  rd2{i}=[rd2{i} syms(y2)];
              else
                  rd2{i}=[rd2{i} syms{y2}];
              end
           else
               if n==1
                   rd2{i}=[rd2{i} syms(y2) '¡ú'];
               else
                   rd2{i}=[rd2{i} syms{y2} '¡ú'];
               end
           end
        end
    end
end
           