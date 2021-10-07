function m_antibody=cloneselection(cloneantibody,m_antibody)
antibodynum=size(m_antibody,2);
num=0;
for i=1:antibodynum
   n=floor(sqrt(antibodynum/i));
   num=num+n;
   if cloneantibody(num-n+1).fitness>m_antibody(i).fitness   %每个克隆子群中的第1个的适应度为最优
        m_antibody(i)=cloneantibody(num-n+1);
   end
end