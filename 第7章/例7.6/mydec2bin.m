function y=mydec2bin(x,n)   %ʮ���ƻ�������
s=dec2bin(x,n);
for i=1:n
    y(i)=str2double(s(i));
end
y=transpose(y);