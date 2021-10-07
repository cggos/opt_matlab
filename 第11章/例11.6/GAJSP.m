function MakeSpan=GAJSP(NIND,MAXGEN,GGAP,P_Cross,P_Mutation,Jm,T)
%====�������====
   %NIND           ����Ⱥ�������������Ŀ
   %MAXGEN         ������Ŵ�����
   %GGAP           ������
   %P_Cross        ���������
   %P_Mutation     ���������
   %Jm             ��������Ŀ�ѡ�������ϣ�Ϊm��n��Ԫ������mΪ�����������򣩣�nΪ������
   %T              ���ӹ�ʱ�����       
%====�������====
   %MakeSpan       ����С������깤ʱ��

gen=0;  %����������
JmNumber=Max_Cell(Jm);               %����Max_Cell�Ӻ��������������
[PNumber,MNumber]=size(Jm);          %PNumberΪ����������MNumberΪ�������
trace=zeros(2,MAXGEN);               %Ѱ�Ž���ĳ�ʼֵ��һ�д�Ÿ��������Ž⣬һ�д�Ÿ�����ľ�ֵ
TotalOP_Number=PNumber*MNumber;      %�����ܸ���
%��ʼ��
Number=zeros(1,PNumber);             %Number���ÿ�������Ĺ�������PNumber��������
for i=1:PNumber
    Number(i)=MNumber;               %MNumber�������
end
%Ⱦɫ�������룺�����Ϊ���㣬��һ���ʾ���򣬵ڶ����ʾ����
Chrom=zeros(NIND,2*TotalOP_Number);  %Ⱦɫ��ĳ���Ϊ2*TotalOP_Number
for j=1:NIND                  %���Ѱ��
    WPNumberTemp=Number;
    for i=1:TotalOP_Number    %�����������
        val=unidrnd(PNumber);
        while WPNumberTemp(val)==0
            val=unidrnd(PNumber);
        end
        %��һ������ʾ����
        Chrom(j,i)=val;
        WPNumberTemp(val)=WPNumberTemp(val)-1;
        %�ڶ�������ʾ����
        Temp=Jm{val,MNumber-WPNumberTemp(val)};
        SizeTemp=length(Temp);
        %������ɹ������
        Chrom(j,i+TotalOP_Number)=unidrnd(SizeTemp);
    end
end
%����Ŀ�꺯��ֵ
[PVal,ObjV,P,S]=cal(Chrom,JmNumber,T,Jm);                %�����Ӻ�������cal�������Ӧ��ֵ
%ѭ��Ѱ�����Ž�
while gen<MAXGEN
    %������Ӧ��ֵ
    FitnV=ranking(ObjV);                                 %�����Ӻ�������ranking�����������
    %ѡ�����
    SelCh=select('rws', Chrom, FitnV, GGAP);             %�����Ӻ�������select����ѡ��������������̶Ĳ��ԣ�rws��
    %�������
    SelCh=across(SelCh,P_Cross,Jm,T);                    %�����Ӻ�������across���н������
    %�������
    SelCh=aberranceJm(SelCh,P_Mutation,Jm,T);            %�����Ӻ�������aberranceJm���б������
    %���¼���Ŀ����Ӧ��ֵ�����壩
    [PVal,ObjVSel,P,S]=cal(SelCh,JmNumber,T,Jm);         %�����Ӻ�������cal�������Ӧ��ֵ
    %���²�������Ⱥ
    [Chrom,ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);    %�����Ӻ�������reins���²�������Ⱥ
    %��������������1
    gen=gen+1;
    %��������ֵ
    trace(1,gen)=min(ObjV);
    trace(2,gen)=mean(ObjV);
    %��¼���ֵ
    if gen==1
        Val1=PVal;
        Val2=P;
        MinVal=min(ObjV);  %��С������̼ӹ�ʱ��
        STemp=S;
    end
    %��¼��С�Ĺ���
    if MinVal>trace(1,gen)
        Val1=PVal;
        Val2=P;
        MinVal=trace(1,gen);
        STemp=S;
    end
end
%��ǰ���ֵ
PVal=Val1;    %����ʱ��
P=Val2;       %���� 
S=STemp;      %���Ȼ��򺬻�������
MakeSpan = min(trace(1,:));      %�����������깤ʱ�䣨����̼ӹ�ʱ�䣩����Сֵ
%��ʾ���Ž�
MP=S(1,PNumber*MNumber+1:PNumber*MNumber*2);  %������������
for i=1:TotalOP_Number  
    val=P(1,i);
    a=(mod(val,100));   %����
    b=((val-a)/100);    %����
    Temp=Jm{b,a};
    mText=Temp(MP(1,i));
    x1=PVal(1,i);       %����Ŀ���ʱ��
    x2=PVal(2,i);       %������깤ʱ��
    y1=mText-0.5;       %0.5Ϊ����ͼ�������Ŀ�ȣ������ֿ��Ը���
    y2=mText;           %���ڻ��ƾ������ĸ߶�
    PlotRec(x1,x2,mText); %����PlotRec�������Ƹ���ͼ����
    hold on;
    fill([x1,x2,x2,x1],[y1,y1,y2,y2],[1-1/b,1/b,b/PNumber]);     %��������
    text((x1+x2)/2,mText-0.25,num2str(P(i)));                    %���ı���ʶ���򣬱�עλ��Ϊ������������
end
function JmNumber = Max_Cell(Jm)
%Jm        ������Ԫ������
%JmNumber  ����������
MMAX=0;
[m,n]=size(Jm);   %Ԫ�����������������������mΪ������������nΪ������
for i=1:m         %���в������ֵ
    for j=1:n     %���в������ֵ
        TMAX=max(Jm{i,j});
        if MMAX<=TMAX;
            MMAX=TMAX;
        end
    end    
end
JmNumber=MMAX;

function [PVal,ObjV,P,S] = cal(Chrom,JmNumber,T,Jm)  % �Ӻ�������2����Ŀ�꺯��ֵ
% �������
   %Chrom       ������Ⱥ  
   %JmNumber    ��������
   %T           �ӹ�ʱ����󡪡�������������ʹ�õ�ʱ�� 
   %Jm          ����Ԫ�����󣨸�����������ʹ�õĻ����� 
% �������
   %PVal        ��ѵ��ȹ���ʱ�䣨������ʼʱ�����깤ʱ�䣩 
   %P           �������ĵ��ȹ��� 
   %ObjV        Ⱥ��ÿ������ĵ��ȹ���ʱ��
   %S           �������ĵ��Ȼ���
%��ʼ��
NIND=size(Chrom,1);
ObjV=zeros(NIND,1);

%���㹤��������������� 
[PNumber,MNumber]=size(Jm);
%��3������Ѱ��
for i=1:NIND
    %ȡһ������
    S=Chrom(i,:);
    %���ݻ��򣬼�����ȹ���
    P=calp(S,PNumber);                     %�����Ӻ�������calp������Ӧ��
    %���ݵ��ȹ��򣬼�������ȹ���ʱ��
    PVal=caltime(S,P,JmNumber,T,Jm);       %�����Ӻ�������caltime        
    %ȡ���ʱ��
    MT=max(PVal);
    TVal=max(MT); 
    %����ʱ��
    ObjV(i,1)=TVal;
    %��ʼ��
    if i==1
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end
    %��¼��С�ĵ��ȹ���ʱ�䡢��ѵ��ȹ���ʱ�� �������ĵ��ȹ���
    if MinVal>ObjV(i,1)
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end
end
%��ѵ��ȹ���ʱ��,�������ĵ��ȹ���
 PVal=Val1;
 P=Val2;
 S=STemp;

function FitnV=ranking(ObjV, RFun, SUBPOP)
   [Nind,Anans]=size(ObjV);
   if nargin<2, RFun=[]; end
   if nargin>1, if isnan(RFun), RFun=[]; end, end
   if numel(RFun)==2,
      if RFun(2)==1, NonLin=1;
      elseif RFun(2) == 0, NonLin = 0;
      else error('Parameter for ranking method must be 0 or 1');
      end
      RFun = RFun(1);
      if isnan(RFun), RFun = 2; end
   elseif numel(RFun) > 2,
      if numel(RFun) ~= Nind, error('ObjV and RFun disagree'); end
   end

   if nargin < 3, SUBPOP = 1; end
   if nargin > 2,
      if isempty(SUBPOP), SUBPOP = 1;
      elseif isnan(SUBPOP), SUBPOP = 1;
      elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar');
      end
   end

   if (Nind/SUBPOP) ~= fix(Nind/SUBPOP), error('ObjV and SUBPOP disagree'); end
   Nind = Nind/SUBPOP;  % Compute number of individuals per subpopulation

% Check ranking function and use default values if necessary
   if isempty(RFun),
      % linear ranking with selective pressure 2
         RFun = 2*(0:Nind-1)'/(Nind-1);
   elseif numel(RFun) == 1
      if NonLin == 1,
         % non-linear ranking
         if RFun(1) < 1
             error('Selective pressure must be greater than 1');
         elseif RFun(1) > Nind-2
             error('Selective pressure too big'); 
         end
         Root1 = roots([RFun(1)-Nind [RFun(1)*ones(1,Nind-1)]]);
         RFun = (abs(Root1(1)) * ones(Nind,1)) .^ [(0:Nind-1)'];
         RFun = RFun / sum(RFun) * Nind;
      else
         % linear ranking with SP between 1 and 2
         if (RFun(1) < 1 || RFun(1) > 2),
            error('Selective pressure for linear ranking must be between 1 and 2');
         end
         RFun = 2-RFun + 2*(RFun-1)*[0:Nind-1]'/(Nind-1);
      end
   end;

   FitnV = [];

% loop over all subpopulations
for irun = 1:SUBPOP,
   % Copy objective values of actual subpopulation
      ObjVSub = ObjV((irun-1)*Nind+1:irun*Nind);
   % Sort does not handle NaN values as required. So, find those...
      NaNix = isnan(ObjVSub);
      Validix = find(~NaNix);
   % ... and sort only numeric values (smaller is better).
      [Anans,ix] = sort(-ObjVSub(Validix));

   % Now build indexing vector assuming NaN are worse than numbers,
   % (including Inf!)...
      ix = [find(NaNix) ; Validix(ix)];
   % ... and obtain a sorted version of ObjV
      Sorted = ObjVSub(ix);

   % Assign fitness according to RFun.
      i = 1;
      FitnVSub = zeros(Nind,1);
      for j = [find(Sorted(1:Nind-1) ~= Sorted(2:Nind)); Nind]',
         FitnVSub(i:j) = sum(RFun(i:j)) * ones(j-i+1,1) / (j-i+1);
         i =j+1;
      end

   % Finally, return unsorted vector.
      [Anans,uix] = sort(ix);
      FitnVSub = FitnVSub(uix);

   % Add FitnVSub to FitnV
      FitnV = [FitnV; FitnVSub];
end
 
function P = calp(S,PNumber)
TotalOP_Number = length(S);   %�����ܸ���
TotalOP_Number = TotalOP_Number/2;
%ȡ�������ȡ�����һ��
S = S(1,1:TotalOP_Number);
%��ʼ��
temp = zeros(1,PNumber);
P = zeros(1,TotalOP_Number);
%�������ɵ��ȹ���
for i = 1:TotalOP_Number
    %�����+1
    temp(S(i)) = temp(S(i))+1;
    P(i) = S(i)*100+temp(S(i));  
end

    function PVal = caltime(S,P,JmNumber,T,Jm)
[PNumber,MNumber] = size(Jm);  %�������������
%ȡ��������ȡ�����һ�루�����һ�룩
M = S(1,PNumber*MNumber+1:PNumber*MNumber*2); 
TotalOP_Number = length(P);        %�����ܸ���
%��ʼ��
TM = zeros(1,JmNumber);
TP = zeros(1,PNumber);
PVal = zeros(2,TotalOP_Number);    %��Ź����ʱ�䣬һ��Ϊ����ʱ�䣬һ��Ϊ�깤ʱ��

%������ȹ���ʱ��
for i = 1: TotalOP_Number 
    % ȡ������
    val = P(1,i);
    a = (mod(val,100));      %����
    b = ((val-a)/100);       %����
    Temp = Jm{b,a};
    m = Temp(M(1,i));
    
    %ȡ�ӹ�ʱ��
    Temp = T{b,a};
    t = Temp(M(1,i));        %����ļӹ�����ʱ��
    
    %ȡ�����ӹ�������Ŀ�ʼʱ���ǰ��һ����������ʱ��
        TMval = TM(1,m);
        TPval = TP(1,b);
    
    %�����ӹ�������Ŀ�ʼʱ�䣺������ڵ���ǰ��һ����������ʱ��ȡ�����ӹ�������Ŀ�ʼʱ��
    if TMval >= TPval 
        val = TMval;
    %����ȡǰ��һ����������ʱ��  
    else
        val = TPval;
    end
    %����ʱ��
    PVal(1,i) = val;    %��ʼʱ��
    PVal(2,i) = val+t;  %���ʱ�䣬���t=0������ζ��ʵ����û�иù���
    
    %��¼���ι���Ļ���ʱ��͹���ʱ��
    TM(1,m) = PVal(2,i);
    TP(1,b) = PVal(2,i); 
end

function SelCh = select(SEL_F,Chrom,FitnV,GGAP,SUBPOP)
   if nargin < 3, error('Not enough input parameter'); end

   % ȷ����Ⱥ�Ĺ�ģ
   [NindCh,Nvar] = size(Chrom);
   [NindF,VarF] = size(FitnV);
   if NindCh ~= NindF, error('Chrom and FitnV disagree'); end
   if VarF ~= 1, error('FitnV must be a column vector'); end
  
   if nargin < 5, SUBPOP = 1; end
   if nargin > 4,
      if isempty(SUBPOP), SUBPOP = 1;
      elseif isnan(SUBPOP), SUBPOP = 1;
      elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar');
      end
   end

   if (NindCh/SUBPOP) ~= fix(NindCh/SUBPOP), error('Chrom and SUBPOP disagree'); end
   Nind = NindCh/SUBPOP;  % ����ÿ������Ⱥ�ĸ�������

   if nargin < 4
       GGAP = 1; 
   end
   if nargin > 3
      if isempty(GGAP)
          GGAP = 1;
      elseif isnan(GGAP)
          GGAP = 1;
      elseif length(GGAP) ~= 1, error('GGAP must be a scalar');
      elseif (GGAP < 0), error('GGAP must be a scalar bigger than 0'); 
      end
   end

%��������Ⱥ�������Ա�ѡ��
   NSel = max(floor(Nind*GGAP+.5),2);

%����Ⱥ��ѡ�����
   SelCh = [];
   for irun = 1:SUBPOP,
      FitnVSub = FitnV((irun-1)*Nind+1:irun*Nind);
      ChrIx = feval(SEL_F,FitnVSub,NSel)+(irun-1)*Nind;
      SelCh = [SelCh; Chrom(ChrIx,:)];
   end

function NewChrom = across(Chrom,P_Cross,Jm,T)
%��ʼ������Ⱥ
[NIND,TotalOP_Number] = size(Chrom);
TotalOP_Number = TotalOP_Number/2;
NewChrom = Chrom;

[PNumber,MNumber] = size(Jm);
Number = zeros(1,PNumber);
for i = 1:PNumber
    Number(i) = 1;
end

%���ѡ�񽻲����(ϴ�ƽ���)
SelNum = randperm(NIND);   

Num = floor(NIND/2);    %������������
for i = 1:2:Num
    if P_Cross > rand;
        Pos = unidrnd(TotalOP_Number);%����λ��
        while Pos == 1
            Pos = unidrnd(TotalOP_Number);
        end
        %ȡ������ĸ��岢��ʼ������
        S1 = Chrom(SelNum(i),1:TotalOP_Number);
        S2 = Chrom(SelNum(i+1),1:TotalOP_Number); 
        S11 = S2;S22 = S1;      %��ʼ���µĸ���     
        %�¸����м�Ƭ�ϵ�COPY      
        S11(1:Pos) = S1(1:Pos);      
        S22(1:Pos) = S2(1:Pos);        
        %�Ƚ�S11���S1,S22���S2�����ȱʧ�Ļ���
        S3 = S11;S4 = S1;
        S5 = S22;S6 = S2;
        for j = 1:TotalOP_Number         
           Pos1 = find(S4 == S3(j),1);
           Pos2 = find(S6 == S5(j),1);
           if Pos1 > 0
               S3(j) = 0;
               S4(Pos1) = 0;
           end                         
           if Pos2 > 0
               S5(j) = 0;
               S6(Pos2) = 0;
           end
        end
        for j = 1:TotalOP_Number          
          if S3(j) ~= 0             %����Ļ���          
            Pos1 = find(S11 == S3(j),1);        
            Pos2 = find(S4,1);      %����ȱʧ�Ļ���
            S11(Pos1) = S4(Pos2);   %��ȱʧ�Ļ����޲�����Ļ���
            S4(Pos2) = 0;       
          end 
          if S5(j) ~= 0              
            Pos1 = find(S22 == S5(j),1); 
            Pos2 = find(S6,1);           
            S22(Pos1) = S6(Pos2);
            S6(Pos2) = 0;
          end
        end

        % ���潻��ǰ�Ļ�������
        S1 = Chrom(SelNum(i),:);
        S2 = Chrom(SelNum(i+1),:); 
       
        for k = 1:TotalOP_Number            
            Pos1 = Find(S11(k),S1);                         %����Find����     
            S11(TotalOP_Number+k) = S1(TotalOP_Number+Pos1);
            S1(Pos1) = 0;
            Pos1 = Find(S22(k),S2);           
            S22(TotalOP_Number+k) = S2(TotalOP_Number+Pos1);
            S2(Pos1) = 0;
        end
        %�����µ���Ⱥ
        NewChrom(SelNum(i),:) = S11;
        NewChrom(SelNum(i+1),:) = S22;
    end
end

function ChromNew = aberranceJm(Chrom,P_Mutation,Jm,T)
%��ʼ��
[NIND,TotalOP_Number] = size(Chrom);
TotalOP_Number = TotalOP_Number/2;
ChromNew = Chrom;
[PNumber,MNumber] = size(Jm);
Number = zeros(1,PNumber);
for i = 1:PNumber
    Number(i) = 1;
end

for i = 1:NIND
    %ȡһ������
    S = Chrom(i,:);
    WPNumberTemp = Number;
    for j = 1:TotalOP_Number
        JMTemp = Jm{S(j), WPNumberTemp(S(j))};
        SizeTemp = length(JMTemp);
        %�Ƿ����
        if P_Mutation > rand;
            %ѡ����������ѡ��
            %S(j+TotalOP_Number) = unidrnd(SizeTemp);
            %ѡ��������ӹ�ʱ���ٵ�ѡ���ʴ�
                if SizeTemp == 1
                    S(j+TotalOP_Number) = 1;
                else
                    S(j+TotalOP_Number) = selectJm(S(j++TotalOP_Number),T{S(j),WPNumberTemp(S(j))});    %����selectJm����
                end
        end
        WPNumberTemp(S(j)) = WPNumberTemp(S(j))+1;
    end
    %���ݷ�����Ⱥ
    ChromNew(i,:)=S;
end

function PlotRec(mPoint1,mPoint2,mText)
%������������ڿ��ƾ�������С������
grid on
vPoint = zeros(4,2);                   %���Ƹ���ͼ���þ���������ʾ��
vPoint(1,:) = [mPoint1,mText-0.5];     %���θ������������ĸ����������
vPoint(2,:) = [mPoint2,mText-0.5];
vPoint(3,:) = [mPoint1,mText];
vPoint(4,:) = [mPoint2,mText];
plot([vPoint(1,1),vPoint(2,1)],[vPoint(1,2),vPoint(2,2)])   %����ˮƽ�ߣ��߶�ΪmText-0.5������ΪvPoint(2,1)��vPoint(1,1)
hold on;
plot([vPoint(1,1),vPoint(3,1)],[vPoint(1,2),vPoint(3,2)])   %���ƴ�ֱ�ߣ��߶�Ϊ0.5��������ΪvPoint(1,1)=vPoint(3,1)
plot([vPoint(2,1),vPoint(4,1)],[vPoint(2,2),vPoint(4,2)])   %���ƴ�ֱ�ߣ��߶�Ϊ0.5��������ΪvPoint(2,1)=vPoint(4,1)
plot([vPoint(3,1),vPoint(4,1)],[vPoint(3,2),vPoint(4,2)])   %����ˮƽ�ߣ��߶�ΪmText������ΪvPoint(4,1)��vPoint(3,1)=vPoint(2,1)��vPoint(1,1)
title('�����ӹ��ĸ���ͼ');
xlabel('�ӹ�ʱ�䣨ʱ�䵥λ��');
ylabel('�ӹ�����');

function NewChrIx = rws(FitnV,Nsel)
% Identify the population size (Nind)
[Nind,Anans] = size(FitnV);
% Perform Stochastic Sampling with Replacement
cumfit  = cumsum(FitnV);
trials = cumfit(Nind) .* rand(Nsel, 1);
Mf = cumfit(:, ones(1, Nsel));
Mt = trials(:, ones(1, Nind))';
[NewChrIx, Anans] = find(Mt < Mf & [ zeros(1, Nsel); Mf(1:Nind-1, :) ] <= Mt);

function SelS = selectJm(S,S_T)

Num = length(S_T);
MaxVal = 0;
for i = 1:Num
  if S_T(i) > MaxVal;
   MaxVal = S_T(i);
  end
end
MaxVal = 2*MaxVal;

for i = 1:Num
 S_T(i) = MaxVal-S_T(i);
end

eVal = 0;
for i = 1:Num
 eVal = eVal+S_T(i);
end

for i = 1:Num
 P(i) = S_T(i)/eVal;
end

for i = 2:Num
  P(i) = P(i)+P(i-1);
end

 num = rand;
 SelS = 1;
 while(num > P(SelS))
     SelS = SelS+1;    
 end 

function [Chrom,ObjVCh] = reins(Chrom,SelCh,SUBPOP,InsOpt,ObjVCh,ObjVSel)
   % Check parameter consistency
   if nargin < 2, error('Not enough input parameter'); end
   if (nargout == 2 && nargin < 6), error('Input parameter missing: ObjVCh and/or ObjVSel'); end

   [NindP, NvarP] = size(Chrom);
   [NindO, NvarO] = size(SelCh);

   if nargin == 2, SUBPOP = 1; end
   if nargin > 2,
      if isempty(SUBPOP), SUBPOP = 1;
      elseif isnan(SUBPOP), SUBPOP = 1;
      elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar'); 
      end
   end

   if (NindP/SUBPOP) ~= fix(NindP/SUBPOP), error('Chrom and SUBPOP disagree'); end
   if (NindO/SUBPOP) ~= fix(NindO/SUBPOP), error('SelCh and SUBPOP disagree'); end
   NIND = NindP/SUBPOP;  % Compute number of individuals per subpopulation
   NSEL = NindO/SUBPOP;  % Compute number of offspring per subpopulation

   IsObjVCh = 0; IsObjVSel = 0;
   if nargin > 4, 
      [mO, nO] = size(ObjVCh);
      if nO ~= 1, error('ObjVCh must be a column vector'); end
      if NindP ~= mO, error('Chrom and ObjVCh disagree'); end
      IsObjVCh = 1;
   end
   if nargin > 5, 
      [mO, nO] = size(ObjVSel);
      if nO ~= 1, error('ObjVSel must be a column vector'); end
      if NindO ~= mO, error('SelCh and ObjVSel disagree'); end
      IsObjVSel = 1;
   end
       
   if nargin < 4, INSR = 1.0; Select = 0; end   
   if nargin >= 4,
      if isempty(InsOpt), INSR = 1.0; Select = 0;   
      elseif isnan(InsOpt), INSR = 1.0; Select = 0;   
      else
         INSR = NaN; Select = NaN;
         if (length(InsOpt) > 2), error('Parameter InsOpt too long'); end
         if (length(InsOpt) >= 1), Select = InsOpt(1); end
         if (length(InsOpt) >= 2), INSR = InsOpt(2); end
         if isnan(Select), Select = 0; end
         if isnan(INSR), INSR =1.0; end
      end
   end
   
   if (INSR < 0 || INSR > 1), error('Parameter for insertion rate must be a scalar in [0, 1]'); end
   if (INSR < 1 && IsObjVSel ~= 1), error('For selection of offspring ObjVSel is needed'); end 
   if (Select ~= 0 && Select ~= 1), error('Parameter for selection method must be 0 or 1'); end
   if (Select == 1 && IsObjVCh == 0), error('ObjVCh for fitness-based exchange needed'); end

   if INSR == 0, return; end
   NIns = min(max(floor(INSR*NSEL+.5),1),NIND);   % Number of offspring to insert   

% perform insertion for each subpopulation
   for irun = 1:SUBPOP,
      % Calculate positions in old subpopulation, where offspring are inserted
         if Select == 1,    % fitness-based reinsertion
            [Dummy, ChIx] = sort(-ObjVCh((irun-1)*NIND+1:irun*NIND));
         else               % uniform reinsertion
            [Dummy, ChIx] = sort(rand(NIND,1));
         end
         PopIx = ChIx((1:NIns)')+ (irun-1)*NIND;
      % Calculate position of Nins-% best offspring
         if (NIns < NSEL),  % select best offspring
            [Dummy,OffIx] = sort(ObjVSel((irun-1)*NSEL+1:irun*NSEL));
         else              
            OffIx = (1:NIns)';
         end
         SelIx = OffIx((1:NIns)')+(irun-1)*NSEL;
      % Insert offspring in subpopulation -> new subpopulation
         Chrom(PopIx,:) = SelCh(SelIx,:);
         if (IsObjVCh == 1 && IsObjVSel == 1), ObjVCh(PopIx) = ObjVSel(SelIx); end
   end
   
function  Pos = Find(FindVal,S)
% S=[1 3 2 3 1 2 1 3 2];
% FindVal=3;
[m,n] = size(S);
Pos = -1;
for i = 1:n 
    if FindVal == S(i)
      Pos = i;
      break;
    end
end