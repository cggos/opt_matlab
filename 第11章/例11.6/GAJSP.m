function MakeSpan=GAJSP(NIND,MAXGEN,GGAP,P_Cross,P_Mutation,Jm,T)
%====输入变量====
   %NIND           ：种群所包含个体的数目
   %MAXGEN         ：最大遗传代数
   %GGAP           ：代沟
   %P_Cross        ：交叉概率
   %P_Mutation     ：变异概率
   %Jm             ：各工序的可选机器集合，为m×n的元胞矩阵，m为机器数（工序），n为工件数
   %T              ：加工时间矩阵       
%====输出变量====
   %MakeSpan       ：最小的最大完工时间

gen=0;  %迭代计数器
JmNumber=Max_Cell(Jm);               %调用Max_Cell子函数求机器的数量
[PNumber,MNumber]=size(Jm);          %PNumber为工件个数，MNumber为工序个数
trace=zeros(2,MAXGEN);               %寻优结果的初始值，一行存放各代的最优解，一行存放各代解的均值
TotalOP_Number=PNumber*MNumber;      %工序总个数
%初始化
Number=zeros(1,PNumber);             %Number存放每个工件的工序数，PNumber工件个数
for i=1:PNumber
    Number(i)=MNumber;               %MNumber工序个数
end
%染色体个体编码：代码分为两层，第一层表示工序，第二层表示机器
Chrom=zeros(NIND,2*TotalOP_Number);  %染色体的长度为2*TotalOP_Number
for j=1:NIND                  %逐代寻优
    WPNumberTemp=Number;
    for i=1:TotalOP_Number    %随机产生工序
        val=unidrnd(PNumber);
        while WPNumberTemp(val)==0
            val=unidrnd(PNumber);
        end
        %第一层代码表示工序
        Chrom(j,i)=val;
        WPNumberTemp(val)=WPNumberTemp(val)-1;
        %第二层代码表示机器
        Temp=Jm{val,MNumber-WPNumberTemp(val)};
        SizeTemp=length(Temp);
        %随机产成工序机器
        Chrom(j,i+TotalOP_Number)=unidrnd(SizeTemp);
    end
end
%计算目标函数值
[PVal,ObjV,P,S]=cal(Chrom,JmNumber,T,Jm);                %调用子函数程序cal求个体适应度值
%循环寻找最优解
while gen<MAXGEN
    %分配适应度值
    FitnV=ranking(ObjV);                                 %调用子函数程序ranking进行排序操作
    %选择操作
    SelCh=select('rws', Chrom, FitnV, GGAP);             %调用子函数程序select进行选择操作，采用轮盘赌策略（rws）
    %交叉操作
    SelCh=across(SelCh,P_Cross,Jm,T);                    %调用子函数程序across进行交叉操作
    %变异操作
    SelCh=aberranceJm(SelCh,P_Mutation,Jm,T);            %调用子函数程序aberranceJm进行变异操作
    %重新计算目标适应度值（个体）
    [PVal,ObjVSel,P,S]=cal(SelCh,JmNumber,T,Jm);         %调用子函数程序cal求个体适应度值
    %重新插入新种群
    [Chrom,ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);    %调用子函数程序reins重新插入新种群
    %迭代计数器增加1
    gen=gen+1;
    %保存最优值
    trace(1,gen)=min(ObjV);
    trace(2,gen)=mean(ObjV);
    %记录最佳值
    if gen==1
        Val1=PVal;
        Val2=P;
        MinVal=min(ObjV);  %最小的最长流程加工时间
        STemp=S;
    end
    %记录最小的工序
    if MinVal>trace(1,gen)
        Val1=PVal;
        Val2=P;
        MinVal=trace(1,gen);
        STemp=S;
    end
end
%当前最佳值
PVal=Val1;    %工序时间
P=Val2;       %工序 
S=STemp;      %调度基因含机器基因
MakeSpan = min(trace(1,:));      %整批工件的完工时间（最长流程加工时间）的最小值
%显示最优解
MP=S(1,PNumber*MNumber+1:PNumber*MNumber*2);  %调出机器排序
for i=1:TotalOP_Number  
    val=P(1,i);
    a=(mod(val,100));   %工序
    b=((val-a)/100);    %工件
    Temp=Jm{b,a};
    mText=Temp(MP(1,i));
    x1=PVal(1,i);       %工序的开工时间
    x2=PVal(2,i);       %工序的完工时间
    y1=mText-0.5;       %0.5为甘特图矩形条的宽度，该数字可以更改
    y2=mText;           %用于绘制矩形条的高度
    PlotRec(x1,x2,mText); %调用PlotRec函数绘制甘特图函数
    hold on;
    fill([x1,x2,x2,x1],[y1,y1,y2,y2],[1-1/b,1/b,b/PNumber]);     %填充矩形条
    text((x1+x2)/2,mText-0.25,num2str(P(i)));                    %用文本标识工序，标注位置为矩形条的中心
end
function JmNumber = Max_Cell(Jm)
%Jm        机器的元胞数组
%JmNumber  机器的数量
MMAX=0;
[m,n]=size(Jm);   %元胞数组的行数与列数，行数m为工件数，列数n为工序数
for i=1:m         %逐行查找最大值
    for j=1:n     %逐列查找最大值
        TMAX=max(Jm{i,j});
        if MMAX<=TMAX;
            MMAX=TMAX;
        end
    end    
end
JmNumber=MMAX;

function [PVal,ObjV,P,S] = cal(Chrom,JmNumber,T,Jm)  % 子函数程序2：求目标函数值
% 输入参数
   %Chrom       基因种群  
   %JmNumber    机器数量
   %T           加工时间矩阵——各工件各工序使用的时间 
   %Jm          工序元胞矩阵（各工件各工序使用的机器） 
% 输出参数
   %PVal        最佳调度工序时间（包括开始时间与完工时间） 
   %P           最佳输出的调度工序 
   %ObjV        群中每个个体的调度工序时间
   %S           最佳输出的调度基因
%初始化
NIND=size(Chrom,1);
ObjV=zeros(NIND,1);

%计算工件个数、工序个数 
[PNumber,MNumber]=size(Jm);
%（3）迭代寻优
for i=1:NIND
    %取一个个体
    S=Chrom(i,:);
    %根据基因，计算调度工序
    P=calp(S,PNumber);                     %调用子函数程序calp计算适应度
    %根据调度工序，计算出调度工序时间
    PVal=caltime(S,P,JmNumber,T,Jm);       %调用子函数程序caltime        
    %取完成时间
    MT=max(PVal);
    TVal=max(MT); 
    %保存时间
    ObjV(i,1)=TVal;
    %初始化
    if i==1
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end
    %记录最小的调度工序时间、最佳调度工序时间 最佳输出的调度工序
    if MinVal>ObjV(i,1)
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end
end
%最佳调度工序时间,最佳输出的调度工序
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
TotalOP_Number = length(S);   %工序总个数
TotalOP_Number = TotalOP_Number/2;
%取工序基因，取基因的一半
S = S(1,1:TotalOP_Number);
%初始化
temp = zeros(1,PNumber);
P = zeros(1,TotalOP_Number);
%解码生成调度工序
for i = 1:TotalOP_Number
    %工序加+1
    temp(S(i)) = temp(S(i))+1;
    P(i) = S(i)*100+temp(S(i));  
end

    function PVal = caltime(S,P,JmNumber,T,Jm)
[PNumber,MNumber] = size(Jm);  %工件、工序个数
%取机器基因，取基因的一半（后面的一半）
M = S(1,PNumber*MNumber+1:PNumber*MNumber*2); 
TotalOP_Number = length(P);        %工序总个数
%初始化
TM = zeros(1,JmNumber);
TP = zeros(1,PNumber);
PVal = zeros(2,TotalOP_Number);    %存放工序的时间，一行为开工时间，一行为完工时间

%计算调度工序时间
for i = 1: TotalOP_Number 
    % 取机器号
    val = P(1,i);
    a = (mod(val,100));      %工序
    b = ((val-a)/100);       %工件
    Temp = Jm{b,a};
    m = Temp(M(1,i));
    
    %取加工时间
    Temp = T{b,a};
    t = Temp(M(1,i));        %工序的加工持续时间
    
    %取机器加工本工序的开始时间和前面一道工序的完成时间
        TMval = TM(1,m);
        TPval = TP(1,b);
    
    %机器加工本工序的开始时间：如果大于等于前面一道工序的完成时，取机器加工本工序的开始时间
    if TMval >= TPval 
        val = TMval;
    %否则，取前面一道工序的完成时间  
    else
        val = TPval;
    end
    %计算时间
    PVal(1,i) = val;    %开始时间
    PVal(2,i) = val+t;  %完成时间，如果t=0，则意味着实际上没有该工序
    
    %记录本次工序的机器时间和工序时间
    TM(1,m) = PVal(2,i);
    TP(1,b) = PVal(2,i); 
end

function SelCh = select(SEL_F,Chrom,FitnV,GGAP,SUBPOP)
   if nargin < 3, error('Not enough input parameter'); end

   % 确定种群的规模
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
   Nind = NindCh/SUBPOP;  % 计算每个子种群的个体数量

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

%计算新种群的数量以备选择
   NSel = max(floor(Nind*GGAP+.5),2);

%从种群中选择个体
   SelCh = [];
   for irun = 1:SUBPOP,
      FitnVSub = FitnV((irun-1)*Nind+1:irun*Nind);
      ChrIx = feval(SEL_F,FitnVSub,NSel)+(irun-1)*Nind;
      SelCh = [SelCh; Chrom(ChrIx,:)];
   end

function NewChrom = across(Chrom,P_Cross,Jm,T)
%初始化新种群
[NIND,TotalOP_Number] = size(Chrom);
TotalOP_Number = TotalOP_Number/2;
NewChrom = Chrom;

[PNumber,MNumber] = size(Jm);
Number = zeros(1,PNumber);
for i = 1:PNumber
    Number(i) = 1;
end

%随机选择交叉个体(洗牌交叉)
SelNum = randperm(NIND);   

Num = floor(NIND/2);    %交叉个体配对数
for i = 1:2:Num
    if P_Cross > rand;
        Pos = unidrnd(TotalOP_Number);%交叉位置
        while Pos == 1
            Pos = unidrnd(TotalOP_Number);
        end
        %取两交叉的个体并初始化个体
        S1 = Chrom(SelNum(i),1:TotalOP_Number);
        S2 = Chrom(SelNum(i+1),1:TotalOP_Number); 
        S11 = S2;S22 = S1;      %初始化新的个体     
        %新个体中间片断的COPY      
        S11(1:Pos) = S1(1:Pos);      
        S22(1:Pos) = S2(1:Pos);        
        %比较S11相对S1,S22相对S2多余和缺失的基因
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
          if S3(j) ~= 0             %多余的基因          
            Pos1 = find(S11 == S3(j),1);        
            Pos2 = find(S4,1);      %查找缺失的基因
            S11(Pos1) = S4(Pos2);   %用缺失的基因修补多余的基因
            S4(Pos2) = 0;       
          end 
          if S5(j) ~= 0              
            Pos1 = find(S22 == S5(j),1); 
            Pos2 = find(S6,1);           
            S22(Pos1) = S6(Pos2);
            S6(Pos2) = 0;
          end
        end

        % 保存交叉前的机器基因
        S1 = Chrom(SelNum(i),:);
        S2 = Chrom(SelNum(i+1),:); 
       
        for k = 1:TotalOP_Number            
            Pos1 = Find(S11(k),S1);                         %调用Find函数     
            S11(TotalOP_Number+k) = S1(TotalOP_Number+Pos1);
            S1(Pos1) = 0;
            Pos1 = Find(S22(k),S2);           
            S22(TotalOP_Number+k) = S2(TotalOP_Number+Pos1);
            S2(Pos1) = 0;
        end
        %生成新的种群
        NewChrom(SelNum(i),:) = S11;
        NewChrom(SelNum(i+1),:) = S22;
    end
end

function ChromNew = aberranceJm(Chrom,P_Mutation,Jm,T)
%初始化
[NIND,TotalOP_Number] = size(Chrom);
TotalOP_Number = TotalOP_Number/2;
ChromNew = Chrom;
[PNumber,MNumber] = size(Jm);
Number = zeros(1,PNumber);
for i = 1:PNumber
    Number(i) = 1;
end

for i = 1:NIND
    %取一个个体
    S = Chrom(i,:);
    WPNumberTemp = Number;
    for j = 1:TotalOP_Number
        JMTemp = Jm{S(j), WPNumberTemp(S(j))};
        SizeTemp = length(JMTemp);
        %是否变异
        if P_Mutation > rand;
            %选择机器（随机选择）
            %S(j+TotalOP_Number) = unidrnd(SizeTemp);
            %选择机器（加工时间少的选择几率大）
                if SizeTemp == 1
                    S(j+TotalOP_Number) = 1;
                else
                    S(j+TotalOP_Number) = selectJm(S(j++TotalOP_Number),T{S(j),WPNumberTemp(S(j))});    %调用selectJm函数
                end
        end
        WPNumberTemp(S(j)) = WPNumberTemp(S(j))+1;
    end
    %数据放入新群
    ChromNew(i,:)=S;
end

function PlotRec(mPoint1,mPoint2,mText)
%输入参数，用于控制矩形条大小的数字
grid on
vPoint = zeros(4,2);                   %绘制甘特图（用矩形条来表示）
vPoint(1,:) = [mPoint1,mText-0.5];     %依次给出矩形条的四个定点的坐标
vPoint(2,:) = [mPoint2,mText-0.5];
vPoint(3,:) = [mPoint1,mText];
vPoint(4,:) = [mPoint2,mText];
plot([vPoint(1,1),vPoint(2,1)],[vPoint(1,2),vPoint(2,2)])   %绘制水平线，高度为mText-0.5，长度为vPoint(2,1)－vPoint(1,1)
hold on;
plot([vPoint(1,1),vPoint(3,1)],[vPoint(1,2),vPoint(3,2)])   %绘制垂直线，高度为0.5，横坐标为vPoint(1,1)=vPoint(3,1)
plot([vPoint(2,1),vPoint(4,1)],[vPoint(2,2),vPoint(4,2)])   %绘制垂直线，高度为0.5，横坐标为vPoint(2,1)=vPoint(4,1)
plot([vPoint(3,1),vPoint(4,1)],[vPoint(3,2),vPoint(4,2)])   %绘制水平线，高度为mText，长度为vPoint(4,1)－vPoint(3,1)=vPoint(2,1)－vPoint(1,1)
title('工件加工的甘特图');
xlabel('加工时间（时间单位）');
ylabel('加工机器');

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