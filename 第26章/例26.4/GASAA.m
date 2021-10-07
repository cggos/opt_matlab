function [JbValue,U_Matrix,A_Matrix,Center] = GASAA(X,CN,T0,Tend,q,sizepop,MAXGEN,pc,pm)
   %X             �����������������
   %CN            ������������  
   %T0            ������
   %Tend          ����ֹ�¶�
   %q             ������ϵ��
   %sizepop       ����Ⱥ��Ŀ
   %MAXGEN        ������������
   %pc            ���������
   %pm            ���������
%====�������====
   %Jbvalue       ��Ŀ�꺯��ֵ
   %U_Matrix      �������Ⱦ���
   %A_Matrix      ������
   %Center        �����ĵ�
m = size(X,2);
%���ĵ㷶Χ[lb;ub]
lb = min(X);
ub = max(X);
options=[3,20,1e-6];
NVAR = m*CN;                   %������ά��
PRECI = 10;                    %�����Ķ�����λ��(Precision of variables)
trace = zeros(NVAR+1,MAXGEN);
FieldD = [rep(PRECI,[1,NVAR]);rep([lb;ub],[1,CN]);rep([1;0;1;1],[1,NVAR])];   %���ƣ�����rep������������gatbx�еĺ�����
Chrom = crtbp(sizepop, NVAR*PRECI);    %������ʼ��Ⱥ������crtbp������������gatbx�еĺ�����
V = bs2rv(Chrom, FieldD);              %ת��������bs2rv������������gatbx�еĺ�����
ObjV = ObjFun(X,CN,V,options);         %�����ܹ�ObjFun���������ʼ��Ⱥ�����Ŀ�꺯��ֵ
T = T0;
while T > Tend     %ģ���˻����ִ�е�����
    gen = 0;                                         %��������
    while gen < MAXGEN                             %�����������е�����
        FitnV = ranking(ObjV);                       %������Ӧ��ֵ������ranking������������gatbx�еĺ�����
        SelCh = select('sus', Chrom, FitnV);         %ѡ�񣬵���select����������sus������������gatbx�еĺ�����
        SelCh = recombin('xovsp', SelCh,pc);         %���飬����recombin����������xovsp������������gatbx�еĺ�����
        SelCh = mut(SelCh,pm);                       %���죬����mut������������gatbx�еĺ�����
        V = bs2rv(SelCh, FieldD);                    %ת��������bs2rv����
        newObjV = ObjFun(X,CN,V,options);            %�����Ӵ�Ŀ�꺯��ֵ
        newChrom = SelCh;

        %�Ƿ��滻�ɸ���
        for i = 1:sizepop
            if ObjV(i) > newObjV(i)
                ObjV(i) = newObjV(i);
                Chrom(i,:) = newChrom(i,:);
            else
                p = rand;
                if p <= exp((newObjV(i)-ObjV(i))/T)
                    ObjV(i) = newObjV(i);
                    Chrom(i,:) = newChrom(i,:);
                end
            end
        end
        gen = gen+1;                                  %������������
        [trace(end,gen),index] = min(ObjV);           %�Ŵ��㷨���ܸ���
        trace(1:NVAR,gen) = V(index,:);
        %fprintf(1,'���ڽ��е����㣺��%d�ε��� \n',gen);
    end
    T = T*q;
    %fprintf('�˻��¶�Ϊ:%1.3f\n',T);
end

[newObjV,center,U] = ObjFun(X,CN,trace(1:NVAR,end)',options);  %������ѳ�ʼ�������ĵ�Ŀ�꺯��ֵ
%��3���鿴������
JbValue = newObjV;        %Ŀ��ֵ
U_Matrix = U{1};          %�����Ⱦ���
Center = center{1};       %���������λ��
figure(1)                 %���������������ͼ�δ���
plot(X(:,1),X(:,2),'o')
hold on
maxU = max(U_Matrix);
Aindex = cell(CN,1);
for i = 1:CN
    Aindex{i} = find(U_Matrix(i,:) == maxU);
end
Newindex = Aindex;
for i = 1:CN
    figure(i+1)
    subplot(1,2,1)
    line(X(Newindex{i},1), X(Newindex{i},2), 'linestyle', 'none','marker', 'x', 'color', 'b');
    hold on
    plot(Center(i,1),Center(i,2),'diamond','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
    title(['GA-SAA����' num2str(i) '�ķֲ�ͼ']);
    xlabel('������');ylabel('������');
    text(Center(i,1),Center(i,2),'����λ��');
    hold off
end
A_Matrix0 = cell(CN,1);
for i = 1:CN
    A_Matrix0{i} = [X(Newindex{i},1), X(Newindex{i}, 2)];
end
A_Matrix = A_Matrix0;

function [Jb,center,U]=ObjFun(X,cn,V,options)
%% ������Ⱥ��ÿ�������Ŀ��ֵ
% ����
%        X����������
%       cn��������
%        V�����еĳ�ʼ�������ľ���
%  options��������ָ����������������Ŀ�꺯������ֹ����
% ���
%    Jb���������Ŀ�����
% center���Ż���ĸ�����ľ�������
%      U�������������Ʒ������
[sizepop,m] = size(V);
ch = m/cn;
Jb = zeros(sizepop,1);
center = cell(sizepop,1);
U = cell(sizepop,1);
for i = 1:sizepop
    v = reshape(V(i,:),cn,ch);
    [Jb(i),center{i},U{i}] = FCMfun(X,cn,v,options);  %�����Ӻ���FCMfun
end

function [obj,center,U]=FCMfun(X,cluster_n,center,options)
%        X����������
%cluster_n��������
%   center����ʼ�������ľ���
%  options��������ָ����������������Ŀ�꺯������ֹ����
% ���
%    obj��Ŀ�����Jbֵ
% center���Ż���ľ�������
%      U�����Ʒ������
X_n = size(X,1);
in_n = size(X,2);
b = options(1);		    % ��Ȩ����
max_iter = options(2);		% ����������
min_impro = options(3);		% �������ε�����С�Ľ��������ж��Ƿ���ǰ��ֹ��
obj_fcn = zeros(max_iter,1);	% ��ʼ��Ŀ��ֵ����
U = initFCM(X,cluster_n,center,b);			% ��ʼ���������ƾ���
% ������ѭ��
for i = 1:max_iter
    [U, center,obj_fcn(i)] = iterateFCM(X,U,cluster_n,b);
    % �˶���ֹ����
    if i > 1
        if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end
    end
end
iter_n = i;	% ��ʵ��������
obj_fcn(iter_n+1:max_iter)=[];
obj=obj_fcn(end);

function U=initFCM(X,cluster_n,center,b)
% ����
%        X����������
%cluster_n��������
%   center����ʼ�������ľ���
%        b��������ָ��
% ���
%      U�����Ʒ������
dist=distfcm(center,X);       % �������������������ĵľ������
% �����µ�U����
tmp=dist.^(-2/(b-1));
U=tmp./(ones(cluster_n,1)*sum(tmp));

function [U_new,center,obj_fcn]=iterateFCM(X,U,cluster_n,b)
% ����
%        X����������
%        U�����Ʒ������
%cluster_n��������
%        b����ָ��
% ���
%obj_fcn����ǰĿ�����Jbֵ
% center���µĵľ�������
%  U_new�����Ʒ������
mf = U.^b;       % ָ���������mf����
center = mf*X./((ones(size(X,2),1)*sum(mf'))'); % �µľ�������
% Ŀ��ֵ
dist = distfcm(center,X);           % �������������������ĵľ������
obj_fcn = sum(sum((dist.^2).*mf));  % Ŀ�꺯��ֵ
% �����µ�U����
tmp = dist.^(-2/(b-1));
U_new = tmp./(ones(cluster_n,1)*sum(tmp));


function MatOut = rep(MatIn,REPN)
% Get size of input matrix
   [N_D,N_L] = size(MatIn);
% Calculate
   Ind_D = rem(0:REPN(1)*N_D-1,N_D) + 1;
   Ind_L = rem(0:REPN(2)*N_L-1,N_L) + 1;

% Create output matrix
   MatOut = MatIn(Ind_D,Ind_L);

function [Chrom, Lind, BaseV] = crtbp(Nind, Lind, Base)
nargs = nargin ;
% Check parameter consistency

if nargs >= 1, [mN, nN] = size(Nind) ; end
if nargs >= 2, [mL, nL] = size(Lind) ; end
if nargs == 3, [mB, nB] = size(Base) ; end

if nN == 2
   if (nargs == 1) 
      Lind = Nind(2) ; Nind = Nind(1) ; BaseV = crtbase(Lind) ;
   elseif (nargs == 2 && nL == 1) 
      BaseV = crtbase(Nind(2),Lind) ; Lind = Nind(2) ; Nind = Nind(1) ; 
   elseif (nargs == 2 && nL > 1) 
      if Lind ~= length(Lind), error('Lind and Base disagree'); end
      BaseV = Lind ; Lind = Nind(2) ; Nind = Nind(1) ; 
   end
elseif nN == 1
   if nargs == 2
      if nL == 1, BaseV = crtbase(Lind) ;
      else
          BaseV = Lind ; Lind = nL ; 
      end
   elseif nargs == 3
      if nB == 1, BaseV = crtbase(Lind,Base) ;               %�����Ӻ���crtbase��������gabtx�е��Ӻ�����
      elseif nB ~= Lind, error('Lind and Base disagree') ; 
      else BaseV = Base ;
      end
   end
else
   error('Input parameters inconsistent') ;
end

% Create a structure of random chromosomes in row wise order, dimensions
% Nind by Lind. The base of each chromosomes loci is given by the value
% of the corresponding element of the row vector base.

Chrom = floor(rand(Nind,Lind).*BaseV(ones(Nind,1),:)) ;

function BaseVec = crtbase(Lind, Base)

[ml LenL] = size(Lind) ;
if nargin < 2 
	Base = 2 * ones(LenL,1) ; % default to base 2
end
[mb LenB] = size(Base) ;
% check parameter consistency
if ml > 1 || mb > 1
	error( 'Lind or Base is not a vector') ;
elseif (LenL > 1 && LenB > 1 && LenL ~= LenB) || (LenL == 1 && LenB > 1 ) 
	error( 'Vector dimensions must agree' ) ;
elseif LenB == 1 && LenL > 1
	Base = Base * ones(LenL,1) ;	
end
BaseVec = [] ;
for i = 1:LenL
	BaseVec = [BaseVec, Base(i)*ones(Lind(i),1)'];
end

function Phen = bs2rv(Chrom,FieldD)
% Identify the population size (Nind)
%      and the chromosome length (Lind)
[Nind,Lind] = size(Chrom);

% Identify the number of decision variables (Nvar)
[seven,Nvar] = size(FieldD);

if seven ~= 7
	error('FieldD must have 7 rows.');
end

% Get substring properties
len = FieldD(1,:);
lb = FieldD(2,:);
ub = FieldD(3,:);
code = ~(~FieldD(4,:));
scale = ~(~FieldD(5,:));
lin = ~(~FieldD(6,:));
uin = ~(~FieldD(7,:));

% Check substring properties for consistency
if sum(len) ~= Lind
	error('Data in FieldD must agree with chromosome length');
end

if ~all(lb(scale).*ub(scale)>0)
	error('Log-scaled variables must not include 0 in their range');
end

% Decode chromosomes
Phen = zeros(Nind,Nvar);

lf = cumsum(len);
li = cumsum([1 len]);
Prec = .5 .^ len;

logsgn = sign(lb(scale));
lb(scale) = log( abs(lb(scale)) );
ub(scale) = log( abs(ub(scale)) );
delta = ub - lb;

Prec = .5 .^ len;
num = (~lin) .* Prec;
den = (lin + uin - 1) .* Prec;

for i = 1:Nvar
    idx = li(i):lf(i);
    if code(i) % Gray decoding
	    Chrom(:,idx)=rem(cumsum(Chrom(:,idx)')',2);
    end
    Phen(:,i) = Chrom(:,idx) * ( (.5).^(1:len(i))');
    Phen(:,i) = lb(i) + delta(i) * (Phen(:,i) + num(i)) ./ (1 - den(i));
end

expand = ones(Nind,1);
if any(scale)
	Phen(:,scale) = logsgn(expand,:) .* exp(Phen(:,scale));
end

function FitnV = ranking(ObjV, RFun, SUBPOP)
% Identify the vector size (Nind)
   [Nind,ans] = size(ObjV);

   if nargin < 2, RFun = []; end
   if nargin > 1, if isnan(RFun), RFun = []; end, end
   if numel(RFun) == 2
      if RFun(2) == 1, NonLin = 1;
      elseif RFun(2) == 0, NonLin = 0;
      else error('Parameter for ranking method must be 0 or 1'); 
      end
      RFun = RFun(1);
      if isnan(RFun), RFun = 2; end
   elseif numel(RFun) > 2
      if numel(RFun) ~= Nind, error('ObjV and RFun disagree'); end
   end

   if nargin < 3, SUBPOP = 1; end
   if nargin > 2
      if isempty(SUBPOP), SUBPOP = 1;
      elseif isnan(SUBPOP), SUBPOP = 1;
      elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar'); 
      end
   end

   if (Nind/SUBPOP) ~= fix(Nind/SUBPOP), error('ObjV and SUBPOP disagree'); end
   Nind = Nind/SUBPOP;  % Compute number of individuals per subpopulation
   
% Check ranking function and use default values if necessary
   if isempty(RFun)
      % linear ranking with selective pressure 2
         RFun = 2*(0:Nind-1)'/(Nind-1);
   elseif numel(RFun) == 1
      if NonLin == 1
         % non-linear ranking
         if RFun(1) < 1, error('Selective pressure must be greater than 1');
         elseif RFun(1) > Nind-2, error('Selective pressure too big'); end
         Root1 = roots([RFun(1)-Nind [RFun(1)*ones(1,Nind-1)]]);
         RFun = (abs(Root1(1)) * ones(Nind,1)) .^ [(0:Nind-1)'];
         RFun = RFun / sum(RFun) * Nind;
      else
         % linear ranking with SP between 1 and 2
         if (RFun(1) < 1 || RFun(1) > 2)
            error('Selective pressure for linear ranking must be between 1 and 2');
         end
         RFun = 2-RFun + 2*(RFun-1)*[0:Nind-1]'/(Nind-1);
      end
   end

   FitnV = [];

% loop over all subpopulations
for irun = 1:SUBPOP
   % Copy objective values of actual subpopulation
      ObjVSub = ObjV((irun-1)*Nind+1:irun*Nind);
   % Sort does not handle NaN values as required. So, find those...
      NaNix = isnan(ObjVSub);
      Validix = find(~NaNix);
   % ... and sort only numeric values (smaller is better).
      [ans,ix] = sort(-ObjVSub(Validix));

   % Now build indexing vector assuming NaN are worse than numbers,
   % (including Inf!)...
      ix = [find(NaNix) ; Validix(ix)];
   % ... and obtain a sorted version of ObjV
      Sorted = ObjVSub(ix);

   % Assign fitness according to RFun.
      i = 1;
      FitnVSub = zeros(Nind,1);
      for j = [find(Sorted(1:Nind-1) ~= Sorted(2:Nind)); Nind]'
         FitnVSub(i:j) = sum(RFun(i:j)) * ones(j-i+1,1) / (j-i+1);
         i =j+1;
      end

   % Finally, return unsorted vector.
      [ans,uix] = sort(ix);
      FitnVSub = FitnVSub(uix);

   % Add FitnVSub to FitnV
      FitnV = [FitnV; FitnVSub];
end

function SelCh = select(SEL_F, Chrom, FitnV, GGAP, SUBPOP)
% Check parameter consistency
   if nargin < 3, error('Not enough input parameter'); end

   % Identify the population size (Nind)
   [NindCh,Nvar] = size(Chrom);
   [NindF,VarF] = size(FitnV);
   if NindCh ~= NindF, error('Chrom and FitnV disagree'); end
   if VarF ~= 1, error('FitnV must be a column vector'); end
  
   if nargin < 5, SUBPOP = 1; end
   if nargin > 4
      if isempty(SUBPOP), SUBPOP = 1;
      elseif isnan(SUBPOP), SUBPOP = 1;
      elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar'); end
   end

   if (NindCh/SUBPOP) ~= fix(NindCh/SUBPOP), error('Chrom and SUBPOP disagree'); end
   Nind = NindCh/SUBPOP;  % Compute number of individuals per subpopulation

   if nargin < 4, GGAP = 1; end
   if nargin > 3
      if isempty(GGAP), GGAP = 1;
      elseif isnan(GGAP), GGAP = 1;
      elseif length(GGAP) ~= 1, error('GGAP must be a scalar');
      elseif (GGAP < 0), error('GGAP must be a scalar bigger than 0');
      end
   end

% Compute number of new individuals (to select)
   NSel=max(floor(Nind*GGAP+.5),2);

% Select individuals from population
   SelCh = [];
   for irun = 1:SUBPOP
      FitnVSub = FitnV((irun-1)*Nind+1:irun*Nind);
      ChrIx = feval(SEL_F, FitnVSub, NSel)+(irun-1)*Nind;
      SelCh = [SelCh; Chrom(ChrIx,:)];
   end

function NewChrIx = sus(FitnV,Nsel)
% Identify the population size (Nind)
   [Nind,ans] = size(FitnV);
% Perform stochastic universal sampling
   cumfit = cumsum(FitnV);
   trials = cumfit(Nind) / Nsel * (rand + (0:Nsel-1)');
   Mf = cumfit(:, ones(1, Nsel));
   Mt = trials(:, ones(1, Nind))';
   [NewChrIx, ans] = find(Mt < Mf & [ zeros(1, Nsel); Mf(1:Nind-1, :) ] <= Mt);

% Shuffle new population
   [ans, shuf] = sort(rand(Nsel, 1));
   NewChrIx = NewChrIx(shuf);

function NewChrom = recombin(REC_F, Chrom, RecOpt, SUBPOP)
% Check parameter consistency
   if nargin < 2, error('Not enough input parameter'); end
   % Identify the population size (Nind)
   [Nind,Nvar] = size(Chrom);
 
   if nargin < 4, SUBPOP = 1; end
   if nargin > 3
      if isempty(SUBPOP), SUBPOP = 1;
      elseif isnan(SUBPOP), SUBPOP = 1;
      elseif length(SUBPOP) ~= 1, error('SUBPOP must be a scalar');
      end
   end

   if (Nind/SUBPOP) ~= fix(Nind/SUBPOP), error('Chrom and SUBPOP disagree'); end
   Nind = Nind/SUBPOP;  % Compute number of individuals per subpopulation

   if nargin < 3, RecOpt = 0.7; end
   if nargin > 2
      if isempty(RecOpt), RecOpt = 0.7;
      elseif isnan(RecOpt), RecOpt = 0.7;
      elseif length(RecOpt) ~= 1, error('RecOpt must be a scalar');
      elseif (RecOpt < 0 || RecOpt > 1), error('RecOpt must be a scalar in [0, 1]'); 
      end
   end
% Select individuals of one subpopulation and call low level function
   NewChrom = [];
   for irun = 1:SUBPOP
      ChromSub = Chrom((irun-1)*Nind+1:irun*Nind,:);  
      NewChromSub = feval(REC_F, ChromSub, RecOpt);
      NewChrom = [NewChrom; NewChromSub];
   end

function NewChrom = xovsp(OldChrom, XOVR)

if nargin < 2, XOVR = NaN; end

% call low level function with appropriate parameters
   NewChrom = xovmp(OldChrom, XOVR, 1, 0);      %�����Ӻ���xovmp��������gabtx�е��Ӻ�����

function NewChrom = xovmp(OldChrom, Px, Npt, Rs)

% Identify the population size (Nind) and the chromosome length (Lind)
[Nind,Lind] = size(OldChrom);
if Lind < 2, NewChrom = OldChrom; return; end
if nargin < 4, Rs = 0; end
if nargin < 3, Npt = 0; Rs = 0; end
if nargin < 2, Px = 0.7; Npt = 0; Rs = 0; end
if isnan(Px), Px = 0.7; end
if isnan(Npt), Npt = 0; end
if isnan(Rs), Rs = 0; end
if isempty(Px), Px = 0.7; end
if isempty(Npt), Npt = 0; end
if isempty(Rs), Rs = 0; end
Xops = floor(Nind/2);
DoCross = rand(Xops,1) < Px;
odd = 1:2:Nind-1;
even = 2:2:Nind;
% Compute the effective length of each chromosome pair
Mask = ~Rs | (OldChrom(odd, :) ~= OldChrom(even, :));
Mask = cumsum(Mask')';
% Compute cross sites for each pair of individuals, according to their
% effective length and Px (two equal cross sites mean no crossover)
xsites(:, 1) = Mask(:, Lind);
if Npt >= 2
        xsites(:, 1) = ceil(xsites(:, 1) .* rand(Xops, 1));
end
xsites(:,2) = rem(xsites + ceil((Mask(:, Lind)-1) .* rand(Xops, 1)) .* DoCross - 1 , Mask(:, Lind) )+1;

% Express cross sites in terms of a 0-1 mask
Mask = (xsites(:,ones(1,Lind)) < Mask) == ...
                        (xsites(:,2*ones(1,Lind)) < Mask);

if ~Npt
        shuff = rand(Lind,Xops);
        [ans,shuff] = sort(shuff);
        for i=1:Xops
          OldChrom(odd(i),:)=OldChrom(odd(i),shuff(:,i));
          OldChrom(even(i),:)=OldChrom(even(i),shuff(:,i));
        end
end

% Perform crossover
NewChrom(odd,:) = (OldChrom(odd,:).* Mask) + (OldChrom(even,:).*(~Mask));
NewChrom(even,:) = (OldChrom(odd,:).*(~Mask)) + (OldChrom(even,:).*Mask);

% If the number of individuals is odd, the last individual cannot be mated
% but must be included in the new population
if rem(Nind,2)
  NewChrom(Nind,:)=OldChrom(Nind,:);
end

if ~Npt
        [ans,unshuff] = sort(shuff);
        for i=1:Xops
          NewChrom(odd(i),:)=NewChrom(odd(i),unshuff(:,i));
          NewChrom(even(i),:)=NewChrom(even(i),unshuff(:,i));
        end
end

function NewChrom = mut(OldChrom,Pm,BaseV)
% get population size (Nind) and chromosome length (Lind)
[Nind, Lind] = size(OldChrom) ;

% check input parameters
if nargin < 2, Pm = 0.7/Lind ; end
if isnan(Pm), Pm = 0.7/Lind; end
if (nargin < 3), BaseV = crtbase(Lind);  end
if (isnan(BaseV)), BaseV = crtbase(Lind);  end
if (isempty(BaseV)), BaseV = crtbase(Lind);  end
if (nargin == 3) && (Lind ~= length(BaseV))
   error('OldChrom and BaseV are incompatible'), end
% create mutation mask matrix
BaseM = BaseV(ones(Nind,1),:) ;
% perform mutation on chromosome structure
NewChrom = rem(OldChrom+(rand(Nind,Lind)<Pm).*ceil(rand(Nind,Lind).*(BaseM-1)),BaseM);

