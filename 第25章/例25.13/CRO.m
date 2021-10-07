function [route,fval]=CRO(varargin)   %化学反应算法
type=varargin{end};
if type==1   %TSP问题
    city=varargin{1};
    popsize=varargin{2};
    iter_max=varargin{3};
    if nargin==4
        molecoll=0.2;alpha=500;beta=10;KElossrate=0.2;
    else
        molecoll=varargin{4};alpha=varargin{5};beta=varargin{6};KElossrate=varargin{7};
    end
    [NC,d]=city2d(city);
elseif type==2   %函数问题
    fun=varargin{1};
    popsize=varargin{2};
    iter_max=varargin{3};
    LB=varargin{4};
    UB=varargin{5};
    if nargin==6
        molecoll=0.2;alpha=500;beta=10;KElossrate=0.2;
    else
        molecoll=varargin{6};alpha=varargin{7};beta=varargin{8};KElossrate=varargin{9};
    end
    NC=size(LB,1);
end  
buffer=0;
for i=1:popsize
   if type==1
      molecule(i).route=randperm(NC);
      molecule(i).PE=value(molecule(i).route,d);
      y(i)=value(molecule(i).route,d);
   elseif type==2
      molecule(i).route=LB'+(UB-LB)'.*rand(1,NC);
      y(i)=fun(molecule(i).route);
      molecule(i).PE=fun(molecule(i).route);
   end
   molecule(i).KE=1000;
   molecule(i).Hitnum=0;  
end
[a,b]=min(y);
minmolecule=molecule(b);
success='false';
for iter=1:iter_max
    if rand>molecoll||popsize==1   %单分子反应
        num=ceil(popsize*rand);
        molecule(num).Hitnum=molecule(num).Hitnum+1;
        if molecule(num).Hitnum-minmolecule.Hitnum>alpha  %分解
            new1=molecule(num);new2=molecule(num);
            if type==1
               n=ceil(NC*rand);
               k=randperm(NC-n);
               aa1=new1.route(n+1:end);
               aa1=aa1(k);
               new1.route(n+1:end)=aa1;
               k=randperm(NC-n);
               aa2=new2.route(n+1:end);
               aa2=aa2(k);
               new2.route(n+1:end)=aa2;
               new1.PE=value(new1.route,d);
               new2.PE=value(new2.route,d);  
            elseif type==2
                for i=1:ceil(NC/2)
                    n1=ceil(NC*rand(1,2));
                    new1.route(n1(1))=new1.route(n1(1))+0.03*randn;
                    new2.route(n1(2))=new2.route(n1(2))+0.03*randn;
                end
                new1.route=boundtest(new1.route,LB,UB);
                new2.route=boundtest(new2.route,LB,UB);
                new1.PE=fun(new1.route);
                new2.PE=fun(new2.route);  
            end
            temp=molecule(num).PE+molecule(num).KE-new1.PE-new2.PE;
            if temp>=0
                temp1=rand;
                new1.KE=temp*temp1;
                new2.KE=temp*(1-temp1);
                success='true';                
            elseif temp+buffer>=0
                m1=rand;m2=rand;m3=rand;m4=rand;
                new1.KE=(temp+buffer)*m1*m2;
                new2.KE=(temp+buffer-new1.KE)*m3*m4;
                buffer=temp+buffer-new1.KE-new2.KE;
                success='true';
            end
            if strcmp(success,'true') 
                new1.Hitnum=0;
                new2.Hitnum=0;
                molecule=[molecule new1 new2];
                molecule=redu(molecule,num,'c');
                popsize=popsize+1;
                if new1.KE<minmolecule.KE
                   minmolecule=new1;
                end
                if new2.KE<minmolecule.KE
                   minmolecule=new2;
                end
                success='false';
            end
            clear k aa1 aa2
        else   % 无效碰撞
            new=molecule(num);
            if type==1
               pos=ceil(NC*rand(1,2));
               while pos(1)==pos(2)
                  pos=ceil(NC*rand(1,2));
               end
               temp2=new.route(pos(1));
               new.route(pos(1))=new.route(pos(2));
               new.route(pos(2))=temp2;
               new.PE=value(new.route,d);
            elseif type==2
                pos=ceil(NC*rand);
                new.route(pos)=new.route(pos)+0.03*randn;
                new.route=boundtest(new.route,LB,UB);
                new.PE=fun(new.route);
            end
            new.Hitnum=0;
            if molecule(num).PE+molecule(num).KE>=new.PE
                r=unifrnd(KElossrate,1);
                new.KE=(molecule(num).PE+molecule(num).KE-new.PE)*r;
                buffer=buffer+(molecule(num).PE+molecule(num).KE-new.PE)*(1-r);
                molecule(num)=new;
                if new.KE<minmolecule.KE
                   minmolecule=new;
                end
            end
        end
    else   %双分子反应
        num=sort(ceil(popsize*rand(1,2)));
        while num(1)==num(2)
            num=ceil(popsize*rand(1,2));
        end
        molecule(num(1)).Hitnum=molecule(num(1)).Hitnum+1;
        molecule(num(2)).Hitnum=molecule(num(2)).Hitnum+1;
        if molecule(num(1)).KE<beta&&molecule(num(2)).KE<beta   %合成反应
           new=molecule(num(1));
           if type==1
              pos=ceil(NC*rand);
              while pos==1||pos==NC
                   pos=ceil(NC*rand);
              end
              temp3=molecule(num(1)).route(1:pos);
              y2=findpos2(molecule(num(2)).route,temp3);          
              temp4=molecule(num(2)).route(sort(y2));
              temp5=redu(molecule(num(2)).route,y2,'c');
              f1=value(temp3,d);
              f2=value(temp4,d);
              f3=value(molecule(num(1)).route(pos+1:end),d);
              f4=value(temp5,d); 
              if f1<f2
                 new.route(1:pos)=temp3;
              else
                 new.route(1:pos)=temp4;
              end
              if f3<f4
                  new.route(pos+1:end)=molecule(num(1)).route(pos+1:end);
              else
                  new.route(pos+1:end)=temp5;
              end
              new.PE=value(new.route,d);
           elseif type==2
               for i=1:NC
                   if rand>0.5
                       new.route(i)=molecule(num(1)).route(i);
                   else
                       new.route(i)=molecule(num(2)).route(i);
                   end
               end
               new.route=boundtest(new.route,LB,UB);
               new.PE=fun(new.route);
           end
           new.KE=molecule(num(1)).PE+molecule(num(2)).PE+molecule(num(1)).KE+molecule(num(2)).KE-new.PE;
           new.Hitnum=0;
           if new.KE>=0
              molecule=[molecule new];
              molecule=redu(molecule,[num(1) num(2)],'c');
              popsize=popsize-1;
              if new.PE<minmolecule.PE
                  minmolecule=new;
              end
           end
        else   %交换
            new1=molecule(num(1));
            new2=molecule(num(2));
            if type==1
               pos=ceil(NC*rand);
               while pos==1||pos==NC
                   pos=ceil(NC*rand);
               end
               [y3,y4]=mycompare(molecule(num(2)).route,new1.route(1:pos));
               new1.route(pos+1:end)=y3;
               [y5,y6]=mycompare(molecule(num(1)).route,new2.route(1:pos));
               new2.route(pos+1:end)=y5;
               new1.PE=value(new1.route,d);
               new2.PE=value(new2.route,d);
            elseif type==2
                for i=1:ceil(NC/2)
                    n1=ceil(NC*rand(1,2));
                    new1.route(n1(1))=new1.route(n1(1))+0.03*randn;
                    new2.route(n1(2))=new2.route(n1(2))+0.03*randn;
                end
                new1.route=boundtest(new1.route,LB,UB);
                new2.route=boundtest(new2.route,LB,UB);
                new1.PE=fun(new1.route);
                new2.PE=fun(new2.route);
            end
            new1.Hitnum=0;
            new2.Hitnum=0;
            temp=molecule(num(1)).PE+molecule(num(2)).PE+molecule(num(1)).KE+molecule(num(2)).KE-(new1.PE+new2.PE);
            if temp>=0
                p=rand;
                new1.KE=temp*p;
                new2.KE=temp*(1-p);
                molecule(num(1))=new1;
                molecule(num(2))=new2;
                if new1.PE<minmolecule.PE
                   minmolecule=new1;
                end
                if new2.PE<minmolecule.PE
                   minmolecule=new2;
                end
            end
            clear y3 y5 
        end
    end
end
route=minmolecule.route;
fval=minmolecule.PE;
if type==1
   TSPplot(city,route);
end
            
            
            
            
        
        
            
                
            
            
            
            
                
            
            
        




