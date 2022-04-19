function [xns,fns,gns]=resortp(x,f,g)
% Re-sort Pareto front
% Unique solutions
[~,uind]=unique(x','rows');
x=x(:,uind);
f=f(:,uind);
g=g(:,uind);
% Eliminate dominated and infeasible solutions
n=size(f,2);
A=zeros(n);
for i=1:n
    fi=f(:,i);
    gi=g(:,i);
    A(i,i)=0;
    for j=1:n
        if j~=i
            fj=f(:,j);
            gj=g(:,j);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            [p_count1,p_count2]=fdominated(fi,gi,fj,gj);
            A(i,j)=p_count1;
            A(j,i)=p_count2;
            %%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=sum(A,1);
Indm=B==0;
Ifea=max(g,[],1)<=0;
if numel(x)>0
    xns=x(:,Indm&Ifea);
else
    xns=[];
end
fns=f(:,Indm & Ifea);
gns=g(:,Indm & Ifea);
infea_ind=max(gns,[],1)>0;
fns(:,infea_ind)=[];
gns(:,infea_ind)=[];
end
%%%%%%%%%%%%%%%%
% Sub-Function %
%%%%%%%%%%%%%%%%
function [p1,p2]=fdominated(f1,g1,f2,g2)
n=length(f1);
mg1=max(g1);
mg2=max(g2);

icount11=0;
icount12=0;
icount21=0;
icount22=0;

if mg1<=0 && mg2<=0
    for i=1:n
        if f1(i) <= f2(i)
            icount11=icount11+1;
        end
        if f1(i) < f2(i)
            icount12=icount12+1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        if f2(i) <= f1(i)
            icount21=icount21+1;
        end
        if f2(i) < f1(i)
            icount22=icount22+1;
        end
    end
    if icount11 == n && icount12 > 0
        p1=1;
    else
        p1=0;
    end
    if icount21 == n && icount22 > 0
        p2=1;
    else
        p2=0;
    end
elseif mg1 <=0 && mg2 > 0
    p1=1;p2=0;
elseif mg2 <=0 && mg1 > 0
    p1=0;p2=1;
else
    if mg1 <= mg2
        p1=1;p2=0;
    else
        p1=0;p2=1;
    end
end
end