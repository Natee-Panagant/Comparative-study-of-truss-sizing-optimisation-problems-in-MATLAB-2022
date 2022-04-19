function rst = MOLA(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
rand('state',sum(100*clock));
%
d=nvar;         %problem dimension
pop=nsol;       %Population
n_iter = nloop; %Max number os iterations/gerations
LB=a';          %lower bounds
UB=b';          %upper bounds
Nr = narchive;  %Maximum number of solutions in PF

ref = 0.4;      %if more than zero, a second LF is created with refinement % the size of the other
Np = 100000;    %Number of Particles (If 3D, better more than 10000)
Rc = 150;       %Creation Radius (if 3D, better be less than 80, untill 150)
S_c = 1;           %Stick Probability: Percentage of particles that can don´t stuck in the
                   %cluster. Between 0 and 1. Near 0 there are more aggregate, the density of
                   %cluster is bigger and difusity is low. Near 1 is the opposite. 
M = 0;             %If M = 0, no lichtenberg figure is created (it is loaded a optimized figure); if 1, a single is created and used in all iterations; If 2, one is created for each iteration.(creating an LF figure takes about 2 min)
ngrid = 30;        %Number of grids in each dimension
IntCon = [];       %Avoid if there are no variables that must be integers. Ex.: IntCon = [1,2];
fnonlin = [];

% Number of objectives
xr=LB+(UB-LB).*rand(1,nvar);
[fr,gr]=feval(fun,xr');
m=numel(fr);
NC=numel(gr);
%
POS=zeros(nsol,d);
POS_fit=zeros(nsol,m);
POS_f=zeros(nsol,m);
POS_g=zeros(nsol,NC);
for i=1:pop
    POS(i,:)=LB+(UB-LB).*rand(1,d);
    [POS_f(i,:),POS_g(i,:)]=feval(fun,POS(i,:)');
    POS_fit(i,1:m)=fpenal00(POS_f(i,:)',POS_g(i,:)');
end
%
PBEST    = POS;
PBEST_fit= POS_fit;
PBEST_f= POS_f;
PBEST_g= POS_g;
%
DOMINATED= checkDomination(POS_fit);
REP.pos  = POS(~DOMINATED,:);
REP.pos_fit = POS_fit(~DOMINATED,:);
REP.pos_f   = POS_f(~DOMINATED,:);
REP.pos_g   = POS_g(~DOMINATED,:);
%
REP      = updateGrid(REP,ngrid);
if M==0 || M==1
LF = LA_figure(d,Np,Rc,S_c,M);
end
%
iter=0;
Neval=0;
rst.Neval_History=[];
while Neval<pop*n_iter
    iter=iter+1;
    x_start = REP.pos(randperm(size(REP.pos,1),1),:);
    if M == 2
        LF = LA_figure(d,Np,Rc,S_c,M);
    end
    scale_factor = 1.2*rand;
    X = LA_points(LF,LB,UB,x_start,scale_factor,d);
    if ref ~=0
        X_local = LA_points(LF,LB*ref,UB*ref,x_start,scale_factor,d);
    end
    for i=1:pop
        if ref ~=0
            pop1 = round((0.4)*pop);
            pop2 = pop-pop1;
            S_global = X(randperm(length(X),pop2),:);
            S_ref = X_local(randperm(length(X),pop1),:);
            POS = [S_global;S_ref];
        else
            POS = X(randperm(length(X),pop),:);
        end
        
    end
    for kk=1:length(POS)
        index1 = find(POS(kk,:) > UB);
        index2 = find(POS(kk,:) < LB);
        POS(kk,index1) = UB(index1);
        POS(kk,index2) = LB(index2);
    end
    %
    for i=1:pop
        [POS_f(i,:),POS_g(i,:)]=feval(fun,POS(i,:)');
        POS_fit(i,1:m)=fpenal00(POS_f(i,:)',POS_g(i,:)');
    end
    Neval=Neval+pop;
    %
    REP = updateRepository(REP,POS,POS_fit,POS_f,POS_g,ngrid);
    if(size(REP.pos,1)>Nr)
        REP = deleteFromRepository(REP,size(REP.pos,1)-Nr,ngrid);
    end
    pos_best = dominates(POS_fit, PBEST_fit);
    best_pos = ~dominates(PBEST_fit, POS_fit);
    best_pos(rand(pop,1)>=0.5) = 0;
    if(sum(pos_best)>1)
        PBEST_fit(pos_best,:) = POS_fit(pos_best,:);
        PBEST(pos_best,:) = POS(pos_best,:);
    end
    if(sum(best_pos)>1)
        PBEST_fit(best_pos,:) = POS_fit(best_pos,:);
        PBEST(best_pos,:) = POS(best_pos,:);
    end

    % Save Results
    ppareto=REP.pos';
    fppareto=REP.pos_fit';
    fpareto=REP.pos_f';
    gpareto=REP.pos_g';
    fea_idx=max(gpareto,[],1)<=0;
    rst.ppareto{iter}=ppareto(:,fea_idx);
    rst.fppareto{iter}=fppareto(:,fea_idx);
    rst.fpareto{iter}=fpareto(:,fea_idx);
    rst.gpareto{iter}=gpareto(:,fea_idx);

    [~,uidx]=unique(rst.ppareto{iter}','rows');
    rst.ppareto{iter}=rst.ppareto{iter}(:,uidx);
    rst.fppareto{iter}=rst.fppareto{iter}(:,uidx);
    rst.fpareto{iter}=rst.fpareto{iter}(:,uidx);
    rst.gpareto{iter}=rst.gpareto{iter}(:,uidx);

    rst.Neval_History=[rst.Neval_History,Neval];
    rst.timestamp=datetime('now');
%     pareto_plot(rst.fppareto{iter},rst.fpareto{iter},rst.gpareto{iter},fun);
    %
end

function REP = updateRepository(REP,POS,POS_fit,POS_f,POS_g,ngrid)
    DOMINATED  = checkDomination(POS_fit);
    REP.pos    = [REP.pos; POS(~DOMINATED,:)];
    REP.pos_fit= [REP.pos_fit; POS_fit(~DOMINATED,:)];
    REP.pos_f  = [REP.pos_f; POS_f(~DOMINATED,:)];
    REP.pos_g = [REP.pos_g; POS_g(~DOMINATED,:)];
    %
    DOMINATED  = checkDomination(REP.pos_fit);
    REP.pos_fit= REP.pos_fit(~DOMINATED,:);
    REP.pos    = REP.pos(~DOMINATED,:);
    REP.pos_f= REP.pos_f(~DOMINATED,:);
    REP.pos_g= REP.pos_g(~DOMINATED,:);
    REP        = updateGrid(REP,ngrid);

function dom_vector = checkDomination(fitness)
    Np = size(fitness,1);
    dom_vector = zeros(Np,1);
    all_perm = nchoosek(1:Np,2);    % Possible permutations
    all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];
    d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
    dominated_particles = unique(all_perm(d==1,2));
    dom_vector(dominated_particles) = 1;

function d = dominates(x,y)
    d = all(x<=y,2) & any(x<y,2);
    
function REP = updateGrid(REP,ngrid)
    ndim = size(REP.pos_fit,2);
    REP.hypercube_limits = zeros(ngrid+1,ndim);
    for dim = 1:1:ndim
        REP.hypercube_limits(:,dim) = linspace(min(REP.pos_fit(:,dim)),max(REP.pos_fit(:,dim)),ngrid+1)';
    end
    npar = size(REP.pos_fit,1);
    REP.grid_idx = zeros(npar,1);
    REP.grid_subidx = zeros(npar,ndim);
    for n = 1:1:npar
        idnames = [];
        for d = 1:1:ndim
            REP.grid_subidx(n,d) = find(REP.pos_fit(n,d)<=REP.hypercube_limits(:,d)',1,'first')-1;
            if(REP.grid_subidx(n,d)==0), REP.grid_subidx(n,d) = 1; end
            idnames = [idnames ',' num2str(REP.grid_subidx(n,d))];
        end
        REP.grid_idx(n) = eval(['sub2ind(ngrid.*ones(1,ndim)' idnames ');']);
    end
    REP.quality = zeros(ngrid,2);
    ids = unique(REP.grid_idx);
    for i = 1:length(ids)
        REP.quality(i,1) = ids(i);  
        REP.quality(i,2) = 10/sum(REP.grid_idx==ids(i));
    end
  
function REP = deleteFromRepository(REP,n_extra,ngrid)
    crowding = zeros(size(REP.pos,1),1);
    for m = 1:1:size(REP.pos_fit,2)
        [m_fit,idx] = sort(REP.pos_fit(:,m),'ascend');
        m_up     = [m_fit(2:end); Inf];
        m_down   = [Inf; m_fit(1:end-1)];
        distance = (m_up-m_down)./(max(m_fit)-min(m_fit));
        [~,idx]  = sort(idx,'ascend');
        crowding = crowding + distance(idx);
    end
    crowding(isnan(crowding)) = Inf;
    [~,del_idx] = sort(crowding,'ascend');
    del_idx = del_idx(1:n_extra);
    REP.pos(del_idx,:) = [];
    REP.pos_fit(del_idx,:) = [];
    REP.pos_f(del_idx,:) = [];
    REP.pos_g(del_idx,:) = [];
    %
    REP = updateGrid(REP,ngrid); 
    
function [X]=LA_points(K,LB,UB,x0,scale_factor,d)
if d ~= 3
    if d < 3 
        p=1; 
    else
        p=3; 
    end
    teta = p*rand;
    for u =1:length(K)
    K(u,1)= K(u,1)*cos(teta)-K(u,2)*sin(teta);
    K(u,2)= K(u,1)*sin(teta)+K(u,2)*cos(teta);
    end
end
if d == 3
    teta = rand; 
    alfa = rand;
    beta = rand;
    for u =1:length(K)
    K(u,1)= K(u,1)*cos(teta)*cos(beta)+K(u,2)*(cos(alfa)*sin(beta)+sin(alfa)*sin(teta)*cos(beta))+K(u,3)*(sin(alfa)*sin(beta)-cos(alfa)*sin(teta)*cos(beta));
    K(u,2)= K(u,1)*(-cos(teta))*sin(beta)+K(u,2)*(cos(alfa)*cos(beta)-sin(alfa)*sin(teta)*sin(beta))+K(u,3)*(sin(alfa)*cos(beta)+cos(alfa)*sin(teta)*sin(beta));
    K(u,3)= K(u,1)*sin(teta)+K(u,2)*(-sin(alfa)*cos(teta))+K(u,3)*cos(alfa)*cos(teta);
    end
end
Xi = zeros(length(K),d);
if d < 4
    for j=1:d
    Xi(:,j) = K(:,j);
    end
end
if d > 3
    for i=1:2:d+1
          gama = rand; 
          Xi(:,i) = K(:,1)*cos(gama)-K(:,2)*sin(gama);
          Xi(:,i+1) = K(:,1)*sin(gama)+K(:,2)*cos(teta);  
    end
    Xi = Xi(:,1:d);
end
for i = 1:d
scale(i) = scale_factor*(UB(i)-LB(i))/(max(Xi(:,i))-min(Xi(:,i)));
Xi(:,i) = scale(i)*Xi(:,i);
end
for i =1:d
    Pcc(i)=(max(Xi(:,i))-min(Xi(:,i)))/2 + min(Xi(:,i));
    Xi(round(length(K)/2),i)= Pcc(i);
    delta(i)=Pcc(i)-x0(i);
end
X=zeros(size(Xi));
for i=1:d
    X(:,i) = Xi(:,i) - delta(i);
end

function [map]=LA_figure(d,Np,Rc,S,M)
if M==1 || M==2
Rk = Rc*1.1;
add = Rk+2;
particle = 0;
stuck = 0;
escape = 0;
die = 0;
walk = 0;
if d<3 || d>3
map = zeros(((Rk+2)*2));
map(add,add)=1;
    X = 0; x = 0;
    Y = 0; y = 0;
    while (( particle >= Np) + (escape))==0 
    particle=particle+1; 
    phi=rand*2*3.14159265359; 
    X=Rc*cos(phi);
    Y=Rc*sin(phi);
    x=round(X);
    y=round(Y);
    stuck = 0; 
    die = 0;
    while ((stuck+die+escape) == 0)
    walk=rand;
    if walk<.25      
        if map(add+x,add+y+1)==0 
        y=y+1;
        end
    elseif walk<.5    
    	if map(add+x+1,add+y)==0 
        x=x+1;
        end
    elseif walk<.75   
        if map(add+x,add+y-1)==0 
        y=y-1;
        end
    else             
        if map(add+x-1,add+y)==0 
        x=x-1;
        end
    end
    if (hypot(x,y)>=Rk) 
    die=1;
    else
    stuck=0;
            if (map(add+x+1,add+y) + map(add+x-1,add+y) + map(add+x,add+y-1) + map(add+x,add+y+1))~=0 
                if (rand<S) 
                	stuck=1;
                end
            end 
          end 
    end 
    if stuck 
        map(add+x,add+y)=1; 
        stuck=0;
        clf
    if ((hypot(x,y)*1.2)>=Rc)
        escape = 1;
      end 
    end 

    end 
    if (escape==1) 
        disp('The cluster has reached the creation radius');
    end
    return
end
if d==3
map = zeros((2*Rk+4),(2*Rk+4),(2*Rk+4));
map(add,add,add)=1;
    X = 0; x = 0;
    Y = 0; y = 0;
    Z = 0; z=  0;
    while (( particle >= Np) + (escape))==0 
    particle=particle+1; 
    alfa=rand*2*pi;
    beta=rand*2*pi;
    X=Rc*sin(alfa)*cos(beta);
    Y=Rc*sin(alfa)*sin(beta);
    Z=Rc*cos(alfa);
    x=round(X);
    y=round(Y);
    z=round(Z);
    stuck = 0; 
    die = 0;
    while ((stuck+die+escape) == 0)
    walk = 1.5*rand;
    if walk<.25       
        if map(add+x,add+y+1, add+z)==0 
        y=y+1;
        end
    elseif walk<.5    
        if map(add+x+1,add+y, add+z)==0 
        x=x+1;
        end
	elseif walk<.75    
        if map(add+x,add+y-1, add+z)==0 
        y=y-1;
        end
    elseif walk<1      
         if map(add+x-1,add+y, add+z)==0 
         x=x-1;
         end
    elseif walk<1.25   
         if map(add+x,add+y, add+z+1)==0 
         z=z+1;
         end
    else              
        if map(add+x, add+y, add+z-1)==0
        z=z-1;
        end
    end
    if (sqrt(abs(x).^2+abs(y).^2+abs(z).^2)>=Rk) 
    die=1;
    else
    stuck=0;
        	if (map(add+x+1,add+y, add+z) + map(add+x-1,add+y, add+z) + map(add+x,add+y+1, add+z) + map(add+x,add+y-1,add+z) + map(add+x,add+y,add+z+1) + map(add+x,add+y,add+z-1))~=0
            	if (rand<S) 
                	stuck=1;
                end
            end 
        end 
    end 
    if stuck 
     map(add+x,add+y,add+z)=1;
     stuck=0;
     clf
    if (sqrt(abs(x).^2+abs(y).^2+abs(z).^2)*1.2>=Rc)
        escape = 1;
    end 
    end 
    end 
    if (escape==1) 
    	disp('The simulation ended before all particle could be tried because boundaries were exceeded');
    end
end
end
if M==0
    if d<3 | d>3
        load('LFND');
        map=LFND;
    else
        load('LF3D');
        map=LF3D;
    end
end      

function z=Fun(fhandle,fnonlin,u)
z=fhandle(u);
z=z+getconstraints(fnonlin,u);

function Z=getconstraints(fnonlin,u)
PEN=10^15;
lam=PEN; lameq=PEN;
Z=0;
[g,geq]=fnonlin(u);
for k=1:length(g),
    Z=Z+ lam*g(k)^2*getH(g(k));
end
for k=1:length(geq),
   Z=Z+lameq*geq(k)^2*geteqH(geq(k));
end

function H=getH(g)
if g<=0, 
    H=0; 
else
    H=1; 
end

function H=geteqH(g)
if g==0,
    H=0;
else
    H=1; 
end
