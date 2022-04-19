%-----------------------------------------------------------------------------------------------------------------%
%  Multi-Objective Jellyfish Search (MOJS) source codes demo version 1.0  Developed in MATLAB R2016a              %
%  Author and programmer:                                                                                         %
%         Professor        Jui-Sheng Chou                                                                         %
%         Ph.D. Candidate  Dinh- Nhat Truong                                                                      %
%  Paper:   Multiobjective optimization inspired by behavior of jellyfish for solving structural design problems, %
%  Chaos, Solitons & Fractals: An interdisciplinary journal of nonlinear science                                  %
%  DOI:  https://doi.org/10.1016/j.chaos.2020.109738                                                              %
%                                     PiM Lab, NTUST, Taipei, Taiwan, March-2020                                  %                                
%-----------------------------------------------------------------------------------------------------------------%
%% Main MOJS optimizer
function rst = MOJS(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
rand('state',sum(100*clock));
% function ARCH = MOJS(params,MultiObj)
% Parameters
Np      = nsol;
Nr      = narchive;
MaxIt   = nloop;
ngrid   = 20;%Default=20
% fun     = fun;
nVar    = nvar;
var_min = a;
var_max = b;

% Number of objectives
xr=a'+(b'-a').*rand(1,nvar);
[fr,gr]=feval(fun,xr');
m=numel(fr);
NC=numel(gr);

% Initialization by Eq. 25
POS=initialchaos(7,Np,nVar,var_max',var_min');
POS_fit=zeros(nsol,m);
POS_f=zeros(nsol,m);
POS_g=zeros(nsol,NC);
for i=1:Np
    [POS_f(i,:),POS_g(i,:)]=feval(fun,POS(i,:)');
    POS_fit(i,1:m)=fpenal00(POS_f(i,:)',POS_g(i,:)');
end

ELI_POS      = POS;
ELI_POS_fit  = POS_fit;
ELI_POS_f    = POS_f;
ELI_POS_g    = POS_g;
DOMINATED    = checkDomination(POS_fit);
ARCH.pos     = POS(~DOMINATED,:);
ARCH.pos_fit = POS_fit(~DOMINATED,:);
ARCH.pos_f = POS_f(~DOMINATED,:);
ARCH.pos_g = POS_g(~DOMINATED,:);
ARCH         = updateGrid(ARCH,ngrid);

%% Main MOJS loop
iter=1;
Neval=0;
rst.Neval_History=[];
while Neval<Np*MaxIt
    % Select leader by Eq. 16
    h = selectLeader(ARCH);
    % Calculate time control by Eq. 15
    Ct=abs((1-iter*((1)/MaxIt))*(2*rand-1));
    if Ct>=0.5
        Meanvl=mean(ELI_POS);
        for i=1:Np
            % The new position is determined by Eq.19 and Eq.20
            POS(i,:) = ELI_POS(i,:) + Levy(nVar).*(ARCH.pos(h,:) - 3*rand([1 nVar]).*Meanvl);
        end
    else
        for i=1:Np
            if rand<(1-Ct)
                % Jellyfish follow type B
                % Determine the direction by Eq. 24
                j=i;
                while j==i
                    j=randperm(Np,1);
                end
                Step = ELI_POS(i,:) - ELI_POS(j,:);
                if dominates(ELI_POS_fit(j,:),ELI_POS_fit(i,:))
                    Step = -Step;
                end
                % The new position is determined by Eq. 22
                POS(i,:) =ARCH.pos(h,:) + rand([1 nVar]).*Step;
            else
                % Jellyfish follow type A
                % The new position is determined by Eq. 21
                POS(i,:)=ARCH.pos(h,:)+Levy(nVar).*(ELI_POS(i,:)-ARCH.pos(h,:));
            end
        end
    end
    %% Update new position by opposition-based jumping using Eq. 26
    if rand <(iter/MaxIt)
        [POS] = OPPOS(POS,var_max,var_min);
    end
    %% Check boundaries
    POS=checksimplebounds(POS,var_min',var_max');
    %% Evaluate the population
    for i=1:Np
        [POS_f(i,:),POS_g(i,:)]=feval(fun,POS(i,:)');
        POS_fit(i,1:m)=fpenal00(POS_f(i,:)',POS_g(i,:)');
    end
    Neval=Neval+Np;
    %
    pos_best = dominates(POS_fit, ELI_POS_fit);
    best_pos = ~dominates(ELI_POS_fit, POS_fit);
    best_pos(rand(Np,1)>=0.5) = 0;
    if(sum(pos_best)>1)
        ELI_POS_fit(pos_best,:) = POS_fit(pos_best,:);
        ELI_POS(pos_best,:) = POS(pos_best,:);
        ELI_POS_f(pos_best,:) = POS_f(pos_best,:);
        ELI_POS_g(pos_best,:) = POS_g(pos_best,:);
    end
    if(sum(best_pos)>1)
        ELI_POS_fit(best_pos,:) = POS_fit(best_pos,:);
        ELI_POS(best_pos,:) = POS(best_pos,:);
        ELI_POS_f(best_pos,:) = POS_f(best_pos,:);
        ELI_POS_g(best_pos,:) = POS_g(best_pos,:);
    end
    %% Update the archive 
    if size(ARCH.pos,1)==1
        ARCH.pos= POS;
        ARCH.pos_fit= POS_fit;
        ARCH.pos_f= POS_f;
        ARCH.pos_g= POS_g;
        ARCH = updateArchive(ARCH,ELI_POS,ELI_POS_fit,ELI_POS_f,ELI_POS_g,ngrid);
    else
        ARCH = updateArchive(ARCH,ELI_POS,ELI_POS_fit,ELI_POS_f,ELI_POS_g,ngrid);
        if size(ARCH.pos,1)==1
            ARCH.pos= ELI_POS;
            ARCH.pos_fit= ELI_POS_fit;
            ARCH.pos_f= ELI_POS_f;
            ARCH.pos_g= ELI_POS_g;
        end
    end
    if(size(ARCH.pos,1)>Nr)
        % Delete the worst members from archive by Eq. 18
        ARCH = deleteFromArchive(ARCH,size(ARCH.pos,1)-Nr,ngrid);
    end
    
%     display(['Iteration #' num2str(it) ' - Archive size: ' num2str(size(ARCH.pos,1))]);

    % Save Results
    ppareto=ARCH.pos';
    fppareto=ARCH.pos_fit';
    fpareto=ARCH.pos_f';
    gpareto=ARCH.pos_g';
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
    
    iter=iter+1;
end
% %% Plotting paretofront
% if(size(ARCH.pos_fit,2)==2)
%     plot(ARCH.pos_fit(:,1),ARCH.pos_fit(:,2),'or'); hold on;
%     grid on; xlabel('f1'); ylabel('f2');
% end
% if(size(ARCH.pos_fit,2)==3)
%     plot3(ARCH.pos_fit(:,1),ARCH.pos_fit(:,2),ARCH.pos_fit(:,3),'or'); hold on;
%     grid on; xlabel('f1'); ylabel('f2'); zlabel('f3');
% end
end

%% This function calucates the leader performance by a roulette wheel selection
% based on the quality of each hypercube
function selected = selectLeader(ARCH)
% Roulette wheel
prob    = cumsum(ARCH.quality(:,2));     % Cumulated probs
sel_hyp = ARCH.quality(find(rand(1,1)*max(prob)<=prob,1,'first'),1); % Selected hypercube
% Select the index leader as a random selection inside that hypercube
idx      = 1:1:length(ARCH.grid_idx);
selected = idx(ARCH.grid_idx==sel_hyp);
selected = selected(randi(length(selected)));
end

%% This function returns 1 if x dominates y and 0 otherwise
function d = dominates(x,y)
d = all(x<=y,2) & any(x<y,2);
end

%% This function checks the domination inside the population.
function domi_vector = checkDomination(fitness)
Np = size(fitness,1);
if Np>2
    domi_vector = zeros(Np,1);
    all_perm = nchoosek(1:Np,2);    % Possible permutations
    all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];
    
    d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
    dominated_particles = unique(all_perm(d==1,2));
    domi_vector(dominated_particles) = 1;
else
    domi_vector=ones(Np,1);
end
end

%% This function updates the archive given a new population 
function ARCH = updateArchive(ARCH,POS,POS_fit,POS_f,POS_g,ngrid)
% Domination between jellyfish
DOMINATED  = checkDomination(POS_fit);
ARCH.pos    = [ARCH.pos; POS(~DOMINATED,:)];
ARCH.pos_fit= [ARCH.pos_fit; POS_fit(~DOMINATED,:)];
ARCH.pos_f= [ARCH.pos_f; POS_f(~DOMINATED,:)];
ARCH.pos_g= [ARCH.pos_g; POS_g(~DOMINATED,:)];
% Domination between nondominated jellyfish and the last archive 
DOMINATED  = checkDomination(ARCH.pos_fit);
ARCH.pos_fit= ARCH.pos_fit(~DOMINATED,:);
ARCH.pos    = ARCH.pos(~DOMINATED,:);
ARCH.pos_f= ARCH.pos_f(~DOMINATED,:);
ARCH.pos_g= ARCH.pos_g(~DOMINATED,:);
% Updating the grid
ARCH        = updateGrid(ARCH,ngrid);
end

%% Function that updates the hypercube grid, the hypercube where belongs
function ARCH = updateGrid(ARCH,ngrid)
% Computing the  hypercube limitation
ndim = size(ARCH.pos_fit,2);
ARCH.hypercube_limits = zeros(ngrid+1,ndim);
for dim = 1:1:ndim
    ARCH.hypercube_limits(:,dim) = linspace(min(ARCH.pos_fit(:,dim)),max(ARCH.pos_fit(:,dim)),ngrid+1)';
end
% Computing where belongs each jellyfish
npar = size(ARCH.pos_fit,1);
ARCH.grid_idx = zeros(npar,1);
ARCH.grid_subidx = zeros(npar,ndim);
for n = 1:1:npar
    idnames = [];
    for d = 1:1:ndim
        ARCH.grid_subidx(n,d) = find(ARCH.pos_fit(n,d)<=ARCH.hypercube_limits(:,d)',1,'first')-1;
        if(ARCH.grid_subidx(n,d)==0), ARCH.grid_subidx(n,d) = 1; end
        idnames = [idnames ',' num2str(ARCH.grid_subidx(n,d))];
    end
    ARCH.grid_idx(n) = eval(['sub2ind(ngrid.*ones(1,ndim)' idnames ');']);
end
% Quality based on the number of jellyfish in each hypercube
ARCH.quality = zeros(ngrid,2);
ids = unique(ARCH.grid_idx);
for i = 1:length(ids)
    ARCH.quality(i,1) = ids(i);                       
    ARCH.quality(i,2) = 10/sum(ARCH.grid_idx==ids(i)); 
end
end

%% This function deletes an excess of jellyfish inside the archive using crowding distances
function ARCH = deleteFromArchive(ARCH,n_extra,ngrid)
% Compute the crowding distances
crowding = zeros(size(ARCH.pos,1),1);
for m = 1:1:size(ARCH.pos_fit,2)
    [m_fit,idx] = sort(ARCH.pos_fit(:,m),'ascend');
    m_up     = [m_fit(2:end); Inf];
    m_down   = [Inf; m_fit(1:end-1)];
    distance = (m_up-m_down)./(max(m_fit)-min(m_fit));
    [~,idx]  = sort(idx,'ascend');
    crowding = crowding + distance(idx);
end
crowding(isnan(crowding)) = Inf;
% This function deletes the extra jellyfish with the smallest crowding distances
[~,del_idx] = sort(crowding,'ascend');
del_idx = del_idx(1:n_extra);
ARCH.pos(del_idx,:) = [];
ARCH.pos_fit(del_idx,:) = [];
ARCH.pos_f(del_idx,:) = [];
ARCH.pos_g(del_idx,:) = [];
ARCH = updateGrid(ARCH,ngrid);
end

%% This function checks the boundary of jellyfish search space
function [POS] = checkBoundaries(POS,var_max,var_min)
% Useful matrices
Np = size(POS,1);
MAXLIM   = repmat(var_max(:)',Np,1);
MINLIM   = repmat(var_min(:)',Np,1);
POS(POS>MAXLIM) = MAXLIM(POS>MAXLIM);
POS(POS<MINLIM) = MINLIM(POS<MINLIM);
end
function POS=checksimplebounds(POS,Lb,Ub)
for i=1:size(POS,1)
    ns_tmp=POS(i,:);
    I=ns_tmp<Lb;
    while sum(I)~=0
        ns_tmp(I)=Ub(I)+(ns_tmp(I)-Lb(I));
        I=ns_tmp<Lb;
    end
    J=ns_tmp>Ub;
    while sum(J)~=0
        ns_tmp(J)=Lb(J)+(ns_tmp(J)-Ub(J));
        J=ns_tmp>Ub;
    end
    POS(i,:)=ns_tmp;
end
Np = size(POS,1);
Ub   = repmat(Ub(:)',Np,1);
Lb   = repmat(Lb(:)',Np,1);
POS(POS>Ub) = Ub(POS>Ub);
POS(POS<Lb) = Lb(POS<Lb);
end

%----------------------------------------------------------------------------------------%
% The equations can be referred to the below paper:                                      %
% Equation 01.                                                                           %
% Paper: Mahdavi S, Rahnamayan S, Deb K. Opposition based learning: A literature review. % 
% Swarm and Evolutionary Computation 2018;39:1-23.                                       %
% DOI:  https://doi.org/10.1016/j.swevo.2017.09.010                                      %
%----------------------------------------------------------------------------------------% 

function [POS] = OPPOS(POS,var_max,var_min)
    Np = size(POS,1);
    MAXLIM   = repmat(var_max(:)',Np,1);
    MINLIM   = repmat(var_min(:)',Np,1);
    POS = (MINLIM+MAXLIM)-POS;
end

%------------------------------------------------------------------------------------------------------------------------%
% The equations can be referred to the below paper:                                                                      %
% Fister I, Perc M, Kamal SM, Fister I. A review of chaos-based firefly algorithms: Perspectives and research challenges % 
% Applied Mathematics and Computation                                                                                    %
% 2015;252:155-65.DOI: https://doi.org/10.1016/j.amc.2014.12.006                                                         %
%------------------------------------------------------------------------------------------------------------------------%
function pop=initialchaos(index,num_pop,nd,Ub,Lb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the population by chaostic map
% This program is written at PiM lab
% pop=initialchaos(index,num_pop,nd,Ub,Lb)
% index: choose the map case
%       1:Chebyshev map
%       2:Circle map
%       3:Gauss/mouse map
%       4:Intermittency map
%       5:Iterative map
%       6:Liebovitch map
%       7:Logistic map
%       8:Piecewise map
%       9:Sine map
%       10:Singer map
%       11:Sinusoidal map
%       12:Tent map
%       13: Kent map
% num_iter: Number of population;
% nd: Number of dimention; e.g: nd=4;
% Ub: Matrix of Upper bound,e.g:[1 1 1 1];
% Lb: Matrix of lower bound,e.g:[-1 0 -2 3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initinationthearms(n,d,Ub,Lb)
%chaos(chaosIndex-1,iteration,max_it,chValue);
if size(Lb,2)==1
    Lb=Lb*ones(1,nd);
    Ub=Ub*ones(1,nd);
end
x(1,:)=rand(1,nd);%0.7;
switch index
%% Chebyshev map
    case 1
        for i=1:(num_pop-1)
            x(i+1,:)=cos(i*acos(x(i,:)));
        end
%% Circle map
    case 2
        a=0.5;
        b=0.2;
        for i=1:(num_pop-1)
            x(i+1,:)=mod(x(i,:)+b-(a/(2*pi))*sin(2*pi*x(i,:)),1);
        end
%% Gauss/mouse map
    case 3
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)==0
                    x(i+1,k)=0;
                else
                    x(i+1,k)=mod(1/x(i,k),1);
                end
            end
        end
%% Intermittency map
    case 4
        P=0.4;
        m=2;
        eps=1e-10;
        c=(1-eps-P)/(P^m);
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)<=P
                    x(i+1,k)=eps+x(i,k)+c*(x(i,k))^m;
                else
                    x(i+1,k)=(x(i,k)-P)/(1-P);
                end
            end
        end
%% Iterative map
    case 5   
        a=0.7;
        for k=1:nd
            for i=1:(num_pop-1)
                x(i+1,k)=sin((a*pi)/x(i,k));
            end
        end
%% Liebovitch map 
    case 6
        
        P1=0.3;
        P2=0.7;
        anpha=P2/P1*(1-(P2-P1));
        beta=(1/(P2-1))*((P2-1)-P1*(P2-P1));
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)<=P1
                    x(i+1,k)=anpha*x(i,k);
                elseif x(i,k)<=P2
                    x(i+1,k)=(P2-x(i,k))/(P2-P1);
                else
                    x(i+1,k)=1-beta*(1-x(i,k));
                end
            end
        end
%% Logistic map
    case 7 
        a=4;
        for i=1:(num_pop-1)
            x(i+1,:)=a*x(i,:).*(1-x(i,:));
        end
%% Piecewise map
    case 8   
        P=0.4;
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)>=0 && x(i,k)<P
                    x(i+1,k)=x(i,k)/P;
                end
                if x(i,k)>=P && x(i,k)<0.5
                    x(i+1,k)=(x(i,k)-P)/(0.5-P);
                end
                if x(i,k)>=0.5 && x(i,k)<1-P
                    x(i+1,k)=(1-P-x(i,k))/(0.5-P);
                end
                if x(i,k)>=1-P && x(i,k)<1
                    x(i+1,k)=(1-x(i,k))/P;
                end
            end
        end
%% Sine map
    case 9      
        for i=1:(num_pop-1)
            x(i+1,:) = sin(pi*x(i,:));
        end
%% Singer map
    case 10       
        u=1.07;
        for k=1:nd
            for i=1:(num_pop-1)
                x(i+1,k) = u*(7.86*x(i,k)-23.31*(x(i,k)^2)+28.75*(x(i,k)^3)-13.302875*(x(i,k)^4));
            end
        end
%% Sinusoidal map
    case 11        
        for k=1:nd
            for i=1:(num_pop-1)
                while mod(x(i,k),1)==0
                    x(i,k)=rand(1,1);
                end
                x(i+1,k) = 2.3*x(i,k)^2*sin(pi*x(i,k));
            end
        end
%% Tent map
    case 12        
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)<0.7
                    x(i+1,k)=x(i,k)/0.7;
                end
                if x(i,k)>=0.7
                    x(i+1,k)=(10/3)*(1-x(i,k));
                end
            end
        end
%% Kent map
    case 13        
        m=0.6;
        for k=1:nd
            for i=1:(num_pop-1)
                if x(i,k)<=m
                    x(i+1,k)=x(i,k)/m;
                else
                    x(i+1,k)=(1-x(i,k))/(1-m);
                end
            end
        end
end
for k=1:nd
    for i=1:num_pop
        pop(i,k)=Lb(k)+x(i,k)*(Ub(k)-Lb(k));
    end
end
end

%--------------------------------------------------------------------------%
% This function can be referred in the below book:                     %
% Levy exponent and coefficient                                            %
% For details, see Chapter 3 of the following book:                        %
% Xin-She Yang, Nature-Inspired Optimization Algorithms, Elsevier, (2014). %
% Chaper 3                                                                 %
%--------------------------------------------------------------------------%
function s=Levy(d) % d is number of dimension
beta=3/2;
% Eq. (3.27)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
% Eq. (3.26)
u=randn(1,d)*sigma;
v=randn(1,d);
% Eq. (3.25)
step=u./abs(v).^(1/beta);
s=0.01*step;
end