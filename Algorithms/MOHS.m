%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified by Natee Panagant           %
% Department of Mechanical Engineering %
% Faculty of Engineering               %
% Khon Kaen University                 %
% Thailand                             %
% natepa@kku.ac.th                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rst=MOHS(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
rand('state',sum(100*clock));
%
nv=nvar;
lb=a';
ub=b';
% Created by Gifari Zulkarnaen, 2020

% This script was created thanks to the author who made this for another
% multi-objective algorithm, here I changed it into the Harmony Search
% algorithm. This script may contain bugs, as I used it for a specific
% project, please use it as a reference.

%clear all;clc;
%{
% Example Multi Objective Function
% (try it by commenting the function lines: Line 1 & Line 107)
fun = @(x) [-10.*(exp(-0.2.*sqrt(x(:,1).^2+x(:,2).^2)) + ,...
                     exp(-0.2.*sqrt(x(:,2).^2+x(:,3).^2))), ...
                     sum(abs(x).^0.8 + 5.*sin(x.^3),2)];
nv = 3;
lb = -5.*ones(1,MultiObj.nv);
ub = 5.*ones(1,MultiObj.nv);

%}
% Parameters
HMS  = nsol;                 % Harmony Memory Size
HMCR = 0.5;                 % Harmony Memory Consideration Rate
PAR  = 0.2;                 % Pitch Adjusting Rate
BW   = [1 0.01];            % Bandwith
ni   = nloop;                % Max total iteration
nr   = narchive;                  % Repository (pareto) size
ns   = 50;


%% Initialization
HMXtemp = repmat((ub-lb),HMS,1).*rand(HMS,nv) + repmat(lb,HMS,1); % HM input
% HMFtemp = fun(HMXtemp); % Harmony memory output

%
% Number of objectives
xr=lb+(ub-lb).*rand(1,nv);
[fr,gr]=feval(fun,xr');
NO=numel(fr);
NC=numel(gr);
%
HMFtemp=zeros(nsol,NO);
F=zeros(nsol,NO);
G=zeros(nsol,NC);
for i=1:HMS
    [F(i,:),G(i,:)] = feval(fun,HMXtemp(i,:)');
    HMFtemp(i,1:NO)=fpenal00(F(i,:)',G(i,:)')';
end
%

dom = checkDomination(HMFtemp);
HS.HMX = HMXtemp(~dom,:);
HS.HMF = HMFtemp(~dom,:);
HMXtemp2 = HS.HMX;
HMFtemp2 = HS.HMF;

%
Fpareto=F(~dom,:);
Gpareto=G(~dom,:);
%

% vol = (1:ns)*100;
% lastVol = 1;
iter = 1;
Neval=0;
rst.Neval_History=[];
%% Main MOHS Loop
% stop = false;
while Neval<nsol*ni

    % Generating new harmony memory
    HMXtemp = zeros(HMS,nv);
    for m = 1:HMS
        for v = 1:nv
            if rand <= HMCR
                % New harmony with HMC
                HMXtemp(m,v) = HMXtemp2(randi(size(HMXtemp2,1)),v);

                % Pitch Adjusment
                if rand <= PAR
                    if v <= nv/2
                        p = [-BW(1) BW(1)];
                    else
                        p = [-BW(2) BW(2)];
                    end
                    %         p = [-BW(2) BW(2)];
                    HMXtemp(m,v) = HMXtemp(m,v) + p(randi(2));
                end
            else
                % New harmony without HMC nor PA
                HMXtemp(m,v) = (ub(v)-lb(v)).*rand(1) + lb(v);
            end
        end
    end

    % Check boundaries
    HMXtemp = checkBoundaries(HMXtemp,lb,ub);

    % New harmony fitness
%     HMFtemp = fun(HMXtemp);
    %
    for i=1:HMS
        [F(i,:),G(i,:)] = feval(fun,HMXtemp(i,:)');
        %
        Neval=Neval+1;
        %
        HMFtemp(i,1:NO)=fpenal00(F(i,:)',G(i,:)')';
    end
    %

    % Update the repository, compared with previous harmony memory
%     HS = updateRepository(HS,HMXtemp,HMFtemp);
    %
    HS.HMX=[HS.HMX;HMXtemp];
    HS.HMF=[HS.HMF;HMFtemp];
    Fpareto=[Fpareto;F];
    Gpareto=[Gpareto;G];
    dom  = checkDomination(HS.HMF);
    HS.HMX=HS.HMX(~dom,:);
    HS.HMF=HS.HMF(~dom,:);
    Fpareto=Fpareto(~dom,:);
    Gpareto=Gpareto(~dom,:);
    %

    if(size(HS.HMX,1)>HMS)
        [HS,del_idx] = deleteFromRepository(HS,size(HS.HMX,1)-nr);
        %
        Fpareto(del_idx,:)=[];
        Gpareto(del_idx,:)=[];
        %
    end

    % Memorize best solutions
    HMXtemp2 = [HS.HMX; HMXtemp2];
    HMFtemp2 = [HS.HMF; HMFtemp2];
    [~,idxTemp2,~] = unique(HMFtemp2(:,1),'stable');
    HMXtemp2 = HMXtemp2(idxTemp2,:);
    HMFtemp2 = HMFtemp2(idxTemp2,:);
    if size(HMXtemp2,1) > HMS
        HMXtemp2 = HMXtemp2(1:HMS,:);
        HMFtemp2 = HMFtemp2(1:HMS,:);
    end

    % Stop?
    % newRelVol = abs(1-measureArea(HS.HMF)/lastVol);
    % lastVol = measureArea(HS.HMF);
    % vol = [vol(2:end) newRelVol];
    % if(iter>ni || mean(vol) < 1e-4), stop = true; end
    % iter = iter + 1;
    % end
    % hold off
    % HS.vol = vol;

    %Save Results
    ppareto=HS.HMX';
    fppareto=HS.HMF';
    fpareto=Fpareto';
    gpareto=Gpareto';
    
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
end
%% ==================  Supplementary Functions ===========================

% Function that measure the hypervolume of front
function vol = measureVolume(frontOut)
refPoint = max(frontOut,[],1)+1;
Np = size(frontOut,1);
vol = 0;
for i=1:Np
    vol = vol + prod(abs(frontOut(i,:)-refPoint));
end
end

% Function that measure the area of front for two objectives problem
function vol = measureArea(frontOut)
refPoint = max(frontOut,[],1);
polyg = [refPoint; frontOut; refPoint];
vol = polyarea(polyg(:,1),polyg(:,2));
end

% Function that measure distance of front
% function dist = distConv(x)
% end

% Function that updates the repository given a new population and its
% fitness
% function HS = updateRepository(HS,HMX,HMF)
% % Domination between individuals
% DOMINATED  = checkDomination(HMF);
% HS.HMX     = [HS.HMX; HMX(~DOMINATED,:)];
% HS.HMF	   = [HS.HMF; HMF(~DOMINATED,:)];
% % Domination between nondominated individuals and the last repository
% DOMINATED  = checkDomination(HS.HMF);
% HS.HMF     = HS.HMF(~DOMINATED,:);
% HS.HMX     = HS.HMX(~DOMINATED,:);
% end

% Function for checking the domination between the population. It
% returns a vector that indicates if each individual is dominated (1) or not
function dom_vector = checkDomination(fitness)
Np = size(fitness,1);
dom_vector = zeros(Np,1);
all_perm = nchoosek(1:Np,2);    % Possible permutations
all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];

d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
dominated_individuals = unique(all_perm(d==1,2));
dom_vector(dominated_individuals) = 1;
end

% Function that returns 1 if x dominates y and 0 otherwise
function d = dominates(x,y)
d = all(x<=y,2) & any(x<y,2);
end

% Function that deletes an excess of solutions inside the repository using
% crowding distances
function [HS,del_idx] = deleteFromRepository(HS,n_extra)
% Compute the crowding distances
crowding = zeros(size(HS.HMX,1),1);
for m = 1:1:size(HS.HMF,2)
    [m_fit,idx] = sort(HS.HMF(:,m),'ascend');
    m_up     = [m_fit(2:end); Inf];
    m_down   = [Inf; m_fit(1:end-1)];
    distance = (m_up-m_down)./(max(m_fit)-min(m_fit));
    [~,idx]  = sort(idx,'ascend');
    crowding = crowding + distance(idx);
end
crowding(isnan(crowding)) = Inf;

% Delete the extra solutions with the smallest crowding distances
[~,del_idx] = sort(crowding,'ascend');
del_idx = del_idx(1:n_extra);
HS.HMX(del_idx,:) = [];
HS.HMF(del_idx,:) = [];
end

% Function that corrects the individuals that exceed the boundaries
function HMX = checkBoundaries(HMX,lb,ub)
% Useful matrices
Np = size(HMX,1);
maxLim   = repmat(ub(:)',Np,1);
minLim   = repmat(lb(:)',Np,1);

% Correct individuals
HMX(HMX>maxLim) = maxLim(HMX>maxLim);
HMX(HMX<minLim) = minLim(HMX<minLim);
end
