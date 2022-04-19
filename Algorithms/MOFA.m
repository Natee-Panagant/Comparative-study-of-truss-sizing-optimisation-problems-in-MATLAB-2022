%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified by Natee Panagant           %
% Department of Mechanical Engineering %
% Faculty of Engineering               %
% Khon Kaen University                 %
% Thailand                             %
% natepa@kku.ac.th                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This demo shows how the multiobjective firefly algorithm (MOFA) works  %
%% The standard firefly algorithm (FA) was developed by Xin-She Yang and  %
%% its Matlab code is alrealy available online at Mathworks, etc.         %
%% This demo focuses on the MOFA with non-dominated sorting               %
%% ---------------------------------------------------------------------- %

%% Programmed by Xin-She Yang in 2011 and 2012 @ Cambridge Univ and
%% National Physics Laboratory, London. Updated and last modified
%% by X S Yang in 2014, and then 2015.
% References:
% (1) Xin-She Yang, Firefly algorithm, stochastic test functions and
%     design optimizaiton, Int. Journal of Bio-Inspired Computation,
%     vol. 2, no. 2, 78-84 (2010).
% (2) Xin-She Yang, Multiobjective firely algorithm for continuous
%     optimizatoin, Engineering with Computers, vol. 29, no. 2,
%     175--184 (2013).
% (3) Xin-She Yang, Nature-Inspired Optimization Algorithms,
%     Elsevier Insight, (2014).  [Book]
% -----------------------------------------------------------------------

function rst=MOFA(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
rand('state',sum(100*clock));
%
n=nsol;               % Population size (number of fireflies)
tMax=nloop;            % Maximum number of iterations
alpha=1.0;              % Randomness strength 0--1 (highly random)
beta0=1.0;              % Attractiveness constant
gamma=0.1;              % Absorption coefficient
theta=10^(-4/tMax);     % The parameter theta can be taken as 0.97 to 0.99
% This is a randomness reduction factor for alpha

% For the ZDT Function #3 with m=2 objectives
% m=2;             % Number of objectives
RnD=zeros(n,2);  % Initilize the rank and distance matrix
% Dimension of the search/independent variables
d=nvar;
Lb=a';   % Lower bounds/limits
Ub=b';   % Upper bounds/limits

% Number of objectives
xr=Lb+(Ub-Lb).*rand(1,nvar);
[fr,gr]=feval(fun,xr');
m=numel(fr);
NC=numel(gr);

% Generating the initial locations of n fireflies
F=zeros(nsol,m);
G=zeros(nsol,NC);
for i=1:n,
    Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
%     f(i,1:m) = obj_funs(Sol(i,:), m);
    [F(i,:),G(i,:)] = feval(fun,Sol(i,:)');
    f(i,1:m)=fpenal00(F(i,:)',G(i,:)')';
end
% Store the fitness or objective values
f_new=f;
%% Sort the initialized population
x=[Sol f];  % combined into a single input
% Non-dominated sorting for the initila population
[Sorted,PFind,frontRanks_Index]=solutions_sorting(x, m,d);
% Decompose into solutions, fitness, rank and distances
Sol=Sorted(:,1:d);  S_new=Sol;       % Record solutions
f=Sorted(:,(d+1):(d+m));  f_new=f;   % Record objectives
RnD=Sorted(:,(d+m+1):end);           % Record ranks
%
x=Sol(frontRanks_Index,:);
F=F(frontRanks_Index,:);
G=G(frontRanks_Index,:);
x1=zeros(nsol,nvar);
F1=zeros(nsol,m);
G1=zeros(nsol,NC);
%
t=1;
Neval=0;
rst.Neval_History=[];
while Neval<n*tMax,          %%%%% start the firely algorithm iterations %%%%%
    alpha=alpha*theta;  % Reduce alpha by a factor 0<theta<1
    scale=abs(Ub-Lb);   % Scale of the optimization problem
    Sol_old=Sol;        % Save the old population
    f_old=f;            % Save the old population objectives
    % Two loops over all the n fireflies
    for i=1:n,
        for j=i:n,
            % Update moves and move to the brighter/more attractive
            % That is, all m objectives [i.e., f(,1:m)] should improve.
            % For example, for m=2, this means that the logical
            % condition (f(j,1)<=f(i,1) & f(j,2) <=f(i,2)) is true.
            if (f_old(j,1:m)<=f_old(i,1:m)),
                r=sqrt(sum((Sol(i,:)-Sol(j,:)).^2));
                beta=beta0*exp(-gamma*r.^2);     % Attractiveness
                steps=alpha.*(rand(1,d)-0.5).*scale;
                % The FA equation for updating position vectors
                % That is, to move firefly i torwards firefly j
                Sol(i,:)=Sol(i,:)+beta*(Sol(j,:)-Sol(i,:))+steps;
                Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
            end
%             f(i,1:m)=obj_funs(Sol(i,1:d),m);
        end % end for j
    end % end for i

    %% Evalute the fitness/function values of the new population
    Sind0=[];
    Sind1=[];
    for j=1:n,
%         f(j, 1:m)=obj_funs(Sol(j,1:d),m);
        x1(j,:)=Sol(j,1:d);
        [F1(j,:),G1(j,:)] = feval(fun,x1(j,:)');
        %
        Neval=Neval+1;
        %
        f_new(j,1:m)=fpenal00(F1(j,:)',G1(j,:)')';
        %
        if (f_old(j,1:m) <= f(j,1:m)),   % if all improve
            Sol(j,1:d) = Sol_old(j,1:d);%Fixed by Natee
            f(j,1:m)=f_old(j,1:m);
            Sind1=[Sind1,i];%select index
        end
        % Update the current best (stored in the first row)
        if (f(j,1:m) <= f_old(1,1:m)),
            Sol(1,1:d) = Sol(j,1:d);
            f(1,:)=f(j,:);
            Sind0=i;%best index
        end
    end % end of for loop j

    %% ! It's very important to combine both populations, otherwise,
    %% the results may look odd and will be very inefficient. !
    %% The combined population consits of both the old and new solutions
    %% So the total size of the combined population for sorting is 2*n
    X(1:n,:)=[Sol f_new];               % Combine new solutions
    X((n+1):(2*n),:)=[Sol_old f_old];   % Combine old solutions
    [Sorted,PFind,frontRanks_Index]=solutions_sorting(X, m, d);

    %Update Pareto
    x(Sind1,:)=x1(Sind1,:);
    F(Sind1,:)=F1(Sind1,:);
    G(Sind1,:)=G1(Sind1,:);
    if numel(Sind0)>0
        x(1,:)=x1(Sind0,:);
        F(1,:)=F1(Sind0,:);
        G(1,:)=G1(Sind0,:);
    end
    x2=[x1;x];
    F2=[F1;F];
    G2=[G1;G];
    x2=x2(frontRanks_Index,:);
    F2=F2(frontRanks_Index,:);
    G2=G2(frontRanks_Index,:);
    if numel(PFind)>nsol
       	CD=Sorted(1:numel(PFind),nvar+2);
        [~,isort]=sort(CD,'descend');
        ppareto=x2(isort(1:nsol),:)';
        fpareto=F2(isort(1:nsol),:)';
        gpareto=G2(isort(1:nsol),:)';
        fppareto=Sorted(isort(1:nsol),nvar+1:nvar+m)';
    else
       	ppareto=x2(1:numel(PFind),:)';
        fpareto=F2(1:numel(PFind),:)';
        gpareto=G2(1:numel(PFind),:)';
        fppareto=Sorted(1:numel(PFind),nvar+1:nvar+m)';
    end
    %

    %% Select n solutions among a combined population of 2*n solutions
    [new_Sol,Sind2]=Select_pop(Sorted, m, d, n);
    % Decompose into solutions, fitness and ranking
    Sol=new_Sol(:,1:d);             % Sorted solutions
    f=new_Sol(:,(d+1):(d+m));       % Sorted objective values
    RnD=new_Sol(:,(d+m+1):end);     % Sorted ranks and distances

    %% Running display at each 100 iterations
%     if ~mod(t,100),
%         disp(strcat('Iterations t=',num2str(t)));
%         plot(f(:, 1), f(:, 2),'ro','MarkerSize',3);
%         axis([0 1 -0.8 1]);
%         xlabel('f_1'); ylabel('f_2');
%         drawnow;
%     end

    %Update Population
    x=x2(Sind2,:);
    F=F2(Sind2,:);
    G=G2(Sind2,:);
    
    %Save Results
    fea_idx=max(gpareto,[],1)<=0;
    rst.ppareto{t}=ppareto(:,fea_idx);
    rst.fppareto{t}=fppareto(:,fea_idx);
    rst.fpareto{t}=fpareto(:,fea_idx);
    rst.gpareto{t}=gpareto(:,fea_idx);

    [~,uidx]=unique(rst.ppareto{t}','rows');
    rst.ppareto{t}=rst.ppareto{t}(:,uidx);
    rst.fppareto{t}=rst.fppareto{t}(:,uidx);
    rst.fpareto{t}=rst.fpareto{t}(:,uidx);
    rst.gpareto{t}=rst.gpareto{t}(:,uidx);
    
    rst.Neval_History=[rst.Neval_History,Neval];
    rst.timestamp=datetime('now');
%     pareto_plot(rst.fppareto{t},rst.fpareto{t},rst.gpareto{t},fun);
    %

    t=t+1;
end % End of t loop (up to tMax) and end of the main FA loop

%% Make sure that new fireflies are within the bounds/limits
function s=simplebounds(s,Lb,Ub)
% Apply the lower bound
ns_tmp=s;
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);
% Apply the upper bounds
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);
% Update this new move
s=ns_tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective functions
function f = obj_funs(x, m)
% Zitzler-Deb-Thiele's funciton No 3 (ZDT function 3)
% m = # of objectives   % d = # of variables/dimensions
d=length(x);  % d=30 for ZDT 3
% First objective f1
f(1) = x(1);
g=1+9/29*sum(x(2:d));
h=1-sqrt(f(1)/g)-f(1)/g*sin(10*pi*f(1));
% Second objective f2
f(2) = g*h;

%%%%%%%%%%%%%%%%%% end of the definitions of obojectives %%%%%%%%%%%%%%%%%%

function [new_Sol,Sind2] = Select_pop(firefly, m, ndim, npop)
% The input population to this part has twice (ntwice) of the needed
% population size (npop). Thus, selection is done based on ranking and
% crowding distances, calculated from the non-dominated sorting
ntwice= size(firefly,1);
% Ranking is stored in column Krank
Krank=m+ndim+1;
% Sort the population of size 2*npop according to their ranks
[~,Index] = sort(firefly(:,Krank));
sorted_firefly=firefly(Index,:);
% Get the maximum rank among the population
RankMax=max(firefly(:,Krank));

%% Main loop for selecting solutions based on ranks and crowding distances
K = 0;  % Initialization for the rank counter
% Loop over all ranks in the population
Sind2=[];
for i =1:RankMax,
    % Obtain the current rank i from sorted solutions
    RankSol = max(find(sorted_firefly(:, Krank) == i));
    % In the new solutions, there can be npop solutions to fill
    if RankSol<npop,
        new_Sol(K+1:RankSol,:)=sorted_firefly(K+1:RankSol,:);
        Sind2=[Sind2,K+1:RankSol];
    end
    % If the population after addition is large than npop, re-arrangement
    % or selection is carried out
    if RankSol>=npop
        % Sort/Select the solutions with the current rank
        candidate_firefly=sorted_firefly(K + 1 : RankSol, :);
        [~,tmp_Rank]=sort(candidate_firefly(:,Krank+1),'descend');
        % Fill the rest (npop-K) fireflies/solutions up to npop solutions
        for j = 1:(npop-K),
            new_Sol(K+j,:)=candidate_firefly(tmp_Rank(j),:);
            Sind2=[Sind2,K+tmp_Rank(j)];
        end
    end
    % Record and update the current rank after adding new solutions
    K = RankSol;
end

function [sorted_x,PFind,frontRanks_Index] = solutions_sorting(x, m, ndim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs and outputs are the extended solutions x with a dimension of    %
%% npop by (ndim+m+2). The objective values are already included in x.    %
% More specifically, the first ndim columns are the actual solutions or
% variable values (1:ndim), followed by the m columns of objective values.
% Then, the next column (i.e.,ndim+m+1) corresponds to the ranks, whereas
% the final column (i.e., ndim+m+2) records the crowd distances.
% ----------------------------------------------------------------------- %

%% Get the parameters from the inputs such as the size of input population
npop=size(x,1);    % Population size
frontRank=1;       % Pareto frontRank (counter) initialization
Rcol=ndim+m+1;     % Store the ranks in the column Rcol=ndim+m+1
% Define the Parato Front as a class (PF) and initilization of xSol
PF(frontRank).R=[];   xSol=[];
%% The main non-dominated sorting starts here             %%%%%%%%%%%%%%%%%
for i = 1:npop,
    % Set the number (initially, 0) of solutions dominating this solution
    xSol(i).n=0;
    % Find all the solutions (that dominated by this solution)
    xSol(i).q=[];
    % Sorting into 3 categories: better (minimization), equal & otherwise
    for j=1:npop,
        % Definte 3 counters for 3 categories
        ns_categ_1=0; ns_categ_2=0; ns_categ_3=0;
        for k=1:m,  % for all m objectives
            % Update the counters for 3 different categories
            if (x(i,ndim+k) < x(j,ndim+k)),      % better/non-dominated
                ns_categ_1=ns_categ_1+1;
            elseif (x(i,ndim+k)==x(j,ndim+k)),   % equal
                ns_categ_2=ns_categ_2+1;
            else                                 % dominated
                ns_categ_3=ns_categ_3+1;
            end
        end % end of k
        % Update the solutions in their class
        if ns_categ_1==0 && ns_categ_2 ~= m
            xSol(i).n=xSol(i).n+1;
        elseif ns_categ_3 == 0 && ns_categ_2 ~= m
            xSol(i).q=[xSol(i).q j];
        end
    end % end of j
    %% Record/Udpate the Pareto Front
    if xSol(i).n==0,
        x(i,Rcol)=1;   % Update the Rank #1 (i.e., the Pareto Front)
        PF(frontRank).R = [PF(frontRank).R i];
    end
end % end of i=1:npop (The first round full rank-sorting process)

% Update the rest frontRanks (close, but not on the Pareto Front)
while ~isempty(PF(frontRank).R),
    nonPF=[];    % Intialization the set
    N=length(PF(frontRank).R);
    for i=1 :N,
        % Get the solution/list
        Sol_tmp_q=xSol(PF(frontRank).R(i)).q;
        % If not empty, update
        if ~isempty(xSol(Sol_tmp_q))
            for j = 1:length(Sol_tmp_q),
                % Get the solutions dominated by the current solution
                Sol_tmp_qj=xSol(PF(frontRank).R(i)).q(j);
                xSol(Sol_tmp_qj).n=xSol(Sol_tmp_qj).n-1;
                if xSol(Sol_tmp_qj).n==0
                    x(Sol_tmp_qj, Rcol)=frontRank + 1;
                    nonPF = [nonPF Sol_tmp_qj];
                end
            end % end of j
        end
    end  % end of i
    frontRank=frontRank+1;
    PF(frontRank).R=nonPF;
end % end of PF(frontRank)

% Now carry out the sorting of ranks and then update
[~,frontRanks_Index]=sort(x(:, Rcol));
Sorted_frontRank=x(frontRanks_Index,:);
PFind=frontRanks_Index(1:numel(PF(1).R));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate the crowding distances for each solution for each frontRank   %
% That is, all the non-domonated solutions on the Pareto Front.  %%%%%%%%%%
Qi=0;      % Initialize a counter
for frontRank=1:(length(PF)-1),
    % Define/initialize a generalized distance matrix
    dc = [];    past_Q=Qi+1;
    for i=1:length(PF(frontRank).R),
        dc(i,:)=Sorted_frontRank(Qi+i,:);
    end
    Qi=Qi+i;
    % Solutions are sorted according to their fitness/objective values
    fobj_sorted=[];
    for i=1:m,
        [~, f_Rank]=sort(dc(:,ndim+i));
        fobj_sorted=dc(f_Rank,:);
        % Find the max and min of the fobj values
        fobj_max=fobj_sorted(length(f_Rank), ndim+i);
        fobj_min=fobj_sorted(1, ndim+i);
        % Calculate the range of the fobj
        f_range=fobj_max-fobj_min;
        % If the solution is at the end/edge, set its distance as infinity
        dc(f_Rank(length(f_Rank)), Rcol+i)=Inf;
        dc(f_Rank(1), Rcol+i) = Inf;
        for j=2:length(f_Rank)-1,
            fobj2=fobj_sorted(j+1,ndim + i);
            fobj1=fobj_sorted(j-1,ndim + i);
            % Check the range or special cases
            if (f_range==0),
                dc(f_Rank(j), Rcol+i)=Inf;
            else
                % Scale the range for distance normalization
                dc(f_Rank(j),Rcol+i)=(fobj2-fobj1)/f_range;
            end
        end % end of j
    end % end of i

    % Calculate and update the crowding distances on the Pareto Front
    dist = []; dist(:,1)=zeros(length(PF(frontRank).R),1);
    for i=1:m,
        dist(:,1)=dist(:,1)+dc(:, Rcol+i);
    end
    % Store the crowding distrance (dc) in the column of Rcol+1=ndim+m+2
    dc(:, Rcol+1)=dist;  dc=dc(:,1:Rcol+1);
    % Update for the output
    xy(past_Q:Qi,:)=dc;
end  % end of all ranks search/update
sorted_x=xy();    % Output the sorted solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of non-dominated sorting %%%%%%%%%%%%%%
