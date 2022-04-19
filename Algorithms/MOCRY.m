%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified by Natee Panagant           %
% Department of Mechanical Engineering %
% Faculty of Engineering               %
% Khon Kaen University                 %
% Thailand                             %
% natepa@kku.ac.th                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%__________________________________________________________________ %
%                          Multi-Objective                          %
%        Crystal Structure Algorithm (CryStAl) (MOCryStAl)          %
%                                                                   %
%                                                                   %
%                  Developed in MATLAB R2021a (MacOs)               %
%                                                                   %
%                      Author and programmer                        %
%                ---------------------------------                  %
%                      Nima Khodadadi (ʘ‿ʘ)                         %
%                       Siamak Talatahari                           %
%                         Mahdi Azizi                               %
%                         Pooya Sareh                               %
%                                                                   %
%                             e-Mail                                %
%                ---------------------------------                  %
%                         inimakhan@me.com                          % 
%                                                                   %
%                            Homepage                               %
%                ---------------------------------                  %
%                    https://nimakhodadadi.com                      %
%                                                                   %
%                                                                   %
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------------------------------------------------------- %

function rst = MOCRY(fun,fout,nloop,nsol,nvar,nbit,narchive,a,b)
rand('state',sum(100*clock));
%
MaxIteation=nloop;
Archive_size=narchive;
Cr_Number=nsol;
Var_Number=nvar;
LB=a';
UB=b';

%% Updating the Size of ProblemParameters
alpha=0.1;  % Grid Inflation Parameter
nGrid=30;   % Number of Grids per each Dimension
beta=4; %=4;    % Leader Selection Pressure Parameter
gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure


if length(LB)==1
    LB=repmat(LB,1,Var_Number);
end
if length(UB)==1
    UB=repmat(UB,1,Var_Number);
end
%% Initialization


Crystal=CreateEmptyParticle(Cr_Number);

% Initializing the Position of first probs
for i=1:Cr_Number
        
    Crystal(i).Velocity=0;
    Crystal(i).Position=zeros(1,Var_Number);
    for j=1:Var_Number
%         Crystal(i).Position(1,j)=unifrnd(LB(j),UB(j),1);
        Crystal(i).Position(1,j)=LB(j)+(UB(j)-LB(j))*rand;
    end
%     Crystal(i).Cost=ObjFuncName(Crystal(i).Position')';
    %Function Evaluation
    [f1,g1]=feval(fun,Crystal(i).Position');
    fp1 = fpenal00(f1,g1);
    Crystal(i).Cost=fp1';
    Crystal(i).f=f1';
    Crystal(i).g=g1';
    %
    Fun_eval(i,:)=norm(Crystal(i).Cost);
    Crystal(i).Best.Position=Crystal(i).Position;
    Crystal(i).Best.Cost=Crystal(i).Cost;
    Crystal(i).Best.f=Crystal(i).f;
    Crystal(i).Best.g=Crystal(i).g;
end

[~,idbest]=min(Fun_eval);
Crb=Crystal(idbest,:);
Crystal=DetermineDominations(Crystal);
Archive=GetNonDominatedParticles(Crystal);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end
% The best Crystal

%% Search Process
Iter=1;
Neval=0;
rst.Neval_History=[];
while Neval<Cr_Number*MaxIteation
    for i=1:Cr_Number
        Leader=SelectLeader(Archive,beta);
        %% Generate New Crystals
        % Main Crystal
        Crmain=Crystal(randperm(Cr_Number,1),:);
        % Random-selected Crystals
        RandNumber=randperm(Cr_Number,1);
        RandSelectCrystal=randperm(Cr_Number,RandNumber);
        CrystalP=vertcat( Crystal.Position );
        % Mean of randomly-selected Crystals
        Fc=mean(CrystalP(RandSelectCrystal,:)).*(length(RandSelectCrystal)~=1)...
            +CrystalP(RandSelectCrystal(1,1),:).*(length(RandSelectCrystal)==1);
        % Random numbers (-1,1)
        r=2*rand-1;       r1=2*rand-1;
        r2=2*rand-1;     r3=2*rand-1;
        % New Crystals
        Crystal(1+4*(i-1)+Cr_Number,:).Position=Leader.Position+r*Crmain.Position;
        Crystal(2+4*(i-1)+Cr_Number,:).Position=Leader.Position+r1*Crmain.Position+r2*Crb.Position;
        Crystal(3+4*(i-1)+Cr_Number,:).Position=Leader.Position+r1*Crmain.Position+r2*Fc;
        Crystal(4+4*(i-1)+Cr_Number,:).Position=Leader.Position+r1*Crmain.Position+r2*Crb.Position+r3*Fc;
        
        for i2=1:4
            % Checking/Updating the boundary limits for Crystals
            Crystal(i2+4*(i-1)+Cr_Number,:).Position=min(max(Crystal(i2+4*(i-1)+Cr_Number).Position,LB),UB);
            
            % Evaluating New Crystals
%             Crystal(i2+Cr_Number,:).Cost=ObjFuncName(Crystal(i2+Cr_Number).Position')';
            %Function Evaluation
            [f0,g0]=feval(fun,Crystal(i2+4*(i-1)+Cr_Number).Position');
            %
            Neval=Neval+1;
            %
            fp0 = fpenal00(f0,g0);
            Crystal(i2+4*(i-1)+Cr_Number,:).Cost=fp0';
            Crystal(i2+4*(i-1)+Cr_Number,:).f=f0';
            Crystal(i2+4*(i-1)+Cr_Number,:).g=g0';
            %
%             Crystal(i2+4*(i-1)+Cr_Number).Best.Position=Crystal(i2+4*(i-1)+Cr_Number).Position;
%             Crystal(i2+4*(i-1)+Cr_Number).Best.Cost=Crystal(i2+4*(i-1)+Cr_Number).Cost;
            % Updating the Crystals
            
        end
    end % End of One Iteration

    [Crystal,level]=DetermineDominations(Crystal);
    non_dominated_Crystal=GetNonDominatedParticles(Crystal);
    
    % Crytal Selection (Added by Natee Panagant: *To maintain population size equal to Cr_Number since there are 4xCr_Number solutions are reproduced each iteration)
    Crystal=selection(Crystal,level,Cr_Number);
    %
    
    Archive=[Archive
        non_dominated_Crystal];
    
    Archive=DetermineDominations(Archive);
    Archive=GetNonDominatedParticles(Archive);
    
    for i=1:numel(Archive)
        [Archive(i).GridIndex Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end
    
    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);
        
        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,alpha);
        
    end
%     disp(['In iteration ' num2str(Iter) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
%     save results

    Archive_costs=GetCosts(Archive);
    
    % The best Crystal
    
    % Save Results
    Nnds=numel(Archive);
    ppareto=[];
    fpareto=[];
    gpareto=[];
    fppareto=[];
    for i=1:numel(Archive)
       ppareto(:,i)=Archive(i).Position';
       fpareto(:,i)=Archive(i).f';
       gpareto(:,i)=Archive(i).g';
       fppareto(:,i)=Archive(i).Cost';
    end
    fea_idx=max(gpareto,[],1)<=0;
    rst.ppareto{Iter}=ppareto(:,fea_idx);
    rst.fpareto{Iter}=fpareto(:,fea_idx);
    rst.gpareto{Iter}=gpareto(:,fea_idx);
    rst.fppareto{Iter}=fppareto(:,fea_idx);

    [~,uidx]=unique(rst.ppareto{Iter}','rows');
    rst.ppareto{Iter}=rst.ppareto{Iter}(:,uidx);
    rst.fppareto{Iter}=rst.fppareto{Iter}(:,uidx);
    rst.fpareto{Iter}=rst.fpareto{Iter}(:,uidx);
    rst.gpareto{Iter}=rst.gpareto{Iter}(:,uidx);

    rst.Neval_History=[rst.Neval_History,Neval];
    rst.timestamp=datetime('now');
%     pareto_plot(rst.fppareto{Iter},rst.fpareto{Iter},rst.gpareto{Iter},fun);
    
    Iter=Iter+1;
end % End of Main Looping
end

%%%%%%%%%%%%%%%
% Sub-Funcion %
%%%%%%%%%%%%%%%
%% Boundary Handling
function x=bound(x,UB,LB)
x(x>UB)=UB(x>UB); x(x<LB)=LB(x<LB);
end
function particle=CreateEmptyParticle(n)
    
    if nargin<1
        n=1;
    end

    empty_particle.Position=[];
    empty_particle.Velocity=[];
    empty_particle.Cost=[];
    empty_particle.Dominated=false;
    empty_particle.Best.Position=[];
    empty_particle.Best.Cost=[];
    empty_particle.GridIndex=[];
    empty_particle.GridSubIndex=[];
    
    particle=repmat(empty_particle,n,1);
    
end

function G=CreateHypercubes(costs,ngrid,alpha)

    nobj=size(costs,1);
    
    empty_grid.Lower=[];
    empty_grid.Upper=[];
    G=repmat(empty_grid,nobj,1);
    
    for j=1:nobj
        
        min_cj=min(costs(j,:));
        max_cj=max(costs(j,:));
        
        dcj=alpha*(max_cj-min_cj);
        
        min_cj=min_cj-dcj;
        max_cj=max_cj+dcj;
        
        gx=linspace(min_cj,max_cj,ngrid-1);
        
        G(j).Lower=[-inf gx];
        G(j).Upper=[gx inf];
        
    end

end

function rep=DeleteFromRep(rep,EXTRA,gamma)

    if nargin<3
        gamma=1;
    end

    for k=1:EXTRA
        [occ_cell_index occ_cell_member_count]=GetOccupiedCells(rep);

        p=occ_cell_member_count.^gamma;
        p=p/sum(p);

        selected_cell_index=occ_cell_index(RouletteWheelSelection(p));

        GridIndices=[rep.GridIndex];

        selected_cell_members=find(GridIndices==selected_cell_index);

        n=numel(selected_cell_members);

        selected_memebr_index=randi([1 n]);

        j=selected_cell_members(selected_memebr_index);
        
        rep=[rep(1:j-1); rep(j+1:end)];
    end
    
end

% function pop=DetermineDominations(pop)
% 
%     npop=numel(pop);
%     
%     for i=1:npop
%         pop(i).Dominated=false;
%         for j=1:i-1
%             if ~pop(j).Dominated
%                 if Dominates(pop(i),pop(j))
%                     pop(j).Dominated=true;
%                 elseif Dominates(pop(j),pop(i))
%                     pop(i).Dominated=true;
%                     break;
%                 end
%             end
%         end
%     end
%     
% end
function [pop,level]=DetermineDominations(pop)%Modified by Natee Panagant
    A=zeros(numel(pop));
    npop=numel(pop);
    for i=1:npop
        pop(i).Dominated=false;
        for j=1:npop
            if i~=j
                if Dominates(pop(j),pop(i))
                    pop(i).Dominated=true;
                    A(i,j)=A(i,j)+1;
                end
            end
        end
    end
    level=sum(A,2);
end
function pop=selection(pop,level,Cr_Number) %Added by Natee Panagant
    sind=[];
    for i=0:max(level)
        sind_i=find(level==i);
        if numel(sind)+numel(sind_i)<=Cr_Number
            sind=[sind;sind_i];
        else
            Ndiff=numel(sind)+numel(sind_i)-Cr_Number;
            dind=randperm(numel(sind_i),Ndiff);
            sind_i(dind)=[];
            sind=[sind; sind_i];
            break;
        end
    end
	pop=pop(sind,1);
end
function dom=Dominates(x,y)

    if isstruct(x)
        x=x.Cost;
    end

    if isstruct(y)
        y=y.Cost;
    end
    
    dom=all(x<=y) && any(x<y);

end

function costs=GetCosts(pop)

    nobj=numel(pop(1).Cost);
    costs=reshape([pop.Cost],nobj,[]);

end

function [Index SubIndex]=GetGridIndex(particle,G)

    c=particle.Cost;
    
    nobj=numel(c);
    ngrid=numel(G(1).Upper);
    
    str=['sub2ind(' mat2str(ones(1,nobj)*ngrid)];

    SubIndex=zeros(1,nobj);
    for j=1:nobj
        
        U=G(j).Upper;
        
        i=find(c(j)<U,1,'first');
        
        SubIndex(j)=i;
        
        str=[str ',' num2str(i)];
    end
    
    str=[str ');'];
    
    Index=eval(str);
    
end

function nd_pop=GetNonDominatedParticles(pop)

    ND=~[pop.Dominated];
    
    nd_pop=pop(ND);

end

function [occ_cell_index occ_cell_member_count]=GetOccupiedCells(pop)

    GridIndices=[pop.GridIndex];
    
    occ_cell_index=unique(GridIndices);
    
    occ_cell_member_count=zeros(size(occ_cell_index));

    m=numel(occ_cell_index);
    for k=1:m
        occ_cell_member_count(k)=sum(GridIndices==occ_cell_index(k));
    end
    
end

function [Archive_X_Chopped, Archive_F_Chopped, Archive_mem_ranks_updated, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize)

for i=1:size(Archive_F,1)-ArchiveMaxSize
    index=RouletteWheelSelection(Archive_mem_ranks);
    
    Archive_X=[Archive_X(1:index-1,:) ; Archive_X(index+1:Archive_member_no,:)];
    Archive_F=[Archive_F(1:index-1,:) ; Archive_F(index+1:Archive_member_no,:)];
    Archive_mem_ranks=[Archive_mem_ranks(1:index-1) Archive_mem_ranks(index+1:Archive_member_no)];
    Archive_member_no=Archive_member_no-1;
end

Archive_X_Chopped=Archive_X;
Archive_F_Chopped=Archive_F;
Archive_mem_ranks_updated=Archive_mem_ranks;
end

function i=RouletteWheelSelection(p)

    r=rand;
    c=cumsum(p);
    i=find(r<=c,1,'first');

end

function rep_h=SelectLeader(rep,beta)
    if nargin<2
        beta=1;
    end

    [occ_cell_index occ_cell_member_count]=GetOccupiedCells(rep);
    
    p=occ_cell_member_count.^(-beta);
    p=p/sum(p);
    
    selected_cell_index=occ_cell_index(RouletteWheelSelection(p));
    
    GridIndices=[rep.GridIndex];
    
    selected_cell_members=find(GridIndices==selected_cell_index);
    
    n=numel(selected_cell_members);
    
    selected_memebr_index=randi([1 n]);
    
    h=selected_cell_members(selected_memebr_index);
    
    rep_h=rep(h);
end