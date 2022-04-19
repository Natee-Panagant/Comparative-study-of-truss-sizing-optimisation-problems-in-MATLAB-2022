clear all; clc;
load(['rst_' num2str(1,'%03.f') '_' num2str(1,'%03.f') '.mat']);
Nrun=30;
Ntest=8;
Nalgo=14;
% Final_Pareto.x=cell(Nalgo,Nrun);
% Final_Pareto.f=cell(Nalgo,Nrun);
% Final_Pareto.g=cell(Nalgo,Nrun);
title_list={'10-bar';
            '25-bar';
            '37-bar';
            '60-bar';
            '72-bar';
            '120-bar';
            '200-bar';
            '942-bar'};
        
for i=1:Ntest
    i
    for j=1:Nalgo
        clear PPareto FPareto GPareto
        filename=['rst_' num2str(i,'%03.f') '_' num2str(j,'%03.f') '.mat'];
        load(filename);
        for k=1:Nrun
            if numel(rst(k).ppareto)~=100
                error('check');
            end
            Final_Pareto.x{j,k}=rst(k).ppareto{100};
            Final_Pareto.f{j,k}=rst(k).fpareto{100};
            Final_Pareto.g{j,k}=rst(k).gpareto{100};
        end
    end
    save(['Final_Pareto_' title_list{i} '.mat'],'Final_Pareto','-v7.3');
end
