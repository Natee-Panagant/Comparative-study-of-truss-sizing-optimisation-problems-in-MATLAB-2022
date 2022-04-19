clearvars; close all; clc;
load rst_001_001.mat
rst1=rst;
load rst_001_007.mat
rst2=rst;
load rst_001_012.mat
rst3=rst;

figure(1); hold on;
plot(rst1(1).fpareto{end}(1,:),rst1(1).fpareto{end}(2,:),'gs');
plot(rst2(1).fpareto{end}(1,:),rst2(1).fpareto{end}(2,:),'bo');
legend('Front no. 1','Front no. 2');

figure(2); hold on;
[ppareto,fpareto,gpareto,A]=ref_front_sorting(rst1(1).ppareto{end},rst1(1).fpareto{end},rst1(1).gpareto{end},rst2(1).ppareto{end},rst2(1).fpareto{end},rst2(1).gpareto{end},[],200);
plot(rst1(1).fpareto{end}(1,:),rst1(1).fpareto{end}(2,:),'gs');
plot(rst2(1).fpareto{end}(1,:),rst2(1).fpareto{end}(2,:),'bo');
plot(fpareto(1,:),fpareto(2,:),'rx');
refpoint(1,1)=max([rst1(1).fpareto{end}(1,:),rst2(1).fpareto{end}(1,:)]);
refpoint(2,1)=max([rst1(1).fpareto{end}(2,:),rst2(1).fpareto{end}(2,:)]);
plot(refpoint(1),refpoint(2),'m^');
legend('Front no. 1','Front no. 2','None Dominated Solutions from Front no. 1,2','RefPoint from Front no. 1,2');

figure(3); hold on;
plot(fpareto(1,:),fpareto(2,:),'ko');
plot(refpoint(1),refpoint(2),'k^');
legend('None Dominated Solutions from Front no. 1,2','RefPoint from Front no. 1,2');

figure(4); hold on;
plot(fpareto(1,:),fpareto(2,:),'ko');
plot(rst3(1).fpareto{end}(1,:),rst3(1).fpareto{end}(2,:),'c^');
plot(refpoint(1),refpoint(2),'k^');
legend('None Dominated Solutions from Front no. 1,2','Front no. 3','RefPoint from Front no. 1,2');

figure(5); hold on;
plot(fpareto(1,:),fpareto(2,:),'ko');
plot(rst3(1).fpareto{end}(1,:),rst3(1).fpareto{end}(2,:),'c^');
[ppareto1,fpareto1,gpareto1,A]=ref_front_sorting(rst3(1).ppareto{end},rst3(1).fpareto{end},rst3(1).gpareto{end},ppareto,fpareto,gpareto,[],200);
plot(fpareto1(1,:),fpareto1(2,:),'rx');

refpoint(1,1)=max([refpoint(1),rst3(1).fpareto{end}(1,:)]);
refpoint(2,1)=max([refpoint(2),rst3(1).fpareto{end}(2,:)]);
plot(refpoint(1),refpoint(2),'m^');
legend('None Dominated Solutions from Front no. 1,2','Front no. 3','None Dominated Solutions from Front no. 1,2,3','RefPoint from Front no. 1,2,3');

figure(6); hold on;
plot(fpareto1(1,:),fpareto1(2,:),'ko');
plot(refpoint(1),refpoint(2),'k^');
legend('None Dominated Solutions from Front no. 1,2,3','RefPoint from Front no. 1,2,3');

for i=1:6
    figure(i);
    saveas(gcf,['ref_method' num2str(i)],'tif');
end