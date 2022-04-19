function pareto_plot(fppareto,fpareto,gpareto,varargin)
    if numel(varargin)>0
        fun=[varargin{1} ' : '];
    else
        fun=[];
    end
    figure(1);clf;
    plot(fppareto(1,:),fppareto(2,:),'ro');
    title([fun 'Gmax = ' num2str(max(max(gpareto)))]);
    pause(0);
end

