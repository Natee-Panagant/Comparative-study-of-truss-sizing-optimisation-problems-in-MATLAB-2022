function fp = fpenal00(f,g)
% Penalty Function
fp=f+(1e5)*max(max(g,[],1),0);
end

