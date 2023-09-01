function p0 = setPrior(phi, s)           

%set the prior distribution

    if strcmp(s.dist, 'norm')
        if strcmp(s.scale, 'lin')
            p0 = normpdf(phi,s.mu,s.std);
        elseif strcmp(s.scale, 'log')
            p0 = normpdf(log10(phi),s.mu,s.std);
        end
    elseif strcmp(s.dist, 'flat')
        p0 = ones(size(phi));
    end
end