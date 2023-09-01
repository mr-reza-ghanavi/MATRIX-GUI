function space = setParSpace(s)

%set parameter space on a linear scale, or a logarithmic scale

    if s.N == 1
        space = s.limits(1);
    else
        if strcmp(s.scale, 'lin')
            space = linspace(s.limits(1),s.limits(2),s.N)';
        elseif strcmp(s.scale, 'log')
            space = logspace(log10(s.limits(1)),log10(s.limits(2)),s.N)';
        end          
    end
end