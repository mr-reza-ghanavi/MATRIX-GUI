% eliminate the 1 or 0 probability to prevent the nummerical errors due to
% taking log of zeros.
function p = prepare_prob(p)
    p = p*(1-1e-10);
    p = p/sum(sum(sum(p)));
end
