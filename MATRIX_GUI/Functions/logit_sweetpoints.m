% Find the sweet points for logistic psychometric function
function swpts = logit_sweetpoints(phi)

%calculate the sweet points for a logit psychometric function, for
%parameter alpha, beta and lambda.

alpha = phi(1);
beta = phi(2);
gamma = phi(3);
lambda = phi(4);

swpts(1) = fminsearch(@(x) betavar_est(x,alpha,beta,gamma,lambda)+(x>=alpha)*1e10,alpha-10);
swpts(3) = fminsearch(@(x) betavar_est(x,alpha,beta,gamma,lambda)+(x<=alpha)*1e10,alpha+10);
swpts(2) = fminsearch(@(x) alphavar_est(x,alpha,beta,gamma,lambda),alpha);
swpts = sort(swpts);

    function sigmaalphasq = alphavar_est(x,alpha,beta,gamma,lambda)

    term1 = exp(2*beta*(alpha-x));
    term2 = (1+exp(beta*(x-alpha))).^2;
    term3 = -gamma+(lambda-1)*exp(beta*(x-alpha));
    term4 = 1-gamma+lambda*exp(beta*(x-alpha));
    term5 = beta^2*(-1+gamma+lambda)^2;

    sigmaalphasq = -term1.*term2.*term3.*term4./term5;
    end

    function sigmabetasq = betavar_est(x,alpha,beta,gamma,lambda)

    term1 = exp(2*beta*(alpha-x));
    term2 = (1+exp(beta*(x-alpha))).^2;
    term3 = -gamma+(lambda-1)*exp(beta*(x-alpha));
    term4 = 1-gamma+lambda*exp(beta*(x-alpha));
    term5 = (x-alpha).^2*(-1+gamma+lambda)^2;

    sigmabetasq = -term1.*term2.*term3.*term4./term5;
    end

end