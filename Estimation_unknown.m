function [result,BIC] = Estimation_unknown(W,nbTerms,a,b)
% W: matrix of size nx1 (where n is the sample size), containing the
% observed covariate
% nbTerms: maximum number of terms to be considered in the Bernstein polynomial
% a,b: lower bound and upper bound 

% bound after transformation 
logA = log(a+b)-log(b);
logB = log(b);

% Objects for the output
result = repmat(-999,nbTerms+2,nbTerms);
BIC = repmat(-999,nbTerms,1);

% starting values for v, the measurement error standard deviation
x01temp = (std(W))*[0.2,0.6,1]; 
[m,n] = meshgrid( x01temp, x01temp' );
[sig0(:,1),sig0(:,2) ] = deal( reshape(m,[],1), reshape(n,[],1) ); 
sig0 = sig0';

for nB = 1:nbTerms % nB = 1+order of Bernstein polynomial

    %constraints for the optimization
    beq = 1;
    A = [];
    b = [];
    Aeq = [ones(1,nB),zeros(1,2)];
    lb = [zeros(nB,1);0.01;0.01];
    ub = [ones(nB,1);std(W)*2;std(W)];
    mynlcon = [];   

    % starting values for the coefficients of the Bernstein polynomial
    xtemp=unifrnd(0,1,nB,9); 
    xtemp2 = xtemp./repmat(sum(xtemp,1),nB,1);
    x01 = [xtemp2;sig0];

    optionsIP = optimoptions(@fmincon,'MaxIter',300,'Algorithm','interior-point','Display','iter','UseParallel',true);
    %optionsAS = optimoptions(@fmincon,'MaxIter',300,'Algorithm','sqp','Display','iter','UseParallel',true);

    % objective function
    g_lik = @(coefsig) likelihood_qua(coefsig,nB-1,W,logA,logB);
    [res2,fval2,exitflag2,output2] = fmincon(g_lik,x01(:,1),A,b,Aeq,beq,lb,ub,mynlcon,optionsIP);
    for j = [3,5,7,9]
        [res22,fval22,exitflag22,output22] = fmincon(g_lik,x01(:,j),A,b,Aeq,beq,lb,ub,mynlcon,optionsIP);
        if(fval22<fval2)
            res2=res22;
            fval2=fval22;
            exitflag2=exitflag22;
        end
    end
    result(:,nB) = [res2;zeros(nbTerms+2-length(res2),1)];
    BIC(nB,:) = 2*fval2+log(length(W))*(nB+2);
    end
end