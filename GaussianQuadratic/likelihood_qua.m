function f = likelihood_qua(coefsig,nB,Y,A,B)
% calculate the likelihood function using gaussian quadrature
% coefsig: initial parameters 
% nB: the order of Bernstein
% Y: observations
% A,B: boundary
    %fy = @(y,theta,sig_mul,sig_add,nB,A,B) integral2(@(z,t) htoInt(z,t,y,sig_mul,sig_add).*densityBernstein(t,nB,theta)./A,0,Inf,B,A+B);
    %fy = @(y,theta,sig_mul,sig_add,nB,A,B) quadgk(@(t) arrayfun(@(t) quadgk(@(z) htoInt(z,t,y,sig_mul,sig_add).*densityBernstein(t,nB,theta)./A,0,Inf),t),B,A+B);
    %fy = @(y,theta,sig_mul,sig_add,nB,A,B) integral(@(z) arrayfun(@(z) integral(@(t) htoInt(z,t,y,sig_mul,sig_add).*densityBernstein((t-B)/A,nB,theta)./A,B,A+B),z),0,Inf);
    fy = @(y,theta,sig_mul,sig_add,nB,A,B) quadgk(@(z) arrayfun(@(z) quadgk(@(t) htoInt(z,t,y,sig_mul,sig_add).*densityBernstein((t-B)/A,nB,theta)./A,B,A+B),z),0,Inf);
    theta = coefsig(1:(nB+1)); % coefficients 
    sig_mul = coefsig(nB+2); % multiplicative std
    sig_add = coefsig(nB+3); % additive std
    %fy_res = arrayfun(@(Y) fy(Y,theta,sig_mul,sig_add,nB,A,B), Y);
    fy_res = zeros(length(Y),1);
    %Y = gpuArray(Y);
    parfor i = 1:length(Y)
        fy_res(i) = fy(Y(i),theta,sig_mul,sig_add,nB,A,B);
        %fy_res(i) = parfeval(backgroundPool,@fy,Y(i),theta,sig_mul,sig_add,nB,A,B);
    end
    f = -sum(log(fy_res)); % to be minimized
    if f == Inf
        f = 99999; 
    end
end
