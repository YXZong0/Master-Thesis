function f = likelihood_qua(coefsig,nB,Y,A,B)
% calculate the likelihood function
% coefsig: initial parameters 
% nB: the order of Bernstein
% Y: observations
% A,B bound
    %h = @(t,y,sig_mul,sig_add) quadgk(@(z) arrayfun(@(t) htoInt(z,t,y,sig_mul,sig_add),t,'UniformOutput',false), 0, 9999);
    %fy = @(y,theta,sig_mul,sig_add,nB,A,B) integral2(@(z,t) htoInt(z,t,y,sig_mul,sig_add).*densityBernstein(t,nB,theta)./A,0,Inf,B,A+B);
    %fy = @(y,theta,sig_mul,sig_add,nB,A,B) quadgk(@(t) arrayfun(@(t) quadgk(@(z) htoInt(z,t,y,sig_mul,sig_add).*densityBernstein(t,nB,theta)./A,0,Inf),t),B,A+B);
    fy = @(y,theta,sig_mul,sig_add,nB,A,B) quadgk(@(z) arrayfun(@(z) quadgk(@(t) htoInt(z,t,y,sig_mul,sig_add).*densityBernstein((t-B)/A,nB,theta)./A,B,A+B),z),0,Inf);
    %fy = @(y,theta,sig_mul,sig_add,nB,A,B) integral(@(z) arrayfun(@(z) integral(@(t) htoInt(z,t,y,sig_mul,sig_add).*densityBernstein((t-B)/A,nB,theta)./A,B,A+B),z),0,Inf);
    theta = coefsig(1:(nB+1));
    sig_mul = coefsig(nB+2);
    sig_add = coefsig(nB+3);
    fy_res = arrayfun(@(Y) fy(Y,theta,sig_mul,sig_add,nB,A,B), Y);
    f = -sum(log(fy_res));
    if f == Inf
        f = 99999; 
    end
    %f=fy_res;
end