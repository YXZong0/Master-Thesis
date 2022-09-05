clc;
clear all;
close all;
sig_mul=0.2;
sig_add=0.3;
a = 2; b = 2; 
n = 300; % sample size
W = (a*betarnd(3,4,n,1)+b).*lognrnd(0,sig_mul,n,1)+normrnd(0,sig_add,n,1); %toy sample
nB = 5; % order
[res1,bic1] = Estimation_unknown(W,nB,a,b); % assume two std unknown 
[res2,bic2] = Estimation_known(W,nB,a,b,sig1,sig2); % assume two std known
