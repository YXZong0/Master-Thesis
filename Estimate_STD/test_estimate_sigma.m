clc;
clear all;
close all;

% pd = makedist('Normal');
% t = truncate(pd,-2,2);
% X = repmat(random(t,300,1)+3,1,5);
% 
% [sig_mul,sig_add,sig_mul_boots,sig_add_boots,sig_mul_CI_low,sig_mul_CI_up,sig_add_CI_low,sig_add_CI_up] = estimate_sigma_uncons(X,0.4,0.6,200);

%sig_mul,sig_add,sig_mul_boots,sig_add_boots,sig_mul_CI_low,sig_mul_CI_up,sig_add_CI_low,sig_add_CI_up] = estimate_sigma_cons(X,0.4,0.6,200);
% 
pd = makedist('Exponential',10);
t = truncate(pd,0,20);
X = repmat(random(t,300,1)+3,1,5);

[sig_mul,sig_add,sig_mul_boots,sig_add_boots,sig_mul_CI_low,sig_mul_CI_up,sig_add_CI_low,sig_add_CI_up] = estimate_sigma_con(X,0.4,0.8,200);
