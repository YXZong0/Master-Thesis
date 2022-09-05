function plotres(a,b,res)
    logA = log(a+b)-log(b);
    logB = log(b);
    x=logB:0.01:(logA+logB);
    plot(x,pdf('beta',(exp(x)-b)/a,3,4)./a.*exp(x))
    hold on;
    nB = length(res);
    k=0:1:(nB-1);
    parAlpha=repmat(k'+1,1,length(x));
    parBeta=repmat(nB-1-k'+1,1,length(x));
    %x_hat = reshape(res2(1:nB),1,length(res2(1:nB)))*(betapdf(repmat((log(x)-b)./a,nB,1),parAlpha,parBeta).*repmat(1./(x*a),nB,1));
    x_hat = reshape(res(1:nB),1,length(res(1:nB)))*(betapdf(repmat((x-logB)./logA,nB,1),parAlpha,parBeta).*repmat(1./(logA),nB,1));
    hold on;
    plot(x,x_hat)
end