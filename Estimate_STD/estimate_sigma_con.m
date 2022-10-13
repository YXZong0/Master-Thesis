function [sig_mul,sig_add,sig_mul_boots,sig_add_boots,sig_mul_CI_low,sig_mul_CI_up,sig_add_CI_low,sig_add_CI_up] = estimate_sigma_con(X,I,J,n)
    beq = [];
    A = [];
    b = [];
    Aeq = [];
    lb = [1,0];
    
    sig_mul = zeros(I*10+1,J*10+1);
    sig_add = zeros(I*10+1,J*10+1);
    sig_mul_boots = zeros(I*10+1,J*10+1,n);
    sig_add_boots = zeros(I*10+1,J*10+1,n);
    sig_mul_CI_low = zeros(I*10+1,J*10+1);
    sig_mul_CI_up = zeros(I*10+1,J*10+1);
    sig_add_CI_low = zeros(I*10+1,J*10+1);
    sig_add_CI_up = zeros(I*10+1,J*10+1);
    

    for i = 0:0.1:I 
        for j = 0:0.1:J

            W = X.*lognrnd(0,i,size(X,1),size(X,2))+normrnd(0,j,size(X,1),size(X,2));
            Y = mean(W.^2,2);
            ub = [];
            X0 = (-(sum(W.^2,2)-sum(W,2).^2)/20);
            Data_ori = [Y,X0,ones(length(X0),1)];
            
            [res,fval,exitflag,output] = lsqlin(Data_ori(:,2:3),Data_ori(:,1),A,b,Aeq,beq,lb,ub);
            %[res,fval,exitflag,output] = lsqlin(Data_ori(:,2:3)./std(Data_ori(:,2)),Data_ori(:,1)./std(Data_ori(:,2)),A,b,Aeq,beq,lb,ub);
            sig_mul(round(i*10)+1,round(j*10)+1) = sqrt(log(res(1)));
            sig_add(round(i*10)+1,round(j*10)+1) = sqrt(res(2));

            for k=1:n
                Data_boots = datasample(Data_ori,size(Data_ori,1),1);
                [res2,fval2,exitflag2,output2] = lsqlin(Data_boots(:,2:3),Data_boots(:,1),A,b,Aeq,beq,lb,ub);
                %[res2,resnorm,residual,exitflag,output,lambda] = lsqlin(Data_boots(:,2:3)./sqrt(Data_boots(:,2)),Data_boots(:,1)./sqrt(Data_boots(:,2)),A,b,Aeq,beq,lb,ub);
                sig_mul_boots(round(i*10)+1,round(j*10)+1,k) = sqrt(log(res2(1)));
                sig_add_boots(round(i*10)+1,round(j*10)+1,k) = sqrt(res2(2));
            end
        end
    end
            sig_mul_CI_low = 2*sig_mul - quantile(sig_mul_boots,0.975,3);
            sig_mul_CI_up = 2*sig_mul - quantile(sig_mul_boots,0.025,3);
            sig_add_CI_low = 2*sig_add - quantile(sig_add_boots,0.975,3);
            sig_add_CI_up = 2*sig_add - quantile(sig_add_boots,0.025,3);
end
