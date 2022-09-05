function f = htoInt(z,t,y,sig_mul,sig_add)
    % h function
    %f = (lognpdf(z,t,sig_mul).*normpdf(y,z,sig_add));
    f = normpdf((log(z)-t)./sig_mul).*(normpdf((y-z)./sig_add))./(z*sig_add*sig_mul);
end

