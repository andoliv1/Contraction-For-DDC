function [u] = compute_feedback(rho_1,Q,vars,x,G)
 % u = integral(@(s)
 % -(1/2)*replace(rho_1,vars,(1-s)*x)*G'*inv_Q*(x),0,1,'ArrayValued',true);
 X = 0:0.01:1;
 Y = zeros(size(X,1));
 counter = 1;
 for i = 0:0.01:1
    Y(counter) = -(1/2)*replace(rho_1,vars,(1-i)*x)*G'*inv(replace(Q,vars,(1-i)*x))*(x);
    counter = counter +1;
 end
 u = trapz(X,Y);
end