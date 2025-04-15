function [u] = compute_feedback(rho_1,inv_Q,vars,x,G)
 u = integral(@(s) -(1/2)*replace(rho_1,vars,(1-s)*x)*G'*inv_Q*(x),0,1,'ArrayValued',true);
end