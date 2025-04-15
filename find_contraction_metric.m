function [M,C,x,y,rho_1] = find_contraction_metric(degree,G)
x = sdpvar(3,1);
y = sdpvar(3,1);

% Define the differential equations
% dxdt = [3*x(2) - 3/2*(x(1)^2) - 1/2*(x(1)^3); -x(1) - x(2)];
% dxdt = [-x(1) + x(1)*x(2); -x(2)];
% dxdt = [4 - 0.5*x(1) + x(2) - (1-x(2))*x(1); -x(1) + (1-x(2))*x(1)];
% dxdt = [x(2); -(x(1)^2 + 100)*x(2) + x(1)];

dxdt = [-x(1) + x(3); x(1)^2 - x(2) - 2*x(1)*x(3) + x(3); -x(2)];

%Finds a contraction Metric M and R that satisifes contraction condition
v = monolist(x',degree,0);
C = sdpvar(length(v)*3,length(v)*3);
M = kron(v,eye(3))'*C*kron(v,eye(3));


v_rho = monolist(x',2,0);
W_rho = sdpvar(length(v_rho),length(v_rho));
rho_1 = v_rho'*W_rho*v_rho;

Z = jacobian(dxdt,x)*M + M*(jacobian(dxdt,x))' - jacobian(M,x(1))*dxdt(1) ... 
 - jacobian(M,x(2))*dxdt(2) - jacobian(M,x(3))*dxdt(3) + 0.05*M - rho_1*G*G';
c = coefficients(jacobian(M,x(1))*G(1,1) + jacobian(M,x(2))*G(2,1) + jacobian(M,x(3))*G(3,1),x);

F = [
    C >= 0,...
    sum(diag(C)) >= 0.5,...
    c == 0,...
    sum(abs(W_rho(:))) <= 10,...
    sos(-y'*Z*y),...
    ];

options = sdpsettings('solver','mosek','verbose',1); % Below is nonlinear metric

solvesos(F,[],options,[C(:);W_rho(:)]);

for i = 1:size(C,1)
    for j = 1:size(C,2)
        if (abs(value(C(i,j))) < 1e-5)
            C(i,j) = 0;
        end
    end
end
M = kron(v,eye(3))'*value(C)*kron(v,eye(3));

[c,p] = coefficients(rho_1,x); 
for i = 1:size(c,1)
    if(abs(value(c(i))) < 1e-5)
        c(i) = 0;
    end
end
rho_1 = p'*value(c);

end