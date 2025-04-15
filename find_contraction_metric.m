function [M,C,x,y,sigma] = find_contraction_metric(degree,G)
x = sdpvar(2,1);
y = sdpvar(2,1);

% Define the differential equations
% dxdt = [3*x(2) - 3/2*(x(1)^2) - 1/2*(x(1)^3); -x(1) - x(2)];
% dxdt = [-x(1) + x(1)*x(2); -x(2)];
% dxdt = [4 - 0.5*x(1) + x(2) - (1-x(2))*x(1); -x(1) + (1-x(2))*x(1)];
% dxdt = [x(2); -(x(1)^2 + 100)*x(2) + x(1)];
% dxdt = [-x(1) + x(3); x(1)^2 - x(2) - 2*x(1)*x(3) + x(3); -x(2)];
dxdt = [-x(2) - 3/2*(x(1)^2) - 1/2*(x(1)^3);3*x(1) - x(2)];

%Finds a contraction Metric M and R that satisifes contraction condition
v = monolist(x',degree,0);
C = sdpvar(length(v)*2,length(v)*2);
M = kron(v,eye(2))'*C*kron(v,eye(2));

rho_1 = polynomial(x,degree*2,0);
rho_1 = 0;

v_ball = monolist([x(:)]',2,0);
L = sdpvar(length(v_ball),length(v_ball));
sigma = v_ball'*L*v_ball;
% sigma = 0;

% Z = jacobian(dxdt,x)*M + M*(jacobian(dxdt,x))' - jacobian(M,x(1))*dxdt(1) ... 
%  - jacobian(M,x(2))*dxdt(2) + 0.05*M - rho_1*G*G';

Z = jacobian(dxdt,x)'*M + M*(jacobian(dxdt,x)) + jacobian(M,x(1))*dxdt(1) ... 
 + jacobian(M,x(2))*dxdt(2) + 0.05*M - rho_1*G*G';

F = [
    C >= 0,...
    L >= 0,...
    sum(diag(L)) >= 0.5,...
    sum(diag(C)) >= 0.5,...
    % sos(-y'*(Z)*y + y'*y*sigma*(4-x'*x)),...
    sos(-y'*(Z)*y),...
    ];

options = sdpsettings('solver','mosek','verbose',1); % Below is nonlinear metric

% solvesos(F,[],options,[C(:);L(:)]);
% solvesos(F,[],options,[C(:);coefficients(rho_1,x)]);
obj = sum(abs(C));
% solvesos(F,[],options,[C(:)]);
optimize(F,obj,options);

for i = 1:size(C,1)
    for j = 1:size(C,2)
        if (abs(value(C(i,j))) < 1e-5)
            C(i,j) = 0;
        end
    end
end
M = kron(v,eye(2))'*value(C)*kron(v,eye(2));

% for i = 1:size(L,1)
%     for j = 1:size(L,2)
%         if (abs(value(L(i,j))) < 1e-5)
%             L(i,j) = 0;
%         end
%     end
% end
% sigma = v_ball'*value(L)*v_ball;

end