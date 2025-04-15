function [x,y,mu,rho_1,Q,W,c,rho_Y] = ddc_poly(phi_data,T,d,degree,degree_rho,G)
% Define sdpvariables
x = sdpvar(3,1);
y = sdpvar(3,1);
allvars = [x(:);y(:)];
% Define the contraction metric
v = monolist(x',degree,0);
W = sdpvar(length(v)*size(x,1),length(v)*size(x,1));
Q = kron(v,eye(size(x,1)))'*W*kron(v,eye(size(x,1)));
% define polynomial to \|y\|_2 = 1
v_Y = monolist(allvars, 2, 0);
W_Y = sdpvar(length(v_Y),length(v_Y));
rho_Y = v_Y'*W_Y*v_Y;
% feedback gain
v_rho = monolist(x',degree_rho,0);
W_rho = sdpvar(length(v_rho),length(v_rho));
rho_1 = v_rho'*W_rho*v_rho;
% Set up vector of monomials
Phi = [x(1);x(2);x(3);x(1)*x(3);x(1)^2];
% Define the array and scalar to store the constraints
mu = sdpvar(2*size(x,1)*T,1);
sum_d = 0;
% define monomials that will define the constraints
d_f = 2;
p_x = monolist(x,d_f + 2*degree, 0);
p_y = monolist(y,2,1);

% define variable for formatting data
index = 1;

for i = 1:2*size(x,1):2*size(x,1)*T
    C_1 = sdpvar(size(p_x,1),size(p_y,1));
    C_2 = sdpvar(size(p_x,1),size(p_y,1));
    C_3 = sdpvar(size(p_x,1),size(p_y,1));
    C_4 = sdpvar(size(p_x,1),size(p_y,1));
    C_5 = sdpvar(size(p_x,1),size(p_y,1));
    C_6 = sdpvar(size(p_x,1),size(p_y,1));

    [mu(i)] = p_x'*C_1*p_y;
    [mu(i+1)] = p_x'*C_2*p_y;
    [mu(i+2)] = p_x'*C_3*p_y;
    [mu(i+3)] = p_x'*C_4*p_y;
    [mu(i+4)] = p_x'*C_5*p_y;
    [mu(i+5)] = p_x'*C_6*p_y;
    sum_d = sum_d + mu(i)*d(i) + mu(i+1)*d(i+1) + mu(i+2)*d(i+2) + mu(i+3)*d(i+3) + mu(i+4)*d(i+4) + mu(i+5)*d(i+5);
    [phi_mat_pos_1,phi_mat_pos_2,phi_mat_pos_3,phi_mat_neg_1,phi_mat_neg_2,phi_mat_neg_3] = formatData(phi_data(index,:)');
    if i == 1
        sum_phi = mu(i)*phi_mat_pos_1 + mu(i+1)*phi_mat_pos_2 + mu(i+2)*phi_mat_pos_3 + mu(i+3)*phi_mat_neg_1 + mu(i+4)*phi_mat_neg_2 + mu(i+5)*phi_mat_neg_3;
    else
        sum_phi = sum_phi + mu(i)*phi_mat_pos_1 + mu(i+1)*phi_mat_pos_2 + mu(i+2)*phi_mat_pos_3 + mu(i+3)*phi_mat_neg_1 + mu(i+4)*phi_mat_neg_2 + mu(i+5)*phi_mat_neg_3;
    end
    index = index + 1;
end

% the + 2 below is related to the number of states
tr_phi = zeros(size(Phi,1)*size(x,1),1);
I = eye(size(x,1)); 
for i=1:size(x,1)
    e_i = I(:,i);
    tr = trace(jacobian(Q,x(i))*(y*y'));
    tr_phi = tr_phi + tr*(kron(e_i',Phi'))';
end

% Get the coefficients of the constraints, so that we can enforce those to
% be zero
mat_phi_constraint = -2*jacobian(Phi,x)*Q*(y*y') + sum_phi;

c1 = coefficients(mat_phi_constraint(:) + tr_phi, allvars);
c2 = coefficients(jacobian(Q,x(1))*G(1,1) + jacobian(Q,x(2))*G(2,1) + jacobian(Q,x(3))*G(3,1),x);
options = sdpsettings('solver','mosek','verbose',1,'MSK_DPAR_DATA_TOL_X',1e-4, 'MSK_DPAR_INTPNT_CO_TOL_PFEAS', 1e-4,'MSK_DPAR_DATA_TOL_QIJ', 1e-4); % Below is nonlinear metric

F = [
    W >= 0,...
    sum(diag(W)) >= 0.1,...  
    sum(abs(W_rho(:))) <= 10,...
    sos(-y'*0.05*Q*y + rho_1*y'*G*G'*y - sum_d + rho_Y*(1-y'*y)),...
    sos(mu),...
    c1 == 0,...
    c2 == 0,...
    ];

solvesos(F,[],options,[coefficients(mu,allvars); W_rho(:); W(:); W_Y(:)]);

% throw out W's useless coefficients
for i = 1:size(W,1)
    for j = 1:size(W,2)
        if (abs(value(W(i,j))) < 1e-4)
            W(i,j) = 0;
        end
    end
end
W = value(W);
Q = kron(v,eye(size(x,1)))'*W*kron(v,eye(size(x,1)));

[c,p] = coefficients(rho_1,x);
for i = 1:size(c,1)
    if(abs(value(c(i))) < 1e-4)
        c(i) = 0;
    end
end
rho_1 = p'*value(c);

end