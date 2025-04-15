function [phi_data,x_dot,d] = data_matrix(T,eps,dimensions)
% build according to eq.(14)
global x_glo
x1 = sym('x1', 'real');           % X1 is used as data
x2 = sym('x2', 'real');
x3 = sym('x3', 'real');
phi = [x1;x2;x3;x1*x3;x1^2];
s_f = size(phi,1);       % size of monomials above
phi_fun = matlabFunction(phi);
phi_data = zeros(T, s_f);
for i = 1:T
    phi_data(i, :) = feval(phi_fun, x_glo(i, 2),x_glo(i, 3), x_glo(i, 4));
end
x_dot = x_glo(1:T, 2+dimensions:end);
d = x_dot(1,:)' + ones(size(x_dot,2),1)*eps;
d = vertcat(d, -x_dot(1,:)' + ones(size(x_dot,2),1)*eps);
for i = 2:T
    d = vertcat(d, x_dot(i,:)' + ones(size(x_dot,2),1)*eps);
    d = vertcat(d, -x_dot(i,:)' + ones(size(x_dot,2),1)*eps );
end
end