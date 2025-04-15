function [phi_data,x_dot,d] = data_matrix(T,epsilon)
% build according to eq.(14)
global x_glo
global name_of_system
x1 = sym('x1', 'real');           % X1 is used as data
x2 = sym('x2', 'real');
% phi = monomials(vars, 1:d_f);     % monomial vector of f
if(strcmp(name_of_system,'unstable-linear'))
    phi = [x1;x2];
else
    phi = [x1; x2; x1^2; x1^3];
end
% phi = [x1; x2; x1*x2];
% phi = [x1; x2];
s_f = size(phi,1);       % size of monomials above
phi_fun = matlabFunction(phi);
phi_data = zeros(T, s_f);
for i = 1:T
    phi_data(i, :) = feval(phi_fun, x_glo(i, 2), x_glo(i, 3));
end

x_dot = x_glo(1:T, 4:5);
d = x_dot(1,:)' + ones(size(x_dot,2),1)*epsilon;
d = vertcat(d, -x_dot(1,:)' + ones(size(x_dot,2),1)*epsilon);
for i = 2:T
    d = vertcat(d, x_dot(i,:)' + ones(size(x_dot,2),1)*epsilon);
    d = vertcat(d, -x_dot(i,:)' + ones(size(x_dot,2),1)*epsilon);
end
end