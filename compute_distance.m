function [d] = compute_distance(Q,x,p,j)
 % u = integral(@(s)
 X = 0:0.01:1;
 Y = zeros(size(X,1));
 counter = 1;
 for i = 0:0.01:1
    Y(counter) = sqrt(p'*inv(replace(Q,x,(1-i)*p))*p);
    counter = counter +1;
 end
 'done'
 j
 d = trapz(X,Y);
end