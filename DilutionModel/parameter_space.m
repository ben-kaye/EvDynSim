% parameter space / contour plot:

R = [ 0.1 0.5 1 2 ];
alpha = 1/13.3;
beta = 0.8;
s0 = 0.7;
m = 100;

% f = @(a, b) (s0 - m*a)/(a/alpha - 1) - (R*(s0 - m*a) + s0 - m*b)/(R*(1 - a/alpha) + (1 - b/alpha));


N = 20;



% for i = 1:size(x2,2)
%     f_k = @(a) f(a, x2(i));
%     x1(i) = fsolve(f_k, 0.01);
% end


x2 = linspace(0, 6e-2, N);
x1 = zeros(length(R),length(x2));

for ii = 1:length(R)
    R_k = R(ii);
    f = (s0 - m*a)/(a/alpha - beta) - (R_k*(s0 - m*a) + s0 - m*b)/(R_k*(beta - a/alpha) + (beta - b/alpha)) == 0;

    for i = 1:N
        f_temp = subs(f, b, x2(i));
        z = solve(f_temp);
        x1(ii,i) = z(1);
    end
end

plot(x2, x1);
xlabel('$\chi_2$', 'Interpreter','latex', 'FontSize', 16)
ylabel('$\chi_1$', 'Interpreter','latex','FontSize', 16)


labels = cell(1,length(R));
for o = 1:length(R)
    labels{o} = sprintf('R=%.3f',R(o));
end
legend(labels) 
