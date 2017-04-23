% Create exact function and its derivative
N_exact = 301; % number of sample points for exact function
x_exact = linspace(0, 6*pi, N_exact);
f_exact = sin(x_exact).*exp(-0.3*x_exact);
f_derivative_exact = cos(x_exact).*exp(-0.3*x_exact) - 0.3*sin(x_exact).*exp(-0.3*x_exact);

% plot exact function
figure (1);
plot (x_exact, f_exact, 'k-', 'linewidth', 1.5);
axis ([0 6*pi -1 1]); grid on;
xlabel ('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 16);

% create exact function for pi/5 sampling period
% and its finite difference derivatives
N_a = 31; % number of points for pi/5 sampling period
x_a = linspace(0, 6*pi, N_a);
f_a = sin(x_a).*exp(-0.3*x_a);
f_derivative_a = cos(x_a).*exp(-0.3*x_a)-0.3*sin(x_a).*exp(-0.3*x_a);

dx_a = pi/5;
f_derivative_forward_a = zeros(1, N_a);
f_derivative_backward_a = zeros(1, N_a);
f_derivative_central_a = zeros(1, N_a);
f_derivative_forward_a(1:N_a-1) = (f_a(2:N_a) - f_a(1:N_a-1))/dx_a;
f_derivative_backward_a(2:N_a) = (f_a(2:N_a) - f_a(1:N_a-1))/dx_a;
f_derivative_central_a(2:N_a-1) = (f_a(3:N_a) - f_a(1:N_a-2))/(2*dx_a);

% create exact function for pi/10 sampling period
% and its finit difference derivatives
N_b = 61; % number of points for pi/10 sampling period
x_b = linspace(0, 6*pi, N_b);
f_b = sin(x_b).*exp(-0.3*x_b);
f_derivative_b = cos(x_b).*exp(-0.3*x_b)-0.3*sin(x_b).*exp(-0.3*x_b);

dx_b = pi/10;
f_derivative_forward_b = zeros(1, N_b);
f_derivative_backward_b = zeros(1, N_b);
f_derivative_central_b = zeros(1, N_b);
f_derivative_forward_b(1:N_b-1) = (f_b(2:N_b) - f_b(1:N_b - 1))/dx_b;
f_derivative_backward_b(2:N_b) = (f_b(2:N_b) - f_b(1:N_b - 1))/dx_b;
f_derivative_central_b(2:N_b-1) = (f_b(3:N_b) - f_b(1:N_b - 2))/(2*dx_b);

% plot exact derivatives of the function and its finit difference
% derivatives using pi/5 sampling period
figure(2);
plot(x_exact, f_derivative_exact, 'k', ...
    x_a(1:N_a-1), f_derivative_forward_a(1:N_a-1), 'b-', ...
    x_a(2:N_a), f_derivative_backward_a(2:N_a), 'r-', ...
    x_a(2:N_a-1), f_derivative_central_a(2:N_a-1), ':ms', ...
    'MarkerSize', 4, 'linewidth', 1.5);
set (gca, 'FontSize', 12, 'fontweight', 'demi');
axis([0 6*pi -1 1]);
grid on;
legend('exact', 'forward difference', ...
    'backward difference', 'central difference');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$f''(x)$', 'Interpreter', 'latex', 'FontSize', 16);
text(pi, 0.6, '$\Delta_x = \pi/5$', 'Interpreter', ...
    'latex', 'fontsize', 16, 'BackgroundColor', 'w', 'EdgeColor', 'k');

% plot error for finite difference derivatives
% using pi/5 sampling period
error_forward_a = f_derivative_a - f_derivative_forward_a;
error_backward_a = f_derivative_a - f_derivative_backward_a;
error_central_a = f_derivative_a - f_derivative_central_a;

figure(3);
plot(x_a(1:N_a-1), error_forward_a(1:N_a-1), 'b-',...
    x_a(2:N_a), error_backward_a(2:N_a), 'r-',...
    x_a(2:N_a-1), error_central_a(2:N_a-1), ':ms', ...
    'MarkerSize', 4, 'linewidth', 1.5);

set(gca, 'FontSize', 12, 'fontweight', 'demi');
axis([0 6*pi -0.2 0.2]);
grid on;
legend('forward difference', 'backward difference', ...
    'central difference');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('error $[f''(x)]$', 'Interpreter', 'latex', 'FontSize', 16);
text(pi, 0.15, '$\Delta_x = \pi/5$', 'Interpreter', 'latex', ...
    'fontsize', 16, 'BackgroundColor', 'w', 'EdgeColor', 'k');

% plot error for finite difference dirivatives
% using pi/10 sampling period
error_forward_b = f_derivative_b - f_derivative_forward_b;
error_backward_b = f_derivative_b - f_derivative_backward_b;
error_central_b = f_derivative_b - f_derivative_central_b;

figure(4);
plot(x_b(1:N_b-1), error_forward_b(1:N_b-1), 'b-', ...
    x_b(2:N_b), error_backward_b(2:N_b), 'r-', ...
    x_b(2:N_b-1), error_central_b(2:N_b-1), ':ms', ...
    'MarkerSize', 4, 'linewidth', 1.5);
set(gca, 'FontSize', 12, 'fontweight', 'demi');
axis([0 6*pi -0.2 0.2]);
grid on;
legend('forward difference', 'backward difference', 'central difference');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('error $[f''(x)]$', 'Interpreter', 'latex', 'FontSize', 16);
text(pi, 0.15, '$\Delta_x = \pi/10$', 'Interpreter', 'latex', ...
    'FontSize', 16, 'BackgroundColor', 'w', 'EdgeColor', 'k');