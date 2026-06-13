%% Ερώτημα (α) – Γραμμικοποίηση μέσω ανάδρασης (σωστός νόμος)
% Σύστημα:
%   x1_dot = x2 + theta*g(x1)
%   x2_dot = u
% με g(x1) = x1^3, θ* = 2.

clear; clc; close all;

%% Παράμετροι συστήματος και ελεγκτή
theta = 2;      % θ*
k1 = 2;         % π.χ. πόλοι στο -1 και -2 -> k1=2, k2=3
k2 = 3;

Tfinal = 10;    % χρόνος προσομοίωσης [s]

% g(x1) και g'(x1)
g  = @(x1) x1.^3;
gp = @(x1) 3*x1.^2;

%% Αρχικές συνθήκες (όπως δίνει η εκφώνηση)
x0_list = [ 0.5   0;    % 1: (0.5,0)
           -2    20;    % 2: (-2,20)
            2    20;    % 3: (2,20)
           -1    10];   % 4: (-1,10)

nIC = size(x0_list,1);
colors = lines(nIC);

%% -------- Figure 1: x1, x2 για μικρές αρχικές --------
figure(1);
subplot(2,1,1); hold on; grid on;
ylabel('x_1(t)');
title('Καταστάσεις x_1(t), x_2(t) για (0.5,0) - (-2,20) - (-1,10) μέθοδος (α)');

subplot(2,1,2); hold on; grid on;
ylabel('x_2(t)'); xlabel('t [s]');

%% -------- Figure 2: x1, x2 για (2,20) --------
figure(2);
subplot(2,1,1); hold on; grid on;
ylabel('x_1(t)');
title('Καταστάσεις x_1(t), x_2(t) για αρχική (2,20) μέθοδος (α)');

subplot(2,1,2); hold on; grid on;
ylabel('x_2(t)'); xlabel('t [s]');

%% -------- Figure 3: u(t) για μικρές αρχικές --------
figure(3); hold on; grid on;
xlabel('t [s]'); ylabel('u(t)');
title('Σήμα ελέγχου u(t) για (0.5,0) - (-2,20) - (-1,10) μέθοδος (α)');

%% -------- Figure 4: u(t) για (2,20) --------
figure(4); hold on; grid on;
xlabel('t [s]'); ylabel('u(t)');
title('Σήμα ελέγχου u(t) για αρχική (2,20) μέθοδος (α)');

%% Προσομοίωση
for i = 1:nIC
    x0 = x0_list(i,:).';

    [t,x] = ode45(@(t,x) dyn_cl(t,x,theta,k1,k2,g,gp), [0 Tfinal], x0);
    u = control_law(x,theta,k1,k2,g,gp);

    if i ~= 3
        % Μικρές αρχικές: Figures 1 & 3
        figure(1);
        subplot(2,1,1);
        plot(t, x(:,1), 'LineWidth', 1.5, 'Color', colors(i,:));
        subplot(2,1,2);
        plot(t, x(:,2), 'LineWidth', 1.5, 'Color', colors(i,:));

        figure(3);
        plot(t, u, 'LineWidth', 1.5, 'Color', colors(i,:));
    else
        % Μεγάλη αρχική (2,20): Figures 2 & 4
        figure(2);
        subplot(2,1,1);
        plot(t, x(:,1), 'LineWidth', 1.5, 'Color', colors(i,:));
        subplot(2,1,2);
        plot(t, x(:,2), 'LineWidth', 1.5, 'Color', colors(i,:));

        figure(4);
        plot(t, u, 'LineWidth', 1.5, 'Color', colors(i,:));
    end
end

figure(1);
subplot(2,1,1);
legend('(0.5,0)','(-2,20)','(-1,10)','Location','best');
subplot(2,1,2);
legend('(0.5,0)','(-2,20)','(-1,10)','Location','best');

figure(2);
subplot(2,1,1);
legend('(2,20)','Location','best');
subplot(2,1,2);
legend('(2,20)','Location','best');

figure(3);
legend('(0.5,0)','(-2,20)','(-1,10)','Location','best');
ylim([-50 50]);

figure(4);
legend('(2,20)','Location','best');

%% ----- Τοπικές συναρτήσεις -----

function dx = dyn_cl(~,x,theta,k1,k2,g,gp)
    x1 = x(1);
    x2 = x(2);
    gx  = g(x1);
    gpx = gp(x1);

    % z2 = x2 + theta*g(x1)
    z2 = x2 + theta*gx;

    % Σωστός νόμος ελέγχου:
    % u = -k1*z1 - k2*z2 - theta*g'(x1)*z2
    u = -k1*x1 - k2*z2 - theta*gpx*z2;

    dx1 = x2 + theta*gx;
    dx2 = u;
    dx  = [dx1; dx2];
end

function u = control_law(x,theta,k1,k2,g,gp)
    x1 = x(:,1);
    x2 = x(:,2);
    gx  = g(x1);
    gpx = gp(x1);
    z2  = x2 + theta*gx;
    u   = -k1.*x1 - k2.*z2 - theta.*gpx.*z2;
end
