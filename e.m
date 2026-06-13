%% Ερώτημα (ε) - Προσαρμοστικό backstepping
% Σύστημα:
%   x1_dot = x2 + theta_star * g(x1)
%   x2_dot = u
% με g(x1) = x1^3, άγνωστη παράμετρο θ*, και adaptive backstepping ελεγκτή.

clear; clc; close all;

%% Πραγματική παράμετρος & συνάρτηση g
theta_star = 2;             % θ* στο φυτό

g  = @(x1) x1.^3;
gp = @(x1) 3*x1.^2;

%% Κέρδη backstepping
k1 = 2;      % >0
k2 = 3;      % >0
gamma = 1;   % κέρδος προσαρμογής γ > 0

%% Αρχικές συνθήκες από το (ς)
x0_list = [ 0.5   0;
           -2    20;
            2    20;
           -1    10];

Tfinal = 10;  % διάρκεια προσομοίωσης

%% Προσομοιώσεις για όλες τις αρχικές συνθήκες
for i = 1:size(x0_list,1)

    x10 = x0_list(i,1);
    x20 = x0_list(i,2);

    % Αρχική εκτίμηση παραμέτρου (μπορείς να την αλλάξεις)
    theta_hat0 = 0;

    % Κατάσταση: [x1; x2; theta_hat]
    z0 = [x10; x20; theta_hat0];

    % Προσομοίωση με ode45
    ode_fun = @(t,z) closedloop_eps(t, z, theta_star, g, gp, k1, k2, gamma);
    [t, z] = ode45(ode_fun, [0 Tfinal], z0);

    x1 = z(:,1);
    x2 = z(:,2);
    theta_hat = z(:,3);

    % Υπολογισμός σφαλμάτων z1,z2 και ελέγχου u(t) για plotting
    z1 = x1;
    u  = zeros(size(t));
    for k = 1:length(t)
        x1k = x1(k); x2k = x2(k);
        thk = theta_hat(k);

        gk  = g(x1k);
        gpk = gp(x1k);

        % Backstepping variables
        z1k = x1k;
        z2k = x2k + k1*x1k + thk*gk;
        Bk  = k1 + thk*gpk;

        % Έλεγχος (σύμφωνα με τη θεωρία του (ε)):
        u(k) = -(Bk + gamma*Bk*gk^2 + k2) * z2k ...
               -(1 - Bk*k1 + gamma*gk^2) * z1k;
    end

    %% Γραφήματα
    figure;

    % x1, x2
    subplot(3,1,1);
    plot(t, x1, 'LineWidth', 1.5); hold on;
    plot(t, x2, 'LineWidth', 1.5);
    grid on;
    xlabel('t (s)');
    ylabel('x_1, x_2');
    legend('x_1','x_2');
    title(sprintf('(ε) Adaptive backstepping: x_1(0)=%.1f, x_2(0)=%.1f', x10, x20));

    % Εκτίμηση παραμέτρου
    subplot(3,1,2);
    plot(t, theta_hat, 'LineWidth', 1.5); hold on;
    yline(theta_star,'--','\theta^*');
    grid on;
    xlabel('t (s)');
    ylabel('hat {\theta}(t)');
    legend('hat {\theta}','\theta^*');

    % Σήμα ελέγχου
    subplot(3,1,3);
    plot(t, u, 'LineWidth', 1.5);
    grid on;
    xlabel('t (s)');
    ylabel('u(t)');
    title('Σήμα ελέγχου');

end

%% ====== Τοπική συνάρτηση δυναμικών κλειστού βρόχου για (ε) ======
function dz = closedloop_eps(~, z, theta_star, g, gp, k1, k2, gamma)
% z = [x1; x2; theta_hat]

x1 = z(1);
x2 = z(2);
theta_hat = z(3);

% Φυτό
g_val  = g(x1);
gp_val = gp(x1);

x1_dot = x2 + theta_star * g_val;   % \dot x1 = x2 + θ* g(x1)

% Backstepping variables
z1 = x1;
z2 = x2 + k1*x1 + theta_hat*g_val;
B  = k1 + theta_hat*gp_val;

% Νόμος προσαρμογής:
theta_hat_dot = gamma * g_val * (z1 + B*z2);

% Νόμος ελέγχου:
u = -(B + gamma*B*g_val^2 + k2) * z2 ...
    -(1 - B*k1 + gamma*g_val^2) * z1;

% Δυναμική x2
x2_dot = u;

% Τελικό διάνυσμα παραγώγων
dz = [x1_dot;
      x2_dot;
      theta_hat_dot];
end
