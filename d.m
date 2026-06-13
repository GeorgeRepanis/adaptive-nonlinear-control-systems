%% Ερώτημα (δ) - Προσαρμοστική γραμμικοποίηση μέσω ανάδρασης
% Σύστημα:
%   x1_dot = x2 + θ* g(x1)
%   x2_dot = u
% με g(x1) = x1^3 και θ* = 2

clear; clc; close all;

%% Πραγματική παράμετρος & συνάρτηση g
theta_star = 2;

g  = @(x1) x1.^3;
gp = @(x1) 3*x1.^2;

%% Κέρδη επιθυμητής 2ης τάξης δυναμικής για x1
k1 = 2;          % >0
k2 = 3;          % >0

% Πίνακες A_c, B για την κανονική μορφή χ = [x1; x1_dot]
Ac = [0 1;
     -k1 -k2];
B  = [0; 1];

% Λύση εξίσωσης Lyapunov: A_c^T P + P A_c = -Q
Q = eye(2);
P = lyap(Ac', Q);     % P συμμετρικός, θετικά ορισμένος

% Πίνακας κερδών προσαρμογής Γ (2x2)
Gamma = diag([2 2]);

%% Πακετάρισμα παραμέτρων σε struct
params.theta_star = theta_star;
params.k1 = k1;
params.k2 = k2;
params.P  = P;
params.Gamma = Gamma;
params.B  = B;
params.g  = g;
params.gp = gp;

%% Αρχικές συνθήκες από το (ς)
x0_list = [ 0.5   0;
           -2    20;
            2    20;
           -1    10];

Tfinal = 10;  % διάρκεια προσομοίωσης (s)

%% Προσομοιώσεις για όλες τις αρχικές συνθήκες
for i = 1:size(x0_list,1)

    x10 = x0_list(i,1);
    x20 = x0_list(i,2);

    % Αρχικές τιμές εκτιμήσεων παραμέτρων (μπορείς να τις αλλάξεις)
    theta_hat0 = [0; 0];

    % Συνολικό αρχικό διάνυσμα καταστάσεων: [x1; x2; theta_hat1; theta_hat2]
    z0 = [x10; x20; theta_hat0];

    % Ορισμός συνάρτησης δυναμικών κλειστού βρόχου
    ode_fun = @(t,z) closedloop_delta(t, z, params);

    % Προσομοίωση
    [t, z] = ode45(ode_fun, [0 Tfinal], z0);

    x1 = z(:,1);
    x2 = z(:,2);
    th1_hat = z(:,3);
    th2_hat = z(:,4);

    % Υπολογισμός σήματος ελέγχου u(t) για τα γραφήματα
    u = zeros(size(t));
    for k = 1:length(t)
        x1k = x1(k); x2k = x2(k);
        gk  = g(x1k);
        gpk = gp(x1k);

        % Διανυσμα φ(x)
        phi1 = k2*gk + gpk*x2k;
        phi2 = gpk*gk;
        phi  = [phi1; phi2];

        theta_hat_k = [th1_hat(k); th2_hat(k)];

        % Νόμος ελέγχου:
        % u = -k1 x1 - k2 x2 - θ̂^T φ(x)
        u(k) = -k1*x1k - k2*x2k - theta_hat_k.'*phi;
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
    title(sprintf('(δ) Προσαρμοστική γραμμικοποίηση: x_1(0)=%.1f, x_2(0)=%.1f', x10, x20));

    % Εκτιμήσεις παραμέτρων
    subplot(3,1,2);
    plot(t, th1_hat, 'LineWidth', 1.5); hold on;
    plot(t, th2_hat, 'LineWidth', 1.5);
    yline(theta_star,'--','\theta^*');
    yline(theta_star^2,'--','(\theta^*)^2');
    grid on;
    xlabel('t (s)');
    ylabel('hat {\theta}_1, hat {\theta}_2');
    legend('hat {\theta}_1','hat{\theta}_2','\theta^*','(\theta^*)^2');

    % Σήμα ελέγχου
    subplot(3,1,3);
    plot(t, u, 'LineWidth', 1.5);
    grid on;
    xlabel('t (s)');
    ylabel('u(t)');
    title('Σήμα ελέγχου');

end

%% ====== Τοπική συνάρτηση δυναμικών κλειστού βρόχου ======
function dz = closedloop_delta(~, z, params)
% z = [x1; x2; theta_hat1; theta_hat2]

x1 = z(1);
x2 = z(2);
theta_hat = z(3:4);

theta_star = params.theta_star;
k1   = params.k1;
k2   = params.k2;
P    = params.P;
Gamma = params.Gamma;
B    = params.B;
g    = params.g;
gp   = params.gp;

% Πραγματική δυναμική συστήματος
g_val  = g(x1);
gp_val = gp(x1);

% \dot x1 = x2 + θ* g(x1)
x1_dot = x2 + theta_star * g_val;

% Διανυσμα φ(x)
phi1 = k2*g_val + gp_val*x2;
phi2 = gp_val*g_val;
phi  = [phi1; phi2];

% Νόμος ελέγχου u(t) = -k1 x1 - k2 x2 - θ̂^T φ(x)
u = -k1*x1 - k2*x2 - theta_hat.'*phi;

% \dot x2 = u
x2_dot = u;

% Κανονική μορφή χ = [x1; xdot1]
chi = [x1; x1_dot];

% Scalar s = B^T P χ
s = (B.' * P * chi);

% Νόμος προσαρμογής: \dot{θ̂} = Γ φ(x) (B^T P χ)
theta_hat_dot = Gamma * phi * s;

% Τελικό διάνυσμα παραγώγων
dz = [x1_dot;
      x2_dot;
      theta_hat_dot(:)];

end
