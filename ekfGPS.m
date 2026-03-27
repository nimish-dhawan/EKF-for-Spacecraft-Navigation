% Nimish Dhawan
% March 2, 2026
% EXTENDED KALMAN FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates the state vector of a satellite in ECI along with
% a clock bias. Filter expects an array of psuedorange measurements of n
% number of satellites. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y   [m x 1] : Vector of pseudorange measurements for the visible satellites at the current epoch (m = nVis)
% x0  [6 x 1] : Initial position and velocity (m, m/s)
% t0  [1 x 1] : Initial epoch time (s)
% P0  [7 x 7] : Prior state covariance matrix
% muE [1 x 1] : Earth's gravitational parameter (m^3/s^2)
% wE  [1 x 1] : Earth's rotation rate (rad/s)
% Q   [7 x 7] : Process noise covariance matrix
% t   [1 x 1] : Time elapsed since initial epoch (s)
% clkGps_m [m x 1] : Satellite clock corrections (seconds) for the visible satellites at this epoch
% c0   [1 x 1] : Speed of light (m/s)
% r_sv [3 x m] : Satellite position vectors in ECEF frame (m)
% R    [m x m] : Measurement noise covariance matrix
% theta0 [1 x 1] : Greenwich sidereal angle at initial epoch (rad)
% nVis [1 x 1] : Number of visible satellites 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% est        struct : Structure storing EKF outputs
% est.x_est [7 x 1] : Posterior state estimate
% est.P     [7 x 7] : Posterior covariance
% est.res   [m x 1] : Measurement residual vector

function [est] = ekfGPS(y,x0,t0,P0,muE,wE,Q,t,...
                        clkGps_m,c0,r_sv,R,theta0,nVis)

%% Initializing state and covariance --------------------------------------
persistent x 
persistent P 
persistent t_prev 

if isempty(x)
    x = [x0;t0];
    P = P0;
    t_prev = t;
end

% Filter output structure
est = struct;

% Forcing column vector
y = y(:);

% Internal clock
dt     = t-t_prev;
t_prev = t;

%% Propogation ------------------------------------------------------------
% State transition matrix
r = x(1:3); r_norm = norm(r);
A = muE*(3*(r*r')/r_norm^5 - eye(3)/r_norm^3);

F = [zeros(3,3)  eye(3,3)    zeros(3,1);
     A           zeros(3,3)  zeros(3,1);
     zeros(1,3)  zeros(1,3)  0         ];

Fdt = F*dt;

Phi = eye(7,7) + Fdt + Fdt^2/2 + Fdt^3/6;

% Priori estimates
% x_minus = Phi*x;
x_minus = rk4(x, dt, muE);
P_minus = Phi*P*Phi' + Q;

%% Correction -------------------------------------------------------------
% Measurement model
C_FI = rotate_z(wrapTo2Pi(theta0+wE*t));

% Checking valid GPS measurements
vis = length(y);

% Allocating memory
h = zeros(nVis,1);
H = zeros(nVis,7);

% Measurement vector and Jacobian vector
for k = 1:vis
    if ~isfinite(y(k)) 
        continue;
    end
    d = C_FI*x_minus(1:3)-r_sv(:,k); rho = norm(d);
    h(k)   = rho + c0*x_minus(7) - clkGps_m(:,k);
    H(k,:) = [d.'/rho*C_FI, zeros(1,3), c0];
end

% Innovation
valid = isfinite(y);
y = y(valid);
res = y-h;

% Kalman gain calculation
S = H*P_minus*H' + R;
K = P_minus*H' * inv(S);

% Postiori estimates
x_plus = x_minus + K*res;
P_plus = (eye(7)-K*H)*P_minus*(eye(7)-K*H)' + K*R*K';

% Updating memory
P = P_plus;
x = x_plus;

%% Storing output values 
est.x_est = x(1:6);
est.bias  = x(7);
est.res   = res;
est.x_minus = x_minus;
est.P_est = P;

end

%% Alternative dynamics
function x_minus = rk4(x, dt, muE)
    r0 = x(1:3);
    v0 = x(4:6);
    b0 = x(7);

    f = @(r,v) (-muE * r / norm(r)^3);  

    k1_r = v0;
    k1_v = f(r0, v0);

    k2_r = v0 + 0.5*dt*k1_v;
    k2_v = f(r0 + 0.5*dt*k1_r, v0 + 0.5*dt*k1_v);

    k3_r = v0 + 0.5*dt*k2_v;
    k3_v = f(r0 + 0.5*dt*k2_r, v0 + 0.5*dt*k2_v);

    k4_r = v0 + dt*k3_v;
    k4_v = f(r0 + dt*k3_r, v0 + dt*k3_v);

    r1 = r0 + (dt/6)*(k1_r + 2*k2_r + 2*k3_r + k4_r);
    v1 = v0 + (dt/6)*(k1_v + 2*k2_v + 2*k3_v + k4_v);

    x_minus = [r1; v1; b0];
end