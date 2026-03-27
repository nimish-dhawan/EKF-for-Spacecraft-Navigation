% Nimish Dhawan
% March 2, 2026
% EXTENDED KALMAN FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y   [6x1]: GPS measurements
% x0  [6x1]: Initial conditions
% P0  [6x6]: Initial covariance matrix
% muE [1x1]: Standard gravitational parameter
% wE  [1x1]: Earth's rotation
% Qk  [6x6]: Process noise
% t   [1x1]: Time elapsed since epoch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% est      struct : Structure storing EKF outputs
% est.x_est [6x1] : Filter state estimates

function [est] = ekf(y,x0,P0,muE,wE,Qk,t)

% Initializing state and covariance
persistent x 
persistent P

if isempty(x)
    x = x0;
    P = 1e-06*P0;
end

% Filter output structure
est = struct;

% Filter parameters
dt = 0.01;
Q  = Qk;
R  = 1-09*eye(6,6);

% Extracting position
r = x(1:3); r_norm = norm(r);

% Propogation -------------------------------------------------------------
% State transition matrix
Phi21 = muE*(3*(r*r')/r_norm^5 - eye(3)/r_norm^3);
Phi   = [eye(3), eye(3)*dt; 
         Phi21 , eye(3)];

% Priori estimates
x_minus = Phi*x;
P_minus = Phi*P*Phi' + Q;

% Correction --------------------------------------------------------------
% Measurement model
C_FI = rotate_z(wE*t);
wE   = [0;0;wE];
h    = [C_FI*x_minus(1:3); -skew(wE)*C_FI*x_minus(1:3) + C_FI*x_minus(4:6)];

% Residual
res = y-h;

% Measurement Jacobian
H = [C_FI, zeros(3,3); -skew(wE)*C_FI, C_FI];

% Kalman gain
W = H*P_minus*H' + R;
K = P_minus*H'*inv(W);

% Postiori estimates
x_est = x_minus + K*res;
P     = (eye(6)-K*H)*P_minus*(eye(6)-K*H)' + K*R*K';
x     = x_est;

% Storing values 
est.x_est = x_est;
est.res   = res;

end