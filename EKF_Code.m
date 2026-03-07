% ===========================
% EKF: Target Satellite OD
% ===========================

%% Setup
clc;
clear all;

% Physical constants
c0   = 299792458.0;          % Speed of light [m/s]
fL1  = 1575.42e6;            % L1 frequency [Hz] 
fL2  = 1227.60e6;            % L2 frequency [Hz] 
lamL1 = c0/fL1;              % L1 wavelength [m]
lamL2 = c0/fL2;              % L2 wavelength [m]

muE   = 3.986004419e14;      % Earth gravitational parameter [m^3/s^2]
gmE   = 3.986005e14;         % Alternative GM value [m^3/s^2] 
wE    = 7292115.1467e-11;    % Earth rotation rate [rad/s] (for Earth rotation correction if needed)

% Time parameters
dtSec = 10;                  % Time step between measurement epochs [s]
nEpoch = 11;                 % Total number of measurement epochs
tEpoch = 0:dtSec:(nEpoch-1)*dtSec;  % Time vector for each epoch [s], t=0 at first epoch

% Interpolation parameters for GPS satellite positions/clocks from SP3 file
interpOrdPos = 11;           % Polynomial order for position interpolation
interpOrdClk = 3;            % Polynomial order for clock interpolation
tSp3 = 60 * (-78:15:73);     % SP3 time grid relative to t=0 [seconds] (15-min intervals)

% Earth parameters
rEarth_km = 6371;            % Mean Earth radius [km] (for altitude calculation)

% Receiver clock bias for each epoch [seconds]
% These values represent the target satellite receiver clock offset from GPS time
% Estimated from prior processing or single-point positioning
del_t = [0.0031, -0.0011, -0.0012, -0.0012, -0.0012, -0.0012, -0.0012,  -0.0012, -0.0012, -0.0012, -0.0012];  % [s] for 11 epochs


fprintf('========================================\n');
fprintf('Running EKF for Target Orbit Determination\n');
fprintf('========================================\n\n');

%% GPS SP3 load
fprintf('Reading GPS satellite data from SP3 file...\n');
[xSp3_km, ySp3_km, zSp3_km, clkSp3_s] = get_clock('GPS_Sat_Readings.sp3');
fprintf('  --> GPS data loaded successfully.\n\n');

%% Target measurements (pseudorange) + interpolate GPS sats to epoch
fprintf('Reading Target GPS observations...\n');

prnMat = nan(nEpoch,8);
rhoMat = nan(nEpoch,8);

data = rinexread('TargetSat_Orbit_Measurement.dat');
tt = data.GPS;

% Group by receiver epoch time
[G, trGps] = findgroups(tt.Time);
nEpoch = numel(trGps);  % Forcing measured Epochs

% Receiver time in seconds since first epoch 
t0_dt = trGps(1);
tRx_s = seconds(trGps - t0_dt);

JD = juliandate(t0_dt);                
Du = JD - 2451545.0;           

% IAU 2000 ERA Formula; DOI: 10.1051/0004-6361:20030817
revs = 0.7790572732640 + 1.00273781191135448*Du; 

theta0 = mod(revs, 1) * 2*pi;    % radians in [0,2pi)

% Allocate memory
xGps_m   = nan(nEpoch, 8);
yGps_m   = nan(nEpoch, 8);
zGps_m   = nan(nEpoch, 8);
clkGps_s = nan(nEpoch, 8);   
rhoGps_m = nan(nEpoch, 8);   
nVis     = zeros(nEpoch,1);

% Interpolatingn GPS
for e = 1:nEpoch
    idx = (G == e);

    prn_list = tt.SatelliteID(idx);
    rho_list = tt.C1C(idx);

    nVis(e) = numel(prn_list);
    prnMat(e,1:nVis(e)) = prn_list';
    rhoMat(e,1:nVis(e)) = rho_list';

    for k = 1:nVis(e)
        prn = prn_list(k);
        rho = rho_list(k);

        % Receiver time for this epoch
        t_rx = tRx_s(e);

        % Transmit time (seconds since first receiver epoch)
        t_tx = t_rx - rho/c0;

        % Interpolate SP3 
        dt_sv = interpolate_gps(tSp3, clkSp3_s(:,prn), t_tx, 3);          % [s]
        x_sv  = interpolate_gps(tSp3, xSp3_km(:,prn)*1e3, t_tx, 11);      % [m]
        y_sv  = interpolate_gps(tSp3, ySp3_km(:,prn)*1e3, t_tx, 11);      % [m]
        z_sv  = interpolate_gps(tSp3, zSp3_km(:,prn)*1e3, t_tx, 11);      % [m]

        % Storing in arrays
        rhoGps_m(e, k)   = rho;
        clkGps_s(e, k)   = dt_sv;
        xGps_m(e, k)     = x_sv;
        yGps_m(e, k)     = y_sv;
        zGps_m(e, k)     = z_sv;
    end
end

clkGps_m = clkGps_s * c0;

%% Load precise target orbit (truth data)
fprintf('Reading precise Target Sat orbit...\n');
fidTruth = fopen('TargetSat_Orbit_Truth.dat','r');

for k = 1:22
    fgets(fidTruth);   
end

xTruth_km = zeros(1,11);
yTruth_km = zeros(1,11);
zTruth_km = zeros(1,11);

for e = 1:11
    fgets(fidTruth);
    ln = fgets(fidTruth);
    xTruth_km(e) = read_string(ln(5:18));
    yTruth_km(e) = read_string(ln(19:32));
    zTruth_km(e) = read_string(ln(33:46));

    fgets(fidTruth);
end
fclose(fidTruth);

rTruth_km = sqrt(xTruth_km.^2 + yTruth_km.^2 + zTruth_km.^2);
altTruth_km = rTruth_km - rEarth_km;
altTruth_m = altTruth_km*1e03;

fprintf('  --> Precise orbit loaded successfully.\n\n');

%% EKF init 
xg  = xSp3_km; yg = ySp3_km; zg = zSp3_km; cc = clkSp3_s;
xgi = xGps_m;  ygi = yGps_m; zgi = zGps_m; cci = clkGps_s; ccf = clkGps_m;

xp = xTruth_km; yp = yTruth_km; zp = zTruth_km; rp = rTruth_km;

fprintf('Initializing Extended Kalman Filter...\n');

x0 = [ 1.9393e6;
      -5.8377e6;
       2.9332e6;
       1.0010e3;
      -3.1000e3;
      -6.9000e3];

P0 = diag([10e03^2*ones(1,3), 100^2*ones(1,3), 1e-04]);       

qPos = (400)^2;
qVel = 66.74;
qBias= 1e-02^2;
Qk   = diag([qPos qPos qPos qVel qVel qVel qBias]);

sigRho = 5;

dtRx_s = 0.000123456;
P0(7,7)= dtRx_s^2;

xHist     = zeros(7, nEpoch);
xMinus    = zeros(7, nEpoch);
PHist     = zeros(7, 7, nEpoch);
rEst_km   = zeros(1, nEpoch);
altEst_km = zeros(1, nEpoch);

res  = nan(nEpoch,8);
r_sv = zeros(nEpoch, 8);
NIS  = nan(nEpoch, 8);

fprintf('  --> EKF initialized\n\n');

%% EKF loop
fprintf('Running EKF for %d epochs...\n', nEpoch);
fprintf('------------------------------------------\n');
clear ekfGPS
for e = 1:nEpoch
    fprintf('Epoch %2d/%2d: ', e, nEpoch);
    t = tRx_s(e);                   
    R = 1*sigRho^2*eye(nVis(e));                % Measurement covariance matrix
    C_FI = rotate_z(wrapTo2Pi(theta0 + wE*t));  % ECI to ECEF rotation matrix

    % Running the EKF
    est = ekfGPS(rhoGps_m(e,:),x0,del_t(1),P0,muE,wE,Qk,t,clkGps_m(e,:),...
                 c0,[xgi(e,:);ygi(e,:);zgi(e,:)],R,theta0,nVis(e));

    xHist(:,e)       = [est.x_est*1e-03; est.bias];   % Storing estimated states [m, m/s, s]
    xHist(1:3,e)     = C_FI*xHist(1:3,e);             % Rotating to ECEF [m]
    rEst_km(e)       = sqrt(xHist(1,e)^2 + xHist(2,e)^2 + xHist(3,e)^2);  % Orbit distance [km]
    altEst_km(e)     = rEst_km(e) - rEarth_km;        % Altitude [km]
    res(e,1:nVis(e)) = est.res';                      % Innovation
    xMinus(:,e)      = est.x_minus;                   % Priori estimates
    PHist(:,:,e)     = est.P_est;                     % State covariance
end

fprintf('------------------------------------------\n');
fprintf('EKF processing complete\n\n');

%% Results
fprintf('Computing results and errors...\n');
Plots
