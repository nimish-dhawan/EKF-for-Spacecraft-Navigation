function [x_coord,y_coord,z_coord,clock_correction] = get_clock(orbfil)
% read clock correction from an IGS SP3 file
%
% INPUTS
% orbfil    file with sp3 data of GPS satellite orbits
% prn       PRN of satellite for which you want coordinates

% OUTPUTS
% x_coords  x-coordinate in meters
% y_coords  y-coordinate in meters
% z_coords  z-coordinate in meters
% clock_correction in seconds

% time is in seconds from the beginning of the day
% the matrix clk contains the clock corrections at the position
% corresponding to sats

% open file
obs = fopen(orbfil,'r');
%
% tell user if file was successfully opened or not
%
if obs ~= -1
    disp('Orbit file successfully opened')
else
    disp('Unable to open orbit file')
    return
end

for ii = 1:22
    line = fgets(obs);
end

% read the date information
line = fgets(obs);
Y = str2num(line(4:7));
M = str2num(line(10));
D = str2num(line(12:13));
H = str2num(line(16));
MI = str2num(line(19));
S = str2num(line(22:31));

% Loop over the epochs
for ii = 1:11
    while 1
        line = fgets(obs);
        prn = read_string(line(3:4));
        x_coord(ii,prn) = read_string(line(5:18));
        y_coord(ii,prn) = read_string(line(19:32));
        z_coord(ii,prn) = read_string(line(33:46));
        clock_correction(ii,prn) = 1e-6*read_string(line(47:60)); % file contains microseconds
        if prn == 32
            break
        end
    end
    line = fgets(obs);
end
fclose(obs);
