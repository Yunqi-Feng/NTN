function gasAtt = gasAttenuation676(scenario)
%GASATTENUATION676 Gaseous (oxygen + water-vapour) attenuation [dB], ITU-R P.676.
%
%   gasAtt = openntn.gasAttenuation676(scenario)
%
% Implements the compact P.676 model used by TR 38.811 Section 6.6.4: it computes the
% zenith-path specific attenuations and equivalent heights for water vapour and dry
% air separately, sums them, and projects onto the slant path via 1/sin(elevation).
% Valid roughly for 1-350 GHz; mirrors utils.compute_pathloss_gas in the Python OpenNTN.
%
% Atmospheric inputs are taken from scenario.environment (temperature, pressure,
% water-vapour density); see openntn.defaultEnvironment.
%
% Output
%   gasAtt : scalar slant-path gaseous attenuation [dB], applied to all links.

% Environment and frequency.
env = scenario.environment;
T = env.temperature;              % temperature [K]
p = env.atmosphericPressure;      % surface pressure [hPa]
f = scenario.carrierFrequency/1e9;% carrier frequency [GHz]
pw = env.waterVaporDensity;       % water-vapour density [g/m^3]
rt = 288/T;                       % normalized inverse temperature
rp = p/1013;                      % normalized pressure

% --- Water-vapour specific attenuation Yw [dB/km] (P.676 line series) ---
% Sum of the principal H2O resonance lines with temperature/pressure-scaled strengths.
n1 = 0.955*rp*(rt^0.68) + 0.006*pw;
n2 = 0.735*rp*(rt^0.5) + 0.0353*(rt^4)*pw;
yw = ((3.98*n1*exp(2.23*(1-rt))) / ((f-22.235)^2 + 9.42*n1^2) * (1+((f-22)/(f+22))^2) + ...
      (11.96*n1*exp(0.7*(1-rt))) / ((f-183.31)^2 + 11.14*n1^2) + ...
      (0.081*n1*exp(6.44*(1-rt))) / ((f-321.226)^2 + 6.29*n1^2) + ...
      (3.66*n1*exp(1.6*(1-rt))) / ((f-325.153)^2 + 9.22*n1^2) + ...
      (25.37*n1*exp(1.09*(1-rt))) / ((f-380)^2) + ...
      (17.4*n1*exp(1.46*(1-rt))) / ((f-448)^2) + ...
      (844.6*n1*exp(0.17*(1-rt))) / ((f-557)^2) * (1+((f-557)/(f+557))^2) + ...
      (290*n1*exp(0.41*(1-rt))) / ((f-752)^2) * (1+((f-752)/(f+752))^2) + ...
      (83328*n2*exp(0.99*(1-rt))) / ((f-1780)^2) * (1+((f-1780)/(f+1780))^2)) * ...
      (f^2 * rt^2.5 * (pw*10^(1-5)));

% Equivalent height for water vapour hw [km], then the zenith water-vapour
% attenuation Aw = Yw*hw [dB].
conw = 1.013/(1+exp((0-8.1)*(rp-0.57)));
hw = 1.66*(1 + (1.39*conw)/((f-22.235)^2+2.56*conw) + ...
              (3.37*conw)/((f-183.31)^2+4.69*conw) + ...
              (1.5*conw)/((f-325.1)^2+2.89*conw));
aw = yw*hw;

% --- Dry-air (oxygen) specific attenuation Yo [dB/km] ---
% Empirical coefficients ee1..ee7 capture the temperature/pressure dependence used in
% the frequency-band branches below.
ee1 = rp^0.0717 * rt^(-1.8132) * exp(0.0156*(1-rp) - 1.6515*(1-rt));
ee2 = rp^0.5146 * rt^(-4.6368) * exp(-0.1921*(1-rp) - 5.7416*(1-rt));
ee3 = rp^0.3414 * rt^(-6.5851) * exp(0.2130*(1-rp) - 8.5854*(1-rt));
ee4 = rp^(-0.0112) * rt^0.0092 * exp(-0.1033*(1-rp) - 0.0009*(1-rt));
ee5 = rp^0.2705 * rt^(-2.7192) * exp(-0.3016*(1-rp) - 4.1033*(1-rt));
ee6 = rp^0.2445 * rt^(-5.9191) * exp(0.0422*(1-rp) - 8.0719*(1-rt));
ee7 = rp^(-0.1833) * rt^6.5589 * exp(-0.2402*(1-rp) + 6.131*(1-rt));

% The dry-air attenuation uses a different closed form in each frequency band around
% the 60 GHz oxygen complex; above 120 GHz it is taken as 0 (outside the NTN bands).
if f <= 54
    yo = (((7.2*rt^2.8)/(f^2 + 0.34*rp^2*rt^1.6)) + ...
          ((0.62*ee3)/((54-f)^(1.16*ee1) + 0.83*ee2))) * f^2 * rp^2 * 1e-3;
elseif f <= 66
    % ITU-R P.676 dry-air interpolation across the 60 GHz O2 complex.
    % The exponent n switches at 60 GHz (matches utils.compute_pathloss_gas).
    if f <= 60
        n = 0;
    else
        n = -15;
    end
    g54 = 2.136 * rp^1.4975 * rt^(-1.5852) * exp(-2.5196*(1-rt));
    g57 = 9.984 * rp^0.9313 * rt^2.6732 * exp(0.8563*(1-rt));
    g60 = 15.42 * rp^0.8595 * rt^3.6178 * exp(1.1521*(1-rt));
    g63 = 10.63 * rp^0.9298 * rt^2.3284 * exp(0.6287*(1-rt));
    g66 = 1.944 * rp^1.6657 * rt^(-3.3714) * exp(-4.1643*(1-rt));
    yo = exp(((54^-n)*log(g54)*(f-57)*(f-60)*(f-63)*(f-66)/1944) - ...
             ((57^-n)*log(g57)*(f-54)*(f-60)*(f-63)*(f-66)/486) + ...
             ((60^-n)*log(g60)*(f-54)*(f-57)*(f-63)*(f-66)/324) - ...
             ((63^-n)*log(g63)*(f-54)*(f-57)*(f-60)*(f-66)/486) + ...
             ((66^-n)*log(g66)*(f-54)*(f-57)*(f-60)*(f-63)/1944)) * f^n;
elseif f <= 120
    yo = 3.02e-4*rt^3.5 + (0.283*rt^3.8)/((f-118.75)^2 + 2.91*rp^2*rt^1.6) + ...
         (0.502*ee6*(1 - 0.0163*ee7*(f-66))) / ((f-66)^(1.4346*ee4 + 1.15*ee5));
else
    yo = 0;
end

% Equivalent height for dry air ho [km] (t1..t3 are its frequency correction terms),
% then the zenith dry-air attenuation Ao = Yo*ho [dB].
t1 = (4.64/(1+0.066*rp^(-2.3))) * exp(-((f-59.7)/(2.87+12.4*exp(-7.9*rp)))^2);
t2 = (0.14*exp(2.12*rp))/((f-118.75)^2 + 0.031*exp(2.2*rp));
t3 = (0.0114/(1+0.14*rp^(-2.6))) * f * ...
     ((-0.0247 + 0.0001*f + 1.61e-6*f^2)/(1 - 0.0169*f + 4.1e-5*f^2 + 3.2e-7*f^3));
ho = (6.1/(1+0.17*rp^(-1.1))) * (1+t1+t2+t3);
ao = yo*ho;

% Total zenith attenuation (dry air + water vapour) projected onto the slant path.
gasAtt = abs((ao + aw) / sind(scenario.elevationAngle));
end
