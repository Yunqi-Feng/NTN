function scenario = setTopology(scenario, utLoc, bsLoc, varargin)
%SETTOPOLOGY Attach UT/BS positions, derive geometry, and sample LOS states.
%
%   scenario = openntn.setTopology(scenario, utLoc, bsLoc, Name, Value, ...)
%
% Computes everything that depends on geometry and then evaluates the path loss,
% mirroring SystemLevelScenario.set_topology in the Python reference.
%
% Required inputs
%   utLoc : [numUT x 3] user-terminal positions [m] ([x y z]).
%   bsLoc : [numBS x 3] satellite/HAPS positions [m]; the z-coordinate is the orbital
%           height and drives the curved-Earth slant range and Doppler.
%
% Name-Value options
%   "UTOrientations"/"BSOrientations" : [N x 3] array orientations [rad] (stored for
%                                       downstream use; not needed by the SISO model).
%   "UTVelocities" : [numUT x 3] velocity vectors [m/s] (user-motion Doppler).
%   "Indoor"       : [numUT x 1] logical indoor flags (forced NLOS; default all false).
%   "LOS"          : [] to sample LOS per the TR 38.811 elevation-dependent probability,
%                    or a logical scalar/[numBS x numUT] array to force the state.
%   "Environment"  : atmospheric/earth-station struct (default openntn.defaultEnvironment).
%
% Populates scenario.topology (distances, LOS angles, LOS states, ZoD offset, ...) and
% scenario.pathloss. The 3D distance uses the spherical-Earth slant range of TR 38.811
% Eq. 6.6-3 rather than a flat-Earth Euclidean distance.

p = inputParser;
addParameter(p, "UTOrientations", zeros(size(utLoc)), @isnumeric);
addParameter(p, "BSOrientations", zeros(size(bsLoc)), @isnumeric);
addParameter(p, "UTVelocities", zeros(size(utLoc)), @isnumeric);
addParameter(p, "Indoor", false(size(utLoc,1), 1), @(x)islogical(x) || isnumeric(x));
addParameter(p, "LOS", [], @(x)islogical(x) || isnumeric(x) || isempty(x));
addParameter(p, "Environment", scenario.environment, @isstruct);
parse(p, varargin{:});

utLoc = double(utLoc);
bsLoc = double(bsLoc);
numUT = size(utLoc, 1);
numBS = size(bsLoc, 1);

if size(utLoc, 2) ~= 3 || size(bsLoc, 2) ~= 3
    error("OpenNTN:InvalidTopology", "utLoc and bsLoc must have three columns [x y z].");
end

top = struct();
top.utLoc = utLoc;
top.bsLoc = bsLoc;
top.utOrientations = double(p.Results.UTOrientations);
top.bsOrientations = double(p.Results.BSOrientations);
top.utVelocities = double(p.Results.UTVelocities);
top.indoor = logical(p.Results.Indoor(:)).';

deltaX = reshape(utLoc(:,1), 1, []) - reshape(bsLoc(:,1), [], 1);
deltaY = reshape(utLoc(:,2), 1, []) - reshape(bsLoc(:,2), [], 1);
deltaZ = reshape(utLoc(:,3), 1, []) - reshape(bsLoc(:,3), [], 1);

top.distance2D = hypot(deltaX, deltaY);
top.distance3D = slantRange(scenario.elevationAngle, bsLoc(:,3), numUT);
top.losAOD = wrap360(rad2deg(atan2(deltaY, deltaX)));
top.losAOA = wrap360(top.losAOD + 180);
top.losZOD = wrap360(rad2deg(atan2(top.distance2D, deltaZ)));
top.losZOA = wrap360(top.losZOD - 180);
top.matrixUTDistance2D = squareformDistance(utLoc(:,1:2));

% Zenith-angle-of-departure offset, TR 38.901 Eq. 7.5-17 / Table 7.5-7.
% Computed exactly as in the Python reference (UrbanScenario._compute_lsp_log_mean_std):
% it is 0 for LOS links and a small NLOS correction otherwise. NOTE: the reference
% evaluates this with atan() (radians) and adds it to ZoD angles expressed in degrees;
% we replicate that behaviour verbatim so MATLAB results match the reference numerically.
% The magnitude is negligible for satellite NTN (ZSD spread is forced to 0), but the term
% is retained for parity and for future HAPS scenarios.
zodOffset = atan((35 - 3.5)./max(top.distance2D, eps)) ...
          - atan((35 - 1.5)./max(top.distance2D, eps));

angle = round(scenario.elevationAngle/10)*10;
losProbability = scenario.paramsLOS.(sprintf("LoS_p_%d", angle));
if isempty(p.Results.LOS)
    los = rand(numBS, numUT) < losProbability;
else
    los = logical(p.Results.LOS);
    if isscalar(los)
        los = repmat(los, numBS, numUT);
    end
end
top.los = los & ~repmat(top.indoor, numBS, 1);
top.losProbability = losProbability + zeros(numBS, numUT);

% ZoD offset is zero on LOS links (3GPP TR 38.901, Table 7.5-7).
zodOffset(top.los) = 0;
top.zodOffset = zodOffset;

scenario.environment = p.Results.Environment;
scenario.topology = top;
scenario.pathloss = openntn.pathloss(scenario);
end

function d3 = slantRange(elevationAngle, heights, numUT)
earthRadius = 6371000;
el = deg2rad(elevationAngle);
heights = heights(:);
d = sqrt(earthRadius^2*sin(el)^2 + heights.^2 + 2*heights*earthRadius) ...
    - earthRadius*sin(el);
d3 = repmat(d, 1, numUT);
end

function d = squareformDistance(xy)
n = size(xy, 1);
d = zeros(n, n);
for i = 1:n
    for j = 1:n
        d(i,j) = norm(xy(i,:) - xy(j,:));
    end
end
end

function y = wrap360(x)
y = mod(x, 360);
end
