function [scenario, utLoc, bsLoc] = generateSingleSectorTopology(scenario, numUT, varargin)
%GENERATESINGLESECTORTOPOLOGY Randomly drop UTs in a sector and attach the topology.
%
%   [scenario, utLoc, bsLoc] = openntn.generateSingleSectorTopology(scenario, numUT, Name, Value, ...)
%
% MATLAB analogue of the Python OpenNTN gen_single_sector_topology /
% drop_uts_in_sector utilities. A single satellite (BS) is placed at
% [0 0 BSHeight]; numUT ground users are dropped around the sub-satellite ground
% point so that their elevation angle matches scenario.elevationAngle, then
% scattered uniformly within +/- ISD/2 in x and y. The resulting topology is
% attached via openntn.setTopology so the returned scenario is ready for
% sampleLSP / sampleRays / generateChannel.
%
% Inputs
%   scenario : struct from openntn.createScenario.
%   numUT    : number of user terminals to drop.
%
% Name-Value options
%   "BSHeight"    : satellite height [m]. Default 600e3 (LEO, per TR 38.811 Table 4.5-1).
%   "ISD"         : inter-site distance [m]. Default depends on scenario type
%                   (dense_urban 200, urban 500, sub_urban 5000), per the reference.
%   "UTHeight"    : user height [m]. Default 1.5.
%   "MaxVelocity" : max user speed [m/s] (direction uniform in azimuth). Default 0.
%   "LOS"         : [] to sample LOS (default), or logical scalar/array to force it.
%   "Seed"        : optional RNG seed for reproducibility.
%
% Outputs
%   scenario : input scenario with .topology and .pathloss populated.
%   utLoc    : [numUT x 3] sampled user positions [m].
%   bsLoc    : [1 x 3] satellite position [m].

p = inputParser;
addParameter(p, "BSHeight", 600e3, @(x)isnumeric(x) && isscalar(x));
addParameter(p, "ISD", [], @(x)isnumeric(x));
addParameter(p, "UTHeight", 1.5, @(x)isnumeric(x) && isscalar(x));
addParameter(p, "MaxVelocity", 0.0, @(x)isnumeric(x) && isscalar(x));
addParameter(p, "LOS", [], @(x)islogical(x) || isnumeric(x) || isempty(x));
addParameter(p, "Seed", [], @(x)isempty(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});

if ~isempty(p.Results.Seed)
    rng(p.Results.Seed);
end

bsHeight = p.Results.BSHeight;

% Default inter-site distance per 3GPP TR 38.811 system-level scenarios.
isd = p.Results.ISD;
if isempty(isd)
    switch scenario.type
        case "dense_urban", isd = 200;
        case "urban",       isd = 500;
        case "sub_urban",   isd = 5000;
        otherwise,          isd = 500;
    end
end

el = deg2rad(scenario.elevationAngle);

% Slant range to the satellite for the target elevation, and the ground-range of
% the sub-satellite point that yields that elevation (flat-earth approximation used
% by the reference's drop_uts_in_sector for placing the sector centre).
actualBsUtDistance = bsHeight / sin(el);
distanceCenterToUT = sqrt(max(actualBsUtDistance^2 - bsHeight^2, 0));

% Sector-centre ground point (random split of the ground range across x and y).
xBase = rand(1) * distanceCenterToUT;
yBase = sqrt(max(actualBsUtDistance^2 - bsHeight^2 - xBase^2, 0));

% Scatter the UTs around the sector centre.
xDis = (rand(numUT, 1) - 0.5) * isd;
yDis = (rand(numUT, 1) - 0.5) * isd;
utLoc = [xBase + xDis, yBase + yDis, repmat(p.Results.UTHeight, numUT, 1)];

bsLoc = [0, 0, bsHeight];

% Random velocity vectors in the horizontal plane.
velAngle = (rand(numUT, 1) * 2 - 1) * pi;
velNorm  = rand(numUT, 1) * p.Results.MaxVelocity;
utVelocities = [velNorm .* cos(velAngle), velNorm .* sin(velAngle), zeros(numUT, 1)];

scenario = openntn.setTopology(scenario, utLoc, bsLoc, ...
    "UTVelocities", utVelocities, ...
    "Indoor", false(numUT, 1), ...
    "LOS", p.Results.LOS);
end
