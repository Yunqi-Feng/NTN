function rays = sampleRays(scenario, lsp)
%SAMPLERAYS Sample small-scale cluster/ray parameters, TR 38.901 steps 5-9.
%
%   rays = openntn.sampleRays(scenario, lsp)
%
% Given a topology-bearing scenario and one realization of large-scale parameters
% (from openntn.sampleLSP), this draws the geometry-based stochastic small-scale
% parameters for every BS-UT link, following 3GPP TR 38.901 Section 7.5 as
% parametrized for NTN by TR 38.811 Section 6. It mirrors RaysGenerator in the
% Python OpenNTN reference (rays.py).
%
% The per-link procedure covers:
%   * step 5 - cluster delays (exponential draw, scaled by the delay spread, with an
%              extra LOS compression based on the Rician K-factor);
%   * step 6 - cluster powers (exponential power-delay profile + log-normal per-cluster
%              shadowing, normalized to sum to one);
%   * step 7 - per-ray azimuth (AoA/AoD) and zenith (ZoA/ZoD) angles, via the
%              clusterAngles and zenithAngles helpers below;
%   * step 9 - cross-polarization power ratios (XPR), log-normal per ray.
% (Step 8, random coupling of rays within a cluster, is a reshuffle that does not
% change the per-cluster statistics and is not applied in this compact port.)
%
% Inputs
%   scenario : struct from openntn.setTopology (carries topology, LOS states, tables).
%   lsp      : struct from openntn.sampleLSP (delaySpread, asa/asd/zsa/zsd, kFactor, logMean).
%
% Output (rays struct), each angle field sized [numBS x numUT x numClusters x numRays]:
%   delays : [numBS x numUT x numClusters] cluster delays [s], sorted, first = 0.
%   powers : [numBS x numUT x numClusters] normalized cluster powers (sum to 1).
%   aoa/aod/zoa/zod : per-ray arrival/departure azimuth/zenith angles [rad].
%   xpr    : per-ray cross-polarization power ratios (linear).

% LOS state per BS-UT link, and the link-count dimensions.
los = scenario.topology.los;
numBS = size(los, 1);
numUT = size(los, 2);

% Number of clusters is state-dependent (TR 38.811 tables). The tensors are sized to
% the larger of the LOS/NLOS counts; per link, surplus clusters are masked out below.
nLOS = openntn.param(scenario, "numClusters", true);
nNLOS = openntn.param(scenario, "numClusters", false);
numClusters = max(nLOS, nNLOS);
numRays = 20;

% Intra-cluster ray angular offsets, TR 38.901 Table 7.5-3 (degrees, scaled by the
% per-cluster spread later). NOTE: the 16th entry is kept at -0.1481 to match the
% Python/Sionna reference; the strict 3GPP table value is -1.1481. See README_MATLAB.md.
rayOffsets = [0.0447 -0.0447 0.1413 -0.1413 0.2492 -0.2492 0.3715 -0.3715 ...
              0.5129 -0.5129 0.6797 -0.6797 0.8844 -0.8844 1.1481 -0.1481 ...
              1.5195 -1.5195 2.1551 -2.1551];

% Preallocate the per-link, per-cluster (and per-ray, for angles/XPR) outputs.
delays = zeros(numBS, numUT, numClusters);
powers = zeros(numBS, numUT, numClusters);
aoa = zeros(numBS, numUT, numClusters, numRays);
aod = zeros(numBS, numUT, numClusters, numRays);
zoa = zeros(numBS, numUT, numClusters, numRays);
zod = zeros(numBS, numUT, numClusters, numRays);
xpr = zeros(numBS, numUT, numClusters, numRays);

% Each BS-UT link is sampled independently from its own LOS/NLOS parameter set.
for ibs = 1:numBS
    for iut = 1:numUT
        isLOS = los(ibs,iut);

        % Active-cluster count for this link state; the mask zeros surplus clusters
        % (mask=1 marks clusters that should not exist for this link).
        n = nNLOS;
        if isLOS
            n = nLOS;
        end
        mask = [zeros(1,n), ones(1,numClusters-n)];

        % Per-link parameters: delay-scaling rTau, per-cluster shadowing std zeta [dB],
        % delay spread ds [s], and Rician K-factor k (linear).
        rTau = openntn.param(scenario, "rTau", isLOS);
        zeta = openntn.param(scenario, "zeta", isLOS);
        ds = lsp.delaySpread(ibs,iut);
        k = lsp.kFactor(ibs,iut);

        % --- Step 5: cluster delays (TR 38.901 Eq. 7.5-1..4) ---
        % Exponential draw scaled by rTau*ds; masked clusters are pushed to a huge
        % delay (=1 s) so that, after subtracting the minimum and sorting, they land
        % last and their (zeroed) power does not perturb the active clusters.
        u = max(rand(1,numClusters), 1e-6);
        unscaled = -rTau * ds .* log(u);
        unscaled = unscaled .* (1-mask) + mask;
        unscaled = sort(unscaled - min(unscaled));
        % LOS links compress the delays by a K-factor-dependent factor (Eq. 7.5-3).
        kdb = 10*log10(max(k, eps));
        scale = 0.7705 - 0.0433*kdb + 0.0002*kdb^2 + 0.000017*kdb^3;
        if isLOS
            delays(ibs,iut,:) = unscaled ./ scale;
        else
            delays(ibs,iut,:) = unscaled;
        end

        % --- Step 6: cluster powers (TR 38.901 Eq. 7.5-5..6) ---
        % Exponential power-delay profile times a log-normal per-cluster shadowing
        % term (z in dB), then mask surplus clusters and normalize to unit sum.
        z = zeta .* randn(1,numClusters);
        pwr = exp(-unscaled .* (rTau - 1) ./ (rTau * max(ds, eps))) .* 10.^(-z/10);
        pwr = pwr .* (1-mask);
        pwr = pwr ./ max(sum(pwr), eps);
        powers(ibs,iut,:) = pwr;

        % Power weights used for angle generation. On LOS links the specular component
        % is folded in via the K-factor so the first cluster dominates (Eq. 7.5-8).
        pAngles = pwr;
        if isLOS
            pAngles = pwr/(k+1);
            pAngles(1) = pAngles(1) + k/(k+1);
        end

        % --- Step 7: per-ray azimuth and zenith angles ---
        % AoA/AoD use the inverse-Gaussian mapping (clusterAngles); ZoA/ZoD use the
        % inverse-Laplacian mapping (zenithAngles). The isDeparture flag selects the
        % departure-side LOS angle/cluster spread; ZoD additionally takes the ZoD offset.
        aoa(ibs,iut,:,:) = clusterAngles(scenario, lsp.asa(ibs,iut), k, pAngles, ...
            scenario.topology.losAOA(ibs,iut), openntn.param(scenario, "cASA", isLOS), ...
            openntn.param(scenario, "CPhiNLoS", isLOS), rayOffsets, false, isLOS);
        aod(ibs,iut,:,:) = clusterAngles(scenario, lsp.asd(ibs,iut), k, pAngles, ...
            scenario.topology.losAOD(ibs,iut), openntn.param(scenario, "cASD", isLOS), ...
            openntn.param(scenario, "CPhiNLoS", isLOS), rayOffsets, true, isLOS);
        zoa(ibs,iut,:,:) = zenithAngles(scenario, lsp.zsa(ibs,iut), k, pAngles, ...
            scenario.topology.losZOA(ibs,iut), openntn.param(scenario, "cZSA", isLOS), ...
            openntn.param(scenario, "CThetaNLoS", isLOS), rayOffsets, false, isLOS, 0);
        zod(ibs,iut,:,:) = zenithAngles(scenario, lsp.zsd(ibs,iut), k, pAngles, ...
            scenario.topology.losZOD(ibs,iut), max((3/8)*10^lsp.logMean(ibs,iut,7), 0), ...
            openntn.param(scenario, "CThetaNLoS", isLOS), rayOffsets, true, isLOS, ...
            scenario.topology.zodOffset(ibs,iut));

        % --- Step 9: cross-polarization power ratios (TR 38.901 Eq. 7.5-21) ---
        % XPR is log-normal: drawn in dB as N(muXPR, sigmaXPR^2), then to linear.
        muXPR = openntn.param(scenario, "muXPR", isLOS);
        sigmaXPR = openntn.param(scenario, "sigmaXPR", isLOS);
        xpr(ibs,iut,:,:) = 10.^((muXPR + sigmaXPR.*randn(numClusters,numRays))/10);
    end
end

% Assemble the output. Angles are converted from degrees (used internally and in the
% 3GPP tables) to radians, matching the convention of the Python reference's Rays.
rays = struct();
rays.delays = delays;
rays.powers = powers;
rays.aoa = deg2rad(aoa);
rays.aod = deg2rad(aod);
rays.zoa = deg2rad(zoa);
rays.zod = deg2rad(zod);
rays.xpr = xpr;
end

function angles = clusterAngles(scenario, spread, k, pwr, losAngle, clusterSpread, cNLOS, rayOffsets, isDeparture, isLOS)
%CLUSTERANGLES Sample per-ray azimuth angles (AoA/AoD), TR 38.901 step 7 (Eq. 7.5-9..13).
%
% The cluster centres follow an inverse-Gaussian mapping of the normalized cluster
% powers; rays are then spread around each centre by the fixed offsets of Table 7.5-3.
% Mirrors RaysGenerator._azimuth_angles in the Python OpenNTN reference, including the
% LOS re-centring of the first cluster onto the geometric LOS angle.
%
% Inputs
%   spread        : azimuth angle spread (ASA for arrival, ASD for departure) [deg].
%   k             : Rician K-factor (linear); shapes the LOS scaling constant.
%   pwr           : per-cluster power weights used for the angle mapping.
%   losAngle      : geometric LOS azimuth for this link (los_aoa or los_aod) [deg].
%   clusterSpread : intra-cluster angular spread (cASA/cASD) [deg].
%   cNLOS         : NLOS C-phi constant from the tables (CPhiNLoS).
%   isDeparture   : true for AoD (departure), false for AoA (arrival).
%   isLOS         : LOS state of the link.

% Satellite departures have no azimuth spread at the transmitter (point source);
% the bs_loc>=160 km guard matches the reference.
if isDeparture && scenario.topology.bsLoc(1,3) >= 160000
    spread = 0;
end

% C-phi scaling constant: NLOS value from the tables, with a K-factor-dependent
% correction on LOS links (Eq. 7.5-10).
kdb = 10*log10(max(k, eps));
cLOS = cNLOS*(1.1035 - 0.028*kdb - 0.002*kdb^2 + 0.0001*kdb^3);
if isLOS
    c = cLOS;
else
    c = cNLOS;
end

% Inverse-Gaussian angle magnitude per cluster (Eq. 7.5-9); the normalized power
% ratio is floored at 1e-6 before the log, as in the reference.
z = max(pwr./max(pwr), 1e-6);
prime = (2*spread/1.4) .* sqrt(-log(z)) ./ c;

% Random +/-1 sign and a small Gaussian perturbation per cluster (Eq. 7.5-11), then
% offset by the geometric LOS angle.
sgn = 2*randi([0 1], size(prime)) - 1;
randomComp = (spread/7) .* randn(size(prime));
center = sgn.*prime + randomComp + losAngle;
% On LOS links, re-centre so the first (strongest) cluster lands exactly on losAngle.
if isLOS
    center = center - (sgn(1)*prime(1) + randomComp(1));
end

% Expand cluster centres into rays and wrap to the (-180,180) azimuth range.
angles = center(:) + clusterSpread .* rayOffsets(:).';
angles = mod(angles, 360);
angles(angles > 180) = angles(angles > 180) - 360;
end

function angles = zenithAngles(scenario, spread, k, pwr, losAngle, clusterSpread, cNLOS, rayOffsets, isDeparture, isLOS, zodOffset)
%ZENITHANGLES Sample per-ray zenith angles (ZoA/ZoD), TR 38.901 step 7 (Eq. 7.5-14..20).
%
% The cluster centres follow an inverse-Laplacian mapping of the normalized cluster
% powers; rays are then spread around each centre by the fixed offsets of Table 7.5-3.
% Mirrors RaysGenerator._zenith_angles in the Python OpenNTN reference, including the
% LOS re-centring of the first cluster onto the geometric LOS angle.
%
%   zodOffset : ZoD offset (Eq. 7.5-17), applied only to NLOS departure angles. Pass 0
%               for arrival angles, where it is unused.

% Satellite departures have no angular spread at the transmitter (the satellite is a
% point source seen from the ground), matching the reference bs_loc>=160 km guard.
if isDeparture && scenario.topology.bsLoc(1,3) >= 160000
    spread = 0;
end

kdb = 10*log10(max(k, eps));
cLOS = cNLOS*(1.3086 + 0.0339*kdb - 0.0077*kdb^2 + 0.0002*kdb^3);
if isLOS
    c = cLOS;
else
    c = cNLOS;
end

% Inverse-Laplacian angle magnitude per cluster (Eq. 7.5-16). Clip the normalized
% power ratio to [1e-6, 1] exactly as the reference does before the log.
z = min(max(pwr./max(pwr), 1e-6), 1.0);
prime = -spread .* log(z) ./ c;

sgn = 2*randi([0 1], size(prime)) - 1;       % random +/-1 sign per cluster
randomComp = (spread/7) .* randn(size(prime)); % N(0,(spread/7)^2) perturbation
base = sgn.*prime + randomComp;

% Cluster-centre offset depends on the link state (TR 38.901 Eq. 7.5-18/7.5-19):
%   LOS  -> re-centre so the first (strongest) cluster lands exactly on the LOS angle.
%   NLOS -> add the geometric LOS angle, plus the ZoD offset for departure angles only.
if isLOS
    additional = -(sgn(1)*prime(1) + randomComp(1) - losAngle);
elseif isDeparture
    additional = losAngle + zodOffset;
else
    additional = losAngle;
end
center = base + additional;

% Expand cluster centres into rays and wrap to the (0,180) zenith range.
angles = center(:) + clusterSpread .* rayOffsets(:).';
angles = mod(angles, 360);
angles(angles > 180) = 360 - angles(angles > 180);
end
