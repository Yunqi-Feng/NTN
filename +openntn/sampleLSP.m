function lsp = sampleLSP(scenario)
%SAMPLELSP Sample correlated large-scale parameters (LSPs), TR 38.901 step 4.
%
%   lsp = openntn.sampleLSP(scenario)
%
% Draws the seven log-domain LSPs in the fixed order
%   [DS, ASD, ASA, SF, K, ZSA, ZSD]
% from a multivariate normal whose 7x7 cross-correlation matrix is built from the
% TR 38.811 corr* table entries, then maps them to the linear domain. This mirrors
% LSPGenerator.__call__ in the Python reference.
%
% Implementation notes
%   * Cross-LSP correlation is applied per BS-UT link via a Cholesky factor; a
%     nearest-SPD projection guards against tables that are not positive definite.
%   * Spatial correlation *across* UTs (the second Cholesky in the reference) is not
%     applied here — single-link/independent-UT studies are the intended use. See
%     README_MATLAB.md "Scope and simplifications".
%   * ASA/ASD are capped at 104 deg and ZSA/ZSD at 52 deg per the spec.
%   * K and SF table values are given in dB, hence the /10 scaling on their log means/stds.
%
% Output struct fields: delaySpread [s], asd/asa/zsa/zsd [deg], shadowFading (linear),
% kFactor (linear), and logMean/logStd ([numBS x numUT x 7]) for downstream ray sampling.

if isempty(scenario.topology)
    error("OpenNTN:MissingTopology", "Call openntn.setTopology before sampleLSP.");
end

los = scenario.topology.los;
mu = zeros([size(los), 7]);
sigma = zeros([size(los), 7]);

mu(:,:,1) = openntn.param(scenario, "muDS", los);
mu(:,:,2) = finiteOr(openntn.param(scenario, "muASD", los), -100);
mu(:,:,3) = openntn.param(scenario, "muASA", los);
mu(:,:,4) = zeros(size(los));
mu(:,:,5) = openntn.param(scenario, "muK", los)/10;
mu(:,:,6) = openntn.param(scenario, "muZSA", los);
mu(:,:,7) = finiteOr(openntn.param(scenario, "muZSD", los), -100);

sigma(:,:,1) = openntn.param(scenario, "sigmaDS", los);
sigma(:,:,2) = openntn.param(scenario, "sigmaASD", los);
sigma(:,:,3) = openntn.param(scenario, "sigmaASA", los);
sigma(:,:,4) = openntn.param(scenario, "sigmaSF", los)/10;
sigma(:,:,5) = openntn.param(scenario, "sigmaK", los)/10;
sigma(:,:,6) = openntn.param(scenario, "sigmaZSA", los);
sigma(:,:,7) = openntn.param(scenario, "sigmaZSD", los);

z = randn(size(mu));
for ibs = 1:size(los, 1)
    for iut = 1:size(los, 2)
        corr = crossCorrelationMatrix(scenario, los(ibs,iut));
        root = chol(nearestSPD(corr), 'lower');
        z(ibs,iut,:) = reshape(root * reshape(z(ibs,iut,:), [], 1), 1, 1, []);
    end
end

logLSP = mu + sigma .* z;
linear = 10.^logLSP;

lsp = struct();
lsp.delaySpread = linear(:,:,1);
lsp.asd = min(linear(:,:,2), 104);
lsp.asa = min(linear(:,:,3), 104);
lsp.shadowFading = linear(:,:,4);
lsp.kFactor = linear(:,:,5);
lsp.zsa = min(linear(:,:,6), 52);
lsp.zsd = min(linear(:,:,7), 52);
lsp.logMean = mu;
lsp.logStd = sigma;
end

function y = finiteOr(x, replacement)
y = x;
y(~isfinite(y)) = replacement;
end

function corr = crossCorrelationMatrix(scenario, losMask)
names = ["corrASDvsDS", "corrASAvsDS", "corrASAvsSF", "corrASDvsSF", ...
         "corrDSvsSF", "corrASDvsASA", "corrASDvsK", "corrASAvsK", ...
         "corrDSvsK", "corrSFvsK", "corrZSDvsSF", "corrZSAvsSF", ...
         "corrZSDvsK", "corrZSAvsK", "corrZSDvsDS", "corrZSAvsDS", ...
         "corrZSDvsASD", "corrZSAvsASD", "corrZSDvsASA", "corrZSAvsASA", ...
         "corrZSDvsZSA"];
pairs = [1 2; 1 3; 4 3; 4 2; 4 1; 2 3; 2 5; 3 5; 1 5; 4 5; ...
         4 7; 4 6; 7 5; 6 5; 7 1; 6 1; 7 2; 6 2; 7 3; 6 3; 6 7];
corr = eye(7);
for k = 1:numel(names)
    v = openntn.param(scenario, names(k), losMask);
    corr(pairs(k,1), pairs(k,2)) = v;
    corr(pairs(k,2), pairs(k,1)) = v;
end
end

function a = nearestSPD(a)
a = (a + a')/2;
[v,d] = eig(a);
d = max(diag(d), 1e-8);
a = v*diag(d)*v';
a = (a + a')/2;
end
