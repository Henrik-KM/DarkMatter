% EXAMPLE_MINIMAL Run small analytic examples without saved workspaces.
% From the repository root: run('examples/minimal_example.m')
% This example does not load MAT files or write output files.

exampleDirectory = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(exampleDirectory);
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(repositoryRoot);

v0 = 230e5;
vE = 240e5;
vesc = 600e5;

velocityIntegral = Eta(1000, 1, 1e7, v0, vE, vesc);
assert(isscalar(velocityIntegral) && isfinite(velocityIntegral) && velocityIntegral >= 0);

energyBins = 1:500;
chargeBins = 1:500;
siliconBackground = background(energyBins, 0.1, chargeBins, true);
assert(isequal(size(siliconBackground), [1, 500]));
assert(all(isfinite(siliconBackground)) && all(siliconBackground >= 0));

fprintf('Minimal example completed. Eta = %.12g\n', velocityIntegral);
fprintf('Silicon background: %d finite, nonnegative bins.\n', numel(siliconBackground));
clear pathCleanup;
