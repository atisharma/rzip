function [Gc] = choosecoils(G,coils)
% Choose which inputs will be used in control of plant

% E.g. choose all E, F coils
Gc = G(:,coils);