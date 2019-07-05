function val = thickness(x)
%THICKNESS Calculate thickness of simulated sheets.
%
%   Example:
%{
      randnet = forming([21 21 30 40], [1 2 3 4], 200, 200, 'nfib', 1500);
      thickness(randnet)
%}

%   Reference:
%     Hellen, E.K.O., Ketoja, J.A., Niskanen, K.J., and Alava, M.J. (2002).
%     Diffusion through fibre networks.
%     J. Pulp Pap. Sci. 28(2), 56-62.

%   Copyright 2008, University of Coimbra
%   Copyright 2018, Eduardo L. T. Conceicao
%   Available under the GPL-3


Nz = size(x.web, 1);
thickness = zeros(length(x.sheet_top(:)), 1);
idx = x.sheet_top(:) < Nz+1;
thickness(idx) = x.sheet_bottom(idx) - x.sheet_top(idx) + 1;
val.effective = mean(thickness);
val.apparent = prctile(thickness, 80);

end