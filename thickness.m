function val = thickness(x)
%THICKNESS Calculate thickness of simulated sheets.
%
%   Example:
%{
      randnet = forming([21 21 30 40], [1 2 3 4], 200, 200, 'nfib', 1500);
      thickness(randnet)
%}

%   References:
%     [1] Hellen, E.K.O., Ketoja, J.A., Niskanen, K.J., and Alava, M.J. (2002).
%         Diffusion through fibre networks.
%         J. Pulp Pap. Sci. 28(2), 56-62.
%
%     [2] Charry, E. Machado, Neumann, M., Lahti, J., Schennach, R.,
%         Schmidt, V. and Zojer, K. (2018).
%         Pore space extraction and characterization of sack paper using u-CT.
%         Journal of Microscopy 272(1), 35-46.

%   Copyright 2008, University of Coimbra
%   Copyright 2018, 2020, Eduardo L. T. Conceicao
%   Available under the GPL-3


% NOTE: in the matrix coordinate system origin is at the upper
% left corner. The physical z-axis is going down, not up.

Nz = size(x.web, 1);
thickness = zeros(length(x.sheet_top(:)), 1);
idx = x.sheet_top(:) < Nz+1;
thickness(idx) = x.sheet_bottom(idx) - x.sheet_top(idx) + 1;
val.effective = mean(thickness);
val.apparent = prctile(thickness, 80);
val.box = prctile(x.sheet_bottom, 95, 'all') - prctile(x.sheet_top, 5, 'all') + 1;

end