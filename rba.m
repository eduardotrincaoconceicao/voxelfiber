function val = rba(x)
%RBA Calculate relative bonded area (RBA) for 3D fiber web.
%
%   Example:
%{
      randnet = forming([21 21 30 40], [1 2 3 4], 200, 200, nfib = 1500);
      rba(randnet)
%}

%   Copyright 2008, 2009, University of Coimbra
%   Available under the GPL-3


[~, n, p] = size(x.top);

% Insure that fibers at the top and bottom paper surfaces are not
% identified as adjacent.
top = cat( 1, x.top, false(1, n, p) );
bottom = cat( 1, x.bottom, false(1, n, p) );

% Count each pair of adjacent fiber surface cells.
num_cell_contact_area = sum( bottom(:) & circshift(top(:), -1) );

% Calculate RBA
num_cell_total_area = sum(top(:));
val = num_cell_contact_area/num_cell_total_area;

end