function val = numcrossings(x)
%NUMCROSSINGS Calculate number of crossings (contacts) between fibers.
%
%   Example:
%{
      randnet = forming([21 21 30 40], [1 2 3 4], 200, 200, 'nfib', 1500);
      number_crossings_per_fiber = numcrossings(randnet);
      ksdensity(number_crossings_per_fiber)
%}

%   Copyright 2010, University of Coimbra
%   Available under the GPL-3


[~, n, p] = size(x.top);

% Insure that fibers at the top and bottom paper surfaces are not
% identified as adjacent.
top = cat( 1, x.top, false(1, n, p) );
bottom = cat( 1, x.bottom, false(1, n, p) );

% Locate crossing voxels
bottom = bottom & circshift(top, -1); % bottom surface / fiber on top
bottom(end, :, :) = [];
top = circshift(bottom, 1); % top surface / underneath fiber

% Make sure that multiple crossing voxels between the same two fibers are
% counted as only one point.
crossings = unique([x.web(bottom(:)) x.web(top(:))], 'rows');

% Count the number of crossings for each fiber.
val = arrayfun( @(x) sum( crossings(:) == x ), 1:x.nfib );

end