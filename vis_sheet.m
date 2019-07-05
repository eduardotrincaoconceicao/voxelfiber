function vis_sheet(x, ratio)
%VIS_SHEET Volume visualization of simulated sheet.
%
%   Example:
%{
      randnet = forming([21 21 30 40], [1 2 3 4], 200, 200, 'nfib', 1500);
      vis_sheet(randnet, [1 4 1])
%}

%   Copyright 2008, University of Coimbra
%   Copyright 2018, Eduardo L. T. Conceicao
%   Available under the GPL-3


% Check for a legal number of input arguments
narginchk(2, 2)
% Check input parameters
validateattributes(ratio, {'numeric'}, ...
    {'size', [1 3], 'integer', 'positive'}, mfilename, 'ratio')

figure('Units', 'normalized', ...
       'OuterPosition', [0 0 1/3 1/3], ...
       'Name', 'Default three-dimensional view', ...
       'NumberTitle', 'off', ...
       'Toolbar', 'none');
cameratoolbar
view(3)
p1 = patch(isosurface(x.web, 0), 'FaceColor', 'blue', 'EdgeColor', 'none');
reducepatch(p1, 0.5)
isonormals(x.web, p1)
daspect(ratio)
axis tight vis3d
material dull
camlight left; lighting flat

figure('Units', 'normalized', ...
       'OuterPosition', [1/3 0 1/3 1/3], ...
       'Name', 'Side view', ...
       'NumberTitle', 'off', ...
       'Toolbar', 'none');
cameratoolbar
view([0 -90])
p2 = copyobj(p1, gca);
isonormals(x.web, p2)
daspect(ratio)
axis tight vis3d
material dull
camlight left; lighting flat

figure('Units', 'normalized', ...
       'OuterPosition', [0 0.39 1/3 1/3], ...
       'Name', 'Top view', ...
       'NumberTitle', 'off', ...
       'Toolbar', 'none');
cameratoolbar
view([0 0])
p2 = copyobj(p1, gca);
isonormals(x.web, p2)
daspect(ratio)
axis tight vis3d
material dull
camlight right; lighting flat

end