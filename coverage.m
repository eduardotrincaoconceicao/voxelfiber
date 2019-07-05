function val = coverage(x)
%COVERAGE Calculate number of fibres covering a point.
%
%   Example:
%{
      randnet = forming([21 21 30 40], [1 2 3 4], 200, 200, 'nfib', 1500);
      c = coverage(randnet);
      subplot(1,2,1), imagesc(c.point), axis image, colormap bone
      subplot(1,2,2), stem3(c.point), view(-25, 30)
%}

%   Copyright 2008, 2011, University of Coimbra
%   Available under the GPL-3


val.point = squeeze(sum(x.top, 1));
val.mean = mean(val.point(:));

end