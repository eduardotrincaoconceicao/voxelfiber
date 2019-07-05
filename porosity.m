function val = porosity(x)
%POROSITY Calculate porosity for 3D fiber web.
%
%   Example:
%{
      randnet = forming([21 21 30 40], [1 2 3 4], 200, 200, 'nfib', 1500);
      porosity(randnet)
%}

%   Copyright 2008, 2009, University of Coimbra
%   Copyright 2018, Eduardo L. T. Conceicao
%   Available under the GPL-3


Nz = size(x.web, 1);
void_above_top_surface = sum(x.sheet_top(:) - 1);
void_between_net_substrate = sum(max(Nz - x.sheet_bottom(:), 0));
val.interfiber = 1 - nnz(x.web)/( numel(x.web) - ...
                                  void_above_top_surface - ...
                                  void_between_net_substrate );
val.intrafiber = sum( x.lumen(:)/nnz(x.web) );

end