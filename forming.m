function out = forming(fib_length, flex, Nx, Ny, options)
%FORMING Simulation of 3D planar random fiber network using the KCL-PAKKA model for fiber flexibility.
%
%   Examples:
%{
      forming([21 21 30 40], [1 2 3 4], 200, 200, ...
              nfib = 1500, trace = true, delay_trace = .01);
      forming([21 21 30 40], [1 1 1 1], 100, 100, ...
              nfib = 100, horzspan = [1 2 4 6], trace = true);
      forming([21 21 30 40], [1 2 3 4], 100, 100, ...
              nfib = 100, fib_width = [1 2 3 4], trace = true);

      randnet = forming([21 21 30 40], [1 2 3 4], 200, 200, nfib = 1500);
      figure
      subplot(2,2,1), spy(randnet.surface)
      subplot(2,2,2), imagesc(randnet.surface), axis image
      subplot(2,2,[3,4]), surf(randnet.surface), daspect([1 1 4])

      tic
      randnet = forming(100, 1, 400, 400, ...
                        acceptanceprob = 0.01, alpha = 0.94, ...
                        stop_criterion = 'Grammage', grammage = 40, ...
                        mass = 'coarseness', coarseness = 0.21e-3, ...
                        delxy = 5e-6, fib_width = 6, ...
                        wall_thick = 2, lumen_thick = 0);
      toc
%}

%   References:
%     [1] Niskanen, K.J. and Alava, M.J. (1994).
%         Planar random networks with flexible fibers.
%         Phys. Rev. Lett. 73(25), 3475-3478.
%
%     [2] Niskanen, Kaarlo, Nilsen, Niko, Hellen, Erkki, and Alava, Mikko (1997).
%         KCL-PAKKA: Simulation of the 3D structure of paper.
%         In: C.F. Baker (Ed.), The Fundamentals of Papermaking Materials,
%         Proceedings of the 11th Fundamental Research Symposium,
%         Cambridge, UK, 22-26 September 1997, Pira International, pp. 1273-1292.
%
%     [3] Nilsen, Niko, Zabihian, Mari, and Niskanen, Kaarlo (1998).
%         KCL-PAKKA: a tool for simulating paper properties.
%         Tappi J. 81(5), 163-166.

%   Copyright 2008-2011, University of Coimbra
%   Copyright 2014, 2018, 2019, 2023, Eduardo L. T. Conceicao
%   Available under the GPL-3


    arguments
        fib_length
        flex
        Nx (1, 1) {mustBeNonempty, mustBeInteger, mustBePositive}
        Ny (1, 1) {mustBeNonempty, mustBeInteger, mustBePositive}
        options.acceptanceprob (1, 1) ...
            {mustBeNonempty, mustBeNumeric, ...
             mustBeInRange(options.acceptanceprob, 0, 1)} = 1
        options.alpha (1, 1) ...
            {mustBeNonempty, mustBeNumeric, ...
             mustBeInRange(options.alpha, 0, 1)} = 1
        options.stop_criterion ...
            {mustBeTextScalar, ...
             mustBeMember(options.stop_criterion, ["NumFibers", "Grammage"])} ...
            = "NumFibers"
        options.nfib = []
        options.grammage = []
        options.mass ...
            {mustBeTextScalar, ...
             mustBeMember(options.mass, ["coarseness", "wall_density"])} ...
            = "coarseness"
        options.coarseness = []
        options.wall_density = []
        options.delxy = []
        options.delz = []
        options.fib_width = repelem(1, 1, length(fib_length))
        options.wall_thick = repelem(1, 1, length(fib_length))
        options.lumen_thick = repelem(2, 1, length(fib_length))
        options.horzspan = repelem(1, 1, length(fib_length))
        options.fraction = repelem(1/length(fib_length), 1, length(fib_length))
        options.choose_species = []
        options.blocksize (1, 3) ...
            {mustBeNonempty, mustBeInteger, mustBePositive} = [10000 10 10000]
        options.trace ...
            {mustBeA(options.trace, "logical")} = false
        options.delay_trace (1, 1) ...
            {mustBeNonempty, mustBeNumeric, ...
             mustBeNonnegative, mustBeFinite} = 0.3
    end

    function draw_outofplane_grid()
        [zPos, lPos] = find(bend_plane);
        [m, n] = size(bend_plane);

        figure(hndlFig_vertical)
        clf
        plot(lPos, zPos, ...
             'LineStyle', 'none', ...
             'Marker', 'x', ...
             'Color', 'black')
        axis ij image
        grid on
        set(gca, ...
            'XLim', [0.5 (n+0.5)], ...
            'YLim', [0.5 (m+0.5)], ...
            'XTick', 0.5:(n+0.5), ...
            'YTick', 0.5:(m+0.5), ...
            'XTickLabel', '', ...
            'YTickLabel', '', ...
            'TickLength', [0 0], ...
            'GridLineStyle', '-')
        hold on
    end

    function [line_handle, currentColor] = draw_fibre_initPos()
        currentColor = col{mod( cur_species, length(col) ) + 1};
        n = size(bend_plane, 2);
        line_handle = plot(1:n, repmat(startInd, 1, n), ...
                           'LineStyle', 'none', ...
                           'Marker', 'square', ...
                           'MarkerEdgeColor', currentColor, ...
                           'MarkerFaceColor', currentColor);
        pause(delay_trace)
    end

    % NOTE: in the matrix coordinate system origin is at the upper
    % left corner. The physical z-axis is going down, not up.
    % Thus elevation is expressed as depth.

    function [bend_plane, bend_slice, bheight, net_idx] = get_outofplane_slice()
        % This function extracts an out-of-plane slice from the 3D network
        z_extent = firstVPos:lastVPos;
        bheight = length(z_extent);
        X = repmat(reshape(fib_cell_xcoord, [1 blength bwidth]), [bheight 1 1]);
        Y = repmat(reshape(fib_cell_ycoord, [1 blength bwidth]), [bheight 1 1]);
        Z = repmat(z_extent', [1 blength bwidth]);
%         net_idx = sub2ind(net_dims, Z, X, Y);
        net_idx = NzNx*(Y-1) + Nz*(X-1) + Z;
        bend_slice = net(net_idx);
        bend_plane = sum(bend_slice, 3)/bwidth;
    end

    function species = choose_species_default(fraction, N)
        % Choice of fiber species.
        species = randsample(nspecies, N, true, fraction);
    end

    function val = rule_numfibers()
        val = fibCount <= nfib;
    end

    function val = rule_grammage()
        val = cur_grammage < grammage;
    end

    function val = complete_numfibers()
        val = (double(fibCount-1))/nfib;
    end

    function val = complete_grammage()
        val = cur_grammage/grammage;
    end

    function val = fiber_areal_mass_wall_density()
        val = projected_area*cur_fib_width*cur_wall_thick*multiplier;
    end

    function val = fiber_areal_mass_coarseness()
        if ismultispecies
            cur_coarseness = coarseness(cur_species);
        else
            cur_coarseness = coarseness;
        end
        if cur_fib_width > 1
            val = cur_coarseness * numel(in_bounds(:))/cur_fib_width * multiplier;
        else
            val = cur_coarseness * numel(in_bounds(:)) * multiplier;
        end
    end


% Check for a legal number of input arguments
narginchk(6, 42)

% Extract arguments
acceptanceprob = options.acceptanceprob;
alpha = options.alpha;
stop_criterion = options.stop_criterion;
nfib = options.nfib;
grammage = options.grammage; % g/m2
mass = options.mass;
coarseness = options.coarseness; % g/m
wall_density = options.wall_density; % g/m3
delxy = options.delxy; % m
delz = options.delz; % m
fib_width = options.fib_width;
wall_thick = options.wall_thick;
lumen_thick = options.lumen_thick;
horzspan = options.horzspan;
fraction = options.fraction;
if isempty(options.choose_species)
    choose_species = @choose_species_default;
end
blocksize = options.blocksize;
% blocksize(1) = 1e4;
% if acceptanceprob >= 0.1
%     blocksize(2) = 10;
% else
%     blocksize(2) = 1;
% end
% blocksize(3) = 1e4;
trace = options.trace;
delay_trace = options.delay_trace;

% Check input parameters
validateattributes(fib_length, {'numeric'}, ...
    {'row', 'integer', 'positive'}, mfilename, 'fib_length')
nspecies = length(fib_length);
validateattributes(flex, {'numeric'}, ...
    {'size', [1 nspecies], 'integer', 'nonnegative'}, mfilename, 'flex')
switch stop_criterion
    case 'NumFibers'
        validateattributes(nfib, {'numeric'}, ...
            {'nonempty', 'scalar', 'integer', 'positive', '<=', intmax('int32')}, ...
            mfilename, 'nfib')
    case 'Grammage'
        validateattributes(grammage, {'numeric'}, ...
            {'nonempty', 'scalar', 'positive', 'finite'}, mfilename, 'grammage')
        switch mass
            case 'coarseness'
                validateattributes(coarseness, {'numeric'}, ...
                    {'nonempty', 'size', [1 nspecies], 'positive', 'finite'}, ...
                    mfilename, 'coarseness')
                validateattributes(delxy, {'numeric'}, ...
                    {'nonempty', 'scalar', 'positive', 'finite'}, ...
                    mfilename, 'delxy')
            case 'wall_density'
                validateattributes(wall_density, {'numeric'}, ...
                    {'nonempty', 'size', [1 nspecies], 'positive', 'finite'}, ...
                    mfilename, 'wall_density')
                validateattributes(delz, {'numeric'}, ...
                    {'nonempty', 'scalar', 'positive', 'finite'}, ...
                    mfilename, 'delz')
        end
end
validateattributes(fib_width, {'numeric'}, ...
    {'size', [1 nspecies], 'integer', 'positive'}, mfilename, 'fib_width')
validateattributes(wall_thick, {'numeric'}, ...
    {'size', [1 nspecies], 'integer', 'positive'}, mfilename, 'wall_thick')
validateattributes(lumen_thick, {'numeric'}, ...
    {'size', [1 nspecies], 'integer', 'nonnegative'}, mfilename, 'lumen_thick')
validateattributes(horzspan, {'numeric'}, ...
    {'size', [1 nspecies], 'integer', 'positive'}, mfilename, 'horzspan')
validateattributes(fraction, {'numeric'}, ...
    {'size', [1 nspecies], 'nonnegative', 'finite'}, mfilename, 'fraction')
validateattributes(sum(fraction), {'numeric'}, ...
    {'scalar', '>=', 1, '<=', 1}, mfilename, 'sum(fraction)')
validateattributes(choose_species, {'function_handle'}, {'nonempty'}, ...
                   mfilename, 'choose_species')

% NOTE: in-plane cells must be square because we locate discrete positions
% along the fiber by using a raster line-drawing algorithm.

% Initialize
fib_thick = lumen_thick + 2*wall_thick;
XMIN = 1;
XMAX = Nx;
YMIN = 1;
YMAX = Ny;
spacing = max(fib_width) - 1;
Nx = Nx + spacing; % Add east and
Ny = Ny + spacing; % north margins
Nz = 10*max(fib_thick);
% net_dims = [Nz Nx Ny];
% sheet_dims = [Nx Ny];
NzNx = Nz*Nx;
sheet_area = XMAX*YMAX;
SUBSTRATE = Nz+1;
ismultispecies = nspecies > 1;
deposition_rule_adjustment = Nz*(1-alpha);

switch stop_criterion
    case 'Grammage'
        switch mass
            case 'coarseness'
                multiplier = 1/(delxy*sheet_area);
                fiber_areal_mass = @() fiber_areal_mass_coarseness();
            case 'wall_density'
                multiplier = 2*wall_density*delz/sheet_area;
                fiber_areal_mass = @() fiber_areal_mass_wall_density;
        end
    otherwise
        fiber_areal_mass = @() [];
end

Lx = XMAX;
Ly = YMAX;
scales = [(XMAX-1)/Lx (YMAX-1)/Ly]; %XY scaling factors for axes

% Location in the fiber
TOP    = int32(4);
LUMEN  = int32(3);
WALL   = int32(2);
BOTTOM = int32(1);

% States of the cellular automaton
MOVING = int32(-2);
IDLE   = int32(-1);
EMPTY  = int32(0);

%Preallocate
net = repmat(EMPTY, [Nz Nx Ny]); % 3D fiber network
sheet_top_surface = repmat(int32(SUBSTRATE), Nx, Ny); % elevation grid
wire = true(XMAX, YMAX);

if ~ismultispecies % there is only one species of fibers
    cur_species = 1;
    cur_fib_width = fib_width;
    cur_fib_thick = fib_thick;
    cur_wall_thick = wall_thick;
    cur_lumen_thick = lumen_thick;
    cur_horzspan = horzspan;
    cur_flex = flex;
    cur_flexL = zeros(max(Nx,Ny), 1);
    cur_flexR = cur_flexL;
    cur_flexL(1:horzspan:max(Nx,Ny)) = flex;
    cur_flexR(horzspan:horzspan:max(Nx,Ny)) = flex;
else
    flexL = zeros(max(Nx,Ny), nspecies);
    flexR = flexL;
    for k = 1:nspecies
        flexL(1:horzspan(k):max(Nx,Ny), k) = flex(k);
        flexR(horzspan(k):horzspan(k):max(Nx,Ny), k) = flex(k);
    end
end
speciesSequence = []; % variable to hold the species ordering
last_deposition_event_accepted = true;
offset = 0;
species_idx = Inf;
deposition_idx = Inf;
fibCount = int32(1);
cur_grammage = 0;
switch stop_criterion
    case 'NumFibers'
        converge = @() rule_numfibers();
    case 'Grammage'
        converge = @() rule_grammage();
end
if trace
    hndl_waitbar = waitbar(0, '', 'Name', 'Forming sheet...', ...
        'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
    setappdata(hndl_waitbar, 'canceling', 0);
    switch stop_criterion
        case 'NumFibers'
            complete = @() complete_numfibers();
        case 'Grammage'
            complete = @() complete_grammage();
    end
end

while converge() % Throw one fiber at a time
    if last_deposition_event_accepted
        if species_idx > blocksize(3)
            if ismultispecies
                species = int32(choose_species(fraction, blocksize(3)));
                halffib_length = poissrnd(fib_length(species))/2;
                speciesSequence = [speciesSequence; species];
            else
                halffib_length = poissrnd(fib_length, blocksize(3), 1)/2;
            end
            % Fill in elements corresponding to illegal values.
            halffib_length(halffib_length <= 0) = 1/2;

            species_idx = 1;
        end

        if ismultispecies
            cur_species = species(species_idx);
            cur_fib_width = fib_width(cur_species);
            cur_fib_thick = fib_thick(cur_species);
            cur_wall_thick = wall_thick(cur_species);
            cur_lumen_thick = lumen_thick(cur_species);
        end

        cur_halffib_length = halffib_length(species_idx);
    end

    % Generating fiber deposition on the continuous in-plane direction.
    if deposition_idx > blocksize(1)
        % The center of each fiber is positioned onto a rectangle
        % [0,Lx] x [0,Ly].
        midpoint = unifrnd(zeros(blocksize(1), 2), ...
                           repmat([Lx Ly], blocksize(1), 1), ...
                           blocksize(1), 2);
        orientation = unifrnd(-1, 1, blocksize(1), 2);
        orientation = orientation./hypot( orientation(:, 1), orientation(:, 2) );
        deposition_idx = 1;
    end
    mov = cur_halffib_length * orientation(deposition_idx, :);

    % Map forth between the continuous and discrete world so that the
    % center is placed within the boundaries of the discrete network area
    % [1,Nx] x [1,Ny].
    %
    %   1   2   3   4    <-- discrete coords
    %   o   o   o   o
    % -----------------
    % |   |   |   |   |  <-- continuous coords
    % 0   1   2   3   4
    %
    % Heckbert, Paul S. (1990). What are the coordinates of a pixel?
    % In: Andrew S. Glassner (Ed.), Graphics Gems, Academic Press
    % Professional, Cambridge, MA, pp. 246-248.
    head = floor( (midpoint(deposition_idx, :) - mov) .* scales ) + 1;
    tail = floor( (midpoint(deposition_idx, :) + mov) .* scales ) + 1;

    deposition_idx = deposition_idx + 1;

    % Placing fiber at cells on the discretized in-plane direction.
    % NOTE: Have to eliminate first or last cells along the fiber since
    % each sample point should be treated as representing a finite area
    % (rasterization rule 2 in p. 61 of Watt, Alan and Policarpo, Fabio,
    % "The Computer Image", Addison-Wesley, 1998).
    [fib_cell_xcoord, fib_cell_ycoord] = ...
        bresline(head, tail, cur_fib_width, XMIN, XMAX, YMIN, YMAX);

    [blength, bwidth] = size(fib_cell_xcoord);
    in_bounds = fib_cell_xcoord <= XMAX & fib_cell_ycoord <= YMAX;
    projected_area = sum(in_bounds(:));

%     inplane_idx = sub2ind(sheet_dims, fib_cell_xcoord, fib_cell_ycoord);
    inplane_idx = Nx*(fib_cell_ycoord-1) + fib_cell_xcoord;
    zPos = double( sheet_top_surface(inplane_idx) );

    % Initial positioning of the fiber bottom cell layer.
    if all(zPos(:) == SUBSTRATE)
        overallStartInd = Nz;
        firstVPos = overallStartInd - (cur_fib_thick-1);
        lastVPos = Nz;
        [bend_plane, bend_slice, bheight, net_idx] = get_outofplane_slice();

        if trace
            if fibCount == 1
                col = {'yellow', 'blue', 'magenta', 'red', 'green', 'cyan'};
                hndlFig_vertical = ...
                    figure('Units', 'normalized', ...
                           'OuterPosition', [0.63 0.6 0.37 0.4], ...
                           'Name', 'Bending mechanism of fibers', ...
                           'NumberTitle', 'off', ...
                           'Toolbar', 'none');
            end
            draw_outofplane_grid()
        end

        last_deposition_event_accepted = true;

        startInd = cur_fib_thick;
        bend_plane(startInd, :) = IDLE; % stopped fiber at lattice substrate

        if trace
            [hndlLine_vertical, cur_col] = draw_fibre_initPos();
        end

    else
        % Particle deposition rule
        % in p. 333 of Provatas, N. and Uesaka, T. (2003).
        % Modelling paper structure and paper-press interactions.
        % J. Pulp Pap. Sci. 29(10), 332-340.

        % Rejection model
        % in p. 214 of Provatas, N., Haataja, M., Asikainen, J., Majaniemi,
        % S., Alava, M., and Ala-Nissila T. (2000).
        % Fiber deposition models in two and three spatial dimensions.
        % Colloids Surf. A 165, 209-229.
        rejection_model_rule = rand > acceptanceprob;
        if last_deposition_event_accepted && rejection_model_rule
            avg_thickness = sum(sheet_top_surface(wire))/sheet_area;
        end
        if rejection_model_rule && sum(zPos(in_bounds))/projected_area < ...
                alpha*avg_thickness + deposition_rule_adjustment
            last_deposition_event_accepted = false;
            continue
        else
            last_deposition_event_accepted = true;
        end

        overallStartInd = min(zPos(:)) - 1;

        % Determines if there is enough vertical space for the fiber.
        % If not, or if it is exactly the required one, add space for
        % blocksize(2) more fibers. The equality eliminates the need to
        % include an extra layer on top of bend_plane when cur_flex equals
        % cur_fib_thick for the bending simulation.
        if overallStartInd <= cur_fib_thick
            Gap = blocksize(2)*cur_fib_thick;
            net = cat(1, repmat(EMPTY, [Gap Nx Ny]), net);
            % Recompute the parameters
            Nz = Nz + Gap;
            SUBSTRATE = SUBSTRATE + Gap;
            overallStartInd = overallStartInd + Gap;
            zPos = zPos + Gap;
            sheet_top_surface = sheet_top_surface + int32(Gap);
            if exist('avg_thickness', 'var') == 1
                avg_thickness = avg_thickness + Gap;
            end
            deposition_rule_adjustment = Nz*(1-alpha);
%             net_dims(1) = Nz;
            NzNx = Nz*Nx;
        end

        if ismultispecies
            cur_horzspan = horzspan(cur_species);
            cur_flex = flex(cur_species);
            cur_flexL = flexL(:, cur_species);
            cur_flexR = flexR(:, cur_species);
        end

        % Locate the top and bottom limits of the out-of-plane slice to be
        % extracted from the 3D volume.
        firstVPos = overallStartInd - cur_fib_thick;
        touch_websurface = repmat((1:blength+2)', 1, bwidth);
        touch_websurface = touch_websurface([true(1, bwidth); ...
                                            min(zPos(:)) == zPos; ...
                                            true(1, bwidth)]);
        len = diff(touch_websurface(:)) - 1;
        if length(len) > 3*bwidth - 1
            len(2:end-1) = ceil( len(2:end-1)/2 );
        end
        lastVPos = min(max(zPos(:)) - 1, overallStartInd + max(len)*cur_flex);

        [bend_plane, bend_slice, bheight, net_idx] = get_outofplane_slice();

        if trace
            draw_outofplane_grid();
        end

        startInd = overallStartInd - (firstVPos-1);
        bend_plane(startInd, :) = MOVING; % runnable fiber

        if trace
            [hndlLine_vertical, cur_col] = draw_fibre_initPos();
        end

        % Add edges to complete the neighborhood of first and last cell
        % columns.
        bzero = zeros(bheight, 1, 'int32');
        bend_plane = [bzero bend_plane bzero];

        % Fiber bending simulated using a Cellular Automata framework.
        stop = false(1, blength+2);
        run = [false true(1, blength) false];
        jj = 2:blength+1;
        for ii = startInd:(bheight-1)
            Lidx = sparse(max(ii-cur_flexL(jj-1), 1), jj-1, true, bheight, blength+2);
            Ridx = sparse(max(ii-cur_flexR(jj-1), 1), jj+1, true, bheight, blength+2);
            for k = 1:cur_horzspan
                stop(jj) = bend_plane(Lidx) == IDLE | ...
                           (bend_plane(ii+1, jj) ~= EMPTY)' | ...
                           bend_plane(Ridx) == IDLE;
                bend_plane(ii, run & stop) = IDLE;
                run = run & ~stop;
            end
            if all(~run), break, end % fiber bending stops
            bend_plane(ii, run)   = EMPTY; % move one cell down
            bend_plane(ii+1, run) = MOVING;

            if trace
                [zPos, lPos] = find( bend_plane(:, 2:end-1) == IDLE | ...
                                     bend_plane(:, 2:end-1) == MOVING );
                set(hndlLine_vertical, 'XData', lPos, 'YData', zPos-offset)
                pause(delay_trace)
            end

        end
        bend_plane( end, bend_plane(end, :) == MOVING ) = IDLE;

        % Back to the original grid
        bend_plane(:, [1 end]) = [];
    end

    % Tagging occupied positions according to fiber ID and layer position.
    [zPos, lPos] = find(bend_plane == IDLE);

    lPos = repmat(lPos, [1 bwidth]);
    wPos = repmat(1:bwidth, [blength 1]);
    zPos = repmat(zPos, [1 bwidth]);
%     bottom_idx = sub2ind([bheight blength bwidth], zPos, lPos, wPos);
    bottom_idx = bheight*blength*(wPos-1) + bheight*(lPos-1) + zPos;
    bottom_idx = bottom_idx(:);
    wall_idx   = bottom_idx - ...
        [1:(cur_wall_thick-1) (cur_wall_thick+cur_lumen_thick):(cur_fib_thick-2)];
    lumen_idx  = bottom_idx - ...
        ( cur_wall_thick:(cur_wall_thick+cur_lumen_thick-1) );
    top_idx    = bottom_idx - (cur_fib_thick-1);

    bend_slice(bottom_idx) = complex(fibCount, BOTTOM); % bottom layer
    % Fill in the remaining cell layers
    bend_slice(wall_idx)   = complex(fibCount, WALL);   % wall
    bend_slice(lumen_idx)  = complex(fibCount, LUMEN);  % lumen
    bend_slice(top_idx)    = complex(fibCount, TOP);    % top layer

    % Update elevation grid
    zPos = zPos - (cur_fib_thick-1) + (firstVPos-1);
    sheet_top_surface(inplane_idx) = int32(zPos);

    % Update 3D fiber network
    net(net_idx) = bend_slice;

    if trace
        if fibCount == 1
            lastNz = Nz;
            hndlFig_web = ...
                figure('Units', 'normalized', ...
                       'OuterPosition', [0.63 0.17 0.37 0.4], ...
                       'Name', 'View of the fiber network', ...
                       'NumberTitle', 'off', ...
                       'Toolbar', 'none');
        end
        if fibCount == 1 || Nz > lastNz
            if Nz > lastNz
                warning('Drawing new 3-D plot will erase the previous fibers.');
                lastNz = Nz;
            end
            figure(hndlFig_web)
            clf
            plot3(fib_cell_xcoord, ...
                  fib_cell_ycoord, ...
                  zPos, ...
                  'Color', cur_col);
            axis image
            box on
            set(gca, ...
                'XLim', [0.5 (Nx+0.5)], ...
                'YLim', [0.5 (Ny+0.5)], ...
                'ZLim', [0.5 (Nz+0.5)], ...
                'ZDir', 'reverse')
            hold on
            pause(delay_trace)
        else
            plot3(hndlFig_web.CurrentAxes, ...
                  fib_cell_xcoord, ...
                  fib_cell_ycoord, ...
                  zPos, ...
                 'Color', cur_col)
            pause(delay_trace)
        end
    end

    species_idx = species_idx + 1;
    fibCount = fibCount + 1;
    cur_grammage = cur_grammage + fiber_areal_mass();

    if trace
        if getappdata(hndl_waitbar, 'canceling')
            break
        end
        waitbar(complete())
    end
end

if trace
    delete(hndl_waitbar)
end

% Remove top empty space and lateral margins
ZMAX = min(sheet_top_surface(:));
Nx = Nx - spacing;
Ny = Ny - spacing;
Gap = ZMAX - 1;

net = net(ZMAX:end, 1:Nx, 1:Ny);
sheet_top_surface = sheet_top_surface(1:Nx, 1:Ny);
if ismultispecies
    speciesSequence(fibCount:end) = []; % Make sure vector has appropriate size
end

sheet_top_surface = sheet_top_surface - Gap;
SUBSTRATE = SUBSTRATE - Gap;

sheet_bottom_surface = repmat(SUBSTRATE, Nx, Ny);
for ii = 1:Nx
    for jj = 1:Ny
        if sheet_top_surface(ii, jj) < SUBSTRATE
            sheet_bottom_surface(ii, jj) = find(net(:, ii, jj), 1, 'last');
        end
    end
end

out = struct('web',          real(net), ...
             'surface',      SUBSTRATE - sheet_top_surface, ...
             'sheet_top',    sheet_top_surface, ...
             'sheet_bottom', sheet_bottom_surface, ...
             'grammage', cur_grammage, ...
             'nfib',     fibCount - 1, ...
             'order',    speciesSequence, ...
             'top',    imag(net) == TOP, ...
             'lumen',  imag(net) == LUMEN, ...
             'bottom', imag(net) == BOTTOM);

end