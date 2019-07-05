function [xcoord, ycoord] = bresline(first, last, width, XMIN, XMAX, YMIN, YMAX)
%BRESLINE Draws 2-D line from (x,y) to (xend-1,yend-1) using Bresenham's algorithm.
%   This implementation is a port of the source code found in Agoston (2005).

%   Source:
%      Agoston, Max K. (2005).
%      Computer Graphics and Geometric Modeling: Implementation &
%      Algorithms, Springer-Verlag, London, pp. 38-41.

%   Copyright 2008, 2011, University of Coimbra
%   Copyright 2018, Eduardo L. T. Conceicao
%   Available under the GPL-3


x = first(1);
y = first(2);
xend = last(1);
yend = last(2);
dx = xend - x;
dy = yend - y;
sgnx = sign(dx);
sgny = sign(dy);
ax = abs(dx);
ay = abs(dy);
two_ax = 2*ax;
two_ay = 2*ay;

if ax > ay
    % x increases faster than y
    d = two_ay - ax;

    xcoord = zeros(ax, 1);
    ycoord = xcoord;
    for i = 1:ax
        xcoord(i) = x;
        ycoord(i) = y;
        while d >= 0
            y = y + sgny;
            d = d - two_ax;
        end
        x = x + sgnx;
        d = d + two_ay;
    end

    % Clipping (scissoring) to restrict the line segment to the area
    % [XMIN,XMAX] x [YMIN,YMAX].
    clipped = xcoord < XMIN | xcoord > XMAX | ...
              ycoord < YMIN | ycoord > YMAX;
    xcoord(clipped) = [];
    ycoord(clipped) = [];

    if width > 1
        xcoord = repmat(xcoord, 1, width);
        ycoord = ycoord + (0:(width-1));
    end
else
    % y increases faster than x
    d = two_ax - ay;

    xcoord = zeros(ay, 1);
    ycoord = xcoord;
    for i = 1:ay
        xcoord(i) = x;
        ycoord(i) = y;
        while d >= 0
            x = x + sgnx;
            d = d - two_ay;
        end
        y = y + sgny;
        d = d + two_ax;
    end

    % Clipping (scissoring) to restrict the line segment to the area
    % [XMIN,XMAX] x [YMIN,YMAX].
    clipped = xcoord < XMIN | xcoord > XMAX | ...
              ycoord < YMIN | ycoord > YMAX;
    xcoord(clipped) = [];
    ycoord(clipped) = [];

    if width > 1
        ycoord = repmat(ycoord, 1, width);
        xcoord = xcoord + (0:(width-1));
    end
end

end