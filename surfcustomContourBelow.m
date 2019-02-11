function h = surfcustomContourAbove(varargin)
    %SURFC  Combination surf/contour plot.
    %   SURFC(...) is the same as SURF(...) except that a contour plot
    %   is drawn beneath the surface.
    %
    %   See also SURF, SHADING.
    
    %   Clay M. Thompson 4-10-91
    %   Copyright 1984-2010 The MathWorks, Inc.
    %   $Revision: 5.11.4.5 $  $Date: 2011/04/04 22:59:19 $
    
    % First we check which HG plotting API should be used.
    if ishg2parent( varargin{:} )
        [~, cax, args] = parseplotapi(varargin{:},'-mfilename',mfilename);
        if nargout > 0
            h = surfcHGUsingMATLABClasses(cax, args{:});
        else
            surfcHGUsingMATLABClasses(cax, args{:});
        end
    else
        [cax, args] = axescheck(varargin{:});
        [reg, prop] = parseparams(args);
        nargs = length(reg);
        
        error(nargchk(1, 4, nargs, 'struct'));
        if rem(length(prop), 2) ~= 0
            error(message('MATLAB:surfc:InvalidPVPair'))
        end
        
        if nargs == 1  % Generate x, y matrices for surface z.
            z = reg{1};
            [m, n] = size(z);
            [x, y] = meshgrid(1 : n, 1 : m);
        elseif nargs == 2
            z = reg{1};
            [m, n] = size(z);
            [x, y] = meshgrid(1 : n, 1 : m);
        elseif nargs == 3
            [x, y, z] = deal(reg{1 : 3});
        elseif nargs == 4
            [x, y, z] = deal(reg{1 : 3});
        end
        
        if min(size(z)) == 1
            error(message('MATLAB:surfc:NonMatrixInput'));
        end
        
        % Determine state of system
        if isempty(cax)
            cax = gca;
        end
        next = lower(get(cax, 'NextPlot'));
        hold_state = ishold(cax);
        
        % Plot surface
        hs = surf(cax, args{:});
        
        hold(cax, 'on');
        
        a = get(cax, 'ZLim');
        
        % Always put contour below the plot.
        zpos = -30;
        
        % Get D contour data
        [~, hh] = contour3(cax, x, y, z);
        
        % size zpos to match the data
        for i = 1 : length(hh)
            zz = get(hh(i), 'ZData');
            set(hh(i), 'ZData', zpos * ones(size(zz)));
        end
        
        if ~hold_state
            set(cax, 'NextPlot', next);
        end
        
        if nargout > 0
            h = [hs; hh(:)];
        end
    end
end

