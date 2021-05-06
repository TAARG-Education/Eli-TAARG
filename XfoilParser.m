function [ flag, Cl, Cd ] = XfoilParser( airfoil, alfa, Re, varargin)

Nv = length( varargin );
if Nv > 0
    if ischar( varargin{1} ) || isstring( varargin{1} )
        plotFlag = true;
    end
else
    plotFlag = false;
end



% Variable initialization
Cl = nan;
Cd = nan;
flag = 0;

span = 2;
step = .5;
alfaVet = alfa - span : step : alfa + span;
Np = 200;                                  % number of panel

[ ClVet , CdVet, alfarVet ] = CdCl_xfoil( airfoil, Np, Re, alfaVet, 'plot' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two angles have to converge at least
Ns = length( alfarVet );
if Ns < 2
    error( "Poor Xfoil convergence. Less then 2 alfa convergence." + ...
        "Impossible interpolation.");
else
    % Extrapolation check
    if ( alfa < alfarVet( 1 ) || alfa > alfarVet( end ) )
        flag = 2;
    else
        flag = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%f
% Interpolation to extract the desired alfa-Cl-Cd value
Cl = spline( alfarVet, ClVet, alfa );
Cd = spline( alfarVet, CdVet, alfa );


if plotFlag
    f1 = figure( 1 );
    movegui( f1, 'northeast' );
    plot( alfarVet, ClVet, alfa, Cl, 'o' );
    xlabel('Angle of attach \alpha');
    ylabel('Lift coefficient C_{l}');
end


end

