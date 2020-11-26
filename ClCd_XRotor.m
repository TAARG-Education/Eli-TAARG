function [Cl_vec, Cd, Cd_bp] = ClCd_XRotor(Cd_min, dCd_dCl2, Cl_Cd_min, Re_ref, Re_inf, f, Cl_max, Cl_min, Cl_bp)
    
    Cl_vec = Cl_min:.01:Cl_max;

    for i = 1:length(Cl_vec)
        Cl = Cl_vec(i);
        Cd(i) = (Cd_min + dCd_dCl2*(Cl_Cd_min - Cl)^2)*(Re_inf/Re_ref)^f;
    end

    % plot
    plot(Cd, Cl_vec);
    
    %% evaluation of Cd values @ requested Cl breakpoints 
    % note that Cl_bp can be both a single value or a vector
    % definition of an anonymous function through the handle (@) symbol
    Cd_interp = @(Cl_interp) interp1(Cl_vec, Cd, Cl_interp, 'pchip');
    % anonymous fcn used to compute requested Cd
    Cd_bp = Cd_interp(Cl_bp);


end