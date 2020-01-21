function mR2eff = model_MQRD_CRJ( vcpmg, TCPMG, params )
% Calculate the R2Eff(vCPMG) RD curve using MQ Carver-Richards equations and exchange parameters
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2008/07/01 Start code
% 2011/01/13 Change dwC -> dwX, vectorize some calculations
% 2011/04/22 Rename to model_MQRD_CRJ
%
% This returns the model function for the MQ R2eff(Vcpmg) given a set of parameters
%  D.M. Korzhnev, K. Kloiber, V. Kanelis, V. Tugarinov, L.E. Kay, Probing slow dynamics in high
%  molecular weight proteins by methyl-trosy nmr spectroscopy: application to a 723-residue enzyme,
%  J. Am. Chem. Soc. (2004) 3964-3973.
%
% INPUT  
%   vcpmg   = CPMG frequency array
%   TCPMG       = Total CPMG time
%   params       = Parameter array (dwH,dwX,Pa,Kex,R2mq)
%
% OUTPUT
%  mR2eff   = Model R2Eff(vCPMG) function at each input point vCPMG
 
    % Number of observations
    N = length(vcpmg);    

    %vcpmg = vcpmg';    
    
    % TODO % When using small values of TCPMG (e.g., 5 ms) then RD curve drops
    % below R2mq (I think this shouldn't happen, nor should it be sensitive
    % to TCPMG!)    

    % Set values of the parameters to vary for best fit
    dw_h	= params(1);
	dw_x	= params(2);
	pa      = params(3);
	pb      = 1 - params(3);
	kex     = params(4);
	R2mq	= params(5);

    % These arrays accompany R2eff values
    delta = 1 ./ (4*vcpmg);
    n = TCPMG ./ (4*delta);
    
    % (3.10) scalars
    z_plus  = dw_h - dw_x + sqrt(-1)*kex;
    z_minus = dw_h - dw_x - sqrt(-1)*kex;

    d_plus  = dw_h + dw_x + sqrt(-1)*kex;
    d_minus = dw_h + dw_x - sqrt(-1)*kex;
    
        
    % (3.9) vector
    mZ = -sqrt(-1)*kex*sqrt(pa*pb) / (d_minus*z_minus) * ...
        ( d_minus - 2*dw_x*sin(d_minus.*delta) ./ sin((d_minus+z_minus).*delta) );

    % (3.8) vector
    mD = sqrt(-1)*kex*sqrt(pa*pb) / (d_plus*z_plus) * ...
        ( z_plus + 2*dw_x*sin(z_plus.*delta) ./ sin((d_plus+z_plus).*delta) ) ;

    % (3.7) - vector
    Q = real( ones(1,N) - mD.^2 + mD.*mZ - mZ.^2 + 0.5*(mD + mZ)*sqrt(pb/pa) );
    
    %mZ
    %mD
    %Q
    
    
    % (3.6) - scalar
    zeta = -2*dw_x*(sqrt(-1)*dw_h + (pa-pb)*kex);
    
    % (3.5) - scalar
    Psi = (sqrt(-1)*dw_h + (pa-pb)*kex)^2 - dw_x^2 + 4*pa*pb*kex*kex;
    
    % (3.4) - vectors
    eta_plus = sqrt(2)*delta*sqrt( sqrt(Psi^2 + zeta^2) + Psi );
    eta_minus= sqrt(2)*delta*sqrt( sqrt(Psi^2 + zeta^2) - Psi );
    
    % (3.3) - scalars    
    D_plus = 0.5*( (Psi + 2*dw_x^2) / sqrt(Psi^2 + zeta^2) + 1 );
    D_minus= 0.5*( (Psi + 2*dw_x^2) / sqrt(Psi^2 + zeta^2) - 1 );

        
    % (3.2) - vector
    lambda1 = R2mq + 0.5*(kex - (1./(2.*delta)) .* ...
           acosh( D_plus*cosh(eta_plus) - D_minus*cos(eta_minus) ) );
        
    % (3.1) - vector
    mR2eff = real(lambda1) - log(Q)./(4*n.*delta);
           
    %lambda1;
    %mR2eff;    
    
end
