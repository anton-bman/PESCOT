function [omega0,r_pseudo,phi_pseudo] = pseudotrue_omega0(omega_vec,r_vec,phi_vec,t,harm_order_vec,omega_int,nbr_of_omegas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the pseudo-true fundamental frequency, as well as the
% pseudo-true amplitudes and phases or the harmonics of the harmonic
% approximation.
%
%
% "Defining Fundamental Frequency for Almost Harmonic Signals", Elvander
% and Jakobsson, IEEE Transaction on Signal Processing vol 68, 2020.
%
% DOI: 10.1109/TSP.2020.3035466
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT
% omega_vec             -       vector of frequencies.
% r_vec                 -       vector of magnitudes, i.e., the absolute
%                               values of the complex amplitudes.
% phi_vec               -       vector of initial phases.
% t                     -       vector of sampling times.
% 
% INPUT (optional)
% harm_order_vec        -       vector of harmonic orders. If left empty,
%                               it is assumed that the orders are 1:K,
%                               where K are the number of elements of
%                               omega_vec.
% omega_int             -       search interval for the 
%                               pseudo-true fundamental frequency.
% nbr_of_omegas         -       initial number of grid points for the
%                               search of the pseudo-true fundamental 
%                               frequency.
%
% OUTPUT
% omega0                -       pseudo-true fundamental frequency.
% r_pseudo              -       vector of pseudo-true magnitudes.
% phi_pseudo            -       vector of pseudo-true initial phases.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7 || isempty(nbr_of_omegas)
    nbr_of_omegas = 200;
end
if nargin < 6 || isempty(omega_int)
    omega_int = mean(omega_vec./harm_order_vec)*[0.7,1.5];
end
if nargin<5 || isempty(harm_order_vec)
    harm_order_vec = (1:length(omega_vec))';
end

x = exp(1i*t*omega_vec')*(r_vec.*exp(1i*phi_vec));
[omega0,alpha_vec] = estimate_omega0(x,t,harm_order_vec,omega_int,nbr_of_omegas);

r_pseudo = abs(alpha_vec);
phi_pseudo = angle(alpha_vec);

end
