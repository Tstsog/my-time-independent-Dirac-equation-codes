% This matlab code solves the time-independent Dirac equation for
% hydrogenic (ns, n=1,2,3...) atom Ref.[1], in which the pseudospectral method is used [2].  
%
% The atomic unit (au) is used in the calculation.
% Note that to get more accurate results, one would increase numerical parameters (N, a, b, etc., )
%
% Refs: [1] W. Greiner, "Relativistic quantum mechanics", Springer-Verlag, Berlin  (1990);
%       [2] Tsogbayar (PhD) thesis, York University, Canada (2014), available at https://yorkspace.library.yorku.ca/xmlui/handle/10315/28160
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% Dec 20, 2024 & University of North Dakota 
%%%
%
function [] = time_independent_dirac_hydrogenic_atom_ns_0p5_states
clc; clear;
format long
%
clear; clc;
Z = 10.00;   % charge of nucleus; you may change it
N = 512.;   % number of grid point; you may change it
rb = 0.0;   % beginning point of coordinate r 
re = 80.;   % ending point of coordnate r; you may change it
%
a = rb;
b = re;
%
kappa = -1.;         % for ns_{1/2} states (n=1,2,3...)
sp_light = 137.;      % speed of light
alpha = 1./sp_light;  % fine-structure constant 
%
[r,wr,D] = legDC2(N,a,b);
wr = wr(2:N);
r = r(2:N);
D1 = (2/(b-a))*D;
%
D1 = D1(2:N,2:N);
%
h11 = -diag(Z./r);
h12 = sp_light*kappa.*diag(1./r) - sp_light*D1;
%
h21 = sp_light*kappa.*diag(1./r) + sp_light*D1;
h22 = -2.*sp_light^2*eye(N-1) - diag(Z./r);
%
H_ham = [h11, h12;
         h21, h22];
%
%%% Lannczos begins
rle = -30.0; % -2, Z = 2; -815 for Z = 40.;
opts.p = 64.;
[Vecf,Enf] = eigs ( H_ham, 5, rle, opts );
Enf = diag(Enf);
[foo, ij] = sort(Enf);
Enf = Enf(ij)
%%% ends

Vecf = Vecf(:,ij);  % The unnormalized eigenfunctions
V1f = Vecf(:,1);    % The unnormalized eigenfunction for the ground state, 
                  % here 1 is for the ground state, 2,3... are corresponded 
                  % for the excited states
%                            
V1f = reshape(V1f,[N-1,2]);
V1f = V1f(:,:); % n = -1, if 
%
wf = V1f; 
wf1 = wf(:,1); % larger (l) component
wf2 = wf(:,2); % smaller (s) component 
%
n_c = sum(wr.*wf1.*wf1 + wr.*wf2.*wf2 ); % normalization 
wf1 = (1./sqrt(n_c)).*wf1;
wf2 = (1./sqrt(n_c)).*wf2;
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% exact energy values
n = (1:1:5);
jj = 0.5; % 
En_exact = -sp_light^2*(Z*alpha)^2.*( (1./(2.*n.^2)) + ((Z*alpha)^2./(2.*n.^3)).*((1/(jj+0.5)) - (3./(4.*n))) )
%

%%%%%
% Z = 2
% En_exact = [-2.000106558687197  -0.500033299589749  -0.222234062076355  -0.125005411183334 -0.080002898396292];
% Enf      = [-2.000106570043467  -0.500033303315963  -0.222234063275764  -0.125005411685265 -0.080002898647383];

%%%%%
% Z = 10
% En_exact = [-50.066599179498098 -12.520812243593157  -5.562955464388677 -3.128381989583888 -2.001811497682348];
% Enf      = [-50.066777163048748 -12.520870655060254  -5.562974263870252 -3.128389853518296 -2.001815426646981];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytic radial wave functions
n = 1.; % for 1s_{1/2} states
n_prime = n - abs(kappa);
En_1 = (1. + (alpha*Z/(n-abs(kappa)+sqrt(kappa)^2-alpha^2*Z^2))^2)^(-1/2);
%%%
gam_s = sqrt(kappa^2 - Z^2.*alpha^2);
lambda = sp_light.*sqrt(1-En_1^2);
aa1_l = (1 + En_1)*gamma(2*gam_s+n_prime+1);
aa1_s = (1 - En_1)*gamma(2*gam_s+n_prime+1);
%
aa2 = 4*((n_prime+gam_s)/En_1)*(((n_prime+gam_s)/En_1) - kappa)*factorial(n_prime);
aa3_l = sqrt(aa1_l/aa2);
aa3_s = sqrt(aa1_s/aa2);
%
g_l_component = ((2.*lambda)^(3/2)/(gamma(2*gam_s+1))).*aa3_l.*(2.*lambda.*r).^(gam_s-1).*exp(-lambda.*r).*(((n_prime+gam_s)/En_1)-kappa);
f_s_component = -((2.*lambda)^(3/2)/(gamma(2*gam_s+1))).*aa3_s.*(2.*lambda.*r).^(gam_s-1).*exp(-lambda.*r).*(((n_prime+gam_s)/En_1)-kappa);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure(1)
%
hold on
plot(r, wf1, 'b', 'LineWidth',1.5)
plot(r, wf2, 'g-', 'LineWidth',1.5)
%%%
plot(r(1:10:length(r)), r(1:10:length(r)).*g_l_component(1:10:length(r)), 'bo', 'MarkerSize',6, 'LineWidth',1.5)
plot(r(1:10:length(r)), r(1:10:length(r)).*f_s_component(1:10:length(r)), 'go', 'MarkerSize',6, 'LineWidth',1.5)
%%%
hold off
xlabel('$r\,(au)$', 'Interpreter','latex')
ylabel('Radial wave function')
set(gca,'FontSize',20) 
axis([0 1.5 -0.30 2.5])  
box on


%%%
return
end



        function [xi,w,D]=legDC2(N,a,b)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % legDc.m
            %
            % Computes the Legendre differentiation matrix with collocation at the
            % Legendre-Gauss-Lobatto nodes.
            %
            % Reference:
            %   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
            %   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
            %
            % Written by Greg von Winckel - 05/26/2004
            % Contact: gregvw@chtm.unm.edu
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % Truncation + 1
            N1=N+1;
            
            % CGL nodes
            xc=cos(pi*(0:N)/N)';
            
            % Uniform nodes
            xu=linspace(-1,1,N1)';
            
            % Make a close first guess to reduce iterations
            if N<3
                x=xc;
            else
                x=xc+sin(pi*xu)./(4*N);
            end
            
            % The Legendre Vandermonde Matrix
            P=zeros(N1,N1);
            
            % Compute P_(N) using the recursion relation
            % Compute its first and second derivatives and
            % update x using the Newton-Raphson method.
            
            xold=2;
            while max(abs(x-xold))>eps
                
                xold=x;
                
                P(:,1)=1;    P(:,2)=x;
                
                for k=2:N
                    P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
                end
                
                x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
            end
            
            X=repmat(x,1,N1);
            Xdiff=X-X'+eye(N1);
            
            L=repmat(P(:,N1),1,N1);
            L(1:(N1+1):N1*N1)=1;
            D=(L./(Xdiff.*L'));
            D(1:(N1+1):N1*N1)=0;
            D(1)=(N1*N)/4;
            D(N1*N1)=-(N1*N)/4;
            
            % Linear map from[-1,1] to [a,b]
            xi=(a*(1-x)+b*(1+x))/2;        % added by Tsogbayar Tsednee
            
            % Compute the weights
            w=(b-a)./(N*N1*P(:,N1).^2);    % added by Tsogbayar Tsednee
            
        end
