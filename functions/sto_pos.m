%   Stochastic subspace identification (Algorithm 3)
%
%           [A,K,C,R] = sto_pos(y,i);
% 
%   Inputs:
%           y: matrix of measured outputs
%           i: number of block rows in Hankel matrices 
%              (i * #outputs) is the max. order that can be estimated 
%              Typically: i = 2 * (max order)/(#outputs)
%           
%   Outputs:
%           A,K,C,R: stochastic state space system
%           
%                 x_{k+1} = A x_k + K e_k
%                   y_k   = C x_k + e_k
%                cov(e_k) = R
%                
%   Optional:
%
%           [A,K,C,R,AUX,G,L0,ss] = sto_pos(y,i,n,AUX,W,sil);
%   
%           n:    optional order estimate (default [])
%           AUX:  optional auxilary variable to increase speed (default [])
%           W:    optional weighting flag
%                    CVA: canonical variate analysis (default)
%                    PC:  principal components
%                    UPC: unweighted principal components
%           G,L0: covariance model
%           ss:   column vector with singular values
%           sil:  when equal to 1 no text output is generated
%           
%   Note:        
%           The computed system is ALWAYS positive real
%
%   Example:
%   
%           [A,K,C,R,AUX] = sto_pos(y,10,2);
%           for k=3:6
%              [A,K,C,R] =sto_pos(y,10,k,AUX);
%           end
%           
%   Reference:
%   
%           Subspace Identification for Linear Systems
%           Theory - Implementation - Applications
%           Peter Van Overschee / Bart De Moor
%           Kluwer Academic Publishers, 1996, Page 90 (Fig 3.13): sto_pos.m algorithm
%           
%   Copyright:
%   
%           Peter Van Overschee, December 1995
%           peter.vanoverschee@esat.kuleuven.ac.be
%   
%   Revised:
%           Andrea saettone, January 2005
%           Giuseppe Abbiati, August 2009
%           Marica Pecorelli, Ottobre 2016

function [A,K,C,Ro,AUX,G,L0,ss] = sto_pos(y,i,n,AUXin,W,sil)
% nargin = number of function input arguments
if (nargin < 6)         % controlla l'output                                                   
    sil = 0;
end
if (nargin < 2)
    error('sto_pos needs at least two arguments');
end
if (nargin < 3)
    n = [];
end
if (nargin < 4)
    AUXin = [];
end
if (nargin < 5)
    W = [];
end
if isempty(W)
    W = 'CVA';
end

% Verifico che i segnali siano sulle righe
[l,ny] = size(y);
if (ny < l)
    y = y';
    [l,ny] = size(y); % ny = campioni del segnale acquisito
end

if (i < 0)
    error('Number of block rows should be positive');
end
if (l < 0)
    error('Need a non-empty output vector');
end
if ((ny-2*i+1) < (2*l*i))
    error('Not enough data points');
end
Wn = 0;

if (length(W) == 3) 
  if ((W == 'CVA') | (W == 'cva') | (W == 'Cva'))
      Wn = 1;
  end 
  if ((W == 'UPC') | (W == 'upc') | (W == 'Upc'))
      Wn = 3;
  end
end    
if (length(W) == 2) 
  if ((W == 'PC') | (W == 'pc') | (W == 'Pc'))
      Wn = 2;
  end 
end
if (Wn == 0)
    error('W should be CVA, PC or UPC');
end
W = Wn;

% Determine the number of columns in Hankel matrices
j = ny-2*i+1;

%=========================================================================
% Check compatibility of AUXin CHE NON FUNZIONA PER MOTIVI APPARENTEMENTE
% SENZA SENSO PER CUI SOSTITUITO DA
%[AUXin,Wflag] = chkaux(AUXin,i,[],y(1,1),2,W,sil); 
AUXin=[];
Wflag=1;
%=========================================================================

% Compute the R factor
if isempty(AUXin)
  Y = blkhank(y/sqrt(j),2*i,j); 		% Output block Hankel
  
  disp('Computing ... R factor');
  R = triu(qr(Y'))'; 			% R factor
  R = R(1:2*i*l,1:2*i*l); 		% Truncate
%   save marica R
%   pause
  clear Y
else
  R = AUXin(2:2*i*l+1,1:2*i*l);
  bb = 2*i*l+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                  BEGIN ALGORITHM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% **************************************
%               STEP 1 
% **************************************

% First compute the orthogonal projections Ob and Obm
if isempty(AUXin)
  Ob  = R(l*i+1:2*l*i,1:l*i);
else
  Ob  = AUXin(bb+1:bb+l*i,1:l*i);
  bb  = bb+l*i;
end

% **************************************
%               STEP 2 
% **************************************

% Compute the SVD
if isempty(AUXin) | (Wflag == 1)
  disp('Computing ... SVD');
  % Compute the matrix WOW we want to take an svd off
  % W == 1 (CVA), W == 2 (PC), W == 3 (UPC)
  if (W == 1)
    W1i = triu(qr(R(l*i+1:2*l*i,1:2*l*i)'));
    W1i = W1i(1:l*i,1:l*i)';
    WOW = W1i\Ob;
  end
  if (W == 2)
    WOW = R(l*i+1:2*l*i,1:l*i)*R(1:l*i,1:l*i)';
  end
  if (W == 3)
    WOW = Ob;
  end
  [U,S,V] = svd(WOW);
  if (W == 1)
      U = W1i*U;
  end
  ss = diag(S);
  clear V S WOW
else
  U = AUXin(bb+1:bb+l*i,1:l*i);
  ss = AUXin(bb+1:bb+l*i,l*i+1);
end

% **************************************
%               STEP 3 
% **************************************

% Determine the order from the singular values
if isempty(n)
  n = 0;
  while (n < 1) | (n > l*i-1)
    n = input('System order ?');
    if isempty(n)
        n = -1;
    end
  end
end

U1 = U(:,1:n); 				% Determine U1

% **************************************
%               STEP 4 
% **************************************

% Determine gam and gamm
gam  = U1*diag(sqrt(ss(1:n)));
gamm = U1(1:l*(i-1),:)*diag(sqrt(ss(1:n)));
% And their pseudo inverses
gam_inv  = pinv(gam);
gamm_inv = pinv(gamm);
clear gam gamm

% **************************************
%               STEP 5 
% **************************************

% Determine the states Xi and Xip
Xi  = gam_inv  * Ob;
Xip = gamm_inv * R(l*(i+1)+1:2*l*i,1:l*(i+1));
clear gamm_inv

% **************************************
%               STEP 2a 
% **************************************

% Determine the state matrices A and C
disp(['Computing ... System matrices A,C (Order ',num2str(n),')']);

Rhs = [Xi,zeros(n,l)];                          % Right hand side
Lhs = [Xip;R(l*i+1:l*(i+1),1:l*(i+1))];         % Left hand side
sol = Lhs/Rhs;                                  % Solve least squares

% Extract the system matrices
A = sol(1:n,1:n);
C = sol(n+1:n+l,1:n);


% **************************************
%               STEP 3a 
% **************************************

disp(['Computing ... System matrices G,L0 (Order ',num2str(n),')']); 
% Determine the residuals
res = Lhs - sol*Rhs;                                                % Residuals
cov = res*res';                                                     % Covariance
Qs = cov(1:n,1:n);
Ss = cov(1:n,n+1:n+l);
Rs = cov(n+1:n+l,n+1:n+l); 

% **************************************
%               STEP 4b 
% **************************************

sig = dlyap(A,Qs);
G = A*sig*C' + Ss;
L0 = C*sig*C' + Rs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                  END ALGORITHM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine K and Ro
disp('Computing ... Riccati solution');
[K,Ro] = gl2kr(A,G,C,L0);

% Make AUX when needed
if nargout > 4
  AUX = zeros((4*l)*i+1,2*l*i);
  info = [2,i,0,y(1,1),W]; % in/out - i - u(1,1) - y(1,1) - W
  AUX(1,1:5) = info;
  bb = 1;
  AUX(bb+1:bb+2*l*i,1:2*l*i) = R;
  bb = bb+2*l*i;
  AUX(bb+1:bb+l*i,1:l*i) = Ob;
  bb = bb+l*i;
  AUX(bb+1:bb+l*i,1:l*i) = U;
  AUX(bb+1:bb+l*i,l*i+1) = ss;
end






