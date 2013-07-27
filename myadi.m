% Alternative Direction Implicit (ADI) method implementation for solving
% parabolic partial differential equations, with Neumann boundary
% conditions:
% u_t = \Nabla(a* \Nabla(u) )+ c*u + f 
%          or
% u_t = \Nabla(a* \Nabla(u) )+ Cg*\Nabla[u \Nabla(g)] + 
%     +   Cphi*\Nabla[u \Nabla(phi)]c*u + f 
% function myadi(u, a, c, f, h_t, Cg, g, Cphi, phi)
%
%
% Input Arguments:
%	u:      [NxN]: the input function. u(x,y,n)
%	a:  	[1x1] or [NxN]: the diffusion coefficient.
%	C:  	[NxN]. c(x,y,n+1/2)
%	f:    [NxN]: f(x,y,n+1/2):=0.5*[ f(x,y,n)+f(x,y,n+1) ]
%	Cg:  	[1x1] or [NxN]: Cg(x,y,n+1/2)
%	g:  	[NxN]: g(x,y,n+1/2)
%	Cphi: [1x1] or [NxN]: Cphi(x,y,n+1/2)
%	phi:  [NxN]: phi(x,y,n+1/2)
%  h_t:  [1x1]: time step.
% Output Arguments:
%	u:	[NxN]: the function u(x,y,n+1)
%
% Author: Giorgos Grekas (grekas.g@gmail.com)

function u= myadi(u, a, C, f, h_t,  Cg, g, Cphi, phi)


N= size(u,1);
h= 1/(N-1);
h2=h*h;

rhs= zeros( size(u) );%

k_t= round( 1/h_t ) +1;


f= 0.5*f;
% tic;
% if( isscalar(a) && nargin == 5 )
%     diag_y_xSweep= k_t + 0.25*C  - a/h2;
%     diag_x_xSweep= k_t - 0.25*C  + a/h2;
% %    [diag_y_xSweep, diag_x_xSweep] = diags_calc(a, C, k_t, N);
%    
% 
%    diag_y_ySweep = diag_x_xSweep;
%    diag_x_ySweep = diag_y_xSweep;
%    
%    
%    y_sub_diagonal_x_sweep= 0.5/h2 * a;
%    y_hyp_diagonal_x_sweep= 0.5/h2 * a;
%    x_sub_diagonal_x_sweep= -0.5/h2 *a;
%    x_hyp_diagonal_x_sweep= -0.5/h2 *a;
% else
   if( nargin == 9)
      g_xx = derivative_xx(g);
      phi_xx = derivative_xx(phi);
      g_yy = derivative_yy(g);
      phi_yy = derivative_yy(phi);
      
      temp_x = derivative_x(a) + Cg.*derivative_x(g) + Cphi.*derivative_x(phi);
      temp_y = derivative_y(a) + Cg.*derivative_y(g) + Cphi.*derivative_y(phi);
      
      diag_temp_y = 0.25*C  - a/h2 + 0.5*(Cg.*g_yy + Cphi.*phi_yy);
      diag_temp_x = 0.25*C  - a/h2 + 0.5*(Cg.*g_xx + Cphi.*phi_xx);
      
   
   diag_y_xSweep2 = k_t + diag_temp_y;
   diag_x_xSweep2 = k_t - diag_temp_x;
   
   diag_y_ySweep2 = k_t - diag_temp_y;
   diag_x_ySweep2 = k_t + diag_temp_x;
   

   y_sub_diagonal_x_sweep2= 0.5/h2 * (a - 0.25*temp_y);
   y_hyp_diagonal_x_sweep2= 0.5/h2 * (a + 0.25*temp_y);
   x_hyp_diagonal_x_sweep2= -0.5/h2 *( a +0.25*temp_x);
   x_sub_diagonal_x_sweep2= -0.5/h2 *( a -0.25*temp_x);

   [diag_y_xSweep, diag_x_xSweep, diag_y_ySweep, diag_x_ySweep,...
         y_hyp_diagonal_x_sweep, y_sub_diagonal_x_sweep, x_hyp_diagonal_x_sweep,...
         x_sub_diagonal_x_sweep] = diags_calc(a, C, k_t, N, Cg, g, Cphi, phi);

   else
      [diag_y_xSweep, diag_x_xSweep, diag_y_ySweep, diag_x_ySweep,...
         y_hyp_diagonal_x_sweep, y_sub_diagonal_x_sweep, x_hyp_diagonal_x_sweep,...
         x_sub_diagonal_x_sweep] = diags_calc(a, C, k_t, N);

   end
u= x_sweep(u, y_sub_diagonal_x_sweep, diag_y_xSweep, y_hyp_diagonal_x_sweep,...
   x_sub_diagonal_x_sweep, diag_x_xSweep, x_hyp_diagonal_x_sweep, f, rhs);


u = myTranspose(u);
x_sub_diagonal_x_sweep = my_minusTranspose(x_sub_diagonal_x_sweep);
diag_x_ySweep = myTranspose(diag_x_ySweep); 
x_hyp_diagonal_x_sweep = my_minusTranspose(x_hyp_diagonal_x_sweep);
y_sub_diagonal_x_sweep = my_minusTranspose(y_sub_diagonal_x_sweep);
diag_y_ySweep = myTranspose(diag_y_ySweep);
y_hyp_diagonal_x_sweep = my_minusTranspose(y_hyp_diagonal_x_sweep);
f = myTranspose(f);
rhs = myTranspose(rhs);

u= x_sweep(u, x_sub_diagonal_x_sweep, diag_x_ySweep, x_hyp_diagonal_x_sweep,...
   y_sub_diagonal_x_sweep, diag_y_ySweep, y_hyp_diagonal_x_sweep, f, rhs);

u = myTranspose(u);

return;



function u= x_sweep(u, y_sub_diag, diag1, y_hyp_diag, x_sub_diag, diag2,...
   x_hyp_diag, f, rhs)
N= size(u,1);

rhs =calculate_rhs(rhs, u, y_sub_diag, diag1, y_hyp_diag, f);


if ( isscalar(x_sub_diag) )
   sub_diag= zeros(N,1);
   sub_diag(1:N-1)= x_sub_diag;
   sub_diag(N)= 2*x_sub_diag;
   
   hyp_diag= zeros(N,1);
   hyp_diag(1)= 2*x_hyp_diag;
   hyp_diag(2:N-1)= x_hyp_diag;

   u=TDMAsolver(u, sub_diag, diag2, hyp_diag, rhs);
%    u=TDMAsolver_par(u, sub_diag, diag2, hyp_diag, rhs);
else
   x_sub_diag(N,:) = x_sub_diag(N,:) + x_hyp_diag(N,:);
   x_hyp_diag(1,:) = x_sub_diag(1,:) + x_hyp_diag(1,:);
   u=TDMAsolver(u, x_sub_diag, diag2, x_hyp_diag, rhs);
%    u=TDMAsolver_par(u, x_sub_diag, diag2, x_hyp_diag, rhs);
end

return;


function df = derivative_xx(f)
% N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(2:end-1,:)= ( f(3:end,:) -2*f(2:end-1,:)+ f(1:end-2,:) ); 
df(1,:)= 2*( f(2,:) -f(1,:) );
df(end,:)= 2*( f(end-1,:) -f(end,:) );

% df = (N-1)^2*df;
return;


function df = derivative_yy(f)
% N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(:,2:end-1)= f(:,3:end) -2*f(:,2:end-1)+ f(:,1:end-2); 
df(:,1)= 2*( f(:,2) -f(:,1) );
df(:,end)= 2*( f(:,end-1) -f(:,end) ) ;

% df = (N-1)^2*df;
return;


function df = derivative_x(f)
N= size(f,1); %h = 1/(N-1)


df= zeros(size(f));

df(2:end-1,:,:)= (f(3:end,:,:) - f(1:end-2,:,:)); 
return;

function df = derivative_y(f)
% N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(:,2:end-1)= (f(:,3:end) - f(:,1:end-2)); 
return;

