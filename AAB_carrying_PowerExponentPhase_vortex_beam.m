function Eout = AAB_carrying_PowerExponentPhase_vortex_beam(r0,w,n,l,Nx,Lx,Ny,Ly)
% Introduce 突然聚焦的带幂指数相位的艾里涡旋光
% generate Spiral autofocusing Airy beams carrying power-exponent-phase vortices
% 论文 https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-22-7-7598&id=282225
% writen by Ziwei Ye

%   r0   参数
%   w    光腰参数
%   n    幂指数参数
%   l    拓扑荷
%   Lx   x方向长度
%   Nx   x方向分成多少个点
%   Ly   y方向长度
%   Ny   y方向分成多少个点

    dx=Lx/Nx;
    x=-Lx/2:dx:Lx/2-dx;
    dy = Ly/Ny;
    y = -Ly/2:dy:Ly/2-dy;
    
    [xx,yy]=meshgrid(x,y);
    [theta,rho]=cart2pol(xx,yy);
    
    Eout = airy((r0-rho)./w) .* exp((r0 - rho)./w) .* exp(1i * 2 * pi * l * (theta/2/pi).^n);
end






