function Eout = AAB_carrying_PowerExponentPhase_vortex_beam(r0,w,n,l,Nx,Lx,Ny,Ly)
% Introduce ͻȻ�۽��Ĵ���ָ����λ�İ���������
% generate Spiral autofocusing Airy beams carrying power-exponent-phase vortices
% writen by Ziwei Ye

%   r0   ����
%   w    ��������
%   n    ��ָ������
%   l    ���˺�
%   Lx   x���򳤶�
%   Nx   x����ֳɶ��ٸ���
%   Ly   y���򳤶�
%   Ny   y����ֳɶ��ٸ���

    dx=Lx/Nx;
    x=-Lx/2:dx:Lx/2-dx;
    dy = Ly/Ny;
    y = -Ly/2:dy:Ly/2-dy;
    
    [xx,yy]=meshgrid(x,y);
    [theta,rho]=cart2pol(xx,yy);
    
    Eout = airy((r0-rho)./w) .* exp((r0 - rho)./w) .* exp(1i * 2 * pi * l * (theta/2/pi).^n);
end






