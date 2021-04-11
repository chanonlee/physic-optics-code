function[Eout] = Fraunhofer(Lx0,Nx0,Ly0,Ny0,E,lambda,distance,Lx,Nx,Ly,Ny)
% ���źͷ��������
% Introduce
% Fraunhofer diffraction integral
% writen by Ziwei Ye

%   Lx0      ������x���򳤶�
%   Nx0      �����x����ֳɶ��ٸ���
%   Ly0      ������y���򳤶�
%   Ny0      �����y����ֳɶ��ٸ���
%   E        ����ⳡ������
%   lambda   ����
%   distance �������� 
%   Lx       x���򳤶�
%   Nx       x����ֳɶ��ٸ���
%   Ly       y���򳤶�
%   Ny       y����ֳɶ��ٸ���

%==================================================
%                   �������ݴ���
%==================================================
dx0 = Lx0/Nx0;
x0 = -Lx0/2:dx0:Lx0/2-dx0;
    
dy0 = Ly0 / Ny0;
y0 = -Ly0/2:dy0:Ly0/2 - dy0;
    
dx = Lx/Nx;
x = -Lx/2:dx:Lx/2-dx;
        
dy = Ly / Ny;
y = -Ly/2:dy:Ly/2 - dy;
    
k = 2 * pi / lambda;
Eout = zeros(Lx,Ly);
%==================================================
%                   �������ݴ���
%==================================================

%==================================================
%                   ѭ���������
%==================================================
% �Գ���ⳡ��ѭ��
for a = 1:Lx
    for b = 1:Ly
        
        temp = 0;
        xishu = exp(1i * k * distance) / (1i * lambda * distance) * exp(1i * k /(2 * distance) * (x(a)^2 + y(b)^2));
        for aa = 1:Lx0
            for bb = 1:Ly0
                temp = temp + E(aa,bb) * exp(-1i * 2 * pi * (x(a)/(lambda * distance) * x0(aa) + y(b)/(lambda * distance) * y0(bb))) * dx0 * dy0;
            end
        end
        Eout(a,b) = xishu * temp;
        
    end
end
%==================================================
%                   ѭ���������
%==================================================

