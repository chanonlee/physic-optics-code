function[Eout] = Fraunhofer(Lx0,Nx0,Ly0,Ny0,E,lambda,distance,Lx,Nx,Ly,Ny)
% 夫琅和费衍射积分
% Introduce
% Fraunhofer diffraction integral
% writen by Ziwei Ye

%   Lx0      入射光的x方向长度
%   Nx0      入射光x方向分成多少个点
%   Ly0      入射光的y方向长度
%   Ny0      入射光y方向分成多少个点
%   E        入射光场，复数
%   lambda   波长
%   distance 传播距离 
%   Lx       x方向长度
%   Nx       x方向分成多少个点
%   Ly       y方向长度
%   Ny       y方向分成多少个点

%==================================================
%                   输入数据处理
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
%                   输入数据处理
%==================================================

%==================================================
%                   循环计算积分
%==================================================
% 对出射光场的循环
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
%                   循环计算积分
%==================================================

