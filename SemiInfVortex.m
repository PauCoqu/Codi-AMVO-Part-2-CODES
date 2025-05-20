function [V_inf] = SemiInfVortex(r1,Ur)

Ur1=r1/norm(r1);
V_inf=1/(4*pi)*(1-dot(Ur,Ur1))/(dot(cross(Ur,r1),cross(Ur,r1)))*cross(Ur,r1);
end