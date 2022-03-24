function A = GetTheodorsenMatrix(param,k)

b   = param.b;
a   = param.a;

if k==0
 Ck=1;
else
 H0=besselh(0,2,k);
 H1=besselh(1,2,k);
 Ck=H1/(H1+1i*H0);
end

Ph = -2*b^2*(-k^2*pi*(1/b)^2)                                 -2*pi*2*b*Ck*1i*k/b;
Pa = -2*b^2*(pi*1i*k/b+pi*b*a*k^2*(1/b)^2)                    -2*pi*2*b*Ck*(1+b*(1/2-a)*1i*k/b);
Mh = -2*b^2*(k^2*a*pi*b*(1/b)^2)                              +2*pi*2*b*Ck*b*(1/2+a)*1i*k/b;
Ma = -2*b^2*(pi*(1/2-a)*b*1i*k/b-pi*b^2*(1/8+a^2)*k^2*(1/b)^2)+2*pi*2*b*Ck*b*(1/2+a)*(1+b*(1/2-a)*1i*k/b);

A = [Ph Pa;Mh Ma];