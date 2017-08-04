function [perpTotal, parrTotal] = totalPowerWebb(n1Real, n2Real, n3Real, zReal, lambdaReal)
global c;
hellow;
k1Real = 2*pi*n1Real/lambdaReal;
const = c*k1Real^4/2/n1Real^3;
perp = @(v, k1, z, n1, n2, n3, lambda) v./sqrt(1-v.^2).*(v.^2).*(1+freselWebbrp(n1, n2, n3, lambda, v, z).*exp(1i*2*k1.*z.*sqrt(1-v.^2))); %u^2 is missing
parr = @(v, k1, z, n1, n2, n3, lambda) v./sqrt(1-v.^2)*(1/2.*(1+freselWebbrs(n1, n2, n3, lambda, v, z)*exp(1i*2*k1.*z.*(1-v.^2).^(1/2))) + 1/2*(1-v.^2).*(1+freselWebbrp(n1, n2, n3, lambda, v, z)*exp(1i*2*k1.*z.*(1-v.^2).^(1/2)))); 
perpReal = @(v) perp(v,k1Real, zReal, n1Real, n2Real, n3Real, lambdaReal); 
parrReal = @(v) parr(v,k1Real, zReal, n1Real, n2Real, n3Real, lambdaReal); 
perp1 = quadc(perpReal, -1, 1, 1e-3, []);
perp2 = integral(perpReal, 1, inf);
parr1 = quadc(parrReal, -1, 1, 1e-3, []);
parr2 = integral(parrReal, 1, inf);
parrTotal = const*(real(0.5*parr1 + parr2));
perpTotal = const*(real(0.5*perp1 + perp2));
end



function rp= freselWebbrp(n1, n2, n3, lambda, v, t)
k1 = 2*pi*n1/lambda;
k2 = 2*pi*n2/lambda;
k3 = 2*pi*n3/lambda;
p = v.*k1;
q1 = sqrt(k1^2 - p.^2);
q2 = sqrt(k2^2 - p.^2);
q3 = sqrt(k3^2 - p.^2);

r12p = (q1.*n2.^2- q2.*n1.^2)/(q1.*n2.^2 + q2.*n1.^2);
r12s = (q1-q2)/(q1+q2);
r23p = (q2.*n3.^2- q3.*n2.^2)/(q2.*n3.^2 + q3.*n2.^2);
r23s = (q2-q3)/(q2+q3);

rs = (r12s+ r23s.*exp(2*1i*k1.*t.*sqrt(n2.^2/n1.^2-v.^2)))/(1+r12s.*r23s.*exp(2*1i*k1.*t.*sqrt(n2.^2/n1.^2-v.^2)));
rp = (r12p+ r23p.*exp(2*1i*k1.*t.*sqrt(n2.^2/n1.^2-v.^2)))/(1+r12p.*r23p.*exp(2*1i*k1.*t.*sqrt(n2.^2/n1.^2-v.^2)));
end



function rs= freselWebbrs(n1, n2, n3, lambda, v, t)
k1 = 2*pi*n1/lambda;
k2 = 2*pi*n2/lambda;
k3 = 2*pi*n3/lambda;
p = v.*k1;
q1 = sqrt(k1^2 - p.^2);
q2 = sqrt(k2^2 - p.^2);
q3 = sqrt(k3^2 - p.^2);

r12s = (q1-q2)/(q1+q2);
r23s = (q2-q3)/(q2+q3);

rs = (r12s+ r23s.*exp(2*1i*k1.*t.*sqrt(n2.^2/n1.^2-v.^2)))/(1+r12s.*r23s.*exp(2*1i*k1.*t.*sqrt(n2.^2/n1.^2-v.^2)));
end

