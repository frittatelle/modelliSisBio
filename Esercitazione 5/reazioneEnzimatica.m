function [dy] = reazioneEnzimatica(y,k1,k_1,k2,k_2)

s = y(1);
e = y(2);
c = y(3);
p = y(4);

dy(1) = -k1 * e*s + k_1 * c;
dy(2) = -k1 * e*s - k_2 * e*p + (k_1+k2) * c;
dy(3) = k1 * e*s - (k_1+k2) * c + k_2 *e*p;
dy(4) = k2*c - k_2 *e*p;

dy = dy';

end

