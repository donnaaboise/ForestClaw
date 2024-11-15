function [xp,yp,zp] = mapc2m(xc,yc)

parms = read_vars();

maplist = {'pillowsphere5','pillow'};

map = maplist{parms.mapping+1};

alpha = parms.alpha;

switch map
    case 'pillowsphere5'
        % Map to [-1,1] x [-1,1]
        [xp,yp,~] = mapc2m_pillowdisk5(xc,yc,alpha);
        r2 = min(xp.^2 + yp.^2,1);
        zp = sqrt(1 - r2);
    case 'pillow'
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);

end        

        
end
