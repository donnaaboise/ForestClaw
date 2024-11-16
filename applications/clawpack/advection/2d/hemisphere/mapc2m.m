function [xp,yp,zp] = mapc2m(xc,yc)

parms = read_vars();

maplist = {'pillowsphere5','pillow'};

map = maplist{parms.mapping+1};

alpha = parms.alpha;

switch map
    case 'pillowsphere5'

        [xp,yp,zp] = mapc2m_fivepatch(xc,yc,alpha);
        xc = (xp + 1)/2;
        yc = (yp + 1)/2;
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);

        zp = abs(zp);  
        
    case 'pillow'
        [xp,yp,zp] = mapc2m_pillowsphere(xc,yc);

end        

        
end
