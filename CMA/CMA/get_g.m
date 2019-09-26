
function gn=get_g(sn,Rp,p)
     gn=sn/(abs(sn)+0.000001)*(abs(sn)+Rp*abs(sn)^(p-1)-abs(sn)^(2*p-1));
     if p==2
        gn=sn*(Rp-abs(sn)^2);
     end
end


