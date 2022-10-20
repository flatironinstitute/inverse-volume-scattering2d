function [C] = LOCAL_construct_cont(tt,flag_geom,len2)

C             = zeros(6,length(tt));


if(strcmp(flag_geom,'square'))
for jcount = 1:length(tt)
    if tt(jcount)<0.25
    % build south edge
    C(1,jcount) = 4*tt(jcount)-0.5;
    C(2,jcount) = 4*1;
    C(3,jcount) = 0;
    C(4,jcount) = -0.5;
    C(5,jcount) = 0;
    C(6,jcount) = 0;
    elseif (tt(jcount)>0.25 && tt(jcount)<0.5)
    % east edge
    C(1,jcount) = 0.5;
    C(2,jcount) = 0;
    C(3,jcount) = 0;
    C(4,jcount) = 4*(tt(jcount)-0.25)-0.5;
    C(5,jcount) = 4*1;
    C(6,jcount) = 0;
    elseif (tt(jcount)>0.5 && tt(jcount)<0.75)
    C(1,jcount) = 4*(0.75-tt(jcount))-0.5;
    C(2,jcount) = -4*1;
    C(3,jcount) = 0;
    C(4,jcount) = 0.5;
    C(5,jcount) = 0;
    C(6,jcount) = 0;
    else
    C(1,jcount) = -0.5;
    C(2,jcount) = 0;
    C(3,jcount) = 0;
    C(4,jcount) = 4*(1-tt(jcount))-0.5;
    C(5,jcount) = -4*1;
    C(6,jcount) = 0;
    end
end

elseif (strcmp(flag_geom,'bsquare'))
for jcount = 1:length(tt)
    if tt(jcount)<0.25
    % build south edge
    C(1,jcount) = 8*tt(jcount)-1.0;
    C(2,jcount) = 8*1;
    C(3,jcount) = 0;
    C(4,jcount) = -1.0;
    C(5,jcount) = 0;
    C(6,jcount) = 0;
    elseif (tt(jcount)>0.25 && tt(jcount)<0.5)
    % east edge
    C(1,jcount) = 1.0;
    C(2,jcount) = 0;
    C(3,jcount) = 0;
    C(4,jcount) = 8*(tt(jcount)-0.25)-1.0;
    C(5,jcount) = 8*1;
    C(6,jcount) = 0;
    elseif (tt(jcount)>0.5 && tt(jcount)<0.75)
    C(1,jcount) = 8*(0.75-tt(jcount))-1.0;
    C(2,jcount) = -8*1;
    C(3,jcount) = 0;
    C(4,jcount) = 1.0;
    C(5,jcount) = 0;
    C(6,jcount) = 0;
    else
    C(1,jcount) = -1.0;
    C(2,jcount) = 0;
    C(3,jcount) = 0;
    C(4,jcount) = 8*(1-tt(jcount))-1.0;
    C(5,jcount) = -8*1;
    C(6,jcount) = 0;
    end
end

elseif (strcmp(flag_geom,'Nsquare'))
for jcount = 1:length(tt)
    if tt(jcount)<0.25
    % build south edge
    C(1,jcount) = 4*len2*tt(jcount)-len2/2;
    C(2,jcount) = 4*len2*1;
    C(3,jcount) = 0;
    C(4,jcount) = -len2/2;
    C(5,jcount) = 0;
    C(6,jcount) = 0;
    elseif (tt(jcount)>0.25 && tt(jcount)<0.5)
    % east edge
    C(1,jcount) = len2/2;
    C(2,jcount) = 0;
    C(3,jcount) = 0;
    C(4,jcount) = 4*len2*(tt(jcount)-0.25)-len2/2;
    C(5,jcount) = 4*len2*1;
    C(6,jcount) = 0;
    elseif (tt(jcount)>0.5 && tt(jcount)<0.75)
    C(1,jcount) = 4*len2*(0.75-tt(jcount))-len2/2;
    C(2,jcount) = -4*len2*1;
    C(3,jcount) = 0;
    C(4,jcount) = len2/2;
    C(5,jcount) = 0;
    C(6,jcount) = 0;
    else
    C(1,jcount) = -len2/2;
    C(2,jcount) = 0;
    C(3,jcount) = 0;
    C(4,jcount) = 4*len2*(1-tt(jcount))-len2/2;
    C(5,jcount) = -4*len2*1;
    C(6,jcount) = 0;
    end
end    

else
  fprintf(1,'This option for the geometry is not implemented.\n');
  keyboard

end

return