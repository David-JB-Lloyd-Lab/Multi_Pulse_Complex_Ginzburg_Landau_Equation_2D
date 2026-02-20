function [value,isterminal,direction] = events_period2(t,y)
x=abs(y(1)-y(4)).*cos(y(3)-y(6));
value=[x;x];
isterminal=[0;1];
direction=[1;-1];
end

