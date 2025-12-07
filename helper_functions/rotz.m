function R = rotz(angle)
% ROTZ Rotation matrix about the z-axis
% angle in radians

c = cos(angle);
s = sin(angle);

R = [ c -s  0;
      s  c  0;
      0  0  1];
end
