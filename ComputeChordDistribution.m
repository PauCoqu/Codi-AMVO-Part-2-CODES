function c = ComputeChordDistribution(cr, ct, b, Xc)
% Calculate chord lengths for each spanwise segment
% c = 2 * (cr - ct) / b * (b / 2 - abs(Xc(2, :))) + ct;
c=cr - (abs(0-Xc(2, :))*2*(cr-ct))/b;
end