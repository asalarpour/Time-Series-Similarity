function lineFrac = TS_freeLine(pt, eps, seg)

%     ************
%
%     Description
%     ----------
%     Return the free space in the segment 'seg', from point 'pt'.
%     This free space is the set of all point in 'seg' whose distance from 'pt' is at most 'eps'.
%     Since 'seg' is a segment, the free space is also a segment.
%     We return a 1 x 2 array with the fraction of the segment 'seg' which are in the free space.
%     If no part of 'seg' are in the free space, return [-1,-1]
% 
%     Parameters
%     ----------
%     param pnt     : 1 x dim, the point centre of the circle
%     param eps     : double, radius of the circle
%     param seg     : 2 x dim , two line segment points
%
%     Returns
%     -------
%     lineFrac      : 1 x dim , fraction of segment which is in the free
%                               space (i.e [0.3,0.7], [0.45,1], ...)
%                               If no part of s are in the free space, return [-1,-1]
%
%     Other m-files required    : TS_euclideanDistance, 
%                                 TS_point2SegmentDistance,
%                                 TS_circleLineIntersection
%     Subfunctions              : none
%     MAT-files required        : none
%     
%     References
%     ----------
%
%     Author
%     ----------
%     Amir Salarpour, Ph.D. Candidate, Artificial Intelligence
%     Bu-Ali Sina University, Hamedan, Iran, Dept. of Computer Engineering
%     email address : amir.salarpour@gmail.com  
%     Website       : http://www.salarpour.com
%     December 2016 : Last revision: 29-Jan-2017
%     
%     ************

px = pt(1);
py = pt(2);
seg1x = seg(1, 1);
seg1y = seg(1, 2);
seg2x = seg(2, 1);
seg2y = seg(2, 2);

if (seg1x == seg2x && seg1y==seg2y)
    if (TS_euclideanDistance(pt, seg(1, :)) > eps)
        lineFrac = [-1, -1];
    else
        lineFrac = [0, 1];
    end
else
    ttmp = TS_point2SegmentDistance(pt, seg(1, :), seg(2, :));
    if (ttmp > eps)
%         fprintf('%f  ', TS_point2SegmentDistance(pt, seg(1, :), seg(2, :)));
%         disp('No Intersection')
        lineFrac = [-1, -1];
    else
        segLen = TS_euclideanDistance(seg(1, :), seg(2, :));
        segLenS = segLen * segLen;
        intersect = TS_circleLineIntersection(pt, seg(1, :), seg(2, :), eps);
        
        if ((intersect(1, 1) ~= intersect(2, 1)) || (intersect(1, 2) ~= intersect(2, 2)))
            
            i1x = intersect(1, 1);
            i1y = intersect(1, 2);
            u1 = (((i1x - seg1x) * (seg2x - seg1x)) + ((i1y - seg1y) * (seg2y - seg1y))) / segLenS;
            
            i2x = intersect(2, 1);
            i2y = intersect(2, 2);
            u2 = (((i2x - seg1x) * (seg2x - seg1x)) + ((i2y - seg1y) * (seg2y - seg1y))) / segLenS;
            ordered_point = sort([0, 1, u1, u2]);
            lineFrac = ordered_point(2:3);
        else
            if (px == seg1x && py == seg1y)
                lineFrac = [0, 0];
            elseif (px == seg2x && py==seg2y)
                lineFrac = [1, 1];
            else
                i1x = intersect(1, 1);
                i1y = intersect(1, 2);
                u1 = (((i1x - seg1x) * (seg2x - seg1x)) + ((i1y - seg1y) * (seg2y - seg1y))) / segLenS;
                if (0 <= u1 && u1 <= 1)
                    lineFrac = [u1, u1];
                else
                    lineFrac = [-1, -1];
                end
            end
        end
    end
end
end

