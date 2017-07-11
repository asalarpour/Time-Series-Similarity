function pt2segDist = TS_point2SegmentDistance(queryPoint, segmentStart, segmentEnd)

%     ************
%
%     Description
%     ----------
%     Compute the shortest distance between queryPoint and segment delimited by segmentStart and segmentEnd
% 
%     Parameters
%     ----------
%     param queryPoint      : 1 x dim, define a point in dim-dimension space
%     param segmentStart    : 1 x dim, define the start point of the line segment
%     param segmentEnd      : 1 x dim, define the end point of the line segment
%
%     Returns
%     -------
%     pt2segDis             :  double, The point to segment distance.
%
%     Other m-files required    : TS_euclideanDistance
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

A = queryPoint(1) - segmentStart(1);
B = queryPoint(2) - segmentStart(2);
C = segmentEnd(1) - segmentStart(1);
D = segmentEnd(2) - segmentStart(2);
dotProd = A * C + B * D;
lenSq = C * C + D * D;

if (lenSq == 0)
    
    pt2segDist = TS_euclideanDistance(queryPoint, segmentStart);
    
else
    param = dotProd / lenSq;
    
    if ( param < 0.00001 || param > 1)
        
        % closest point does not fall within the line segment, take the shorter distance to an endpoint
        distStart = TS_euclideanDistance(queryPoint, segmentStart);
        distEnd = TS_euclideanDistance(queryPoint, segmentEnd);
        
        if ( distStart > distEnd)
            pt2segDist = distEnd;
        else
            pt2segDist = distStart;
        end
    else
        xProj = segmentStart(1) + param * C;
        yProj = segmentStart(2) + param * D;
        pt2segDist = TS_euclideanDistance(queryPoint, [xProj yProj]);
    end
end
% fprintf('%f  ', pt2segDist);
end
        