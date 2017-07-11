function intersect = TS_circleLineIntersection(centerPoint,linePoint1,linePoint2,radius)

%     ************
%
%     Description
%     ----------
%     Find the intersections between the circle with defined radius and centerPoint
%     and the line delimited by points linePoint1 and linePoint2.
%     It is supposed here that the intersection between them exists. If no,  error
% 
%     Parameters
%     ----------
%     param centerPoint : 1 x dim, centre's abscissa and ordinate of the circle
%     param radius      : 1 x 1,   radius of the circle
%     param linePoint1  : 1 x dim, abscissa and ordinate of the first point of the line
%     param linePoint2  : 1 x dim, abscissa and ordinate of the second point of the line
%
%     Returns
%     -------
%     intersect         : 2 x dim, Coordinate of the two intersections.
%     If the two intersections are the same, that's means that the line is a tangent of the circle.
%
%     Other m-files required    : none
%     Subfunctions              : none
%     MAT-files required        : none
%     
%     References
%     ----------
%     link = http://math.stackexchange.com/questions/228841/how-do-i-calculate-the-intersections-of-a-straight-line-and-a-circle
%     http://mathworld.wolfram.com/Circle-LineIntersection.html
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


xCenter = centerPoint(1);
yCenter = centerPoint(2);
lp1x = linePoint1(1);
lp1y = linePoint1(2);
lp2x = linePoint2(1);
lp2y = linePoint2(2);

if (lp2x == lp1x)
    rac = sqrt((radius * radius) - ((lp1x - xCenter) * (lp1x - xCenter)));
    y1 = yCenter + rac;
    y2 = yCenter - rac;
    intersect = [[lp1x,y1];[lp1x,y2]];
else
    m = (lp2y - lp1y) / (lp2x - lp1x);
    b = lp2y - m * lp2x;
    A = m * m + 1;
    B = 2 * (m * b - m * yCenter - xCenter);
    C = (yCenter * yCenter) - (radius * radius) + (xCenter * xCenter) - (2 * b * yCenter) + (b * b);
    delta = B * B - 4 * A * C;
    
    if ( delta <= 0)
        x = -B / (2 * A);
        y = m * x + b;
        intersect = [[x,y];[x,y]];
    elseif( delta > 0)
        sDelta = sqrt(delta);
        x1 = (-B + sDelta) / (2 * A);
        y1 = m * x1 + b;
        x2 = (-B - sDelta) / (2 * A);
        y2 = m * x2 + b;
        intersect = [[x1,y1];[x2,y2]];
    else
        error('The intersection between circle and line is supposed to exist')
        
    end
end

end
        