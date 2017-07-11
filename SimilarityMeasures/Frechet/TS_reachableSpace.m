function[rep, LR, BR] = TS_reachableSpace(LF, BF)

%     ************
%
%     Description
%     ----------
%     Compute all the free space,that are reachable from the origin (freeSpace1{0,0},freeSpace2{0,0})
%     on the boundary of cells in the diagram for polygonal chains ts1 and ts2 and the given free spaces freeSpace1 and freeSpace2
%     freeSpace1{i,j} is the free space, reachable from the origin, of segment [ts1(i, :),ts1(i + 1, :)] from point ts2(j, :)
%     freeSpace1{i,j} is the free space, reachable from the origin, of segment [ts2(j, :),ts2(j + 1, :)] from point ts1(i, :)
% 
%     Parameters
%     ----------
%     LF               : cell, free spaces of segments of ts1 from points of ts2
%     LB               : cell, free spaces of segments of ts2 from points of ts1
%
%     Returns
%     -------
%     rep              : bool, return true if frechet distance is inf to eps
%     LR               : cell, is the free space reachable from the origin of segments of ts1 from points of ts2
%     BR               : cell, is the free space reachable from the origin of segments of ts1 from points of ts1
%
%     Other m-files required    : none
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



m = size(LF, 1) + 1;
n = size(BF, 2) + 1;

BR = zeros(m - 1, n - 1);
LR = zeros(m - 1, n - 1);

if ~(LF{1, 1}(1, 1) == 0 && ...
        BF{1, 1}(1, 1) == 0 && ...
        LF{m - 1, n}(1, 2) == 1 && ...
        BF{m, n - 1}(1, 2) == 1 )
    rep = false;

else
    LR(1, 1) = true;
    BR(1, 1) = true;
    for i = 2: m - 1
        if (LF{i, 1}(1, 1) ~= -1 && ...
                LF{i, 1}(1, 2) ~= -1 && ...
                LF{i - 1, 1}(1, 1) == 0 && ...
                LF{i - 1, 1}(1, 2) == 1)
            LR(i, 1) = true;
        else
            LR(i, 1) = false;
        end
    end
    for j  = 2: n - 1
        if (BF{1, j}(1, 1) ~= -1 && ...
                BF{1, j}(1, 2) ~= -1 && ...
                BF{1, j - 1}(1, 1) == 0 && ...
                BF{1, j - 1}(1, 2) == 1)

            BR(1, j) = true;
        else
            BR(1, j) = false;
        end
    end
    for i = 1: m - 1
        for j = 1: n - 1
            if LR(i, j) || BR(i, j)
                if (LF{i, j + 1}(1, 1) ~= -1 && ...
                        LF{i, j + 1}(1, 2) ~= -1)
                    
                    LR(i, j + 1) = true;
                    
                else
                    
                    LR(i, j + 1) = false;
                    
                end
                if (BF{i + 1, j}(1, 1) ~= -1 && ...
                        BF{i + 1, j}(1, 2) ~= -1)
                    
                    BR(i + 1, j) = true;
                    
                else
                    
                    BR(i + 1, j) = false;
                    
                end
            else
                
                LR(i, j + 1) = false;
                BR(i + 1, j) = false;
            end
        end
    end
    rep = BR(m - 1, n - 1) || LR(m - 1, n - 1);
end

end