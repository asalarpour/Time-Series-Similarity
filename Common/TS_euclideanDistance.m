function euclDist = TS_euclideanDistance(input1, input2)

%     ************
%
%     Description
%     ----------
%     Compute L2-norm between point set input1 and input2.
% 
%     Parameters
%     ----------
%     param input1  :   m x dim, input 1 matrix with the length of m
%     param input2  :   m x dim, input 2 matrix with the length of m
%
%     Returns
%     -------
%     euclDist 		: 	m x 1, matrix of L2-norm between input1 and input2
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
%     December 2016 : Last revision: 27-Jan-2017
%     
%     ************

euclDist = sqrt(sum(((input1 - input2) .^2), 2));