function  [R] = EnsureGategory(S)
% 		"""Built the corresponding relationship between coefficients && its magnitude categories
%
% 		Args
% 			S Given DC DIFF or AC coefficient
%
% 		Returns
% 			SSSS The category SSSS with which size S could be coded
%
% 		Notes
% 				Range of SSSS of DC DIFF is (0, 12) && AC coefficient is (1, 11)
% 		"""
S = abs(S);
R = 0;
while 1
    if S == 0
        break;
    else
        S = bitshift(S, -1);
        R = R+1;
    end
end
end