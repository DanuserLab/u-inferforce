function [idBestMatch]=findBasisClassInTbl(basisClassTbl,basisClassIn,xrange,yrange,meshPtsFwdSol)


%**************************************************************************
% 1) Try to find the basis function in the tablebase:                     *
%**************************************************************************
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
%
% This file is part of TFM_Package.
% 
% TFM_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TFM_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TFM_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
[foundClass]=findBasisClass(basisClassTbl,basisClassIn.neighPos);

idUsable=[];
smplDnst=[];
idBestMatch=[];
% Check in all classes if the range and the meshPtsFwdSol fits:
for idClass=foundClass'
    % read out the important values:
    meshPtsFwdSolTbl=basisClassTbl(idClass).uSol.meshPtsFwdSol;
    xrangeTbl=basisClassTbl(idClass).uSol.xrange;
    yrangeTbl=basisClassTbl(idClass).uSol.yrange;
    
    check1 = (xrangeTbl(1)<=xrange(1) && xrangeTbl(2)>=xrange(2));
    check2 = (yrangeTbl(1)<=yrange(1) && yrangeTbl(2)>=yrange(2));
    check3 = (meshPtsFwdSolTbl>=meshPtsFwdSol);
    
    if check1 && check2 && check3
        % then we can use this tabled solution for this basis class
        idUsable=horzcat(idUsable,idClass);
        
        % store the sample density to select the best one:
        smplDnst=horzcat(smplDnst,meshPtsFwdSolTbl/((xrangeTbl(2)-xrangeTbl(1))*(yrangeTbl(2)-yrangeTbl(1))));
    end
end
if ~isempty(idUsable)
    % The best stored solution is the one with the highes smplDensity:
    [~,posBestMatch]=max(smplDnst);
    idBestMatch=idUsable(posBestMatch);
end