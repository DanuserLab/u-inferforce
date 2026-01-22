function [displField_reg] = perfRotReg(displField,checkTransform)
if nargin<2 || isempty(checkTransform)
    checkTransform=0;
end

for frame=1:length(displField)
%     displField_reg(frame).par   = displField(frame).par;
%
% Copyright (C) 2026, Danuser Lab - UTSouthwestern 
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
    checkVec=(displField(frame).pos(:,1)>50 & displField(frame).pos(:,1)<200 & displField(frame).pos(:,2)>25 & displField(frame).pos(:,2)<700) | (displField(frame).pos(:,1)>1100 & displField(frame).pos(:,2)>25);
    
    X1=displField(frame).pos(checkVec,:);
    X2=X1+displField(frame).vec(checkVec,:);
    % extend to 3D:
    X1(:,3)=0;
    X2(:,3)=0;
    [T R]  = computeICP(X1, X2, 10^5,10^-2);
    
    % transform the whole field using this transformation:
    Xref=displField(frame).pos;
    Xdef=displField(frame).vec+Xref;
    % extend to 3D:
    Xref(:,3)=0;
    Xdef(:,3)=0;
    % apply the transformation:
    Xdef_reg = (R*Xdef'+ repmat(T,1,size(Xdef,1)))';
    
    % store the values:
    displField_reg(frame).pos     = Xref(:,1:2); % This is of course unchanged!
    displField_reg(frame).vec     = Xdef_reg(:,1:2)-Xref(:,1:2);
    displField_reg(frame).par.R   = R;
    displField_reg(frame).par.T   = T;
    
    if checkTransform
        figure()
        quiver(displField_reg(frame).pos(:,1),displField_reg(frame).pos(:,2),displField_reg(frame).vec(:,1),displField_reg(frame).vec(:,2));

        % calculate the residuals of the original displacement field:
        figure()
        subplot(2,2,1)
        quiver(displField(frame).pos(checkVec,1),displField(frame).pos(checkVec,2),displField(frame).vec(checkVec,1),displField(frame).vec(checkVec,2));
        subplot(2,2,2)
        magVec=sqrt(sum(displField(frame).vec(checkVec,:).^2,2));
        hist(magVec,linspace(0,max(magVec),100));
        sum(magVec)

        subplot(2,2,3)
        % apply the transformation:
        X2_reg = (R*X2'+ repmat(T,1,size(X2,1)))';
        % calculate the residuals of the transformed displacement field:
        residualsVec=X2_reg(:,1:2)-X1(:,1:2);
        % position is unchanged of course:
        residualsPos=X1(:,1:2);    
        quiver(residualsPos(:,1),residualsPos(:,2),residualsVec(:,1),residualsVec(:,2));
        subplot(2,2,4)
        magVec_reg=sqrt(sum(residualsVec.^2,2));
        hist(magVec_reg,linspace(0,max(vertcat(magVec,magVec_reg)),100))
        sum(magVec_reg)
    end
end

