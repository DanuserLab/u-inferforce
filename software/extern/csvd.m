function [U,s,V] = csvd(A,tst)
%CSVD Compact singular value decomposition.
%
% s = csvd(A)
% [U,s,V] = csvd(A)
% [U,s,V] = csvd(A,'full')
%
% Computes the compact form of the SVD of A:
%    A = U*diag(s)*V',
% where
%    U  is  m-by-min(m,n)
%    s  is  min(m,n)-by-1
%    V  is  n-by-min(m,n).
%
% If a second argument is present, the full U and V are returned.
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

% Per Christian Hansen, IMM, 06/22/93.

if (nargin==1)
  if (nargout > 1)
    [m,n] = size(A);
    if (m >= n)
      [U,s,V] = svd(full(A),0); s = diag(s);
    else
      [V,s,U] = svd(full(A)',0); s = diag(s);
    end
  else
    U = svd(full(A));
  end
else
  if (nargout > 1)
    [U,s,V] = svd(full(A)); s = diag(s);
  else
    U = svd(full(A));
  end
end