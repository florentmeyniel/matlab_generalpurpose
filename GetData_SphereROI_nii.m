function [y, XYZmm, XYZvx] = GetData_SphereROI_nii(fname, Center_mm, radius)
% Get the data from a nii image, within a spheric ROI.
%
% [y, XYZmm] = GetData_SphereROI_nii(fname, Center_mm, radius)
% Input:
%   fname: full path name of the nii image
%   Center_mm: center of the ROI, in millimeter
%   radius: radius of the ROI, in millimeter
%
% Output:
%   y: data from the ROI
%   XYZmm: coordinate of the ROI voxels, in millimeter
%   XYZvx: indices of voxels in the search volume

% Check input
% =========================================================================
if ~exist(fname)
    ind_comma = strfind(fname, ',');
    if isempty(ind_comma)
        error('cannot find the nii file %s')
    else
    tmpfname = fname(1:ind_comma-1);
    if ~exist(tmpfname)
        error('cannot find the nii file %s')
    end
    end
end

% Extract data
% =========================================================================

% Specify the ROI structure for spm_ROI.m
xY.def  = 'sphere';
xY.xyz  = Center_mm(:);
xY.spec = radius;

% Get indices of voxels of the ROI in the nii image
[~, XYZmm, ind] = spm_ROI(xY, fname);

% get the coordinates of all voxels in the nii image
% (inspired from l. 123 to 138 of spm_ROI)
hdr     = spm_vol(fname);
[R,C,P] = ndgrid(1:hdr.dim(1),1:hdr.dim(2),1:hdr.dim(3));
RCP     = [R(:)';C(:)';P(:)'];

% get the coordinates, in voxels, of the ROI
XYZvx = RCP(:,ind);

% Get data from the ROI
y = spm_get_data(fname, XYZvx);
