function i1 = sf_slice2panel(img, xyzmm, transform, vdims)
% to voxel space of image
vixyz = (transform*img.vol.mat) \ xyzmm;
% raw data 
if mars_struct('isthere', img.vol, 'imgdata')
  V = img.vol.imgdata;
else
  V = img.vol;
end
i1 = spm_sample_vol(V,vixyz(1,:),vixyz(2,:),vixyz(3,:), ...
            [img.hold img.background]);
if mars_struct('isthere', img, 'func')
  eval(img.func);
end
% transpose to reverse X and Y for figure
i1 = reshape(i1, vdims(1:2))';
end