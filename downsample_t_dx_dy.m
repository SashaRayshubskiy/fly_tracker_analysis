function [ t_d, dx_d, dy_d ] = downsample_t_dx_dy( t, dx, dy, downsample )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

t_d = t( 1:downsample:end );

dx_d = squeeze(sum(reshape(dx, [downsample, size(dx,2)/downsample])));
dy_d = squeeze(sum(reshape(dy, [downsample, size(dy,2)/downsample])));

end