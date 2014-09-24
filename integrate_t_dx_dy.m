function [ t_d, dx_d, dy_d ] = integrate_t_dx_dy( t, dx, dy, integrateN )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

t_d = t( 1:integrateN:end );

dx_d = squeeze(sum(reshape(dx, [integrateN, size(dx,2)/integrateN])));
dy_d = squeeze(sum(reshape(dy, [integrateN, size(dy,2)/integrateN])));

end