function R2=matrixR2(A,B);
% function []=matrixR2(A,B);
%
%  Given two matrices, A and B, of the same dimensions, compute an R2-like
%  measure of how similar they are:
%
%      R2=1-RSS/TSS, where RSS = sum of squared deviations (element by element)
%                      and TSS = total sum of squares of each element
%
%    Like SUM((a_ij-b_ij)^2) / SUM(a_ij ^2+b_ij^2)
%
%  Note, we first take out the mean value for each matrix.

a=demean(vector(A)); b=demean(vector(B));  % Stack each matrix into a vector
RSS=sum((a-b).^2);
TSS=sum(a.^2+b.^2);
R2=1-RSS/TSS;