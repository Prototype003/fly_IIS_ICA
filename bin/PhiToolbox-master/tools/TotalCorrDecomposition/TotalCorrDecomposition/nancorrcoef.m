function [ corr,p ] = nancorrcoef(data)
%欠損値(nan)を除いて相関係数を計算する
%function [ corr,p ] = nancorrcoef(data)
nonnan=find(~isnan(sum(data,2)));
[corr,p]=corrcoef(data(nonnan,:));