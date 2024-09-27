
%% HCP结构GND计算，根据官网内容，改写了滑移系定义，以及修正了偶有数据不匹配的错误
% 任何情况下，请勿改动 ebsd('indexed') 部分内容
% 用户可改动部分在滑移系定义那部分
%% 导入EBSD数据
close all;
clc;
clear all;
warning('off')
%
CS = {... 
  'notIndexed',...
  crystalSymmetry('6/mmm', [3.2 3.2 5.2], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Magnesium', 'color', [0.53 0.81 0.98])};
% 画图设定，当前与C5软件一致
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outtoPlane');
% % 画图设定，当前与EDAX软件一致
% setMTEXpref('xAxisDirection','south');
% setMTEXpref('zAxisDirection','outtoPlane');
%% 定义文件名称
% 加载路径
pname = 'D:\MatlabCode\HPC_GND计算（含裂纹）已验证\';   % 此路径使用者需根据自身情况更改
%% 
% 文件名
fname = [pname '3-6-Subset 1.ctf'];    % 此文件名根据自身情况更改
ebsd = EBSD.load(fname,CS,'interface','ctf','convertSpatial2EulerReferenceFrame');


[grains,ebsd.grainId] = calcGrains(ebsd);
grains = smooth(grains);

% % remove very small grains
% ebsd(grains(grains.grainSize<=0.15)) = [];
% 
% % and recompute grains
% [grains,ebsd.grainId] = calcGrains(ebsd('indexed'));
% 
% % smooth the grains a bit
% grains = smooth(grains,4);
%% 定义color key
ipfKey = ipfHSVKey(ebsd('indexed'));
ipfKey.inversePoleFigureDirection = zvector;
% %% 构建晶粒组织
% [grains,ebsd.grainId,ebsd.mis2mean]=calcGrains(ebsd('indexed'),'threshold',[2*degree, 5*degree]);    % 小角度晶界的范围，参考官网Grain boundary部分
% % 去除小晶粒
% ebsd(grains(grains.grainSize<=2)) = [];
% % 重新构建晶粒组织
% [grains,ebsd.grainId,ebsd.mis2mean]=calcGrains(ebsd('indexed'),'threshold',[2*degree, 5*degree]);    % 小角度晶界的范围，参考官网Grain boundary部分
% % 平滑晶界
% grains = smooth(grains,2);

% 数据降噪，但我本人很少用，因为我的原始数据在AztecCrystal中已经是降噪过的了
% ipfKey = axisAngleColorKey(ebsd('indexed'));
% ipfKey.oriRef = grains(ebsd('indexed').grainId).meanOrientation;
% F = halfQuadraticFilter;
% F.alpha=0.25;
% ebsd = smooth(ebsd('indexed'),F,'fill',grains);
% ipfKey.oriRef = grains(ebsd('indexed').grainId).meanOrientation;
% [grains,ebsd.grainId,ebsd.mis2mean]=calcGrains(ebsd('indexed'),'threshold',[2*degree, 20*degree]);    % 小角度晶界的范围，参考官网Grain boundary部分
% grains = smooth(grains,2);

%% EBSD map
h1=figure(1);
plot(ebsd,ebsd.orientations,'micronbar','on')
hold on
plot(grains.boundary('indexed'),'linewidth',2)
% plot(grains.boundary('index'),'linewidth',0.8)
hold off
% %% 如果选择一个稍微小一点的区域，就人工选择一个满意的区域
% 
% % 在大面积EBSD图上选取感兴趣区域，鼠标右键点击矩形区域的两个对角顶点
% figure(h1)
% [x,y]=ginput(2);                           % 获取两个点的像素坐标
% x_min=min(x); x_max=max(x);
% y_min=min(y); y_max=max(y);
% width=x_max-x_min;
% height=y_max-y_min;
% position=[x_min,y_min,width,height];
% rectangle('position',position,'edgecolor','r','linewidth',1)
% hold off
% 
% condition=inpolygon(ebsd,position);
% 
% ebsd=ebsd(condition);                      % 重绘选定的区域

% % 对数据进行降噪
% F = halfQuadraticFilter;
% F.alpha=0.25;
% ebsd = smooth(ebsd,F,'fill',grains);
%% 选定小区域后，再次重新构建晶粒组织
% [grains,ebsd.grainId,ebsd.mis2mean]=calcGrains(ebsd,'threshold',[2*degree, 10*degree],'boundary','tight');    % 晶界角度的范围，参考官网Grain boundary部分
% grains = smooth(grains,2);
