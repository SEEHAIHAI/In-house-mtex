
%% HCP结构GND计算，根据官网内容，改写了滑移系定义，以及修正了偶有数据不匹配的错误
% 任何情况下，请勿改动 ebsd('indexed') 部分内容
% Matlab版本2024a，Mtex版本为5.11.2

%% 不考虑未解析点
cs = ebsd('indexed').CS;
ebsd = ebsd('indexed').gridify;
%% Define the slip system and dislocation system
% 想计算的滑移系，参考附图里的定义，直接更改basal、prismatic和pyramidal即可。
% 基面<a>位错
sSbasal = slipSystem.basal(cs);
sSbasal = sSbasal.symmetrise('antipodal');
% 柱面<a>位错
sSPrismatic = slipSystem.prismaticA(cs);
sSPrismatic = sSPrismatic.symmetrise('antipodal');
% 锥面1阶<c+a>位错
sSPyramidal1 = slipSystem.pyramidalCA(cs);
sSPyramidal1 = sSPyramidal1.symmetrise('antipodal');
% 锥面2阶<c+a>位错
sSPyramidal2 = slipSystem.pyramidal2CA(cs);
sSPyramidal2 = sSPyramidal2.symmetrise('antipodal');

%% 基面<a>位错的位错密度计算
% compute the curvature tensor
kappa = ebsd.curvature;
alpha=kappa.dislocationDensity;
dSBasal = dislocationSystem(sSbasal);

dShcp = dSBasal;
dShcp(dShcp.isEdge).u = 1;
dShcp(dShcp.isScrew).u = 1 - 0.3;
dSRot = ebsd.orientations * dShcp;
[rho,factor] = fitDislocationSystems(kappa,dSRot);   % 会提示Optimal solution found. 然后开始最小化能量，等待即可

alpha = sum(dSRot.tensor .* rho,2);
alpha.opt.unit = '1/um';
kappa=alpha.curvature;

GND_bas=factor*sum(abs(rho .* dSRot.u),2);     
GND_bas=fillmissing(GND_bas,"linear");
GND_bas=filloutliers(GND_bas,"nearest",'median');

% 将取向数据转换为一维数组
% 将取向数据转换为一维数组
orientations_1d = ebsd.orientations(:);

% 找到取向数据中为 NaN 值的位置
nan_mask = isnan(orientations_1d);

% 将对应位置的 GND_bas 值设置为 NaN
GND_bas(nan_mask) = NaN;
% %% 柱面<a>位错的位错密度计算
% % compute the curvature tensor
% kappa = ebsd.curvature;
% alpha=kappa.dislocationDensity;
% dSPris = dislocationSystem(sSPrismatic);
% 
% dShcp = dSPris;
% dShcp(dShcp.isEdge).u = 1;
% dShcp(dShcp.isScrew).u = 1 - 0.3;
% dSRot = ebsd.orientations * dShcp;
% [rho,factor] = fitDislocationSystems(kappa,dSRot);    % 会提示Optimal solution found. 然后开始最小化能量，等待即可
% 
% alpha = sum(dSRot.tensor .* rho,2);
% alpha.opt.unit = '1/um';
% kappa=alpha.curvature;
% 
% GND_pri=factor*sum(abs(rho .* dSRot.u),2);
% GND_pri=fillmissing(GND_pri,"linear");
% GND_pri=filloutliers(GND_pri,"nearest",'median');

%% 1阶锥面<c+a>位错的位错密度计算
% % compute the curvature tensor
% kappa = ebsd.curvature;
% alpha=kappa.dislocationDensity;
% dSPyra = dislocationSystem(sSPyramidal1);
% 
% dShcp = dSPyra;
% dShcp(dShcp.isEdge).u = 1;
% dShcp(dShcp.isScrew).u = 1 - 0.3;
% dSRot = ebsd.orientations * dShcp;
% [rho,factor] = fitDislocationSystems(kappa,dSRot);   % 会提示Optimal solution found. 然后开始最小化能量，等待即可
% 
% alpha = sum(dSRot.tensor .* rho,2);
% alpha.opt.unit = '1/um';
% kappa=alpha.curvature;
% 
% GND_pyr1=factor*sum(abs(rho .* dSRot.u),2);
% GND_pyr1=fillmissing(GND_pyr1,"linear");
% GND_pyr1=filloutliers(GND_pyr1,"nearest",'median');
% %% 2阶锥面<c+a>位错的位错密度计算
% % compute the curvature tensor
% kappa = ebsd.curvature;
% alpha=kappa.dislocationDensity;
% dSPyra = dislocationSystem(sSPyramidal2);
% 
% dShcp = dSPyra;
% dShcp(dShcp.isEdge).u = 1;
% dShcp(dShcp.isScrew).u = 1 - 0.3;
% dSRot = ebsd.orientations * dShcp;
% [rho,factor] = fitDislocationSystems(kappa,dSRot);   % 会提示Optimal solution found. 然后开始最小化能量，等待即可
% 
% alpha = sum(dSRot.tensor .* rho,2);
% alpha.opt.unit = '1/um';
% kappa=alpha.curvature;
% 
% GND_pyr2=factor*sum(abs(rho .* dSRot.u),2);
% GND_pyr2=fillmissing(GND_pyr2,"linear");
% GND_pyr2=filloutliers(GND_pyr2,"nearest",'median');

%%
% GND_total=GND_bas+GND_pri+GND_pyr1+GND_pyr2;
% GND_mean=mean(GND_total(:));
% fprintf('GND的平均值为%e，单位1/m^2\n',GND_mean)
GND_bas_mean=mean(GND_bas(:));
% GND_pri_mean=mean(GND_pri(:));
% GND_pyr1_mean=mean(GND_pyr1(:));
% GND_pyr2_mean=mean(GND_pyr2(:));
fileID = fopen('density.txt', 'w');
fprintf(fileID, 'GND_bas_mean = %f  1/m^2\n', GND_bas_mean);
% fprintf(fileID, 'GND_pri_mean = %f  1/m^2\n', GND_pri_mean);
% fprintf(fileID, 'GND_pyr1_mean = %f 1/m^2\n', GND_pyr1_mean);
% fprintf(fileID, 'GND_pyr2_mean = %f 1/m^2\n', GND_pyr2_mean);
% %%
figure  % 画单个位错类性的GND图
% 也可以把所有位错类性的GND都算出来，然后相加并求出平均值，最后再画图，并且根据各平均值占总平均值百分比，算出贡献比例
plot(ebsd,GND_bas ,'micronbar','off')
mtexColorMap(flipud(parulaColorMap))   % 不喜欢hot的colormap的话，自行更换其他样式。flipud函数是为了保证颜色重的一段对应较大值，颜色轻的对应较小值。据实际情况更改。
hcb = mtexColorbar;
set(gca, 'colorscale', 'log');
% 设置颜色条的刻度为对数刻度
set(gca,'CLim',[2.7e12 1.6e14]);
hold on
plot(grains.boundary,'linewidth',2)
hold off

