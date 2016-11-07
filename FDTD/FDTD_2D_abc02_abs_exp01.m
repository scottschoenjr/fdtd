%% CD Arvanitis BWH, Harvard University 2012
close all;clear all;
%% load date with head properties
Sdir='C:\Users\sschoen3\Dropbox\FDTD-ASA-Scott\data\simulated';
domain=load([Sdir filesep 'head_2d_slc40_thr700_crop_res125.mat']);
%% Some user modifiable parameters
son=14;
bps='n';                            %Choose exitation pulse "y' is from bubble dynamics 'n' is guasian pulse
% trgt=1;                             %Select the location of the point sourec
f0=1.32e6;                        %Mean frequency of the excitation pulse in [Hz]
cw=1479;                            % speed of sound of water at 20 Celcius [m/sec]
rhow=1000;                          %Denstiy of Water [kg/m^3]
w = f0*2*pi;                        %Angular frequency
% lc = c0/f0;                       %Wavelength [m]
intconst=8;                         %data points to collect images for video
dx=domain.res*1e-3;                 %Length of nodes in x direction [m]
%% Initialise the p and u fields and set coefficients + a Courrant compliant dt
rho=flipdim(flipdim(domain.extrhotot,2)',2);
c=flipdim(flipdim(domain.extctot,2)',2);
alpha=flipdim(flipdim(domain.extatttot,2)',2);
alpha(find(alpha==2.9))=2.9;%7.5x is correction to match the shear viscosity of brain tissue (3.5 kg-1sec-1)found in literature (the 2x is for stability)
alpha(find(alpha==0.000288))=0.000288;%30x is correction to mathc the shear viscosity of water (0.9e-4 kg-1sec-1)found in literature
% mi=alpha.*rho.*(c.^3)/w^2;% Shear viscosity
delta=0.0*(f0/1e+6)*(8/3)*alpha.*(c)./(w^2);% convert attenuation to diffusivity [m^2/sec]: 20 times smaller mi was used

xDim =size(rho,2);                %number of nodes in x direction
yDim =size(rho,1);                %number of nodes in y direction

% c=cw*ones(yDim,xDim);
% c(:,xDim/2+xDim/8:xDim/2+xDim/8+xDim/6)=2*cw;
% rho= rhow*ones(yDim,xDim);
% rho(:,xDim/2+xDim/8:xDim/2+xDim/8+xDim/6)=2*rhow;
% delta=(f0/1e+6)*4.6e-15*ones(yDim,xDim);
% mi(:,xDim/2+xDim/8:xDim/2+xDim/8+xDim/6)=5*delta;

ratio=max(domain.ctot(:))/(1.5*f0)/dx

p = zeros(yDim,xDim,3);     %pressure scalar field initialisation
ux = zeros(yDim,xDim-1);    %velocity vector field initialisation
uy = zeros(yDim-1,xDim);    %velocity vector field initialisation

dt = 0.55*dx/(sqrt(ndims(p))*max(c(:))*(1+w^2*max(delta(:))^2/4)); %set time step small enough to satisfy courrant condition
coP = (dt/dx) * rho.* c.^2; %pressure scalar difference co-efficient

stopTime=1.8*xDim*dx/max(c(:)); % Time to run the simulation 1=image size
t=0:dt:stopTime;
%% excitation pulse
% excit_loc=domain.targets/(1e+3*dx);%multiple sources locations
% excit_loc=flipdim(excit_loc,2);

[~, pixx]=min(abs(domain.xx-domain.targets(son,1)));
[~, pixz]=min(abs(domain.zz-domain.targets(son,2)));
excit_loc=[length(domain.xx)-pixx,length(domain.yy)-(pixz+length(domain.yy)-length(domain.zz))];
offset=2.15;
excit_b=0.0466*pulse(t,dt,f0,offset,bps);%4.25* bubbles, 0.0466 gaussian
sum(abs(fft(excit_b))/length(t))
aedata=zeros(yDim,length(t));
%% coefficients for 2nd order ABC
pOldup1=zeros(3,xDim);
pOldup2=zeros(3,xDim);
pOlddown1=zeros(3,xDim);
pOlddown2=zeros(3,xDim);
pOldright1=zeros(yDim,3);
pOldright2=zeros(yDim,3);
pOldleft1=zeros(yDim,3);
pOldleft2=zeros(yDim,3);

cr=2;

temp1=(dt/dx)*(c(1,:)*cr);
temp2=1.0./temp1+2.0+temp1;
abcCoefup=zeros(3,xDim);
abcCoefup(1,:)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefup(2,:)=-2.0*(temp1 - 1.0./temp1)./temp2;
abcCoefup(3,:)= 4.0*(temp1+1.0./temp1)./temp2;

temp1=(dt/dx)*(c(end,:)*cr);
temp2=1.0./temp1+2.0+temp1;
abcCoefdown=zeros(3,xDim);
abcCoefdown(1,:)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefdown(2,:)=-2.0*(temp1-1.0./temp1)./temp2;
abcCoefdown(3,:)= 4.0*(temp1+1.0./temp1)./temp2;

temp1=(dt/dx)*(c(:,end)*cr);
temp2=1.0./temp1+2.0+temp1;
abcCoefright=zeros(yDim,3);
abcCoefright(:,1)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefright(:,2)=-2.0*(temp1 - 1.0./temp1)./temp2;
abcCoefright(:,3)= 4.0*(temp1+1.0./temp1)./temp2;

temp1=(dt/dx)*(c(:,1)*cr);
temp2=1.0./temp1+2.0+temp1;
abcCoefleft=zeros(yDim,3);
abcCoefleft(:,1)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefleft(:,2)=-2.0*(temp1 - 1.0./temp1)./temp2;
abcCoefleft(:,3)= 4.0*(temp1+1.0./temp1)./temp2;

timepnts=1:5*intconst:length(t)-1;
timestmps=zeros(length(timepnts),xDim);
% fig=figure;
t(end)*1e+6
num=0.0;
ampl=1.0;
for n = 1:length(t)-2;
    
    if n<=length(t)/4
        %point source one up in the imge
        p(excit_loc(1,1)-round(ratio*num),excit_loc(1,2),1) = ampl*excit_b(n);%insert hard source
        p(excit_loc(1,1)-round(ratio*num),excit_loc(1,2),2) = ampl*excit_b(n+1);%insert hard source
        
%         p(excit_loc(trgt,1)+round(ratio*num),excit_loc(trgt,2),1) = ampl*excit_b(n);%insert hard source
%         p(excit_loc(trgt,1)+round(ratio*num),excit_loc(trgt,2),2) = ampl*excit_b(n+1);%insert hard source
%         
% %       %point source two down in the image
%         p(excit_loc(trgt,1),excit_loc(trgt,2)-round(ratio*num),1) = ampl*excit_b(n);%insert hard source
%         p(excit_loc(trgt,1),excit_loc(trgt,2)-round(ratio*num),2) = ampl*excit_b(n+1);%insert hard source  
%         
%         p(excit_loc(trgt,1),excit_loc(trgt,2)+round(ratio*num),1) = ampl*excit_b(n);%insert hard source
%         p(excit_loc(trgt,1),excit_loc(trgt,2)+round(ratio*num),2) = ampl*excit_b(n+1);%insert hard source  
    end
    uy = uy-(2*(dt/dx)./(rho(1:end-1,:)+rho(2:end,:))).*...
        ((1+((delta(1:end-1,:)+delta(2:end,:))./(2*dt))).*(p(2:end, :,2)- p(1:end-1,:,2))-((delta(1:end-1,:)+delta(2:end,:))./(2*dt)).*(p(2:end, :,1)- p(1:end-1,:,1)));%calculate velocity in y plane
    ux = ux-(2*(dt/dx)./(rho(:,1:end-1)+rho(:,2:end))).*...
        ((1+((delta(:,1:end-1)+delta(:,2:end))./(2*dt))).*(p(:,2:end,2)- p(:,1:end-1,2))-((delta(:,1:end-1)+delta(:,2:end))./(2*dt)).*(p(:,2:end,1)- p(:,1:end-1,1)));%calculate velocity in y plane
    p(:,:,3) = p(:,:,2) - coP.*([uy; zeros(1,xDim)]- [zeros(1,xDim); uy]+[ux zeros(yDim,1)]-[zeros(yDim,1) ux]);%calculate pressure from xy velocity
    %% 2nd order ABC
    %ABC for uper part of the domain
    p(1,1:end,3)=abcCoefup(1,:).*(p(3,1:end,3)+pOldup2(1,:))...
        +abcCoefup(2,:).*(pOldup1(1,:)+pOldup1(3,:)-p(2,1:end,3)-pOldup2(2,:))...
        +abcCoefup(3,:).*pOldup1(2,:)-pOldup2(3,:);
    
    %ABC for lower part of the domain
    p(end,1:end,3)=abcCoefdown(1,:).*(p(end-2,1:end,3)+pOlddown2(1,:))...
        +abcCoefdown(2,:).*(pOlddown1(1,:)+pOlddown1(3,:)-p(end-1,1:end,3)-pOlddown2(2,:))...
        +abcCoefdown(3,:).*pOlddown1(2,:)-pOlddown2(3,:);
    
    %ABC for right part of the domain
    p(1:end,1,3)=abcCoefright(:,1).*(p(1:end,3,3)+pOldright2(:,1))...
        +abcCoefright(:,2).*(pOldright1(:,1)+pOldright1(:,3)-p(1:end,2,3)-pOldright2(:,2))...
        +abcCoefright(:,3).*pOldright1(:,2)-pOldright2(:,3);
    
    %ABC for left part of the domain
    p(1:end,end,3)=abcCoefright(:,1).*(p(1:end,end-2,3)+pOldleft2(:,1))...
        +abcCoefright(:,2).*(pOldleft1(:,1)+pOldleft1(:,3)-p(1:end,end-1,3)-pOldleft2(:,2))...
        +abcCoefright(:,3).*pOldleft1(:,2)-pOldleft2(:,3);
    %update stored fields
    for kk=1:3
        pOldup2(kk,:)=pOldup1(kk,:);
        pOldup1(kk,:)=p(kk,1:end,3);
        
        pOlddown2(kk,:)=pOlddown1(kk,:);
        pOlddown1(kk,:)=p(end-kk+1,1:end,3);
        
        pOldright2(:,kk)=pOldright1(:,kk);
        pOldright1(:,kk)=p(1:end,kk,3);
        
        pOldleft2(:,kk)=pOldleft1(:,kk);
        pOldleft1(:,kk)=p(1:end,end-kk+1,3);
    end
    p(:,:,1) = p(:,:,2);
    p(:,:,2) = p(:,:,3);
    %% data visualizatioon
    %% data visualizatioon
    if isempty(find(n==1:5*intconst:length(t)-1, 1))
        a=1;
    else
        timestmps((n+5*intconst-1)/(5*intconst),:)=p(excit_loc(1,1),:,3);
    end
    
    %video with sound propagation on the ct images
    imagesc(0:dx:(xDim-1)*dx*1000,0:dx:(yDim-1)*dx*1000,flipdim(p(:,:,2)+c/max(max(c))/580,2))
%     caxis([0.5e-4 2e-3])
      set(gca,'XDir','Reverse')
    colormap(jet)
    axis image
    xlabel('Distance (mm)','FontWeight','bold')
    ylabel('Distance (mm)','FontWeight','bold')
    title(sprintf('Time = %.2f usec',n*dt*1e6),'Color',[0 0 0],'FontSize', 22); drawnow;
    %% Convert imges to tiff for making video
%                 if isempty(find(n==1:5*intconst:length(t)-1, 1))
%                     a=1;
%                 else
%                     currentframe = frame2im(getframe(fig)); % convert fig into image data
%                     %     currentframe = im2bw(currentframe,0.4); % make it binary to save space
%                     if n == 1
%                         imwrite(currentframe, 'attscabs_1bps_trgt01_int02.tif', 'WriteMode', 'overwrite');
%                     else
%                         imwrite(currentframe, 'attscabs_1bps_trgt01_int02.tif', 'WriteMode', 'append');
%                     end
%                     close all;
%                 end
%     
    %%  record and save rf data
    aedata(:,n)= p(:,end-5,2)';
    fprintf('\b\b\b%3.0f',100*n/length(t));
end;

%% waterfall plots at the point source line coordinates
figure
axes('FontSize',16)
simpleWaterfall(timestmps, 1, 1e4) % vertical offset = 1, scale factor = 1.9
xlabel('Space [spatial index]','FontSize',16,'FontWeight','bold')
ylabel('Time [frame number]','FontSize',16,'FontWeight','bold')
axis([0 xDim -1 length(timepnts)]);
% figure
% Waterfall((1:xDim)*dx*1e3,t(timepnts)*1e6,timestmps) % vertical offset = 1, scale factor = 1.9
% colormap(gray)
% view(0,65)
% % axes('FontSize',16)
% xlabel('Space (mm)','FontSize',16,'FontWeight','bold')
% ylabel('Time (usec)','FontSize',16,'FontWeight','bold')
% zlabel('Amplitude','FontSize',16,'FontWeight','bold')
% grid off
% axis([0 (xDim*dx-dx)*1e3 t(1) t(end-1)*1e6 0 max(max(timestmps))]);

%% save the simulation data
% folder1=cd;
% folder1=Sdir;
% folder2=['1320_aeattscabs_2DHead_res' num2str(domain.res*1000) 'um_1gps0075_sontrgt0' num2str(son)]
% savefile =[folder1 '\' folder2 '.mat'];
% save(savefile,'aedata','dx','t','xDim','yDim','domain','excit_loc','timestmps','ratio','num','f0','son');