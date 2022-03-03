%%------NAME Script 1
%---The name of the matlab script is: generate_a_single_crater.m -%%
%---The present script shows how to design a Synthetic Digital Relief Model (SDRM) of a single crater-%%
%---The most important characteristic employed in the design was the crater morphology, -%%
%---Using the  parameters present in the following sentences, a single crater with either a flat bottom, or a central peak cam be drawn numerically -%%
%---The program is based on two previous codes designed for emulating the propagation of a wave in water-%%
%---the first is https://forum.lawebdefisica.com/forum/el-aula/m%C3%A9todos-inform%C3%A1ticos/24017-ecuacion-de-ondas; -%%
%---the second is https://www.fisicalab.com/apartado/ecuaciones-graficas-mas.-%%
%---The input parameters considered in this program such are: wave length (L), number of rows and columns (M, N),-%%
%---period (T), frequency (fr), oscillating time (t), wave speed (v), amplitude (A).-%%
%---In regards amplitude, this parameter is defined as A = C/sqrt(r), where C is a costant dependent on the energy provided by the medium -%%
%---and r is the distance to the focus. Therefore, in order to avoid the amplitude going to infinite, it is necessary to take a value -%%
%---of the amplitude for very close distances. I took C=1, so when the distance to the focus is less than 0.8 than 0.8, we will say that the Amplitude is 1.-%%
%---The wave radius is rr=((X).^2+(Y).^2).^0.5 providing the r value as r=fr*rr -%%
%---Therefore, changing the frequency it is possible to emulate different waves forms; thus is fr is 0.1 a flat bottomed crater is obtained.-%%
%---However, if fr is 2 a crater whith a central peak if as the parameters further presented are changed.-%%
%---Output files. Different possibilities of results are offered. It is possible to define the folder and path where the file will be saved. -%%
%---The first is a txt ascii file, with and without georeferenced information in the header; the second is an excel file.
%% SECTION 1.---DESCRIPTION OF THE VARIABLES INVOLVED--
% Wave number (k = 2*pi/ wavelength)
% The wavenumber is a frequency magnitude that indicates..
% ... the number of times a single wave vibrates in a unit distance.
% fr=1/T where fr is frequency; T is the period
% long, wavelength
% v, is the wave speed
% the frequency has an inverse relationship with landa, wavelength
% f=v/ long
% L, size of the working square matrix, number of rows or columns
%% SECTION 2.---INPUT DATA--close
clear all
clc;
L=26;
M=L;
N=L;
fr= 2;
% long = input('Enter wavelength: ');
long = L;
% T = input('Enter the period of the wave: ');
T = 1;
Tiempo= 1*T; % Time it will be oscillating.
% phi = input('Enter initial phase: ');
phi = 2*pi*fr*T;

%% ---SECTION 2.1.---We create a grid.
x = linspace(-floor(M/2),floor(M/2),M);
y = linspace(-floor(N/2),floor(N/2),N);
[X,Y] = meshgrid(x,y);
% N = length (L);  %Numero de muestras de la señal= L
% Fs = N
% Wave radius
rr=((X).^2+(Y).^2).^0.5;
r=fr*rr;
%% ---SECTION 2.2.---Amplitude.
% the amplitude, A, is defined as A = C/sqrt(r), where C is a
% costant dependent on the energy provided by the medium and r is the
% distance to the focus. Therefore, in order to avoid the amplitude going to
% infinite, it is necessary to take a value of the amplitude for very close distances.
% I took C=1, so when the distance to the focus is less than 0.8
% than 0.8, we will say that the Amplitude is 1.
r0=0.8;
n = length(r);
N = length (L);  %Number of signal samples= L
Fs = N;
Ts = 1/Fs; % Sampling period (time step)
f= Fs.*(0:(L/2)-1)./L; %vector of true frequencies

%% ---SECTION 2.3.---Paing graphically the wave.
%By means of these two loops and the conditional we distinguish between the two cases
% of the distance to the focus to determine the Amplitude of the wave.
for t = 0:0.1:Tiempo
    for i = 1:n
        for j = 1:n
            if (r(i,j))>r0
                % A = 1./(r(i,j).^0.5);
                A=2;
                A1=1;
                A2=0.5;
                Z(i,j) = A.*(cos(2*pi.*((r(i,j)./long)-(t/T)))); % Wave equation.
            end
        end
    end
    surf(X,Y,Z,'Facecolor','blue','Edgecolor','none');
    title('Fig.1, Altimetry'),
    xlabel('x- altimetry'),ylabel('y- altimetry'),
    axis equal;%By means of this command we can equalise the sizes of the axes.
    camlight right; %add Shadow to better visualise the wave
    lighting phong;%add light to better visualise the wave
end

%% ---SECTION 3.---EXPORT DATA--
%
%%BLOCK IN CASE TESTING THE SCRIPT
folder = 'Results';
if ~exist(folder, 'dir')
    mkdir(folder);
end
%---SECTION 3.1.---Create a excel export file of the central row
baseFileName = 'DEMejepru.xlsx';
fullFileName = fullfile(folder, baseFileName);
writematrix(Z, fullFileName,'sheet','WriteMatrix','Range','A13:Z13');%13
%%%---SECTION 3.2.---Create a export file ascii raster without a georeferencing head
% Zx = fopen(sprintf('K:\Luna\manuscrito\DEMejepru.txt', Z), 'wt');
Zx = fopen(sprintf('DEMejepru.txt', Z), 'wt');
% Zx = fopen(sprintf(fullFileName, Z), 'wt');
fprintf(Zx, [repmat('%9.5f\t', 1, size(Z,2)) '\n'], Z');
fclose(Zx);
%---SECTION 3.3.---Create a export file ascii raster with a georeferencing head --
Zx = fopen(sprintf('DEMejepru_concord.txt', Z), 'wt'); 
% % %     the number of columns and row is written, the value of X,
% % %     and in the origin, the size of the cell and the value that is
% % %     to be assigned to cells without data.
fprintf(Zx, 'ncols\t%s\n',num2str(26));
fprintf(Zx, 'nrows\t%s\n',num2str(26));
fprintf(Zx, 'xllcenter\t%s\n',num2str(0));
fprintf(Zx, 'yllcenter\t%s\n',num2str(0));
fprintf(Zx, 'cellsize\t%s\n',num2str(1));
fprintf(Zx, 'NODATA_value\t%s\n',num2str(-9999));
fprintf(Zx, [repmat('%9.5f\t', 1, size(Z,2)) '\n'], Z');
fclose(Zx);
% END OF THE TEST BLOCK______________
%-----END------%
disp('END');
% END