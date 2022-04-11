%{
-----------------------------------------------------------------
                 THE OBJECTIVE FUNCTION CODE OF 
   CONTRAST STRETCHING BASED PANSHARPENING (CSP) IN MATLAB 2022A
-----------------------------------------------------------------

load mydata
rng(100);  % Initialize mergene-twister pseudo random number generator for WDE
% algo_wde(ObjFun,mydata, PATTERNSIZE , DIM , LOW, UP , MAXITERATION , SEED ) ;
algo_wde('ObjFun_CSP',mydata, 5 , 7*mydata.k ,  0.00 , 255 , 2000 , 100 ) ;
[ObjFunVal, PI, solution ] = ObjFun_CSP( outwde.globalminimizer , mydata ) ;
montage({ uint8(mydata.MSI) uint8(PI) },'size',[1 2]), shg

%}                              

function [ ObjFunVal, PI, sol ] = ObjFun_CSP( pattern , mydata )

N = size(pattern,1) ;          % size of pattern-matrix
ObjFunVal = rand(N,1) ;        % Reserved memory for objective function values of patterns to speed-up the Matlab

% Image Data
MSI = double(mydata.MSI) ;     % multispectral image with low spectral resolution
PI = 0 * MSI;                  % Reserved memory for PI to speed-up the Matlab
PAN = double(mydata.PAN) ;     % panchromatic image with high spatial resolution
H = mydata.k ;

% Get monochromatics image band from MSI
[ MSI1, MSI2 , MSI3 ] = deal( MSI(:,:,1), MSI(:,:,2), MSI(:,:,3) ) ;

for i = 1 : N

    % --------------------------------------------------------------------------------------------
    % Sequantial Evaluation of Pattern Vectors
    NODES = sort( reshape( round( pattern(i,:) ) , 7 , H ) , 2 ) ;
    % --------------------------------------------------------------------------------------------
    % Nodes of Linear Contrast Stretching 
    eta1 = NODES ( 1 , : ) ; % nodes for LCS(MSI_BAND1,eta1)
    eta2 = NODES ( 2 , : ) ; % nodes for LCS(MSI_BAND1,eta2)
    eta3 = NODES ( 3 , : ) ; % nodes for LCS(MSI_BAND2,eta3)
    eta4 = NODES ( 4 , : ) ; % nodes for LCS(MSI_BAND2,eta4)
    eta5 = NODES ( 5 , : ) ; % nodes for LCS(MSI_BAND3,eta5)
    eta6 = NODES ( 6 , : ) ; % nodes for LCS(MSI_BAND3,eta6)    
    eta7 = NODES ( 7 , : ) ; % nodes for LCS(PAN,eta7)
    % --------------------------------------------------------------------------------------------
    % Computation of Linear Contrast Stretched Images
    lambda1    = LCS( MSI1 , eta1 ) ;     % 1st version of MSI BAND1
	lambda2    = LCS( MSI1 , eta2 ) ;     % 2nd version of MSI BAND1
	
    lambda3    = LCS( MSI2 , eta3 ) ;     % 1st version of MSI BAND2
    lambda4    = LCS( MSI2 , eta4 ) ;     % 2nd version of MSI BAND2
	
	lambda5    = LCS( MSI3 , eta5 ) ;     % 1st version of MSI BAND3
    lambda6    = LCS( MSI3 , eta6 ) ;     % 2nd version of MSI BAND3
    
    lambda7 = LCS( PAN , eta7 ) ;         % Linear Contrast Stretching of PAN
    % --------------------------------------------------------------------------------------------
    % Generation of initPI,    
    w1  = ( lambda1 ./  lambda2 ) ;
    w2  = ( lambda3 ./  lambda4 ) ;
    w3  = ( lambda5 ./  lambda6 ) ;
    % --------------------------------------------------------------------------------------------
    InitPI = uint8(  cat( 3 , w1 .* lambda7 , w2 .* lambda7 , w3 .* lambda7 )  ) ;
    % --------------------------------------------------------------------------------------------
    % Color Transfer Operation on initPI from MSI to generate PI
    for j = 1:3 
        img = double( InitPI(:,:,j) ); 
        imj = MSI(:,:,j) ; 
        img = ( ( img - mean2(img) ) / std2( img ) ) * std2(imj) + mean2(imj); 
        PI(:,:,j) = uint8(img) ; 
    end
    % --------------------------------------------------------------------------------------------
    % Objective Function Value, Root Mean Squared Error (RMSE)
    ObjFunVal(i) = sqrt( immse( uint8(MSI) , uint8(PI) ) ) ; 
    % --------------------------------------------------------------------------------------------
    sol.lambda1=lambda1; sol.lambda2=lambda2;
    sol.lambda3=lambda3; sol.lambda4=lambda4;
    sol.lambda5=lambda5; sol.lambda6=lambda6;
    sol.C=lambda7;
    sol.InitPI=InitPI;
    sol.PI=PI;
    sol.NODES=NODES;

end






