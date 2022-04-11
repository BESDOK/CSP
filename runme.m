load mydata
rng(100);  % Initialize mergene-twister pseudo random number generator for WDE

% algo_wde(ObjFunction , mydata , populationSize, dimensionOfProblem , low , up , MaxEpk , seed )

algo_wde('ObjFun_CSP',mydata, 5 , 7*mydata.k ,  0.00 , 255.00 , 2000 , 100 ) ;

[ObjFunVal, PI ] = ObjFun_CSP( outwde.globalminimizer , mydata ) ;

montage({ uint8(mydata.MSI) uint8(PI) },'size',[1 2]), shg

