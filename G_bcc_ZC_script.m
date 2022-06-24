function [GmBCC,TCBCC] = GbccPE(XC,T)

%XFE=VARIABLE
XMN=0.01724;
XSI=0.00535;
XMO=0.00157;%0.00157;
%XC=variable;
XFE=1-XMN-XSI-XMO-XC;
XVA=1-XC;
XN=0;


% YJs are site fractions defined as 
%for ferrite q=3 and p=1
p=1;
q=3;


Y1 = XFE / (1 - XC - XN);
Y2 = XMN / (1 - XC - XN);
Y3 = XSI / (1 - XC - XN);
Y4 = XMO / (1 - XC - XN);

YC = (p/q)*(XC/(1 -XC - XN));
YVA = 1 - YC;

Y12 = Y1 - Y2;
Y13 = Y1 - Y3;
Y14 = Y1 - Y4;
Y23 = Y2 - Y3;
Y24 = Y2 - Y4;
Y34 = Y3 - Y4;


P=1E5; %Pa


%FUNCTIONS
%the order (in definition) is imporatant. We onlu need GPFEXXX, but to get
%that we need all the other ones defined. 

%T and P must be defined
R=8.314; %J/Kmol
RT= R*T;


% BFEFCC   298.15 
BFEFCC = +1 + 3.25236341E-11*P+3.36607808E-16*T*P; % 6000 N
% BFEBCC   298.15 
BFEBCC = +1 + 2.80599565E-11*P+3.06481523E-16*T*P; % 6000 N
 
% DFEBCC   298.15 
DFEBCC = +1*log(BFEBCC );  %6000 N
%CFEFCC   298.15 
CFEFCC = +2.62285341E-11+2.71455808E-16*T; %6000 N
% AFEFCC   298.15 
AFEFCC = +7.3097E-05*T;% 6000 N
% DFEFCC   298.15 
DFEFCC = +1*log(BFEFCC ); % 6000 N
% CFEBCC   298.15 
CFEBCC = +2.20949565E-11+2.41329523E-16*T;% 6000 N
%AFEBCC   298.15 
AFEBCC = +2.3987E-05*T+1.2845E-08*T^2; %6000 N

VFEBCC = +7.042095E-06*exp(AFEBCC ); % 6000 N
% EFEBCC   298.15 
EFEBCC = +1*log(CFEBCC ); % 6000 N
% XFEBCC   298.15 
XFEBCC = +1*exp(.7874195*DFEBCC )-1; % 6000 N

%VFEFCC   298.15 
VFEFCC = +6.688726E-06*exp(AFEFCC ); %6000 N
% EFEFCC   298.15 
EFEFCC = +1*log(CFEFCC ); % 6000 N
% XFEFCC   298.15
XFEFCC = +1*exp(.8064454*DFEFCC )-1; %6000 N


% YFEBCC   298.15 
YFEBCC = +VFEBCC *exp(-EFEBCC ); % 6000 N
% ZFEBCC   298.15 
ZFEBCC = +1*log(XFEBCC ); % 6000 N

% YFEFCC   298.15 
YFEFCC = +VFEFCC *exp(-EFEFCC ); % 6000 N
% ZFEFCC   298.15 
ZFEFCC = +1*log(XFEFCC );% 6000 N


%GPFEBCC  298.15
GPFEBCC = +YFEBCC *exp(ZFEBCC ); %6000 N
%GPFEFCC  298.15 
GPFEFCC= +YFEFCC *exp(ZFEFCC ); %6000 N





YLNY= Y1*log(Y1) + Y2*log(Y2)+ Y3*log(Y3) + Y4*log(Y4);
YLNYCVA = YC*log(YC) + YVA*log(YVA);

 


%%%%%%%%%%G(BCC,X,C)%%%%%%%%%%


%  G(BCC_A2,FE:C;0)- 3 H298(GRAPHITE,C;0)-H298(BCC_A2,FE;0) 
%298.15<T< 1811.00: 
GF1= +271170.377 + 711.991*T - 96.4143*T*log(T) - 0.00581442*T^2 - 5.8927E-08 *T^3 + 7765159*T^(-1)-7.929E+08*T^(-2) + 3.6E+10*T^(-3);


%G(BCC_A2,MN:C;0)- 3 H298(GRAPHITE,C;0)-H298(CBCC_A12,MN;0) =
%298.15<T< 1519.00: 
GF2=-50220.603 + 672.249*T - 96.3582*T*log(T) - 0.00876458*T^2 + 7757627*T^(-1) -7.929E+08*T^(-2)+ 3.6E+10*T^(-3);

%G(BCC_A2,SI:C;0)- 3 H298(GRAPHITE,C;0)-H298(DIAMOND_A4,SI;0) =
%298.15<T< 1687.00: 
GF3= +308782.068 + 551.250259*T-95.7317533*T*log(T) - 0.003329804*T^2 - 3.552E-09*T^3 + 7864467*T^(-1) - 7.929E+08*T^(-2)+ 3.6E+10*T^(-3);

%G(BCC_A2,MO:C;0)- 3 H298(GRAPHITE,C;0)-H298(BCC_A2,MO;0) =
%298.15<T< 2896.00: 
GF4=+271148.375+569.1097*T-96.46414*T*log(T) - 0.004860296*T^2 + 5.66283E-07*T^3 + 7753612*T^(-1) -1.30927E-10*T^4 - 7.929E+08*T^(-2) + 3.6E+10*T^(-3);


%L(BCC_A2,FE,MN:C;0) 
LF12= +34052-23.467*T;

%L(BCC_A2,FE,SI:C;0)
LF13= +1000000-100*T;
%L(BCC_A2,FE,SI:C;1) 
LF13o1= -900000;
%L(BCC_A2,FE,MO:C;0) 
LF14= -1250000+667.7*T;

%L(BCC_A2,MN,SI:C;0) 
LF23=  0.0;




GBCC_XC= Y1*GF1 + Y2*GF2 + Y3*GF3 + Y4*GF4 + RT*YLNY + (Y1*Y2)*LF12 +  (Y1*Y3)*(LF13 + Y13*LF13o1) + (Y1*Y4)*LF14;



%GPEBCC

%%%%%%%%%%G(BCC,X,VA)%%%%%%%%%%

%G(BCC_A2,FE:VA;0)-H298(BCC_A2,FE;0) =
% 298.15<T< 1811.00: 
GF1VA = +1225.7+124.134*T-23.5143*T*log(T) - 0.00439752*T^2-5.8927E-08*T^3+77359*T^(-1) + GPFEBCC;
           %GPFEBB function must be defined before!
           
           

%G(BCC_A2,MN:VA;0)-H298(CBCC_A12,MN;0) =
%298.15<T< 1519.00: 
GF2VA= -3235.3+127.85*T-23.7*T*log(T)- 0.00744271*T^2 + 60000*T^(-1);

%G(BCC_A2,SI:VA;0)-H298(DIAMOND_A4,SI;0) =
%298.15<T< 1687.00: 
GF3VA = +38837.391+114.736859*T-22.8317533*T*log(T) - 0.001912904*T^2-3.552E-09*T^3+176667*T^(-1);

%G(BCC_A2,MO:VA;0)-H298(BCC_A2,MO;0) =
%298.15<T< 2896.00: 
GF4VA = -7746.302+131.9197*T-23.56414*T*log(T) - 0.003443396*T^2+5.66283E-07*T^3+65812*T^(-1)-1.30927E-10*T^4;

 %L(BCC_A2,FE,MN:VA;0) = 
 LF12VA = -2759+1.237*T;
 
% L(BCC_A2,FE,SI:VA;0) = 
LF13VA = -153141.13+46.48*T;

% L(BCC_A2,FE,SI:VA;1) = 
LF13o1VA =-92352;

% L(BCC_A2,FE,SI:VA;2) =
LF13o2VA = +62240;


%L(BCC_A2,FE,MO:VA;0) = 
LF14VA = +36818-9.141*T;

%L(BCC_A2,FE,MO:VA;1) = 
LF14o1VA = -362-5.724*T;

% LF4 MN MO

%L(BCC_A2,MN,SI:VA;0) = 
LF23VA = -89620.7+2.94097*T;

%L(BCC_A2,MN,SI:VA;1) =
LF23o1VA = -7500;

%L(BCC_A2,MO,SI:VA;0) = 
LF34VA = -11900.9+.66729*T;

%L(BCC_A2,MO,SI:VA;1) =
LF34o1VA = -78175.9;
  
%L(BCC_A2,FE,MN,SI:VA;0) =
LF123 = -97474;




GBCC_XVA = Y1*GF1VA + Y2*GF2VA + Y3*GF3VA + Y4*GF4VA + RT*YLNY + Y1*Y2*LF12VA + Y1*Y3*(LF13VA + Y13*LF13o1VA + (Y13^2)*LF13o2VA) + Y1*Y4*(LF14VA + Y14*LF14o1VA) + Y2*Y3*(LF23VA + Y23*LF23o1VA) + Y3*Y4*(LF34VA + Y34*LF34o1VA) + Y1*Y2*Y3*LF123;




 %%%%%%%%%% L(BCC, X, C:VA, 0 ) %%%%%%%%%%%
 
 %L(BCC_A2,FE:C,VA;0) = 
 LF1CVA = -190*T;
 
 
 
 LBCC_XCVA = Y1*LF1CVA;
 
 
 %%%%%%%%%% L(BCC, X, C:VA, 1 ) %%%%%%%%%%%

 
 LBCC_XVA1 = 0; 
 
 
 %%%%%%%%%%% BM(BCC, X, C, 0) %%%%%%%%%
 
 %BMAGN(BCC_A2,FE:C;0) 
 BMF1 = +2.22;
 
 
 BMBCC_ZC = Y1*BMF1;
 
 %%%%%%%%%%% BM(BCC, X, VA, 0) %%%%%%%%%
 % BMAGN(BCC_A2,FE:VA;0) = 
 BMF1VA = +2.22;
 %BMAGN(BCC_A2,MN:VA;0) =    298.15<T< 2000.00:
 BMF2VA = -0.27;

 
 BMBCC_ZVA = Y1*BMF1VA + Y2*BMF2VA; 
 
 %%%%%%%%%%% BM(BCC, X, CVA, 0) %%%%%%%%%
 
 
 %%%%%%%%%%% BM(BCC, X, CVA, 1) %%%%%%%%%

 
 
 
%%%%%%%%%%% TC(BCC, X, C, 0) %%%%%%%%%
%TC(BCC_A2,FE:C;0) 
TCF1 = +1043;
  

%TC(BCC_A2,FE,MO:C;0) 
TCF14= +335;
%TC(BCC_A2,FE,MO:C;1) 
TCF14o1= +526;



TCBCC_ZC = Y1*TCF1 + Y1*Y4* (TCF14 + Y14*TCF14o1);

%%%%%%%%%%% TC(BCC, X, VA, 0) %%%%%%%%%
%TC(BCC_A2,FE:VA;0) 
TCF1VA = +1043;
%TC(BCC_A2,FE,MN:VA;0) 
TCF12VA= +123;

%TC(BCC_A2,MN:VA;0) =    298.15<T< 2000.00: 
TCF2VA = -580;

%TC(BCC_A2,FE,MO:VA;0) 
TCF14VA= +335;
%TC(BCC_A2,FE,MO:VA;1) 
TCF14o1VA= +526;
%TC(BCC_A2,FE,SI:VA;0) 
TCF13VA =  0.0;
%TC(BCC_A2,FE,SI:VA;1)
TCF13o1VA= +504;


TCBCC_ZVA = Y1*TCF1VA + Y2*TCF2VA + Y1*Y4* (TCF14VA + Y14*TCF14o1VA)+ Y1*Y3* (TCF13VA + Y13*TCF13o1VA);

%%%%%%%%%%% TC(BCC, X, CVA, 0) %%%%%%%%%

 
%%%%%%%%%%% TC(BCC, X, CVA, 1) %%%%%%%%%



%%%%%%%CALCULATIONS%%%%%%
%%%Curie Temperature%%%%%%

TCBCC = YC*TCBCC_ZC + YVA*TCBCC_ZVA ; 
BMAGBCC = YC*BMBCC_ZC + YVA*BMBCC_ZVA ; 
TAO=T/TCBCC;
if T<TCBCC
    fta= +1-.905299383*TAO^(-1)-.153008346*TAO^3-.00680037095*TAO^9 -.00153008346*TAO^15;
end

if T>TCBCC
    fta=-.0641731208*TAO^(-5)-.00203724193*TAO^(-15) -4.27820805E-04*TAO^(-25);
end

BMAGBCCcor=BMAGBCC/(-1);
MAGBCC= RT*log(BMAGBCC + 1 )* fta;
GmBCC = YC*GBCC_XC + YVA*GBCC_XVA + YC*YVA*LBCC_XCVA + 3*RT*YLNYCVA + MAGBCC; 
 