function [GmFCC,TCFCC] = GfccPE(XC,T)

%%X in mole fractions

%XFE=variable;
 XMN=0.01724;
XSI=0.00535;
XMO=0.00157;
%XC=variable;
XFE=1-XMN-XSI-XMO-XC;
XVA=1-XC;
XN=0;


% YJs are site fractions defined as 
%for ferrite q=3 and p=1
p=1;
q=1;


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

%T=variable; %K
P=1E5; %Pa


%FUNCTIONS
%the order (in definition) is imporatant. We only need GPFEXXX, but to get
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



%%%%%%%%%%%%


YLNY= Y1*log(Y1) + Y2*log(Y2)+ Y3*log(Y3) + Y4*log(Y4);
YLNYCVA = YC*log(YC) + YVA*log(YVA);

%%%%%%%%%%%%%%



%%%%%%%%%%G(fcc,X,C)%%%%%%%%%%
%G(FCC_A1,FE:C;0)-H298(GRAPHITE,C;0)-H298(BCC_A2,FE;0) =
%298.15<T< 1811.00: 
GA1 = +59601.859 + 287.269*T - 48.9643*T*log(T) - 0.00422982*T^2 - 5.8927E-08*T^3 + 2639959*T^(-1)- 2.643E+08*T^(-2) + 1.2E+10*T^(-3);

%G(FCC_A1,MN:C;0)-H298(GRAPHITE,C;0)-H298(CBCC_A12,MN;0) =
%298.15<T< 1519.00: 
GA2 = -24981.721+316.05*T-47.7582*T*log(T) - 0.00781998*T^2 + 2632427*T^(-1)- 2.643E+08*T^(-2)+ 1.2E+10*T^(-3);

%G(FCC_A1,SI:C;0)-H298(GRAPHITE,C;0)-H298(DIAMOND_A4,SI;0) =
%298.15<T< 1687.00:
GA3 = -46041.05 + 346.657259*T - 47.1317533*T*log(T)- 0.002385204*T^2-3.552E-09*T^3+2739267*T^(-1)-2.643E+08*T^(-2) +1.2E+10*T^(-3);

%G(FCC_A1,MO:C;0)-H298(GRAPHITE,C;0)-H298(BCC_A2,MO;0) =
%298.15<T< 2896.00: 
GA4 = -32614.743+294.3497*T-47.86414*T*log(T) - 0.003915696*T^2 + 5.66283E-07*T^3 + 1878412*T^(-1) -1.30927E-10*T^4-2.643E+08*T^(-2)+1.2E+10*T^(-3);


%L(FCC_A1,FE,MN:C;0) 
LA12 = +34052-23.467*T;
% L(FCC_A1,FE,MO:C;0) 
LA14 = +6000;
%L(FCC_A1,FE,SI:C;0) 
LA13 = +226100-34.25*T;
% L(FCC_A1,FE,SI:C;1) 
LA13o1= -202400;
%L(FCC_A1,MN,SI:C;0) 
LA23 =  0.0;


GFCC_XC= Y1*GA1 + Y2*GA2 + Y3*GA3 + Y4*GA4 + RT*YLNY + (Y1*Y2)*LA12 +  (Y1*Y3)*(LA13 + Y13*LA13o1) + (Y1*Y4)*LA14;


%%%%%%%%%%G(FCC,X,VA)%%%%%%%%%%
%G(FCC_A1,FE:VA;0)-H298(BCC_A2,FE;0) =
%298.15<T< 1811.00: 
GA1VA = -236.7+132.416*T - 24.6643*T*log(T) - 0.00375752*T^2 - 5.8927E-08*T^3+77359*T^(-1)+GPFEFCC;
%G(FCC_A1,MN:VA;0)-H298(CBCC_A12,MN;0) =
%298.15<T< 1519.00: 
GA2VA = -3439.3+131.884*T - 24.5177*T*log(T)- 0.006*T^2 + 69600*T^(-1);
%G(FCC_A1,SI:VA;0)-H298(DIAMOND_A4,SI;0) =
%298.15<T< 1687.00: 
GA3VA = +42837.391+115.436859*T-22.8317533*T*log(T) - 0.001912904*T^2-3.552E-09*T^3+176667*T^(-1);
%G(FCC_A1,MO:VA;0)-H298(BCC_A2,MO;0) =
%298.15<T< 2896.00: 
GA4VA = +7453.698+132.5497*T-23.56414*T*log(T)- 0.003443396*T^2 + 5.66283E-07*T^3 + 65812*T^(-1)- 1.30927E-10*T^4;


%L(FCC_A1,FE,MN:VA;0) 
LA12VA = -7762+3.865*T;
%L(FCC_A1,FE,MN:VA;1) 
LA12o1VA = -259;
%L(FCC_A1,FE,MN,SI:VA;0) 
LA123VA= -56655-55.613*T;
%L(FCC_A1,FE,MO:VA;0) 
LA14VA= +28347-17.691*T;
%L(FCC_A1,FE,SI:VA;0) 
LA13VA= -125247.7+41.166*T;
%L(FCC_A1,FE,SI:VA;1) 
LA13o1VA= -142707.6;
%L(FCC_A1,FE,SI:VA;2) 
LA13o2VA= +89907.3;
%L(FCC_A1,MN,SI:VA;0) 
LA23VA= -95600+2.94097*T;
%L(FCC_A1,MN,SI:VA;1) 
LA23o1VA= -7500;



GFCC_XVA = Y1*GA1VA + Y2*GA2VA + Y3*GA3VA + Y4*GA4VA + RT*YLNY + Y1*Y2*(LA12VA + (Y1-Y2)*LA12o1VA) + Y1*Y3*(LA13VA + Y13*LA13o1VA + (Y13^2)*LA13o2VA) + Y1*Y4*LA14VA  + Y2*Y3*(LA23VA + Y23*LA23o1VA)  + Y1*Y2*Y3*LA123VA;


%%%%%%%%%% L(FCC, X, C:VA, 0 ) %%%%%%%%%%%
%L(FCC_A1,MN:C,VA;0) 
LA2CVA = -43433;
%L(FCC_A1,FE:C,VA;0) 
LA1CVA = -34671;
%L(FCC_A1,MO:C,VA;0)
LA4CVA = -41300;
%L(FCC_A1,SI:C,VA;0) 
LA3CVA =  0.0;

LFCC_XCVA = Y1*LA1CVA + Y2*LA2CVA + Y4*LA4CVA ;
  
  
%%%%%%%%%% L(BCC, X, C:VA, 1 ) %%%%%%%%%%%
%nothing was defined here.
 
 
LFCC_XVA1 = 0; 
  
  
%%%%%%%%%%% BM(FCC, X, C, 0) %%%%%%%%%
   
%BMAGN(FCC_A1,FE:C;0) 
BMA1 = -2.1;

BMFCC_ZC = Y1*BMA1;
 
 
%%%%%%%%%%% BMFCCBCC, X, VA, 0) %%%%%%%%%
%BMAGN(FCC_A1,FE:VA;0) 
BMA1VA = -2.1;
%BMAGN(FCC_A1,MN:VA;0) =    298.15<T< 2000.00:
BMA2VA = -1.86;

BMFCC_ZVA = Y1*BMA1VA + Y2*BMA2VA; 


%%%%%%%%%%% BM(FCC, X, CVA, 0) %%%%%%%%%
 
 
%%%%%%%%%%% BM(FCC, X, CVA, 1) %%%%%%%%%

 
 
 
%%%%%%%%%%% TC(FCC, X, C, 0) %%%%%%%%%

%TC(FCC_A1,FE:C;0) 
TCA1 = -201;

TCFCC_ZC = Y1*TCA1;
   

%%%%%%%%%%% TC(FCC, X, VA, 0) %%%%%%%%%

%TC(FCC_A1,FE:VA;0) 
TCA1VA= -201;
%TC(FCC_A1,MN:VA;0) =    298.15<T< 2000.00: 
TCA2VA=-1620;
%TC(FCC_A1,FE,MN:VA;0) 
TCA12VA = -2282;
%TC(FCC_A1,FE,MN:VA;1) 
TCA12o1VA= -2068;

%TC(FCC_A1,FE,MN,SI:VA;0) 
TCA123VA = +13854;

TCFCC_ZVA = Y1*TCA1VA + Y2*TCA2VA + Y1*Y2* (TCA12VA + Y12*TCA12o1VA)+ Y1*Y3*Y2*TCA123VA;


%%%%%%%CALCULATIONS%%%%%%
%%%Curie Temperature%%%%%%

TCFCC = YC*TCFCC_ZC + YVA*TCFCC_ZVA ; 
BMAGFCC = YC*BMFCC_ZC + YVA*BMFCC_ZVA ; 
TAO=T/TCFCC;
if T<TCFCC
    fta= +1-.860338755*TAO^(-1)- 0.17449124*TAO^3-.00775516624*TAO^9 -.0017449124*TAO^15;
end

if T>TCFCC
    fta= -.0426902268*TAO^(-5)-.0013552453*TAO^(-15) -2.84601512E-04*TAO^(-25);
end

BMAGFCCcor=BMAGFCC/(-3);
MAGFCC= RT*log(BMAGFCCcor + 1 )* fta;

GmFCC = YC*GFCC_XC + YVA*GFCC_XVA + YC*YVA*LFCC_XCVA + RT*YLNYCVA + MAGFCC; 
