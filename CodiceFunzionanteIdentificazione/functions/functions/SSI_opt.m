function [IDENT]=SSI_opt(sig,i,fc,nn)
%function [w1,z1,mod1,ph1,IDENT]=SSI(sig,i,fc,nn)

deltat=1/fc;

% ESTRAE FORME MODALI E COMPLESSE CONIUGATE PER CUI VANNO PRESI UNO PER UNO
% W1              = PULSAZIONI
% Z1              = SMORZAMENTI RELATIVI
% MOD1            = MODULI
% PH1             = FASI
%--------------------------------------------------------------------------
% sig                   CANALI DEL SEGNALE DI INGRESSO DA ANALIZZARE
% i                     NUMERO DI BLOCCHI DI RIGHE NELLA MATRICE DI HANKEL
% deltat                TEMPO DI CAMPIONAMENTO
% nn                    NUMERO DI MODI DA ESTRARRE (I VETTORI HANNO
%                       DIMENSIONE DOPPIA POICHE' CONTENGONO PURE I COMPLESSI E CONIUGATI)
%--------------------------------------------------------------------------
%	          Authors:   A.Saettone     January 2005
% Following Author(s):   G.Abbiati      August 2009
%--------------------------------------------------------------------------

[nr,nc]=size(sig);                              % controllo dimensioni del segnale

if nr<nc                                        % i segnali acquisiti devono essere su righe
    sig=sig'; 
end
[nr,nc]=size(sig);

n=2*nn;
[A,K,C,R] = sto_pos_opt(sig,i,n);
[CSI,Z]=eig(A);                                 %Calcolo dei parametri modali

for i=1:1:n
    w(i)=abs(log(Z(i,i)))/deltat;
    z(i)=-real(log(Z(i,i)))/abs(log(Z(i,i)));
end

Forme=C*CSI;                                    %Calcolo forme modali
mod=abs(Forme);                                 %Normalizzazione forme modali 

for i=1:n
    mod(:,i)=mod(:,i)/max(mod(:,i));
end

ph=angle(Forme);                                %le fasi sono in radianti

for i=1:1:nc
    ph1(i,:)=ph(i,:)-ph(1,:);
end

w1=w'/(2*pi);                                      %Passo alle frequenze
z1=z';
mod1=mod';
ph1=ph1';

clear w z mod ph

DUMMY=[w1,z1,mod1.*sign(cos(ph1))];

if exist('DUMMY')==0
    IDENT=[];
else
    for i=2:2:length(DUMMY(:,1));
        IDENT(i/2,:)=DUMMY(i,:);
    end
end
