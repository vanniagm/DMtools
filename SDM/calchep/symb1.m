(*
     ==============================
     *  CalcHEP  3.6.27 *
     ==============================
  process  ~phidm(p1)+~phidm(p2)->b(p3)+b~(p4)
*)

parameters={
 aEWM1 -> 0.00000000000*10^(0)
,Gf -> 0.00000000000*10^(0)
,ymb -> 0.00000000000*10^(0)
,lsk -> 0.00000000000*10^(0)
,MZ -> 0.00000000000*10^(0)
,MB -> 0.00000000000*10^(0)
,MH -> 0.00000000000*10^(0)
,mdm -> 0.00000000000*10^(0)
,Pi -> 0.00000000000*10^(0)
,aEW -> pow[aEWM1,-1]
,MW -> pow[pow[MZ,2]/2.+pow[-(aEW*Pi*pow[2,-0.5]*pow[Gf,-1]*pow[MZ,2])+
 pow[MZ,4]/4.,0.5],0.5]
,sw2 -> 1-pow[MW,2]*pow[MZ,-2]
,EE -> 2*pow[aEW,0.5]*pow[Pi,0.5]
,sw -> pow[sw2,0.5]
,vev -> 2*MW*sw*pow[EE,-1]
,yb -> ymb*pow[2,0.5]*pow[vev,-1]
,x10x0 -> -(yb*pow[2,-0.5])
,WH -> 0.00000000000*10^(0)
           };

substitutions={
 x4x0->lsk*vev
              };

inParticles = {"~phidm","~phidm"}
outParticles = {"b","b~"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p4 = +p1+p2-p3;
p1/: SC[p1,p1] =mdm^2;
p2/: SC[p2,p2] =mdm^2;
p3/: SC[p3,p3] =MB^2;
p2/: SC[p2,p3] = -1*(MB^2-mdm^2-mdm^2-MB^2-2*SC[p1,p2]+2*SC[p1,p3])/2;

initSum;

(*
  Diagram  1 in subprocess 1
*)
totFactor = ((12*x10x0^2*x4x0^2)/(1));
numerator =(SC[p1,p2]+mdm^2-2*MB^2);
denominator =(propDen[-p1-p2,MH,WH]^2);

addToSum;

finishSum;
