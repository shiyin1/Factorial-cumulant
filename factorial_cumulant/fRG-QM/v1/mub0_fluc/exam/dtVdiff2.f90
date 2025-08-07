subroutine dtVdiff2(kk,l,lb,dtVd1l,dtVd1lb)
!This subroutine calculate various differential kernal corresponding k>k_UV

  implicit none

  real(16) kk
  real(16) l,lb,l_con,lb_con
  real(16) k,etaphi,etapsi
  real(16) mf2
  real(16) T,mu
  real(16) xff,xfa
  real(16) zb,zf !distinguish the transverse and longituidanl wave function renormalization
  real(16) nf0,nf1,nf2
  real(16) Fnf0,Fnf1,Fnf2
  external Fnf0,Fnf1,Fnf2
  real(16) Nc,Nf
  parameter(Nc=3.Q+0,Nf=2.Q+0)
  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) v3
  parameter(v3=1.Q+0/(2.Q+0*pi**2))
  real(16) nfd1l,nfd1lb
  real(16) nfd1lf,nfd1lbf,nfd1la,nfd1lba
  real(16) dtVd1l,dtVd1lb

  common /Tmu/ T,mu


  l_con=l
  lb_con=lb

  k=kk

  etaphi=0.Q+0
  etapsi=0.Q+0


  zb=1.Q+0
  zf=1.Q+0

  mf2=0.Q+0

  xff=-mu + (k*Sqrt(1 + mf2))/zf
  xfa=mu + (k*Sqrt(1 + mf2))/zf

  l=l_con
  lb=lb_con

  nf0=Fnf0(xff,T,l,lb)
  nf1=Fnf1(xff,T,l,lb)
  nf2=Fnf2(xff,T,l,lb)

  nfd1l=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
  nfd1lb=-(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))/2.Q+0

  nfd1lf=nfd1l
  nfd1lbf=nfd1lb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l=lb_con
  lb=l_con

  nf0=Fnf0(xfa,T,l,lb)
  nf1=Fnf1(xfa,T,l,lb)
  nf2=Fnf2(xfa,T,l,lb)

  nfd1l=-(nf2*(-1 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))
  nfd1lb=-(nf1*(-2 + 3*nf0 + 3*lb*nf1 + 3*l*nf2))/2.Q+0

  nfd1la=nfd1l
  nfd1lba=nfd1lb

  l=l_con
  lb=lb_con
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dtVd1l=(-2*(1 - etapsi/4.Q+0)*k**4*Nc*Nf*(-nfd1lba - nfd1lf)*v3)/(3.Q+0*Sqrt(1 + mf2)*zf)
  dtVd1lb=(-2*(1 - etapsi/4.Q+0)*k**4*Nc*Nf*(-nfd1la - nfd1lbf)*v3)/(3.Q+0*Sqrt(1 + mf2)*zf)

end
