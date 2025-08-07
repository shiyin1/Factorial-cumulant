subroutine vInf(k_UV,Vinfi)
!calculate the contribution to thermaldynamical potential arising from
!momenta above UV cutoff

  implicit none
  real(16) k_UV,Vinfi
  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) T,mu,muq,muqb
  real(16) l_com,lb_com
  real(16) l,lb,l_con,lb_con
  integer:: npoint, in
  parameter(npoint=512)
  real(16) w(npoint),y(npoint) !Guass integral
  real(16) k,k_infi
  real(16) Nc,Nf
  parameter(Nc=3.Q+0,Nf=2.Q+0)
  real(16) v3
  parameter(v3=1.Q+0/(2.Q+0*pi**2))
  real(16) Fnf0,Fnf1,Fnf2
  external Fnf0,Fnf1,Fnf2
  real(16) nff,nfa
  real(16) nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x
  real(16) xff,xfa
  real(16) muin



  common /Tmu/ T,mu,muq,muqb
  common /polyakov_com/ l_com,lb_com

  muin=550.Q+00/3.Q+00/hc

  l=l_com
  lb=lb_com

  l_con=l
  lb_con=lb

  k_infi=k_UV*10.Q+0

  call gauleg(k_UV, k_infi, y, w, npoint)

  Vinfi=0.Q+0
  do in=1, npoint
    k=y(in)

    xff=-muin-muq  + k
    xfa= muin-muqb + k

    l=l_con
    lb=lb_con

    nf0=Fnf0(xff,T,l,lb)
    nf1=Fnf1(xff,T,l,lb)
    nf2=Fnf2(xff,T,l,lb)

    nff=nf0 + lb*nf1 + l*nf2

    l=lb_con
    lb=l_con

    nf0=Fnf0(xfa,T,l,lb)
    nf1=Fnf1(xfa,T,l,lb)
    nf2=Fnf2(xfa,T,l,lb)

    nfa=nf0 + lb*nf1 + l*nf2

    l=l_con
    lb=lb_con


    Vinfi=Vinfi+w(in)*k**3*(0.Q+0-nfa-nff)
  end do
  Vinfi=Vinfi*v3*4*Nc*Nf/3.Q+0/2.Q+0

end
