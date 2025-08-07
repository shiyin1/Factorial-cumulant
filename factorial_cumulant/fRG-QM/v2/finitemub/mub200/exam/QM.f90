program QM

  implicit none

  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) T,mu,muq,muqb !temperature and chemical potential
  real(16) l,lb !polyakov loop
  real(16) rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,c,kappa
  real(16) sigma_UV,kappa_UV_i,kappa_UV_i_mu,kappa_UV
  real(16) Vtotal0
  real(16) Ti,dT,mu_down,mu_up,mudelta
  integer i,iTmax,j,jmumax,m,mm,i_mm,j1
  parameter(iTmax=300,jmumax=21)
  real(16) pre_res(0:jmumax,0:iTmax),T_res(0:iTmax),mu_res(0:jmumax,0:iTmax),pre_com(jmumax)
  real(16) mu_bound(0:iTmax),mu_bound_low,mu_bound_high
  real(16) T_MeV,muB_MeV
  real(16) mui,muBi,muB,muP,muPb
  integer mmax
  parameter(mmax=10)
!order of chebyshev polynomial
  real(16) dcd0mu(mmax),dcd1mu(mmax),dcd2mu(mmax),dcd3mu(mmax),dcd4mu(mmax),chi(mmax,0:iTmax)
  real(16) factorial
  real(16) fpi_res(0:jmumax,0:iTmax),mPion_res(0:jmumax,0:iTmax),mSigma_res(0:jmumax,0:iTmax),mf_res(0:jmumax,0:iTmax)
  integer iT,iv
  real(16) l_i,lb_i,l_i_mu,lb_i_mu
  real(16) Fnf0,Fnf1,Fnf2
  external Fnf0,Fnf1,Fnf2
  real(16) nfl,nfl0,nfl1,lset
  real(16) chebev,cl,Tc,tc1,t2,cl1
  external chebev


  common /Tmu/ T,mu,muq,muqb
  common /prefit/ pre_com
  common /iTiv/ iT,iv
  common /break/ cl

  dT=1.Q+0

  muBi=200.Q+0/hc
  
  j=j1

  sigma_UV=100.Q+0/hc
  kappa_UV_i=sigma_UV**2/2.Q+0

  l_i=1.Q-10
  lb_i=1.Q-10

  open(unit=11,file='m1.dat')
  read(11,*) j1
  close (11)
  write(*,*)'load OK',j1
  
  j=j1

  mu_up=50.Q+0/hc
  mu_down=-mu_up
  mudelta=cos(pi*(j-0.5Q+00)/real(jmumax,kind=16))*(0.5Q+00*(mu_up-mu_down))+0.5Q+00*(mu_up+mu_down)

  do i=1,iTmax

    T_MeV=i*dT
    T=T_MeV/hc
    
    muB=muBi+mudelta
    mu=muB/3.Q+0
    muB_MeV=muB*hc

    muP=0.Q+0/hc+mudelta
    muPb=0.Q+0/hc
    muq=muP/3.Q+0
    muqb=muPb/3.Q+0

      call selfEQ(kappa_UV_i,l_i,lb_i,kappa_UV,l,lb,rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,c,kappa)
      kappa_UV_i=kappa_UV
      l_i=l
      lb_i=lb

    write(*,"('i=', I4)")i
    write(*,"('muB_MeV=', f15.7, t25, 'T_MeV=', f15.7)")muB_MeV,T_MeV
    write(*,"('fpi=', f15.7, t25, 'mPion=', f15.7)")fpi*hc,mPion*hc
    write(*,"('mf=', f15.7, t25, 'mSigma=', f15.7)")mf*hc,mSigma*hc
    write(*,"('kappa_UV=', e21.14)")kappa_UV

    open(unit=51,file='./buffer/TMeV.dat',position='append')
    write(51,*)T_MeV
    close(51)

    open(unit=51,file='./buffer/mPion.dat',position='append')
    write(51, "(e21.14)")mPion*hc
    close(51)

    open(unit=51,file='./buffer/mSigma.dat',position='append')
    write(51,*)mSigma*hc
    close(51)

    open(unit=51,file='./buffer/mf.dat',position='append')
    write(51,*)mf*hc
    close(51)

    open(unit=51,file='./buffer/Vtotal.dat',position='append')
    write(51,*)Vtotal
    close(51)

    open(unit=51,file='./buffer/fpi.dat',position='append')
    write(51,*)fpi*hc
    close(51)

    open(unit=51,file='./buffer/Zphi.dat',position='append')
    write(51,*)Zphi
    close(51)

    open(unit=51,file='./buffer/h.dat',position='append')
    write(51,*)h
    close(51)

    open(unit=51,file='./buffer/c.dat',position='append')
    write(51,*)cl
    close(51)

    open(unit=51,file='./buffer/l.dat',position='append')
    write(51,*)l
    close(51)

    open(unit=51,file='./buffer/lb.dat',position='append')
    write(51,*)lb
    close(51)

  end do

end






