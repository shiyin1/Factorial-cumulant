program QM

  implicit none

  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) T,mu !temperature and chemical potential
  real(16) l,lb !polyakov loop
  real(16) rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,c,kappa
  real(16) sigma_UV,kappa_UV_i,kappa_UV_i_mu,kappa_UV
  real(16) Vtotal0
  real(16) Ti,dT,mu_down,mu_up
  integer i,iTmax,j,jmumax,m,mm,i_mm,j1
  parameter(iTmax=299,jmumax=21)
  real(16) pre_res(0:jmumax,0:iTmax),T_res(0:iTmax),mu_res(0:jmumax,0:iTmax),pre_com(jmumax)
  real(16) mu_bound(0:iTmax),mu_bound_low,mu_bound_high
  real(16) T_MeV,muB_MeV
  real(16) mui,muBi,muB
  integer mmax
  parameter(mmax=10)
!order of chebyshev polynomial
  real(16) dcd0mu(mmax),dcd1mu(mmax),dcd2mu(mmax),dcd3mu(mmax),dcd4mu(mmax),dcd5mu(mmax)
  real(16) factorial,dcd6mu(mmax),dcd7mu(mmax),dcd8mu(mmax),chi(mmax,0:iTmax)
  real(16) fpi_res(0:jmumax,0:iTmax),mPion_res(0:jmumax,0:iTmax),mSigma_res(0:jmumax,0:iTmax)
  integer iT,iv
  real(16) l_i,lb_i,l_i_mu,lb_i_mu,mf_res(0:jmumax,0:iTmax)
  real(16) Fnf0,Fnf1,Fnf2
  external Fnf0,Fnf1,Fnf2
  real(16) nfl,nfl0,nfl1,lset
  real(16) chebev
  character(14) fileV1
  character(15) fileV2
  external chebev


  common /Tmu/ T,mu
  common /prefit/ pre_com
  common /iTiv/ iT,iv



  Ti=1.Q+0/hc                 !initial temperature
  dT=1.Q+0/hc                 !stepsize of temperature

  do i=0, iTmax
	T_res(i)=Ti+dT*real(i)    !fm**(-1)
  end do



  do j1=1,9
    j=jmumax+1-j1
  write(fileV1,'(a9,i1,a4)')'../data/V',j1,'.dat'
  open(unit=51,file=Trim(fileV1))
  iT=iTmax
  do i=0, iT
     read(51,*)pre_res(j,i) 
  end do
  close(51)
  end do

  do j1=10,21
    j=jmumax+1-j1
  write(fileV2,'(a9,i2,a4)')'../data/V',j1,'.dat'
  open(unit=51,file=Trim(fileV2))
  iT=iTmax
  do i=0, iT
     read(51,*)pre_res(j,i)    
  end do
  close(51)
  end do

  write(*,*)'load OK'

  pre_res=-pre_res   

  do i=0, iTmax

    mu_down=50.Q+0/hc
    mu_up=-50.Q+0/hc

    T=T_res(i)

    do j=1, jmumax
      pre_com(j)=pre_res(j,i)
	  end do
    call chebft(dcd0mu,jmumax,mmax)

    call chder(mu_down,mu_up,dcd0mu,dcd1mu,mmax)
    call chder(mu_down,mu_up,dcd1mu,dcd2mu,mmax)
    call chder(mu_down,mu_up,dcd2mu,dcd3mu,mmax)
    call chder(mu_down,mu_up,dcd3mu,dcd4mu,mmax)
    call chder(mu_down,mu_up,dcd4mu,dcd5mu,mmax)
    call chder(mu_down,mu_up,dcd5mu,dcd6mu,mmax)
    call chder(mu_down,mu_up,dcd6mu,dcd7mu,mmax)
    call chder(mu_down,mu_up,dcd7mu,dcd8mu,mmax)

    chi(1,i)=chebev(mu_down,mu_up,dcd0mu,mmax,0.Q0)/T**4
    chi(2,i)=chebev(mu_down,mu_up,dcd1mu,mmax,0.Q0)/T**3
    chi(3,i)=chebev(mu_down,mu_up,dcd2mu,mmax,0.Q0)/T**2
    chi(4,i)=chebev(mu_down,mu_up,dcd3mu,mmax,0.Q0)/T
    chi(5,i)=chebev(mu_down,mu_up,dcd4mu,mmax,0.Q0)
    chi(6,i)=chebev(mu_down,mu_up,dcd5mu,mmax,0.Q0)*T
    chi(7,i)=chebev(mu_down,mu_up,dcd6mu,mmax,0.Q0)*T**2
    chi(8,i)=chebev(mu_down,mu_up,dcd7mu,mmax,0.Q0)*T**3
    chi(9,i)=chebev(mu_down,mu_up,dcd8mu,mmax,0.Q0)*T**4

  end do

  open(unit=51,file='./buffer/chi0.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(1,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi1.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(2,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi2.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(3,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi3.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(4,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi4.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(5,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi5.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(6,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi6.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(7,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi7.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(8,i)
  end do
  close(51)

  open(unit=51,file='./buffer/chi8.dat')
  do i=0, iTmax
    write(51, "(e21.14)")chi(9,i)
  end do
  close(51)

  open(unit=51,file='./buffer/p.dat')
  do i=0, iTmax
    write(51, "(e21.14)")pre_res(0,i)
  end do
  close(51)

  open(unit=51,file='./buffer/muBi.dat')
  write(51, "(e21.14)")muBi*hc
  close(51)

end






