function Fnb(x,T)          !boson distribution function
  implicit none

  real(16) x, T
  real(16) Fnb
  real(16) over

  over=x/T
  if(over > 60.Q+0)then
    Fnb=0.Q+0
  else
    Fnb=1.Q+0/(exp(over)-1.Q+0)
  end if

  return
end function Fnb

function Fnf0(x,T,l,lb)         !fermion distribution function

  implicit none

  real(16) x,T,l,lb
  real(16) Fnf0
  real(16) over

  over=x/T
  if(over > 60.Q+0)then
    Fnf0=0.Q+0
  else
    Fnf0=1.Q+0/(1.Q+0+3.Q+0*lb*exp(x/T)+3.Q+0*l*exp(2.Q+0*x/T)+exp(3.Q+0*x/T))
  end if

  return
end function Fnf0


function Fnf1(x,T,l,lb)         !fermion distribution function
  implicit none

  real(16) x,T,l,lb
  real(16) Fnf1
  real(16) over

  over=x/T
  if(over > 60.Q+0)then
    Fnf1=0.Q+0
  else
    Fnf1=(2*exp(x/T))/(1.Q+0+3.Q+0*lb*exp(x/T)+3.Q+0*l*exp(2.Q+0*x/T)+exp(3.Q+0*x/T))
  end if

  return
end function Fnf1

function Fnf2(x,T,l,lb)         !fermion distribution function
  implicit none

  real(16) x,T,l,lb
  real(16) Fnf2
  real(16) over

  over=x/T
  if(over > 60.Q+0)then
    Fnf2=0.Q+0
  else
    Fnf2=(exp(2.Q+0*x/T))/(1.Q+0+3.Q+0*lb*exp(x/T)+3.Q+0*l*exp(2.Q+0*x/T)+exp(3.Q+0*x/T))
  end if

  return
end function Fnf2


