program cranck

!Se declaran las variables a usar, recordando siempre que
!La funcion de onda es compleja
implicit none
double precision, allocatable, dimension (:) :: x
!double precision, allocatable, dimension (:) :: pot,kine
double complex, allocatable, dimension(:) ::  wavefun,quanfun, swave,expv,Jcurr
double precision :: Lx, bx, ax, delx, refle, trans, absor,averagex,secmom
double precision :: delt, tf, ti,T,inte,poten,velx,cap,R0,V0,exactsol,x0
double precision :: pot,kine,AvJcurr
double precision :: tau, ome, phi, E0
double complex :: ii = (0.0d0,1.0d0),nux
double complex, allocatable, dimension (:) :: work
double complex, allocatable, dimension(:,:) :: aplus, aminus,c
integer, allocatable, dimension(:):: ipiv
double precision :: pi = 4.0d0*atan(1.0d0)
double complex :: temp, wavefun1,temp1
integer :: nx, xi, j, info,mt,k,nt
!Se dan los valores de la red numerica y el tiempo de difusion
!print*, 'dame (x0)'
!read*, x0
!print*, 'dame la velocidad (velx)'
!read*, velx
!print*, 'dame el tamaño del potencial (R0, V0)'
!read*, R0,V0
!print*, 'Dame los parametros del laser(E0, tau, omega, phi)'
!read*, E0, tau, ome, phi
print*, 'dame el rango de propagacion en x (xi,xf)'
read*, ax, bx
print*, 'dame delta el invervalo de tiempo(tf,ti)'
read*, tf, ti
!Se da el numero de puntos en la red y para la propagacion
print*, 'dame el numero de puntos(x,t)'
read*, nx, mt
!Se abre memoria
allocate(expv(0:nx),x(0:nx),wavefun(0:nx),swave(0:nx),aplus(nx,nx),aminus(nx,nx),ipiv(nx),work(nx))
allocate(quanfun(0:nx),c(nx,nx))
allocate(Jcurr(0:nx))
!Se abre un archivo para graficar la solucion del metodo
open(unit = 11, file = 'fort.11', status = 'unknown')
!Se define el rango espacial y el temporal
!Se define la longitud entre puntos de la red y del tiempo
Lx = bx - ax
T = tf - ti
delx = Lx / nx
delt = T/ mt
nux = ii*delt/(4.0d0*(delx**2))      
!Se manda a imprimir a pantalla el coeficiente de Crank-Nicolson para vericar que el
!programa sea estable
print*, nux


do xi = 1, nx
	do j = 1, nx
		if(xi == j) then
			aplus(xi,j) = 1.0d0 + 2.0d0*nux
			aminus(xi,j) = 1.0d0 - 2.0d0*nux
		else if(xi == j + 1) then 
			aplus(xi,j) = -nux
			aminus(xi,j) = nux
		else if(xi == j - 1) then
			aplus(xi,j) = -nux
			aminus(xi,j) = nux
		else
			aplus(xi,j) = 0.0d0
			aminus(xi,j) = 0.0d0
		end if
	end do
end do

!Se llama a la subrutina que realiza la descomposicion LU
call ZGETRF( nx, nx, aplus, nx, ipiv, INFO )
!Se llama a la subrutina que realiza la inversion de la matriz A+
call ZGETRI( nx,aplus, nx,IPIV, WORK, nx, info)
!Se multiplican las matrices
c = 0.0d0
do xi=1,nx
	do j=1,nx
		do k=1,nx
			c(xi,j) = c(xi,j) + aplus(xi,k)*aminus(k,j)
		end do
	end do
end do
!Se llama a la función de onda incidente a propagar
!y se normaliza
inte = 0.0d0
do xi = 0, nx
	x(xi) = ax + xi*delx
	inte = inte + (abs(wavefun1(x(xi),x0,ti,delt,nt))**2)
end do
inte = inte*delx

do xi = 0, nx
   wavefun(xi) = wavefun1(x(xi),x0,ti,delt,nt)/sqrt(inte)
end do


pot = 0.0d0
kine = 0.0d0
do xi = 1, nx
	pot = pot + wavefun(xi)*poten(x(xi),ti,delt,nt)*wavefun(xi)*delx
	kine = kine -0.5d0*wavefun(xi)*((wavefun(xi+1)-2.0d0*wavefun(xi)+wavefun(xi-1))/delx**2)*delx
end do



do xi = 1, nx
	write(11+1,*) x(xi), (abs(wavefun(xi)))**2, pot,kine,pot+kine !poten(x(xi),ti,delt,nt)!!linea 109  (abs(wavefun(xi)))**2 , poten(x(xi)
end do



!Se realiza la propgacion de la onda incidente
do nt=0,mt-1
	do xi = 1, nx
		swave(xi) = exp(-ii*delt*poten(x(xi),ti,delt,nt))*wavefun(xi)
	end do
	do xi = 1, nx
		temp = (0.0d0,0.0d0)
		do k = 1, nx
			temp = temp + c(xi,k)*swave(k)
		end do
		quanfun(xi) = temp
	end do
	wavefun = quanfun

end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!probability current
!do xi =1,nx
!	Jcurr(xi)=-0.5*ii*( DCONJG(wavefun(xi)) *((wavefun(xi+1)-wavefun(xi))/delx) &
!	-wavefun(xi)*((DCONJG(wavefun(xi+1))-DCONJG(wavefun(xi)))/delx))
!end do

do xi=1,nx
	Jcurr(xi)= DCONJG(wavefun(xi)) *((wavefun(xi+1)-wavefun(xi))/delx)
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!av. position
averagex=0.0d0 
do xi = 1, nx
	averagex = averagex + x(xi)*(abs(wavefun(xi)))**2*delx
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!secondmoment
secmom=0.0d0 
do xi = 1, nx
	averagex = averagex + (x(xi)**2)*(abs(wavefun(xi)))**2*delx
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!average current
AvJcurr=0.0d0
do xi=1,nx
	AvJcurr= AvJcurr + AIMAG(Jcurr(xi))*(abs(wavefun(xi)))**2*delx
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pot = 0.0d0
kine = 0.0d0
do xi = 1, nx
	pot = pot + wavefun(xi)*poten(x(xi),ti,delt,nt)*wavefun(xi)*delx
	kine = kine -0.5d0*wavefun(xi)*((wavefun(xi+1)-2.0d0*wavefun(xi)+wavefun(xi-1))/delx**2)*delx
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do xi = 1, nx
	write(11+nt,*) x(xi), (abs(wavefun(xi)))**2,AIMAG(Jcurr(xi)),averagex,secmom,pot,kine,pot+kine,AvJcurr!, poten(x(xi),ti,delt,nt)!, exactsol(x(xi),tf)
end do

! imprime energias
!do xi = 1, mt
!	write(55+nt,*) x(xi), pot(xi)+kine(xi)
!end do



!
!




!do xi = 1, nx
!	pot = pot + swave(xi)*poten(x(xi),ti,delt,nt)*swave(xi)*delx
!	kine = kine-0.5d0*swave(xi)*((swave(xi+1)-2.0d0*swave(xi)+swave(xi-1))/delx**2)*delx
!end do

!print*, 'potential energy: ',pot
!print*, 'kinetic energy: ', kine
!print*, 'total energy: ', pot + kine







deallocate(expv,x,wavefun,swave,aplus,aminus,ipiv,work)
deallocate(quanfun,c)
close(10)
close(11)
end program








!####################################Potencial con minimo dependiente del tiempo#########################################3
!potencial que solo hace oscilar la posicion del minimo
!double precision function poten(x,ti,delt,nt)
!IMPLICIT NONE
!double precision :: x,ti,delt
!integer :: nt
!poten = 0.5*(x-0.1*sin(2*(ti+ delt*nt)))**2
!return
!end function

!double complex function wavefun1(x,x0,ti,delt,nt)
!double precision :: x,x0,ti,delt
!double precision :: pi = 4.0d0*atan(1.0d0)
!double complex :: ii = (0.0d0,1.0d0)
!integer :: nt
!wavefun1 = exp(-x**2 /2.0d0 )/sqrt(sqrt(pi))
!return
!end function


!double precision function exactsol(x,tf)
!double precision :: x,tf
!double precision :: pi = 4.0d0*atan(1.0d0)
!exactsol = exp(-(x+(0.05*tf*(cos(tf)-sin(tf)/(tf))))**2  )/sqrt(pi)
!return
!end function
!double precision function exactsol(x,tf)
!double precision :: x,tf
!double precision :: pi = 4.0d0*atan(1.0d0)
!exactsol = exp(-(x+(0.25*tf*(cos(tf)-sin(tf)/tf)))**2  )/sqrt(pi)
!return
!end function
!##########################################################################################################################








!Se define la funcion potencial
!double precision function poten(x)
!double precision :: x
!poten = 0.5*(x**2)
!return
!end function







!####################################################Potencial de Dubinko######################################################################
double complex function wavefun1(x,x0,ti,delt,nt)
double precision :: x,x0,ti,delt,xmin
double precision :: pi = 4.0d0*atan(1.0d0),alpha=0.0005d0
double complex :: ii = (0.0d0,1.0d0)
integer :: nt
xmin=-1/(sqrt(2.0d0)*sqrt(sqrt(alpha)))
wavefun1=exp(-(x-xmin)**2 /2.0d0)/sqrt(sqrt(pi))
return
end function

!first excited state
!double complex function wavefun1(x,x0,ti,delt,nt)
!double precision :: x,x0,ti,delt,xmin
!double precision :: pi = 4.0d0*atan(1.0d0),alpha=0.0005d0
!double complex :: ii = (0.0d0,1.0d0)
!integer :: nt
!xmin=0.0d0!-1/(sqrt(2.0d0)*sqrt(sqrt(alpha)))
!wavefun1=sqrt(2.0d0)*x*exp(-(x-xmin)**2 /2.0d0)/sqrt(sqrt(pi))
!return
!end function

!second excited state
!double complex function wavefun1(x,x0,ti,delt,nt)
!double precision :: x,x0,ti,delt,xmin
!double precision :: pi = 4.0d0*atan(1.0d0),alpha=0.0005d0
!double complex :: ii = (0.0d0,1.0d0)
!integer :: nt
!xmin=-1/(sqrt(2.0d0)*sqrt(sqrt(alpha)))
!wavefun1=(1/(2*sqrt(2.0d0)))*(4*(x-xmin)**2 -2)*exp(-(x-xmin)**2 /2.0d0)/sqrt(sqrt(pi))
!return
!end function


double precision function poten(x,ti,delt,nt)
IMPLICIT NONE
double precision :: x,ti,delt,weig,acoef,bcoef
double precision :: OMEGA=1.0d0,alpha=0.0005d0,beta=0.0001d0
integer :: nt
poten = 0.5*(((alpha-beta*cos(2*(ti+ delt*nt)))/(2*sqrt(alpha)))*x**4 - &
	(sqrt((alpha-beta*cos(2*(ti+ delt*nt))))/(2*sqrt(alpha)))*x**2) 
return
end function


!weig=(1.0d0-(beta/(4.0d0*alpha))*cos(2*delt*nt))
!acoef=(alpha-beta*cos(weig*(ti+ delt*nt)))/(2.0d0*sqrt(alpha))
!bcoef=sqrt((alpha-beta*cos(weig*(ti+ delt*nt))))/ (2.0d0*sqrt(alpha))
!poten=0.5*(acoef*x**4-bcoef*x**2)

!wavefun1 = exp(-(x+((sqrt((alpha-beta*cos(OMEGA*(ti+ delt*nt))))/(2*sqrt(alpha)))/(2*((alpha-beta*cos(OMEGA*(ti+ delt*nt)))/ &
!(2*sqrt(alpha)))))**0.5)**2 /2.0d0 )/sqrt(sqrt(pi)) !!!exp(-(x+(x0*(1+beta/(4*alpha)*cos(2*(ti+delt*nt)))))**2 /2)/sqrt(sqrt(pi)) !exp(-x**2 /2.0d0 )/sqrt(sqrt(pi))
!##########################################################################################################################











!####################################Potencial con curvatura dependiente del tiempo#########################################3

!double precision function poten(x,ti,delt,nt)
!IMPLICIT NONE
!double precision :: x,ti,delt,w0=1.0d0,g=0.1d0
!integer :: nt
!poten = 0.5*(w0**2 *(1-g*cos(2*w0*(ti+ delt*nt))))*x**2
!return
!end function


!double complex function wavefun1(x,x0,ti,delt,nt)
!double precision :: x,x0,ti,delt
!double precision :: pi = 4.0d0*atan(1.0d0),w0=1.0d0,eps
!double complex :: ii = (0.0d0,1.0d0)
!integer :: nt
!eps=sqrt(1.0d0/w0)
!wavefun1 = exp(-x**2 /(2.0d0*eps**2)  )/sqrt(sqrt(pi*eps**2))
!return
!end function

!double precision function exactsol(x,tf)
!implicit none
!double precision :: x,tf,g=0.1d0,w0=1.0d0!,w=(1.0d0**2 *(1.0d0-0.1d0*cos(2*tf)))!,w0=1.0d0
!double precision :: pi = 4.0d0*atan(1.0d0),Bt,Azpo,eps
!double complex ::Y,Z
!Y=(0.1332103626068909d0,0.1118914056674523d0) !para t=5 concuerda
!Z=(-0.5948081078849885d0,0.7081384245754117d0) !para t=5 concuerda
!Y=(-0.772651d0,-0.648996d0) !para t=10 concuerda
!Z=(-0.4912172386596277d0,0.5848101209406726d0) !para t=10 concuerda
!Y=(-0.4438689531865810d0,-0.3728322641891675d0) !para t=15 concuerda
!Z=(0.2501251220337612d0,-0.2977821040360947d0) !para t=15  concuerda
!Y=(0.8493813256513792d0,0.7134465263432939d0 ) !para t=25 concuerda
!Z=(0.3266698323813117d0,-0.3889110746680426d0) !para t=25  concuerda
!Y=(-1.298698110344277d0,-1.090854752290632d0) !para t=35 concuerda
!Z=(-0.9588315345952948d0,1.141520169850853d0) !para t=35  concuerda
!Y=(-3.751549061644548,-3.151151980394850d0) !para t=85 
!Z=(-3.08838d0,3.67682d0) !para t=85  
!Y=(-0.7726512531506787d0,-0.6489963176579165d0)!para t=5 con w0=2
!Z=(-0.4912172386596277d0,0.5848101209406726d0)!para t=5 con w0=2
!Bt=(1/sqrt(abs(Y)**2 + w0**2 *abs(Z)**2))
!Azpo=sqrt(cosh(g*w0*tf/2.0d0)/(2.0d0*w0))
!eps=Azpo*dsqrt(2.0d0)
!eps=sqrt(1.0d0/w0)
!exactsol = Bt/eps *exp(-x**2 *Bt**2 /eps**2)/sqrt(pi)
!return
!end function


!############################################################################################################






!double precision function poten(x,R0,V0,velx,E0,tau,ome,phi,nt,ti,delt)
!double precision :: x,R0, V0, velx, E0, tau, ome, phi,ti,delt
!if ( abs(x) <= R0) then
!    poten = V0 - E0*dexp(-(((ti-200.0d0/velx)+delt*nt)/tau)**2)*dcos(ome*((ti-200.0d0/velx)+delt*nt)+phi)*x
!else
!    poten = 0.0d0 - E0*dexp(-(((ti-200.0d0/velx)+delt*nt)/tau)**2)*dcos(ome*((ti-200.0d0/velx)+delt*nt)+phi)*x
!end if
!return
!end function
!Se declara la funcion de onda incidente
!double complex function wavefun1(x,x0,velx)
!double precision :: x,x0,velx
!double precision :: pi = 4.0d0*atan(1.0d0)
!double complex :: ii = (0.0d0,1.0d0)
!wavefun1 = exp(-((x-x0)/20.0d0)**2 + ii*velx*(x-x0))
!return
!end function
