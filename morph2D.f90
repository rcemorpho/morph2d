! *******************************************************************
! 2-Dimensional Morphodynamic Model
! 
! Written by Jongseok Cho and Peter A. Nelson
! Department of Civil and Environmental Engineering
! Colorado State University, Fort Collins, Colorado
! *******************************************************************
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

MODULE precision_kind
IMPLICIT none
integer,	parameter	:: ik = selected_int_kind	(r = 8)
integer,	parameter	:: rk = selected_real_kind	(p = 6)

END MODULE precision_kind
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

MODULE variables
USE precision_kind
IMPLICIT none

! PARAMETERS
integer	(kind=1),	parameter	:: unit_chnl	= 11
integer	(kind=1),	parameter	:: unit_conc	= 12
integer	(kind=1),	parameter	:: unit_var		= 13

integer	(kind=2),	parameter	:: mt	= 3600
integer	(kind=2),	parameter	:: mx	= 1204
integer	(kind=2),	parameter	:: my	= 40
integer	(kind=4),	parameter	:: mxy	= mx*my

real	(kind=rk),	parameter	:: pi	= ATAN(1.0)*4.0_rk
real	(kind=rk),	parameter	:: g	= 9.80665_rk
real	(kind=rk),	parameter	:: pore	= 0.44_rk
real	(kind=rk),	parameter	:: vonk	= 0.408_rk
real	(kind=rk),	parameter	:: visk	= 1.36e-06_rk
real	(kind=rk),	parameter	:: visd	= 1.0e-03_rk

! VARIABLE STATICS
integer	(kind=1)				:: RA,rkn,feed_s
integer	(kind=2)				:: nx,ny,ix,iy,i1,i9,j1,j9,i0,j0,feedxs,feedid,ixn(3),nw
integer	(kind=ik)				:: kcont,ntot,nff,nc_ini,npc,num_errs

real	(kind=4)				:: dx,dy,dt,LC,BC,etaland,etasea,slope,periodh
real	(kind=4)				:: davg,havg,qin,uin,dsmm,dsm
real	(kind=4)				:: rhow,rhos,Rsg,qwtot,scrate,rd,zbini
real	(kind=4)				:: sma,smd,time,qszero,tszero,Ca(2),ar
real	(kind=4)				:: dst,kb,ka,qstot,qs,Qc0,Qs0,thetac0,thetac1,rc(2)
real	(kind=4)				:: Cm0,zbar,zstd,fsigma,flambda,famp,maf

character	(25)				:: nome,idout
character	(8)					:: test,test1

! VARIABLES DYNAMICS
! (mt)
integer	(kind=ik), allocatable :: np		(:)

! (mx,my)
real	(kind=rk), allocatable :: XI		(:,:)
real	(kind=rk), allocatable :: YI		(:,:)

real	(kind=rk), allocatable :: FG0		(:,:)
real	(kind=rk), allocatable :: FG1		(:,:)
real	(kind=rk), allocatable :: FG2		(:,:)
real	(kind=rk), allocatable :: FG3		(:,:)

real	(kind=rk), allocatable :: F1		(:,:)
real	(kind=rk), allocatable :: F2		(:,:)
real	(kind=rk), allocatable :: F3		(:,:)

real	(kind=rk), allocatable :: G1		(:,:)
real	(kind=rk), allocatable :: G2		(:,:)
real	(kind=rk), allocatable :: G3		(:,:)

real	(kind=rk), allocatable :: dddx		(:,:)
real	(kind=rk), allocatable :: dddy		(:,:)
real	(kind=rk), allocatable :: dudx		(:,:)
real	(kind=rk), allocatable :: dudy		(:,:)
real	(kind=rk), allocatable :: dvdx		(:,:)
real	(kind=rk), allocatable :: dvdy		(:,:)
real	(kind=rk), allocatable :: dzdx		(:,:)
real	(kind=rk), allocatable :: dzdy		(:,:)
real	(kind=rk), allocatable :: dzds		(:,:)
real	(kind=rk), allocatable :: dzdn		(:,:)

real	(kind=rk), allocatable :: vis		(:,:)
real	(kind=rk), allocatable :: Dfx		(:,:)
real	(kind=rk), allocatable :: Dfy		(:,:)
real	(kind=rk), allocatable :: Sx		(:,:)
real	(kind=rk), allocatable :: Sy		(:,:)
real	(kind=rk), allocatable :: Sbx		(:,:)
real	(kind=rk), allocatable :: Sby		(:,:)
real	(kind=rk), allocatable :: Sf		(:,:)

real	(kind=rk), allocatable :: yd		(:,:)
real	(kind=rk), allocatable :: hd		(:,:)
real	(kind=rk), allocatable :: zd		(:,:)

real	(kind=rk), allocatable :: Pc		(:,:)
real	(kind=rk), allocatable :: Pe		(:,:)
real	(kind=rk), allocatable :: Cm		(:,:)
real	(kind=rk), allocatable :: cf		(:,:)
real	(kind=rk), allocatable :: cff		(:,:)
real	(kind=rk), allocatable :: cfg		(:,:)
real	(kind=rk), allocatable :: cft		(:,:)
real	(kind=rk), allocatable :: cfs		(:,:)
real	(kind=rk), allocatable :: ks		(:,:)
real	(kind=rk), allocatable :: ksa		(:,:)
real	(kind=rk), allocatable :: ksb		(:,:)
real	(kind=rk), allocatable :: ksf		(:,:)
real	(kind=rk), allocatable :: ksg		(:,:)
real	(kind=rk), allocatable :: kst		(:,:)
real	(kind=rk), allocatable :: ksc		(:,:)
real	(kind=rk), allocatable :: ksn		(:,:)
real	(kind=rk), allocatable :: Cs		(:,:)

real	(kind=rk), allocatable :: hini		(:,:)
real	(kind=rk), allocatable :: dini		(:,:)
real	(kind=rk), allocatable :: zini		(:,:)

real	(kind=rk), allocatable :: fx		(:,:)
real	(kind=rk), allocatable :: fy		(:,:)
real	(kind=rk), allocatable :: fbx		(:,:)
real	(kind=rk), allocatable :: fby		(:,:)

real	(kind=rk), allocatable :: h			(:,:)
real	(kind=rk), allocatable :: d			(:,:)
real	(kind=rk), allocatable :: du		(:,:)
real	(kind=rk), allocatable :: dv		(:,:)
real	(kind=rk), allocatable :: u			(:,:)
real	(kind=rk), allocatable :: v			(:,:)
real	(kind=rk), allocatable :: uv		(:,:)

real	(kind=rk), allocatable :: d1		(:,:)
real	(kind=rk), allocatable :: du1		(:,:)
real	(kind=rk), allocatable :: dv1		(:,:)

real	(kind=rk), allocatable :: d0		(:,:)
real	(kind=rk), allocatable :: du0		(:,:)
real	(kind=rk), allocatable :: dv0		(:,:)

real	(kind=rk), allocatable :: dn		(:,:)
real	(kind=rk), allocatable :: dun		(:,:)
real	(kind=rk), allocatable :: dvn		(:,:)

real	(kind=rk), allocatable :: vc		(:,:)
real	(kind=rk), allocatable :: fp		(:,:)
real	(kind=rk), allocatable :: fm		(:,:)

real	(kind=rk), allocatable :: zra		(:,:)
real	(kind=rk), allocatable :: zavg		(:,:)
real	(kind=rk), allocatable :: zbr		(:,:)
real	(kind=rk), allocatable :: zbr0		(:,:)

real	(kind=rk), allocatable :: vbain		(:,:)
real	(kind=rk), allocatable :: us		(:,:)
real	(kind=rk), allocatable :: usx		(:,:)
real	(kind=rk), allocatable :: usy		(:,:)

real	(kind=rk), allocatable :: z			(:,:)
real	(kind=rk), allocatable :: zb		(:,:)
real	(kind=rk), allocatable :: zba		(:,:)
real	(kind=rk), allocatable :: vb		(:,:)
real	(kind=rk), allocatable :: vba		(:,:)
real	(kind=rk), allocatable :: vbc		(:,:)

real	(kind=rk), allocatable :: zban		(:,:)
real	(kind=rk), allocatable :: vbn		(:,:)
real	(kind=rk), allocatable :: vban		(:,:)

real	(kind=rk), allocatable :: zb0		(:,:)
real	(kind=rk), allocatable :: zba0		(:,:)
real	(kind=rk), allocatable :: vb0		(:,:)
real	(kind=rk), allocatable :: vba0		(:,:)
real	(kind=rk), allocatable :: vbc0		(:,:)

real	(kind=rk), allocatable :: qbc		(:,:)
real	(kind=rk), allocatable :: qb		(:,:)
real	(kind=rk), allocatable :: qbx		(:,:)
real	(kind=rk), allocatable :: qby		(:,:)
real	(kind=rk), allocatable :: qbe		(:,:)
real	(kind=rk), allocatable :: qbn		(:,:)

real	(kind=rk), allocatable :: theta		(:,:)
real	(kind=rk), allocatable :: thetacr	(:,:)
real	(kind=rk), allocatable :: thetacg	(:,:)
real	(kind=rk), allocatable :: phir		(:,:)
real	(kind=rk), allocatable :: phit		(:,:)
real	(kind=rk), allocatable :: rpf		(:,:)

real	(kind=rk), allocatable :: alpha		(:,:)
real	(kind=rk), allocatable :: alphax	(:,:)
real	(kind=rk), allocatable :: alphay	(:,:)
real	(kind=rk), allocatable :: uvx		(:,:)
real	(kind=rk), allocatable :: uvy		(:,:)

real	(kind=rk), allocatable :: M01		(:,:)
real	(kind=rk), allocatable :: M02		(:,:)
real	(kind=rk), allocatable :: M03		(:,:)
real	(kind=rk), allocatable :: M04		(:,:)
real	(kind=rk), allocatable :: M05		(:,:)
real	(kind=rk), allocatable :: vzb		(:,:)

! (mx,my,:)
real	(kind=rk), allocatable :: hL		(:,:,:)
real	(kind=rk), allocatable :: dL		(:,:,:)
real	(kind=rk), allocatable :: duL		(:,:,:)
real	(kind=rk), allocatable :: dvL		(:,:,:)
real	(kind=rk), allocatable :: uL		(:,:,:)
real	(kind=rk), allocatable :: vL		(:,:,:)
real	(kind=rk), allocatable :: zL		(:,:,:)

real	(kind=rk), allocatable :: hR		(:,:,:)
real	(kind=rk), allocatable :: dR		(:,:,:)
real	(kind=rk), allocatable :: duR		(:,:,:)
real	(kind=rk), allocatable :: dvR		(:,:,:)
real	(kind=rk), allocatable :: uR		(:,:,:)
real	(kind=rk), allocatable :: vR		(:,:,:)
real	(kind=rk), allocatable :: zR		(:,:,:)

real	(kind=rk), allocatable :: hE		(:,:,:)
real	(kind=rk), allocatable :: zE		(:,:,:)
real	(kind=rk), allocatable :: dm		(:,:,:)
real	(kind=rk), allocatable :: um		(:,:,:)
real	(kind=rk), allocatable :: vm		(:,:,:)

real	(kind=rk), allocatable :: SW		(:,:,:)
real	(kind=rk), allocatable :: U1		(:,:,:)
real	(kind=rk), allocatable :: U2		(:,:,:)
real	(kind=rk), allocatable :: U3		(:,:,:)
real	(kind=rk), allocatable :: U4		(:,:,:)

contains
!--------------------------------------------------------------------------------------------
	real (kind=rk)	function sgn(r1)
	real (kind=rk),	intent(in)	:: r1
	! sign of a parameter
	if (r1.GT.0.) then
		sgn = 1.
	elseif (r1.LT.0.) then
		sgn = -1.
	else
		sgn = 0.
	endif
	return
	end function sgn
!--------------------------------------------------------------------------------------------
	real (kind=rk)	function minmod1(r1)
	real (kind=rk),	intent(in)	:: r1
	! flux limiter: minmod_1
	minmod1 = max(0.d0,min(r1,1.d0))
	return
	end function minmod1
!--------------------------------------------------------------------------------------------
	real (kind=rk)	function minmod2(r1,r2)
	real (kind=rk), intent(in)	:: r1,r2
	real (kind=rk)				:: lim1
	! flux limiter: minmod_2
	if ((r1*r2).le.0.0) then
		lim1 = 0.0
	else
		if (abs(r1).lt.abs(r2)) then
			lim1 = r1
		elseif (abs(r2).lt.abs(r1)) then
			lim1 = r2
		endif
	endif
	minmod2 = lim1
	return
	end function minmod2
!--------------------------------------------------------------------------------------------
	real (kind=rk)	function superbee(r1)
	real (kind=rk), intent(in)	:: r1
	! flux limiter: superbee
	superbee = max(0.0,min(2.0*r1,1.0),min(r1,2.0))
	return
	end function superbee
!--------------------------------------------------------------------------------------------
	real (kind=rk)	function VanAlbada1(r1)
	real (kind=rk), intent(in)	:: r1
	! flux limiter: van albada_1
	if (r1.le.0.0) then
		VanAlbada1 = 0.0
	else
		VanAlbada1 = (r1**2.+r1)/(1.0+r1**2.)
	endif
	return
	end function VanAlbada1
!--------------------------------------------------------------------------------------------
	real (kind=rk)	function VanAlbada2(r1)
	real (kind=rk), intent(in)	:: r1
	! flux limiter: van albada_2
	VanAlbada2 = (2.0*r1)/(1.0+r1**2.)
	return
	end function VanAlbada2
!--------------------------------------------------------------------------------------------
	real (kind=rk)	function VanLeer(r1)
	real (kind=rk), intent(in)	:: r1
	! flux limiter: van leer
	VanLeer = (r1+abs(r1))/(1.0+abs(r1))
	return
	end function VanLeer
!--------------------------------------------------------------------------------------------
	real (kind=rk)	function Koren(r1)
	real (kind=rk), intent(in)	:: r1
	! flux limiter: Koren
	Koren = max(0.0,min(2.0*r1,(2.0+r1)/3.0,2.0))
	return
	end function Koren
!--------------------------------------------------------------------------------------------
	function roots3(a,Fr1)
	real (kind=rk), intent(in)	:: a(3),Fr1
	real (kind=rk)				:: roots3,c0,c1(2),c2(2),x1(3),r1
	! cubic equation solver
	! x^3 + a(1)x^2 + a(2)x + a(3) = 0
	c1(1) = (a(1)**2.-3.0*a(2))/9.
	c1(2) = (2.*a(1)**3.-9.*a(1)*a(2)+27.*a(3))/54.
	
	if (c1(1)**3 .GT. c1(2)**2) then
		c0 = ACOS(MIN(MAX(c1(2)/c1(1)**1.5,-1.0),1.0))
		x1(1) = -2.*SQRT(c1(1))*COS(c0/3.)-a(1)/3.
		x1(2) = -2.*SQRT(c1(1))*COS((c0+2.*pi)/3.)-a(1)/3.
		x1(3) = -2.*SQRT(c1(1))*COS((c0-2.*pi)/3.)-a(1)/3.
	else
		c2(1) = -sgn(c1(2))*(ABS(c1(2))+SQRT(c1(2)**2.-c1(1)**3.))**(1./3.)
		if (c2(1) .NE. 0.) then
			c2(2) = c1(1)/c2(1)
		else
			c2(2) = 0.
		endif
		x1(1) = SUM(c2)-a(1)/3.
		! x1(2) = hypot(-0.5*(c0+c3)-a(1)/3.,(c2(1)-c2(2))*sqrt(3.0)/2.)
		! x1(3) = hypot(-0.5*(c0+c3)-a(1)/3.,-(c2(1)-c2(2))*sqrt(3.0)/2.)
		x1(2) = -0.5*SUM(c2)-a(1)/3.
		x1(3) = -0.5*SUM(c2)-a(1)/3.
	endif
	
	if (Fr1 .GT. 1.) then
		r1 = MINVAL(x1, MASK = x1 .GT. 0.0)
	else
		r1 = MAXVAL(x1, MASK = x1 .GT. 0.0)
	endif
	
	roots3 = r1
	return
	end function roots3
!--------------------------------------------------------------------------------------------
	real (kind=rk)	function WAF(r1,qL,qR)
	real (kind=rk), intent(in)	:: r1,qL(3),qR(3)
	real (kind=rk)				:: lim1
	! WAF limiter function
	lim1 = 1.0
	if ((qR(2)-qL(2)).ne.0.0) then
		if (r1.gt.0.0) then
			lim1 = (qR(1)-qL(1))/(qR(2)-qL(2))
		elseif (r1.lt.0.0) then
			lim1 = (qR(3)-qL(3))/(qR(2)-qL(2))
		endif
	endif
	if (isnan(lim1)) lim1 = 1.0
	WAF = 1.0-(1.0-abs(r1))*minmod1(lim1)
	return
	end function WAF
!--------------------------------------------------------------------------------------------
	function SLMR(dd,uu,r1)
	real (kind=rk), intent(in)	:: dd(3),uu(3),r1
	real (kind=rk)				:: SLMR(3),ss(3),a1(2),s1
	! middle wave speed
	ss(1) = min(uu(1)-sqrt(g*dd(1)),uu(2)-sqrt(g*dd(2)))
	ss(3) = max(uu(3)+sqrt(g*dd(3)),uu(2)+sqrt(g*dd(2)))
	if (r1==0.0) then
		s1 = ss(1)*dd(3)*(uu(3)-ss(3))-ss(3)*dd(1)*(uu(1)-ss(1))
		ss(2) = s1/(dd(3)*(uu(3)-ss(3))-dd(1)*(uu(1)-ss(1)))
	else
		s1 = dd(3)*uu(3)*(ss(3)-uu(3))-dd(1)*uu(1)*(ss(1)-uu(1))-r1
		ss(2) = s1/(dd(3)*(ss(3)-uu(3))-dd(1)*(ss(1)-uu(1)))
	endif
	SLMR = ss
	return
	end function SLMR
!--------------------------------------------------------------------------------------------
	function fdir(aa)
	real (kind=rk), intent(in)	:: aa(2)
	real (kind=rk)				:: fdir,aa2
	! flow direction determination
	if (aa(1)==0. .AND. aa(2)==0.) then
		fdir = 0.
	else
		aa2 = PRODUCT(aa,MASK = aa.NE.0.)
		if (aa2.GE.0.) then
			fdir = 1.
		elseif (aa2.LT.0.) then
			fdir = -1.
		endif
	endif
	return
	end function fdir
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

END MODULE variables
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE ALLOCATION
USE variables
IMPLICIT none
integer	(kind=ik)	:: iostat
! variable allocation
allocate ( np		(mt), stat = iostat)

allocate ( XI		(mx,my), stat = iostat)
allocate ( YI		(my,my), stat = iostat)

allocate ( FG0		(mx,my), stat = iostat)
allocate ( FG1		(mx,my), stat = iostat)
allocate ( FG2		(mx,my), stat = iostat)
allocate ( FG3		(mx,my), stat = iostat)

allocate ( F1		(mx,my), stat = iostat)
allocate ( F2		(mx,my), stat = iostat)
allocate ( F3		(mx,my), stat = iostat)

allocate ( G1		(mx,my), stat = iostat)
allocate ( G2		(mx,my), stat = iostat)
allocate ( G3		(mx,my), stat = iostat)

allocate ( dddx		(mx,my), stat = iostat)
allocate ( dddy		(mx,my), stat = iostat)
allocate ( dudx		(mx,my), stat = iostat)
allocate ( dudy		(mx,my), stat = iostat)
allocate ( dvdx		(mx,my), stat = iostat)
allocate ( dvdy		(mx,my), stat = iostat)
allocate ( dzdx		(mx,my), stat = iostat)
allocate ( dzdy		(mx,my), stat = iostat)
allocate ( dzds		(mx,my), stat = iostat)
allocate ( dzdn		(mx,my), stat = iostat)

allocate ( vis		(mx,my), stat = iostat)
allocate ( Dfx		(mx,my), stat = iostat)
allocate ( Dfy		(mx,my), stat = iostat)
allocate ( Sx		(mx,my), stat = iostat)
allocate ( Sy		(mx,my), stat = iostat)
allocate ( Sbx		(mx,my), stat = iostat)
allocate ( Sby		(mx,my), stat = iostat)
allocate ( Sf		(mx,my), stat = iostat)

allocate ( yd		(mx,my), stat = iostat)
allocate ( hd		(mx,my), stat = iostat)
allocate ( zd		(mx,my), stat = iostat)

allocate ( Pc		(mx,my), stat = iostat)
allocate ( Pe		(mx,my), stat = iostat)
allocate ( Cm		(mx,my), stat = iostat)
allocate ( cf		(mx,my), stat = iostat)
allocate ( cff		(mx,my), stat = iostat)
allocate ( cfg		(mx,my), stat = iostat)
allocate ( cft		(mx,my), stat = iostat)
allocate ( cfs		(mx,my), stat = iostat)
allocate ( ks		(mx,my), stat = iostat)
allocate ( ksa		(mx,my), stat = iostat)
allocate ( ksb		(mx,my), stat = iostat)
allocate ( ksf		(mx,my), stat = iostat)
allocate ( ksg		(mx,my), stat = iostat)
allocate ( kst		(mx,my), stat = iostat)
allocate ( ksc		(mx,my), stat = iostat)
allocate ( ksn		(mx,my), stat = iostat)
allocate ( Cs		(mx,my), stat = iostat)

allocate ( hini		(mx,my), stat = iostat)
allocate ( dini		(mx,my), stat = iostat)
allocate ( zini		(mx,my), stat = iostat)

allocate ( fx		(mx,my), stat = iostat)
allocate ( fy		(mx,my), stat = iostat)
allocate ( fbx		(mx,my), stat = iostat)
allocate ( fby		(mx,my), stat = iostat)

allocate ( h		(mx,my), stat = iostat)
allocate ( d		(mx,my), stat = iostat)
allocate ( du		(mx,my), stat = iostat)
allocate ( dv		(mx,my), stat = iostat)
allocate ( u		(mx,my), stat = iostat)
allocate ( v		(mx,my), stat = iostat)
allocate ( uv		(mx,my), stat = iostat)

allocate ( d1		(mx,my), stat = iostat)
allocate ( du1		(mx,my), stat = iostat)
allocate ( dv1		(mx,my), stat = iostat)

allocate ( d0		(mx,my), stat = iostat)
allocate ( du0		(mx,my), stat = iostat)
allocate ( dv0		(mx,my), stat = iostat)

allocate ( dn		(mx,my), stat = iostat)
allocate ( dun		(mx,my), stat = iostat)
allocate ( dvn		(mx,my), stat = iostat)

allocate ( vc		(mx,my), stat = iostat)
allocate ( fp		(mx,my), stat = iostat)
allocate ( fm		(mx,my), stat = iostat)

allocate ( zra		(mx,my), stat = iostat)
allocate ( zavg		(mx,my), stat = iostat)
allocate ( zbr		(mx,my), stat = iostat)
allocate ( zbr0		(mx,my), stat = iostat)

allocate ( vbain	(mx,my), stat = iostat)
allocate ( us		(mx,my), stat = iostat)
allocate ( usx		(mx,my), stat = iostat)
allocate ( usy		(mx,my), stat = iostat)

allocate ( z		(mx,my), stat = iostat)
allocate ( zb		(mx,my), stat = iostat)
allocate ( zba		(mx,my), stat = iostat)
allocate ( vb		(mx,my), stat = iostat)
allocate ( vba		(mx,my), stat = iostat)
allocate ( vbc		(mx,my), stat = iostat)

allocate ( zban		(mx,my), stat = iostat)
allocate ( vbn		(mx,my), stat = iostat)
allocate ( vban		(mx,my), stat = iostat)

allocate ( zb0		(mx,my), stat = iostat)
allocate ( zba0		(mx,my), stat = iostat)
allocate ( vb0		(mx,my), stat = iostat)
allocate ( vba0		(mx,my), stat = iostat)
allocate ( vbc0		(mx,my), stat = iostat)

allocate ( qbc		(mx,my), stat = iostat)
allocate ( qb		(mx,my), stat = iostat)
allocate ( qbx		(mx,my), stat = iostat)
allocate ( qby		(mx,my), stat = iostat)
allocate ( qbe		(mx,my), stat = iostat)
allocate ( qbn		(mx,my), stat = iostat)

allocate ( theta	(mx,my), stat = iostat)
allocate ( thetacr	(mx,my), stat = iostat)
allocate ( thetacg	(mx,my), stat = iostat)
allocate ( phir		(mx,my), stat = iostat)
allocate ( phit		(mx,my), stat = iostat)
allocate ( rpf		(mx,my), stat = iostat)

allocate ( alpha	(mx,my), stat = iostat)
allocate ( alphax	(mx,my), stat = iostat)
allocate ( alphay	(mx,my), stat = iostat)
allocate ( uvx		(mx,my), stat = iostat)
allocate ( uvy		(mx,my), stat = iostat)

allocate ( M01		(mx,my), stat = iostat)
allocate ( M02		(mx,my), stat = iostat)
allocate ( M03		(mx,my), stat = iostat)
allocate ( M04		(mx,my), stat = iostat)
allocate ( M05		(mx,my), stat = iostat)
allocate ( vzb		(mx,my), stat = iostat)

allocate ( hL		(mx,my,2), stat = iostat)
allocate ( dL		(mx,my,2), stat = iostat)
allocate ( duL		(mx,my,2), stat = iostat)
allocate ( dvL		(mx,my,2), stat = iostat)
allocate ( uL		(mx,my,2), stat = iostat)
allocate ( vL		(mx,my,2), stat = iostat)
allocate ( zL		(mx,my,2), stat = iostat)

allocate ( hR		(mx,my,2), stat = iostat)
allocate ( dR		(mx,my,2), stat = iostat)
allocate ( duR		(mx,my,2), stat = iostat)
allocate ( dvR		(mx,my,2), stat = iostat)
allocate ( uR		(mx,my,2), stat = iostat)
allocate ( vR		(mx,my,2), stat = iostat)
allocate ( zR		(mx,my,2), stat = iostat)

allocate ( hE		(mx,my,2), stat = iostat)
allocate ( zE		(mx,my,2), stat = iostat)
allocate ( dm		(mx,my,2), stat = iostat)
allocate ( um		(mx,my,2), stat = iostat)
allocate ( vm		(mx,my,2), stat = iostat)

allocate ( SW		(mx,my,3), stat = iostat)
allocate ( U1		(mx,my,3), stat = iostat)
allocate ( U2		(mx,my,3), stat = iostat)
allocate ( U3		(mx,my,3), stat = iostat)
allocate ( U4		(mx,my,3), stat = iostat)

if (iostat/=0) then
	write(*,*) "Error Allocation"
	stop
endif

END SUBROUTINE ALLOCATION
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE INITIALIZATION
USE variables
IMPLICIT none
! variable initialization
np		(1:mt) = 0.0_ik

XI		(1:mx,1:my) = 0.0_rk
YI		(1:my,1:my) = 0.0_rk

FG0		(1:mx,1:my) = 0.0_rk
FG1		(1:mx,1:my) = 0.0_rk
FG2		(1:mx,1:my) = 0.0_rk
FG3		(1:mx,1:my) = 0.0_rk

F1		(1:mx,1:my) = 0.0_rk
F2		(1:mx,1:my) = 0.0_rk
F3		(1:mx,1:my) = 0.0_rk

G1		(1:mx,1:my) = 0.0_rk
G2		(1:mx,1:my) = 0.0_rk
G3		(1:mx,1:my) = 0.0_rk

dddx	(1:mx,1:my) = 0.0_rk
dddy	(1:mx,1:my) = 0.0_rk
dudx	(1:mx,1:my) = 0.0_rk
dudy	(1:mx,1:my) = 0.0_rk
dvdx	(1:mx,1:my) = 0.0_rk
dvdy	(1:mx,1:my) = 0.0_rk
dzdx	(1:mx,1:my) = 0.0_rk
dzdy	(1:mx,1:my) = 0.0_rk
dzds	(1:mx,1:my) = 0.0_rk
dzdn	(1:mx,1:my) = 0.0_rk

vis		(1:mx,1:my) = 0.0_rk
Dfx		(1:mx,1:my) = 0.0_rk
Dfy		(1:mx,1:my) = 0.0_rk
Sx		(1:mx,1:my) = 0.0_rk
Sy		(1:mx,1:my) = 0.0_rk
Sbx		(1:mx,1:my) = 0.0_rk
Sby		(1:mx,1:my) = 0.0_rk
Sf		(1:mx,1:my) = 0.0_rk

yd		(1:mx,1:my) = 0.0_rk
hd		(1:mx,1:my) = 0.0_rk
zd		(1:mx,1:my) = 0.0_rk

Pc		(1:mx,1:my) = 0.0_rk
Pe		(1:mx,1:my) = 0.0_rk
Cm		(1:mx,1:my) = 0.0_rk
cf		(1:mx,1:my) = 0.0_rk
cff		(1:mx,1:my) = 0.0_rk
cfg		(1:mx,1:my) = 0.0_rk
cft		(1:mx,1:my) = 0.0_rk
cfs		(1:mx,1:my) = 0.0_rk
ks		(1:mx,1:my) = 0.0_rk
ksa		(1:mx,1:my) = 0.0_rk
ksb		(1:mx,1:my) = 0.0_rk
ksf		(1:mx,1:my) = 0.0_rk
ksg		(1:mx,1:my) = 0.0_rk
kst		(1:mx,1:my) = 0.0_rk
ksc		(1:mx,1:my) = 0.0_rk
ksn		(1:mx,1:my) = 0.0_rk
Cs		(1:mx,1:my) = 0.0_rk

hini	(1:mx,1:my) = 0.0_rk
dini	(1:mx,1:my) = 0.0_rk
zini	(1:mx,1:my) = 0.0_rk

fx		(1:mx,1:my) = 0.0_rk
fy		(1:mx,1:my) = 0.0_rk
fbx		(1:mx,1:my) = 0.0_rk
fby		(1:mx,1:my) = 0.0_rk

h		(1:mx,1:my) = 0.0_rk
d		(1:mx,1:my) = 0.0_rk
du		(1:mx,1:my) = 0.0_rk
dv		(1:mx,1:my) = 0.0_rk
u		(1:mx,1:my) = 0.0_rk
v		(1:mx,1:my) = 0.0_rk
uv		(1:mx,1:my) = 0.0_rk

d1		(1:mx,1:my) = 0.0_rk
du1		(1:mx,1:my) = 0.0_rk
dv1		(1:mx,1:my) = 0.0_rk

d0		(1:mx,1:my) = 0.0_rk
du0		(1:mx,1:my) = 0.0_rk
dv0		(1:mx,1:my) = 0.0_rk

dn		(1:mx,1:my) = 0.0_rk
dun		(1:mx,1:my) = 0.0_rk
dvn		(1:mx,1:my) = 0.0_rk

vc		(1:mx,1:my) = 0.0_rk
fp		(1:mx,1:my) = 0.0_rk
fm		(1:mx,1:my) = 0.0_rk

zra		(1:mx,1:my) = 0.0_rk
zavg	(1:mx,1:my) = 0.0_rk
zbr		(1:mx,1:my) = 0.0_rk
zbr0	(1:mx,1:my) = 0.0_rk

vbain	(1:mx,1:my) = 0.0_rk
us		(1:mx,1:my) = 0.0_rk
usx		(1:mx,1:my) = 0.0_rk
usy		(1:mx,1:my) = 0.0_rk

z		(1:mx,1:my) = 0.0_rk
zb		(1:mx,1:my) = 0.0_rk
zba		(1:mx,1:my) = 0.0_rk
vb		(1:mx,1:my) = 0.0_rk
vba		(1:mx,1:my) = 0.0_rk
vbc		(1:mx,1:my) = 0.0_rk

zban	(1:mx,1:my) = 0.0_rk
vbn		(1:mx,1:my) = 0.0_rk
vban	(1:mx,1:my) = 0.0_rk

zb0		(1:mx,1:my) = 0.0_rk
zba0	(1:mx,1:my) = 0.0_rk
vb0		(1:mx,1:my) = 0.0_rk
vba0	(1:mx,1:my) = 0.0_rk
vbc0	(1:mx,1:my) = 0.0_rk

qbc		(1:mx,1:my) = 0.0_rk
qb		(1:mx,1:my) = 0.0_rk
qbx		(1:mx,1:my) = 0.0_rk
qby		(1:mx,1:my) = 0.0_rk
qbe		(1:mx,1:my) = 0.0_rk
qbn		(1:mx,1:my) = 0.0_rk

theta	(1:mx,1:my) = 0.0_rk
thetacr	(1:mx,1:my) = 0.0_rk
thetacg	(1:mx,1:my) = 0.0_rk
phir	(1:mx,1:my) = 0.0_rk
phit	(1:mx,1:my) = 0.0_rk
rpf		(1:mx,1:my) = 0.0_rk

alpha	(1:mx,1:my) = 0.0_rk
alphax	(1:mx,1:my) = 0.0_rk
alphay	(1:mx,1:my) = 0.0_rk
uvx		(1:mx,1:my) = 0.0_rk
uvy		(1:mx,1:my) = 0.0_rk

M01		(1:mx,1:my) = 0.0_rk
M02		(1:mx,1:my) = 0.0_rk
M03		(1:mx,1:my) = 0.0_rk
M04		(1:mx,1:my) = 0.0_rk
M05		(1:mx,1:my) = 0.0_rk
vzb		(1:mx,1:my) = 0.0_rk

hL		(1:mx,1:my,1:2) = 0.0_rk
dL		(1:mx,1:my,1:2) = 0.0_rk
duL		(1:mx,1:my,1:2) = 0.0_rk
dvL		(1:mx,1:my,1:2) = 0.0_rk
uL		(1:mx,1:my,1:2) = 0.0_rk
vL		(1:mx,1:my,1:2) = 0.0_rk
zL		(1:mx,1:my,1:2) = 0.0_rk

hR		(1:mx,1:my,1:2) = 0.0_rk
dR		(1:mx,1:my,1:2) = 0.0_rk
duR		(1:mx,1:my,1:2) = 0.0_rk
dvR		(1:mx,1:my,1:2) = 0.0_rk
uR		(1:mx,1:my,1:2) = 0.0_rk
vR		(1:mx,1:my,1:2) = 0.0_rk
zR		(1:mx,1:my,1:2) = 0.0_rk

hE		(1:mx,1:my,1:2) = 0.0_rk
zE		(1:mx,1:my,1:2) = 0.0_rk
dm		(1:mx,1:my,1:2) = 0.0_rk
um		(1:mx,1:my,1:2) = 0.0_rk
vm		(1:mx,1:my,1:2) = 0.0_rk

SW		(1:mx,1:my,1:3) = 0.0_rk
U1		(1:mx,1:my,1:3) = 0.0_rk
U2		(1:mx,1:my,1:3) = 0.0_rk
U3		(1:mx,1:my,1:3) = 0.0_rk
U4		(1:mx,1:my,1:3) = 0.0_rk

END SUBROUTINE INITIALIZATION
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE DEALLOCATION
USE variables
IMPLICIT none
! variable deallocation
deallocate (np		)

deallocate (XI		)
deallocate (YI		)

deallocate (FG0		)
deallocate (FG1		)
deallocate (FG2		)
deallocate (FG3		)

deallocate (F1		)
deallocate (F2		)
deallocate (F3		)

deallocate (G1		)
deallocate (G2		)
deallocate (G3		)

deallocate (dddx	)
deallocate (dddy	)
deallocate (dudx	)
deallocate (dudy	)
deallocate (dvdx	)
deallocate (dvdy	)
deallocate (dzdx	)
deallocate (dzdy	)
deallocate (dzds	)
deallocate (dzdn	)

deallocate (vis		)
deallocate (Dfx		)
deallocate (Dfy		)
deallocate (Sx		)
deallocate (Sy		)
deallocate (Sbx		)
deallocate (Sby		)
deallocate (Sf		)

deallocate (yd		)
deallocate (hd		)
deallocate (zd		)

deallocate (Pc		)
deallocate (Pe		)
deallocate (Cm	 	)
deallocate (cf		)
deallocate (cff		)
deallocate (cfg		)
deallocate (cft		)
deallocate (cfs		)
deallocate (ks		)
deallocate (ksa		)
deallocate (ksb		)
deallocate (ksf		)
deallocate (ksg		)
deallocate (kst		)
deallocate (ksc		)
deallocate (ksn		)
deallocate (Cs		)

deallocate (hini	)
deallocate (dini	)
deallocate (zini	)

deallocate (fx		)
deallocate (fy		)
deallocate (fbx		)
deallocate (fby		)

deallocate (h		)
deallocate (d		)
deallocate (du		)
deallocate (dv		)
deallocate (u		)
deallocate (v		)
deallocate (uv		)

deallocate (d1		)
deallocate (du1		)
deallocate (dv1		)

deallocate (d0		)
deallocate (du0		)
deallocate (dv0		)

deallocate (dn		)
deallocate (dun		)
deallocate (dvn		)

deallocate (vc		)
deallocate (fp		)
deallocate (fm		)

deallocate (zra		)
deallocate (zavg	)
deallocate (zbr		)
deallocate (zbr0	)

deallocate (vbain	)
deallocate (vbc		)
deallocate (us		)
deallocate (usx		)
deallocate (usy		)

deallocate (z		)
deallocate (zb		)
deallocate (zba		)
deallocate (vb		)
deallocate (vba		)
deallocate (Vbc		)

deallocate (zban	)
deallocate (vbn		)
deallocate (vban	)

deallocate (zb0		)
deallocate (zba0	)
deallocate (vb0		)
deallocate (vba0	)
deallocate (vbc0	)

deallocate (qbc		)
deallocate (qb		)
deallocate (qbx		)
deallocate (qby		)
deallocate (qbe		)
deallocate (qbn		)

deallocate (theta	)
deallocate (thetacr	)
deallocate (thetacg	)
deallocate (phir	)
deallocate (phit	)
deallocate (rpf		)

deallocate (alpha	)
deallocate (alphax	)
deallocate (alphay	)
deallocate (uvx		)
deallocate (uvy		)

deallocate (M01		)
deallocate (M02		)
deallocate (M03		)
deallocate (M04		)
deallocate (M05		)
deallocate (vzb		)

deallocate (hL		)
deallocate (dL		)
deallocate (duL		)
deallocate (dvL		)
deallocate (uL		)
deallocate (vL		)
deallocate (zL		)

deallocate (hR		)
deallocate (dR		)
deallocate (duR		)
deallocate (dvR		)
deallocate (uR		)
deallocate (vR		)
deallocate (zR		)

deallocate (hE		)
deallocate (zE		)
deallocate (dm		)
deallocate (um		)
deallocate (vm		)

deallocate (SW		)
deallocate (U1		)
deallocate (U2		)
deallocate (U3		)
deallocate (U4		)

END SUBROUTINE DEALLOCATION
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE INPUT
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j,nplot

! start time for beginning from previous simulation (digits from previous filename)
nc_ini = 0
! seaward for linear slope (m)
etasea = 1.0
! length of the channel (m)
LC = 30.0
! width of the channel (m)
BC = 0.9
! Random Abrasion	(0 = no, 1 = yes)
RA = 1
! grain roughness height coefficient (ksa ~ rd*Ds)
rd = 1.0

! 2-B >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
slope = 0.0200
qwtot = 55.0e-03
kb = 3.4303e-03
dsmm = 7.0
Ca = (/1./vonk,11.0/)

! davg = 0.05
! zbini = 0.000
! qstot = 350.0/2.65e06
! scrate = 1.00

davg = 0.06
zbini = 0.020
qstot = 110.0/2.65e06
scrate = 0.60
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! critical Shields parameters (rc = 0.03**(1/0.6) = 0.0029)
thetac0 = 0.01; thetac1 = 0.60; rc = 0.0029
! number of cells in x- and y-direction
nx = mx-4; ny = my-4
! water and sediment density (kg/m^3)
rhow = 1027.0; rhos = 2650.0
! total simulation duration (hours *changed*)
periodh = 20.
! time step (sec)
dt = 0.008
! number of time steps between printouts (nff *changed*) dt*nff = print output (sec)
nff = 10/dt
! small value for non-zero flow and alluvial depth
sma = 1.0e-04; smd = 1.0e-08
! ==================================================================

! number of cells in x- and y-direction
ix = mx; iy = my
i0 = 2; i1 = 3; i9 = ix-2
j0 = 2; j1 = 3; j9 = iy-2
ixn = (/i1+200,i9-200,i9-i1-399/)

! ------------------------------------------------------------------
! check vectors dimensions
! ------------------------------------------------------------------
if (ix.GT.mx .OR. iy.GT.my) then
	stop
	write(6,*)'increase the size of the vectors mx/my!!!'
endif

! cells along (nx) / across (ny) channel
dx = LC/nx; dy = BC/ny

do i = 1,ix
	XI(i,1) = dx*(real(i)-2.5)
enddo
do j = 1,iy
	YI(1,j) = dy*(real(j)-0.5*(iy+1.0))
enddo
XI = SPREAD(XI(1:ix,1),DIM=2,NCOPIES=iy); YI = SPREAD(YI(1,1:iy),DIM=1,NCOPIES=ix)

! sediment feed distribution pattern
feed_s = 4
fsigma = 80.0; flambda = 0.001; famp = 1.0
feedxs = 13; feedid = 5
nw = 2							! half window size (ws = 2*nw+1)

! sediments parameters
dsm = dsmm*1.0e-3				! median diameter [m]
Rsg = (rhos-rhow)/rhow			! relative sediment and water density ratio  
tszero = g*Rsg*dsm				! dimensionless load discharge
qszero = SQRT(tszero)*dsm

Cm0 = pi*dsm/6./(1.-pore)		! packing concentration Cm
ka = rd*dsm; ksa = ka; ksb = kb	! grain and bedrock roughness height [m]
dst = 0.0204*log(100.*dsm)**2.+0.0220*log(100.*dsm)+0.0709 ! funtion of grain size in kst

! plotting intervals
num_errs = 0
ntot = NINT(periodh/dt)*3600	! number of time steps
nplot = ntot/nff
npc = 1
np = (/1:nplot/)*nff

if (np(nplot).ne.ntot) then
	np(nplot+1) = ntot
endif

END SUBROUTINE input
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SetBed
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j

! upstream end elevation
etaland = LC*slope+etasea

! flat bedrock bed
zbr0(:,1) = etaland-XI(:,1)*slope
zbr0 = spread(zbr0(:,1),2,iy)

! import randomly abraded bedrock
if (RA==1) then
	call BEDROCK
endif

zbr = zbr0+zra
h = zbr0+zbini+davg

if (zbini==0.) then
	zb = 0.	
else
	! intial alluvial bed perturbation
	zb0 = zbini+0.1*zra
	
	! initial upstream perturbation (Lanzoni, 2000)
	zb0(51:60,j1:j9) = zb0(51:60,j1:j9)+0.006*SIN(pi*XI(3:12,j1:j9)/0.5)*SIN(pi*YI(3:12,j1:j9)/BC)
	
	z = zbr0+zb0
	zb = z-zbr
endif

WHERE (zb .LT. 0.) zb = 0.
z = zbr+zb
d = h-z

zd = zra+zb; hd = zd+d
dini = d; hini = hd; zini = zd; havg = SUM(hd(ixn(1):ixn(2),j1:j9))/(ixn(3)*ny)

du = qwtot/BC
dv = 0.0

zba = zb*(1.0-pore)
vba = vb+zba
vbc = 1.0e-15

! set boundary values
! W-E
d(i1-2:i1-1,j1:j9) = spread(d(i1,j1:j9),1,2); d(i9+1:i9+2,j1:j9) = spread(d(i9,j1:j9),1,2)
du(i1-2:i1-1,j1:j9) = spread(du(i1,j1:j9),1,2); du(i9+1:i9+2,j1:j9) = spread(du(i9,j1:j9),1,2)
dv(i1-2:i1-1,j1:j9) = spread(dv(i1,j1:j9),1,2); dv(i9+1:i9+2,j1:j9) = spread(dv(i9,j1:j9),1,2)

vb(i1-2:i1-1,j1:j9) = spread(vb(i1,j1:j9),1,2); vb(i9+1:i9+2,j1:j9) = spread(vb(i9,j1:j9),1,2)
zb(i1-2:i1-1,j1:j9) = spread(zb(i1,j1:j9),1,2); zb(i9+1:i9+2,j1:j9) = spread(zb(i9,j1:j9),1,2)
vba(i1-2:i1-1,j1:j9) = spread(vba(i1,j1:j9),1,2); vba(i9+1:i9+2,j1:j9) = spread(vba(i9,j1:j9),1,2)

! N-S
d(i1:i9,j1-2:j1-1) = spread(d(i1:i9,j1),2,2); d(i1:i9,j9+1:j9+2) = spread(d(i1:i9,j9),2,2)
du(i1:i9,j1-2:j1-1) = spread(du(i1:i9,j1),2,2); du(i1:i9,j9+1:j9+2) = spread(du(i1:i9,j9),2,2)
dv(i1:i9,j1-2:j1-1) = -dv(i1:i9,j1+1:j1:-1); dv(i1:i9,j9+1:j9+2) = -dv(i1:i9,j9:j9-1:-1)

vb(i1:i9,j1-2:j1-1) = spread(vb(i1:i9,j1),2,2); vb(i1:i9,j9+1:j9+2) = spread(vb(i1:i9,j9),2,2)
zb(i1:i9,j1-2:j1-1) = spread(zb(i1:i9,j1),2,2); zb(i1:i9,j9+1:j9+2) = spread(zb(i1:i9,j9),2,2)
vba(i1:i9,j1-2:j1-1) = spread(vba(i1:i9,j1),2,2); vba(i1:i9,j9+1:j9+2) = spread(vba(i1:i9,j9),2,2)

! update U,V,Z,H,Pc
u = du/d
v = dv/d
uv = HYPOT(u,v)
uvx = COS(ATAN2(v,u))
uvy = SIN(ATAN2(v,u))

d0 = d; du0 = du; dv0 = dv; dn = d; dun = du; dvn = dv; vbn = vb; zban = zba; vban = vba

do j = j1,j9
	do i = i1,i9
		if (isnan(d(i,j))) stop '"dini" is a NaN'
		if (isnan(uv(i,j))) stop '"uvini" is a NaN'
	enddo
enddo

do j = 1,iy
	do i = 1,ix
		Pc(i,j) = MAX(0.0,MIN(1.0,zb(i,j)/Cm0)); Pe(i,j) = 1.0-Pc(i,j)
	enddo
enddo

! initial critical Shields number
thetacr = thetac0
thetacg = thetac0

! set u and v to normal conditions based on water discharge
uin = qwtot/(davg*BC)
qin = qwtot/BC
qs = scrate*qstot

! print initial conditions
write(6,'("Slope	: ",f9.5)') slope
write(6,'("d_m	: ",f7.3)') davg
write(6,'("u_m/s	: ",f7.3)') uin
write(6,'("Qw_l/s	: ",f7.3)') uin*davg*BC*1e03
write(6,'("Qc_g/s	: ",f7.3)') qstot*2.65e06
write(6,'("Qs_g/s	: ",f7.3)') qs*2.65e06
write(6,'("qs/qc	: ",f7.3)') scrate
write(6,'("ksa_mm	: ",f7.3)') rd*dsm*1e3
write(6,'("ksb_mm	: ",f7.3)') kb*1e3
write(6,'("zbi_m	: ",f7.3)') zbini

END SUBROUTINE SetBed
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE BEDROCK
USE variables
IMPLICIT none  
integer	(kind=2)	:: i,j
real	(kind=rk)	:: minzra,std,avg
! import bedrock topography elevation

! Read initial topography
OPEN(UNIT=unit_var,FILE="BRtopo.txt",FORM="FORMATTED",STATUS="OLD",ACTION="READ")
do i = 1,ix
	do j = 1,iy
		read(unit_var,*) zra(i,j)
	enddo
enddo
close (unit_var)

! calculate bed topography featuers
minzra = ABS(MINVAL(zra(i1:i9,j1:j9)))
zavg  = zra+minzra
avg = SUM(zavg(i1:i9,j1:j9))
zbar = avg/(nx*ny)	! average bedrock height
std = SUM((zavg(i1:i9,j1:j9)-zbar)**2)
zstd = SQRT(std/((nx*ny)-1))	! bedrock standard deviation

! set boundary values
zra(1:feedxs-1,j1:j9) = SPREAD(zra(feedxs,j1:j9),1,feedxs-1); zra(i9-10:ix,j1:j9) = SPREAD(zra(i9-10,j1:j9),1,13)
zra(1:ix,j1-2:j1-1) = SPREAD(zra(1:ix,j1),2,2); zra(1:ix,j9+1:j9+2) = SPREAD(zra(1:ix,j9),2,2)

end subroutine BEDROCK
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE BoundCond(kBC)
USE variables
IMPLICIT none    
integer	(kind=1),	intent(in)	:: kBC
integer	(kind=2)	:: i,j,k
real	(kind=rk)	:: qw,db,ci,ai,Fri
! set boundary values

if (kBC<=2) then
	! W-E 
	! upstream flow boundary condition
	i = feedxs
	if (kcont .GT. INT(100/dt)) then
		do j = j1,j9
			Fri = SQRT((du(i,j)**2.+dv(i,j)**2.)/g/d(i,j)**3.)
			qw = qin
			ai = du(i,j)/d(i,j)-2.*SQRT(g*d(i,j))
			ci = roots3((/0.5*ai,0.0,-0.5*g*qw/),Fri)
			db = ci**2./g
			d(i-1,j) = db; du(i-1,j) = qw
		enddo
	else
		d(i-1,j1:j9) = d(i,j1:j9)
		du(i-1,j1:j9) = qin
	endif

	do j = j1,j9
		do i = feedxs,2,-1
			d(i-1,j) = d(feedxs-1,j)
			du(i-1,j) = du(feedxs-1,j)
			dv(i-1,j) = dv(i,j)+(dv(i,j)-dv(i+1,j))/15.
		enddo
	enddo

	! downstream flow boundary condition
	do j = j1,j9
		do i = i9-10,ix-1
			hd(i+1,j) = hd(i,j)+(hd(i,j)-hd(i-1,j))/15.
			d(i+1,j) = hd(i,j)-zd(i,j)
			du(i+1,j) = du(i,j)
			dv(i+1,j) = dv(i,j)+(dv(i,j)-dv(i-1,j))/15.
		enddo
	enddo
	
	! N-S
	d(1:ix,j1-2:j1-1) = SPREAD(d(1:ix,j1),2,2); d(1:ix,j9+1:j9+2) = SPREAD(d(1:ix,j9),2,2)
	du(1:ix,j1-2:j1-1) = SPREAD(du(1:ix,j1),2,2); du(1:ix,j9+1:j9+2) = SPREAD(du(1:ix,j9),2,2)
	dv(1:ix,j1-2:j1-1) = -dv(1:ix,j1+1:j1:-1); dv(1:ix,j9+1:j9+2) = -dv(1:ix,j9:j9-1:-1)
	
	u = du/d; v = dv/d; uv = HYPOT(u,v); uvx = COS(ATAN2(v,u)); uvy = SIN(ATAN2(v,u))
elseif (kBC==3) then
	! morphological boundary conditions
	! W-E INLET
	do j = j1,j9
		do i = feedxs,2,-1
			vb(i-1,j) = MAX(0.,vb(i,j)+(vb(i,j)-vb(i+1,j))/15.)
			vba(i-1,j) = MAX(0.,vba(i,j)+(vba(i,j)-vba(i+1,j))/15.)
			zba(i-1,j) = MAX(0.,zba(i,j)+(zba(i,j)-zba(i+1,j))/15.)
		enddo
		do i = i9-10,ix-1
			vb(i+1,j) = MAX(0.,vb(i,j)+(vb(i,j)-vb(i-1,j))/15.)
			vba(i+1,j) = MAX(0.,vba(i,j)+(vba(i,j)-vba(i-1,j))/15.)
			zba(i+1,j) = MAX(0.,zba(i,j)+(zba(i,j)-zba(i-1,j))/15.)
		enddo
	enddo
	
	! N-S
	vb(1:ix,j1-2:j1-1) = SPREAD(vb(1:ix,j1),2,2); vb(1:ix,j9+1:j9+2) = SPREAD(vb(1:ix,j9),2,2)
	zba(1:ix,j1-2:j1-1) = SPREAD(zba(1:ix,j1),2,2); zba(1:ix,j9+1:j9+2) = SPREAD(zba(1:ix,j9),2,2)
	vba(1:ix,j1-2:j1-1) = SPREAD(vba(1:ix,j1),2,2); vba(1:ix,j9+1:j9+2) = SPREAD(vba(1:ix,j9),2,2)
	
	! update Zb,Pc,Z
	zb = zba/(1.-pore)
	Pc = MAX(0.0,MIN(1.0,zb/Cm0)); Pe = 1.0-Pc
	z = zbr+zb; zd = zra+zb
	vban = vba; vbn = vb; zban = zba
endif

h = z+d; hd = zd+d

END SUBROUTINE BoundCond
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE GRAD
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j
! gradient
! dd/dx, dd/dy, du/dx, du/dy, dv/dx, dv/dy, dz/dx, dz/dy

dddx(i1:i9,j1:j9) = 0.5*(d(i1+1:i9+1,j1:j9)-d(i1-1:i9-1,j1:j9))/dx
dddy(i1:i9,j1:j9) = 0.5*(d(i1:i9,j1+1:j9+1)-d(i1:i9,j1-1:j9-1))/dy
dudx(i1:i9,j1:j9) = 0.5*(u(i1+1:i9+1,j1:j9)-u(i1-1:i9-1,j1:j9))/dx
dudy(i1:i9,j1:j9) = 0.5*(u(i1:i9,j1+1:j9+1)-u(i1:i9,j1-1:j9-1))/dy
dvdx(i1:i9,j1:j9) = 0.5*(v(i1+1:i9+1,j1:j9)-v(i1-1:i9-1,j1:j9))/dx
dvdy(i1:i9,j1:j9) = 0.5*(v(i1:i9,j1+1:j9+1)-v(i1:i9,j1-1:j9-1))/dy
dzdx(i1:i9,j1:j9) = 0.5*(zd(i1+1:i9+1,j1:j9)-zd(i1-1:i9-1,j1:j9))/dx
dzdy(i1:i9,j1:j9) = 0.5*(zd(i1:i9,j1+1:j9+1)-zd(i1:i9,j1-1:j9-1))/dy

dudx(i1-2:i1,j1:j9) = 0.; dudx(i9:i9+2,j1:j9) = 0.
dudx(i1:i9,j1-2:j1-1) = SPREAD(dudx(i1:i9,j1),2,2); dudx(i1:i9,j9+1:j9+2) = SPREAD(dudx(i1:i9,j9),2,2)
dudy(i1-2:i1,j1:j9) = SPREAD(dudy(i1+1,j1:j9),1,3); dudy(i9:i9+2,j1:j9) = SPREAD(dudy(i9-1,j1:j9),1,3)
dudy(i1:i9,j1-2:j1-1) = 0.; dudy(i1:i9,j9+1:j9+2) = 0.

dvdx(i1-2:i1,j1:j9) = SPREAD(dvdx(i1+1,j1:j9),1,3); dvdx(i9:i9+2,j1:j9) = SPREAD(dvdx(i9-1,j1:j9),1,3)
dvdx(i1:i9,j1-2:j1-1) = -dvdx(i1:i9,j1+1:j1:-1); dvdx(i1:i9,j9+1:j9+2) = -dvdx(i1:i9,j9:j9-1:-1)
dvdy(i1-2:i1,j1:j9) = SPREAD(dvdy(i1+1,j1:j9),1,3); dvdy(i9:i9+2,j1:j9) = SPREAD(dvdy(i9-1,j1:j9),1,3)
dvdy(i1:i9,j1-2:j1-1) = dvdy(i1:i9,j1+1:j1:-1); dvdy(i1:i9,j9+1:j9+2) = dvdy(i1:i9,j9:j9-1:-1)

dzdx(i1-2:i1,j1:j9) = 0.; dzdx(i9:i9+2,j1:j9) = 0.
dzdx(i1:i9,j1-2:j1-1) = SPREAD(dzdx(i1:i9,j1),2,2); dzdx(i1:i9,j9+1:j9+2) = SPREAD(dzdx(i1:i9,j9),2,2)
dzdy(i1-2:i1,j1:j9) = SPREAD(dzdy(i1+1,j1:j9),1,3); dzdy(i9:i9+2,j1:j9) = SPREAD(dzdy(i9-1,j1:j9),1,3)
dzdy(i1:i9,j1-2:j1-1) = 0.; dzdy(i1:i9,j9+1:j9+2) = 0.

do j = j1,j9
	do i = i1,i9
		if (isnan(du(i,j)) .OR. isnan(dv(i,j))) then
			write(6,*) rkn,i,j,d(i,j),du(i,j),dv(i,j),Cf(i,j),Sf(i,j),uv(i,j)
			pause '"du/dx_du/dy" is a NaN'
		endif
		
		if (isnan(dzdx(i,j)) .OR. isnan(dzdy(i,j))) then
			write(6,*) i,j,dzdx(i,j),dzdy(i,j),zd(i-1:i+1,j-1:j+1)
			pause '"dz/dx" is a NaN'
		endif
	enddo
enddo

END SUBROUTINE GRAD
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SRC
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j
! Streamline Radius of Curvature

do j = j1,j9
	do i = i1,i9
		if (uv(i,j)==0.) then
			Cs(i,j) = 0.
		else
			Cs(i,j) = (u(i,j)**2.*dvdx(i,j)-v(i,j)**2.*dudy(i,j)+u(i,j)*v(i,j)*(dvdy(i,j)-dudx(i,j)))/uv(i,j)**3.
		endif
	enddo
enddo

END SUBROUTINE SRC
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE RoughGRAIN
USE variables
IMPLICIT none
! Roughness (1) skin friction

ksg(i1:i9,j1:j9) = Pc(i1:i9,j1:j9)*ksa(i1:i9,j1:j9)+Pe(i1:i9,j1:j9)*ksb(i1:i9,j1:j9)

ksg(1:feedxs-1,j1:j9) = SPREAD(ksg(feedxs,j1:j9),1,feedxs-1); ksg(i9-9:ix,j1:j9) = SPREAD(ksg(i9-10,j1:j9),1,12)
ksg(i1:i9,j1-2:j1-1) = SPREAD(ksg(i1:i9,j1),2,2); ksg(i1:i9,j9+1:j9+2) = SPREAD(ksg(i1:i9,j9),2,2)

END SUBROUTINE RoughGRAIN
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE RoughSEDTR
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j
real	(kind=rk)	:: kst0
! Roughness (2) near-bed sediment transport
! (Wiberg and Rubin, 1989)

do j = j1,j9
	do i = i1,i9
		if (qbc(i,j).GT.0. .AND. theta(i,j).GT.thetacr(i,j)) then
			kst0 = 0.056*dsm*0.68*theta(i,j)/thetacr(i,j)/(1.+dst*theta(i,j)/thetacr(i,j))
			kst(i,j) = 30.*kst0*Pc(i,j)
		else
			kst(i,j) = 0.0
		endif
		
		if (ISNAN(kst(i,j))) then
			write(6,*) i,j,d(i,j),kst(i,j),theta(i,j),thetacr(i,j),thetacg(i,j)
			pause '"kst" is a NaN'
		endif
	enddo
enddo

kst(1:feedxs-1,j1:j9) = SPREAD(kst(feedxs,j1:j9),1,feedxs-1); kst(i9-9:ix,j1:j9) = SPREAD(kst(i9-10,j1:j9),1,12)
kst(i1:i9,j1-2:j1-1) = SPREAD(kst(i1:i9,j1),2,2); kst(i1:i9,j9+1:j9+2) = SPREAD(kst(i1:i9,j9),2,2)

END SUBROUTINE RoughSEDTR
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE RoughBEDFR
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j,n1,n2,n3,n8(ix),n9(ix)
real	(kind=4)	:: m1,m2,ksf0(ix),ksf1(ix)
! Roughness (3) form drag

n1 = 0; n2 = 0; n3 = 0; n8 = 0; n9 = 0
M01 = 0.; ksf1 = 0.

M01(i1:i9,2) = zd(i1:i9,j1)-zd(i1:i9,j9)
M01(i1-2:i1,2) = M01(i1,2); M01(i9:i9+2,2) = M01(i9,2)

do i = i1,i9
	if (i.LT.nw+i1) then
		n1 = (i+nw)-i1+1
		M01(i,1) = SUM(M01(i1:i+nw,2))/n1
	elseif (i.GT.i9-nw) then
		n1 = i9-(i-nw)+1
		M01(i,1) = SUM(M01(i-nw:i9,2))/n1
	else
		n1 = 2*nw+1
		M01(i,1) = SUM(M01(i-nw:i+nw,2))/n1
	endif
enddo
M01(i1-2:i1,1) = SPREAD(M01(i1,1),1,3); M01(i9:i9+2,1) = SPREAD(M01(i9,1),1,3)

n1 = 1
do i = feedxs,i9-11
	if (M01(i,1).GT.0. .AND. M01(i+1,1).LT.0.) then
		n1 = n1+1; n8(n1) = i
	endif
	if (M01(i,1).LT.0. .AND. M01(i+1,1).GT.0.) then
		n1 = n1+1; n8(n1) = i
	endif
enddo
n8(1) = feedxs-1; n8(n1+1) = i9-11

do i = 1,n1
	n9(i) = MAXLOC(ABS(M01(n8(i)+1:n8(i+1),1)),DIM=1,BACK=.TRUE.)+n8(i)
enddo

do i = 1,n1
	m1 = 0.; n2 = 0
	if (i==1) then
		m1 = ABS(M01(n9(i),1)-M01(n9(i+1),1))
		n2 = ABS(n9(i)-n9(i+1))
		ksf0(i1:n8(i+1)) = 0.25*m1**2./(n2*dx)
	elseif (i==n1) then
		m1 = ABS(M01(n9(i),1)-M01(n9(i-1),1))
		n2 = ABS(n9(i)-n9(i-1))
		ksf0(n8(i)+1:i9) = 0.25*m1**2./(n2*dx)
	else
		m1 = MAXVAL(ABS(M01(n8(i)+1:n8(i+1),1)))
		n2 = ABS(n8(i)-n8(i+1))
		ksf0(n8(i)+1:n8(i+1)) = m1**2./(n2*dx)
	endif
	
	n3 = n3+n2
	if (ISNAN(m1) .OR. n2.LE.0.) then
		write(6,*) i,m1,n2,n8(i:i+1),M01(n8(i),1),M01(n8(i+1),1),ksf0(n8(i)),ksf0(n8(i+1))
		pause '"ksf0" is a NaN'
	endif
enddo

nw = n3/n1/4.
! ar = 0.923 (Grant and Madsen), 0.533 (Raudkivi), 0.267 (Nielsen), 0.3~3 (Soulsby)
do i = i1,i9
	ksf1(i) = 30.*ksf0(i)*0.923
enddo

ksf1(1:feedxs) = ksf1(feedxs); ksf1(i9-10:ix) = ksf1(i9-10)
do i = i1,i9
	if (i.LT.nw+i1) then
		n1 = (i+nw)-i1+1
		M01(i,j1) = SUM(ksf1(i1:i+nw))/n1
	elseif (i.GT.i9-nw) then
		n1 = i9-(i-nw)+1
		M01(i,j1) = SUM(ksf1(i-nw:i9))/n1
	else
		n1 = 2*nw+1
		M01(i,j1) = SUM(ksf1(i-nw:i+nw))/n1
	endif
	
	if (ISNAN(ksf1(i)) .OR. ksf1(i).LT.0.) then
		write(6,*) i,j,d(i,j1),M01(i,j1),zd(i,j1),z(i,j1),ksf1(i),nw,n1
		pause '"ksf1" is a NaN'
	endif
enddo

ksf(i1:i9,j1:j9) = SPREAD(M01(i1:i9,j1),2,ny)
ksf(1:feedxs-1,j1:j9) = SPREAD(ksf(feedxs,j1:j9),1,feedxs-1); ksf(i9-9:ix,j1:j9) = SPREAD(ksf(i9-10,j1:j9),1,12)
ksf(i1:i9,j1-2:j1-1) = SPREAD(ksf(i1:i9,j1),2,2); ksf(i1:i9,j9+1:j9+2) = SPREAD(ksf(i1:i9,j9),2,2)

END SUBROUTINE RoughBEDFR
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE FRICTION
USE variables
IMPLICIT none
! uses log law to calculate roughness, with local roughness proportional to sed concentration
! ks = hydraulic roughness height
! cf = friction roughness
integer	(kind=2)	:: i,j,n1
real	(kind=rk)	:: b1,b2,a3,Rk0,Rkg,kg,k1,k2

call RoughBEDFR
call RoughSEDTR
call RoughGRAIN
call GRAD

do j = j1,j9
	do i = i1,i9
		! total roughness height
		ks(i,j) = MIN(0.5*d(i,j),ksg(i,j)+ksf(i,j)+kst(i,j))
		
		! logarithmic law
		Rk0 = MAX(0.1001,d(i,j)/ks(i,j))
		a3 = (1.0-0.1/Rk0)*Ca(1)
		cf(i,j) = (a3*LOG(Ca(2)*Rk0))**(-2.)
		
		if (rkn==3) then
			! relative roughness to grain size
			k2 = ksg(i,j)+kst(i,j)
			ksn(i,j) = rc(1)*k2+rc(2)*ksf(i,j)
			ksc(i,j) = ksn(i,j)/dsm
			
			! fritction partitioning
			b2 = 1./LOG(Ca(2)*Rk0)
			b1 = b2/(b2+0.5)
			k1 = ksg(i,j)+kst(i,j)
			kg = k1**(1.-b1)*ks(i,j)**b1
			Rkg = MAX(0.1001,d(i,j)/kg)
			cfg(i,j) = (a3*LOG(Ca(2)*Rkg))**(-2.)
		endif
		
		if (isnan(cf(i,j)) .OR. isnan(cfg(i,j))) then
			write(6,*) time,i,j,d(i,j),cf(i,j),cfg(i,j),ksg(i,j),kst(i,j),ksf(i,j),kg,ksc(i,j)
			pause '"cf" is a NaN'
		endif
	enddo
enddo

ks(i1-2:i1,j1:j9) = SPREAD(ks(i1+1,j1:j9),1,3); ks(i9:i9+2,j1:j9) = SPREAD(ks(i9-1,j1:j9),1,3)
ks(i1:i9,j1-2:j1-1) = SPREAD(ks(i1:i9,j1),2,2); ks(i1:i9,j9+1:j9+2) = SPREAD(ks(i1:i9,j9),2,2)
cf(i1-2:i1,j1:j9) = SPREAD(cf(i1+1,j1:j9),1,3); cf(i9:i9+2,j1:j9) = SPREAD(cf(i9-1,j1:j9),1,3)
cf(i1:i9,j1-2:j1-1) = SPREAD(cf(i1:i9,j1),2,2); cf(i1:i9,j9+1:j9+2) = SPREAD(cf(i1:i9,j9),2,2)

ksg(i1-2:i1,j1:j9) = SPREAD(ksg(i1+1,j1:j9),1,3); ksg(i9:i9+2,j1:j9) = SPREAD(ksg(i9-1,j1:j9),1,3)
ksg(i1:i9,j1-2:j1-1) = SPREAD(ksg(i1:i9,j1),2,2); ksg(i1:i9,j9+1:j9+2) = SPREAD(ksg(i1:i9,j9),2,2)
cfg(i1-2:i1,j1:j9) = SPREAD(cfg(i1+1,j1:j9),1,3); cfg(i9:i9+2,j1:j9) = SPREAD(cfg(i9-1,j1:j9),1,3)
cfg(i1:i9,j1-2:j1-1) = SPREAD(cfg(i1:i9,j1),2,2); cfg(i1:i9,j9+1:j9+2) = SPREAD(cfg(i1:i9,j9),2,2)

END SUBROUTINE FRICTION
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE DIFFUSION
USE variables
IMPLICIT none    
integer	(kind=2)	:: i,j
real	(kind=rk)	:: vim,ss,lmh
real	(kind=rk)	:: dndx,dndy,ddudx,ddudy,ddvdx,ddvdy,d2dudx2,d2dvdx2
real	(kind=rk)	:: d2dudxdy,d2dvdxdy,d2dudy2,d2dvdy2
! diffusive term calculation

do j = j1,j9
	do i = i1,i9
		vim = uv(i,j)*SQRT(cf(i,j))*d(i,j)*vonk/6.
		ss = 2.*dudx(i,j)**2.+2.*dvdy(i,j)**2.+(dvdx(i,j)+dudy(i,j))**2.
		lmh = 4./15.*vonk*d(i,j)
		vis(i,j) = d(i,j)*(SQRT(lmh**4.*ss+vim**2.)+visk)
	enddo
enddo

vis(i1-2:i1-1,j1:j9) = SPREAD(vis(i1,j1:j9),1,2); vis(i9+1:i9+2,j1:j9) = SPREAD(vis(i9,j1:j9),1,2)
vis(i1:i9,j1-2:j1-1) = SPREAD(vis(i1:i9,j1),2,2); vis(i1:i9,j9+1:j9+2) = SPREAD(vis(i1:i9,j9),2,2)

do j = j1,j9
	do i = i1,i9
		if (d(i,j).GT.sma) then
			dndx = 0.5*(vis(i+1,j)-vis(i-1,j))/dx
			dndy = 0.5*(vis(i,j+1)-vis(i,j-1))/dy
			ddudx = 0.5*(du(i+1,j)-du(i-1,j))/dx
			ddudy = 0.5*(du(i,j+1)-du(i,j-1))/dy
			ddvdx = 0.5*(dv(i+1,j)-dv(i-1,j))/dx
			ddvdy = 0.5*(dv(i,j+1)-dv(i,j-1))/dy
			d2dudx2 = (du(i-1,j)-2.*du(i,j)+du(i+1,j))/dx**2.
			d2dudy2 = (du(i,j-1)-2.*du(i,j)+du(i,j+1))/dx**2.
			d2dvdx2 = (dv(i-1,j)-2.*dv(i,j)+dv(i+1,j))/dx**2.
			d2dvdy2 = (dv(i,j-1)-2.*dv(i,j)+dv(i,j+1))/dx**2.
			d2dudxdy = (du(i+1,j+1)-du(i+1,j-1)-du(i-1,j+1)+du(i-1,j-1))/(4.*dx*dy)
			d2dvdxdy = (dv(i+1,j+1)-dv(i+1,j-1)-dv(i-1,j+1)+dv(i-1,j-1))/(4.*dx*dy)
			
			Dfx(i,j) = 2.*(dndx*ddudx+vis(i,j)*d2dudx2)+dndy*(ddudy+ddvdx)+vis(i,j)*(d2dudy2+d2dvdxdy)
			Dfy(i,j) = dndx*(ddudy+ddvdx)+vis(i,j)*(d2dvdx2+d2dudxdy)+2.*(dndy*ddvdy+vis(i,j)*d2dvdy2)
		else
			Dfx(i,j) = 0.
			Dfy(i,j) = 0.
		endif
		
		if (isnan(Dfx(i,j)) .OR. isnan(Dfy(i,j))) then
			write(6,*) i,j,Dfx(i,j),Dfy(i,j),du(i-1:i+1,j-1:j+1),dv(i-1:i+1,j-1:j+1)
			pause '"Dfx & Dfy" is a NaN'
		endif
    enddo
enddo

END SUBROUTINE DIFFUSION
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SOURCE
USE variables
IMPLICIT none    
integer	(kind=2)	:: i,j
! source term bed elevation gradient calculation

do j = j1,j9
	do i = i1,i9
		if (d(i,j).GT.sma) then
			Sbx(i,j) = -g*0.5*(hL(i,j,1)+hR(i-1,j,1))*(zE(i,j,1)-zE(i-1,j,1))/dx
			Sby(i,j) = -g*0.5*(hL(i,j,2)+hR(i,j-1,2))*(zE(i,j,2)-zE(i,j-1,2))/dy
		else
			Sbx(i,j) = 0.
			Sby(i,j) = 0.
		endif
		
		if (isnan(Sbx(i,j)) .OR. isnan(Sby(i,j))) then
			write(6,*) i,j,Sbx(i,j),Sby(i,j),hL(i,j,1),hR(i-1,j,1),zE(i,j,1),zE(i-1,j,1)
			pause '"Sbx & Sby" is a NaN'
		endif
	enddo
enddo

end subroutine SOURCE
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE TVD(kW)
USE variables
IMPLICIT none
integer	(kind=1),	intent(in)	:: kW
integer	(kind=2)	:: i,j
real	(kind=rk)	:: r0
! TVD (Total Variation Diminishing)

fp = 0.; fm = 0.
if (kW==1) then
	do j = j1,j9
		do i = i1,i9-1
			r0 = (vc(i+1,j)-vc(i,j))/(vc(i,j)-vc(i-1,j))
			fm(i,j) = vc(i,j)+0.5*minmod1(r0)*(vc(i,j)-vc(i-1,j))
			r0 = (vc(i+2,j)-vc(i+1,j))/(vc(i+1,j)-vc(i,j))
			fp(i,j) = vc(i+1,j)-0.5*minmod1(r0)*(vc(i+1,j)-vc(i,j))
		enddo
	enddo
elseif (kW==2) then
	do j = 3,iy-3
		do 	i = i1,i9
			r0 = (vc(i,j+1)-vc(i,j))/(vc(i,j)-vc(i,j-1))
			fm(i,j) = vc(i,j)+0.5*minmod1(r0)*(vc(i,j)-vc(i,j-1))
			r0 = (vc(i,j+2)-vc(i,j+1))/(vc(i,j+1)-vc(i,j))
			fp(i,j) = vc(i,j+1)-0.5*minmod1(r0)*(vc(i,j+1)-vc(i,j))
		enddo
	enddo
endif

END SUBROUTINE TVD
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SPLIT
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j,k,n
! flux splitting

do n = 1,4
	vc = 0.
	if (n==1) then
		vc = h
	elseif (n==2) then
		vc = d
	elseif (n==3) then
		vc = du
	elseif (n==4) then
		vc = dv
	endif

	do k = 1,2
		call TVD(k)
		
		if (k==1) then
			fm(1:3,j1:j9) = vc(1:3,j1:j9)
			fp(1:3,j1:j9) = vc(2:4,j1:j9)
			fm(ix-3:ix-1,j1:j9) = vc(ix-3:ix-1,j1:j9)
			fp(ix-3:ix-1,j1:j9) = vc(ix-2:ix,j1:j9)
		elseif (k==2) then
			fm(i1:i9,1:2) = vc(i1:i9,1:2)
			fp(i1:i9,1:2) = vc(i1:i9,2:3)
			fm(i1:i9,iy-2:iy-1) = vc(i1:i9,iy-2:iy-1)
			fp(i1:i9,iy-2:iy-1) = vc(i1:i9,iy-1:iy)
		endif
		
		if (n==1) then
			hL(:,:,k) = fm; hR(:,:,k) = fp
		elseif (n==2) then
			dL(:,:,k) = fm; dR(:,:,k) = fp
		elseif (n==3) then
			duL(:,:,k) = fm; duR(:,:,k) = fp
		elseif (n==4) then
			dvL(:,:,k) = fm; dvR(:,:,k) = fp
		endif
	enddo
enddo

dvL(i1:i9,1:2,2) = 0.; dvL(i1:i9,j9:j9+1,2) = 0.

uL = duL/dL; uR = duR/dR
vL = dvL/dL; vR = dvR/dR
zL = hL-dL; zR = hR-dR

END SUBROUTINE SPLIT
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE RECON
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j
! flux reconstruction

do j = j1,j9
	do i = 1,ix-1
		zE(i,j,1) = MAX(zL(i,j,1),zR(i,j,1))
	enddo
enddo

do j = 2,iy-2
	do i = i1,i9
		zE(i,j,2) = MAX(zL(i,j,2),zR(i,j,2))
	enddo
enddo
zE(i1:i9,2,2) = z(i1:i9,3); zE(i1:i9,j9,2) = z(i1:i9,j9)

! wetting/drying
do j = j1,j9
	do i = 3,ix-3
		if (dL(i,j,1).LE.sma .AND. dR(i,j,1).GT.sma) then
			if (zL(i,j,1) .GE. hR(i,j,1)) hL(i,j,1) = hR(i,j,1)
		endif
		if (dL(i,j,1).GT.sma .AND. dR(i,j,1).LE.sma) then
			if (zR(i,j,1) .GE. hL(i,j,1)) hR(i,j,1) = hL(i,j,1)
		endif
	enddo
enddo
do j = 3,iy-3
	do i = i1,i9
		if (dL(i,j,2).LE.sma .AND. dR(i,j,2).GT.sma) then
			if (zL(i,j,2) .GE. hR(i,j,2)) hL(i,j,2) = hR(i,j,2)
		endif
		if (dL(i,j,2).GT.sma .AND. dR(i,j,2).LE.sma) then
			if (zR(i,j,2) .GE. hL(i,j,2)) hR(i,j,2) = hL(i,j,2)
		endif
	enddo
enddo

do j = j1,j9
	do i = 1,ix-1
		dL(i,j,1) = MAX(0.,hL(i,j,1)-zE(i,j,1))
		dR(i,j,1) = MAX(0.,hR(i,j,1)-zE(i,j,1))
	enddo
enddo

do j = 2,iy-2
	do i = i1,i9
		dL(i,j,2) = MAX(0.,hL(i,j,2)-zE(i,j,2))
		dR(i,j,2) = MAX(0.,hR(i,j,2)-zE(i,j,2))
	enddo
enddo

hL = dL+zE; hR = dR+zE
duL = dL*uL; duR = dR*uR
dvL = dL*vL; dvR = dR*vR

! Roe average
! Song et al. (2011); Liang (2011)
um = (uR*SQRT(dR)+uL*SQRT(dL))/(SQRT(dR)+SQRT(dL))
vm = (vR*SQRT(dR)+vL*SQRT(dL))/(SQRT(dR)+SQRT(dL))
dm = 0.5*(dR+dL)
! Guan et al. (2013)
! um = (uR*SQRT(dR)+uL*SQRT(dL))/(SQRT(dR)+SQRT(dL))
! vm = (vR*SQRT(dR)+vL*SQRT(dL))/(SQRT(dR)+SQRT(dL))
! dm = SQRT(dR*dL)
! Liang and Borthwick (2011); Loukili and Soulaimani (2007); Kim and Lee (2012)
! um = 0.5*(uL+uR)+SQRT(g*dL)-SQRT(g*dR)
! vm = 0.5*(vL+vR)+SQRT(g*dL)-SQRT(g*dR)
! dm = 1.0/g*(0.5*(SQRT(g*dL)+SQRT(g*dR))+0.25*(uL-uR))**2.

END SUBROUTINE RECON
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE HLLC_CONV(kC)
USE variables
IMPLICIT none
integer	(kind=1),	intent(in)	:: kC
integer	(kind=2)	:: i,j,k,n
real	(kind=rk)	:: P0(3,4),SI(3)
! HLLC conservative variables
! UL U*L U*R UR

k = kC
U1 = 0.0; U2 = 0.0; U3 = 0.0; U4 = 0.0; SW = 0.0

! F component
if (k==1) then
	do j = j1,j9
		do i = 2,ix-2
			P0 = 0.; SI = 0.
			P0(2,1) = 0.5*g*(dL(i,j,k)**2.-zE(i,j,k)**2.)
			P0(2,4) = 0.5*g*(dR(i,j,k)**2.-zE(i,j,k)**2.)
			
			! Estimate wave speeds: SI = (/SL, SM, SR/)
			SI = SLMR((/dL(i,j,k),dm(i,j,k),dR(i,j,k)/),(/uL(i,j,k),um(i,j,k),uR(i,j,k)/),P0(2,4)-P0(2,1))
			! SI = SLMR((/dL(i,j,k),dm(i,j,k),dR(i,j,k)/),(/uL(i,j,k),um(i,j,k),uR(i,j,k)/),0.0)
			SW(i,j,1:3) = SI
			
			U1(i,j,1:3) = (/dL(i,j,k),duL(i,j,k),dvL(i,j,k)/)
			U4(i,j,1:3) = (/dR(i,j,k),duR(i,j,k),dvR(i,j,k)/)
			U2(i,j,1:3) = dL(i,j,k)*(SI(1)-uL(i,j,k))/(SI(1)-SI(2))*(/1.0,SI(2),vL(i,j,k)/)
			U3(i,j,1:3) = dR(i,j,k)*(SI(3)-uR(i,j,k))/(SI(3)-SI(2))*(/1.0,SI(2),vR(i,j,k)/)
		enddo
	enddo

! G component
elseif (k==2) then
	do j = 3,iy-3
		do i = i1,i9
			P0 = 0.0; SI = 0.
			P0(3,1) = 0.5*g*(dL(i,j,k)**2.-zE(i,j,k)**2.)
			P0(3,4) = 0.5*g*(dR(i,j,k)**2.-zE(i,j,k)**2.)
			
			! Estimate wave speeds: SI = (/SL, SM, SR/)
			SI = SLMR((/dL(i,j,k),dm(i,j,k),dR(i,j,k)/),(/vL(i,j,k),vm(i,j,k),vR(i,j,k)/),P0(3,4)-P0(3,1))
			! SI = SLMR((/dL(i,j,k),dm(i,j,k),dR(i,j,k)/),(/vL(i,j,k),vm(i,j,k),vR(i,j,k)/),0.0)
			SW(i,j,1:3) = SI
			
			U1(i,j,1:3) = (/dL(i,j,k),duL(i,j,k),dvL(i,j,k)/)
			U4(i,j,1:3) = (/dR(i,j,k),duR(i,j,k),dvR(i,j,k)/)
			U2(i,j,1:3) = dL(i,j,k)*(SI(1)-vL(i,j,k))/(SI(1)-SI(2))*(/1.0,uL(i,j,k),SI(2)/)
			U3(i,j,1:3) = dR(i,j,k)*(SI(3)-vR(i,j,k))/(SI(3)-SI(2))*(/1.0,uR(i,j,k),SI(2)/)
		enddo
	enddo
	
	do i = i1,i9
		U1(i,2,1:3) = (/dL(i,2,k),duL(i,2,k),dvL(i,2,k)/)
		U4(i,2,1:3) = (/dR(i,2,k),duR(i,2,k),dvR(i,2,k)/)
		U1(i,j9,1:3) = (/dL(i,j9,k),duL(i,j9,k),dvL(i,j9,k)/)
		U4(i,j9,1:3) = (/dR(i,j9,k),duR(i,j9,k),dvR(i,j9,k)/)
	enddo
	
	U2(i1:i9,2.,1:3) = U1(i1:i9,2.,1:3)
	U3(i1:i9,2.,1:3) = U4(i1:i9,2.,1:3)
	U2(i1:i9,j9,1:3) = U1(i1:i9,j9,1:3)
	U3(i1:i9,j9,1:3) = U4(i1:i9,j9,1:3)
endif

END SUBROUTINE HLLC_CONV
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE HLLC_FLUX(kF)
USE variables
IMPLICIT none
integer	(kind=1),	intent(in)	:: kF
integer	(kind=2)	:: i,j,k,n
real	(kind=rk)	:: U0(3,4),M0(3,4),P0(3,4),M1(3,3),SI(3),ck(3),Ak(3)
! HLLC flux with WAF

k = kF
! F component
if (k==1) then
	do j = j1,j9
		do i = 2,ix-2
			U0 = 0.; M0 = 0.; P0 = 0.; SI = SW(i,j,:)
			
			U0(1:3,1) = U1(i,j,1:3)
			U0(1:3,2) = U2(i,j,1:3)
			U0(1:3,3) = U3(i,j,1:3)
			U0(1:3,4) = U4(i,j,1:3)
			P0(2,1:4) = 0.5*g*(U0(1,1:4)**2.-zE(i,j,k)**2.)
			
			M0(1:3,1) = (/duL(i,j,k),duL(i,j,k)**2./dL(i,j,k),duL(i,j,k)*dvL(i,j,k)/dL(i,j,k)/)
			M0(1:3,4) = (/duR(i,j,k),duR(i,j,k)**2./dR(i,j,k),duR(i,j,k)*dvR(i,j,k)/dR(i,j,k)/)
			M0(1:3,2) = M0(1:3,1)+SI(1)*(U0(1:3,2)-U0(1:3,1))
			M0(1:3,3) = M0(1:3,4)-SI(3)*(U0(1:3,4)-U0(1:3,3))
			M0 = M0+P0
			
			ck = SI*dt/dx
			Ak(1) = WAF(ck(1),U1(i-1:i+1,j,1),U2(i-1:i+1,j,1))
			Ak(2) = WAF(ck(2),U2(i-1:i+1,j,3)/U2(i-1:i+1,j,1),U3(i-1:i+1,j,3)/U3(i-1:i+1,j,1))
			Ak(3) = WAF(ck(3),U3(i-1:i+1,j,1),U4(i-1:i+1,j,1))
			do n = 1,3
				M1(1:3,n) = sgn(ck(n))*Ak(n)*(M0(1:3,n+1)-M0(1:3,n))
			enddo
			
			F1(i,j) = 0.5*(M0(1,1)+M0(1,4))-0.5*sum(M1(1,1:3))
			F2(i,j) = 0.5*(M0(2,1)+M0(2,4))-0.5*sum(M1(2,1:3))
			F3(i,j) = 0.5*(M0(3,1)+M0(3,4))-0.5*sum(M1(3,1:3))
		enddo
	enddo

! G component
elseif (k==2) then
	do j = 3,iy-3
		do i = i1,i9
			U0 = 0.; M0 = 0.; P0 = 0.; SI = SW(i,j,:)
			
			U0(1:3,1) = U1(i,j,1:3)
			U0(1:3,2) = U2(i,j,1:3)
			U0(1:3,3) = U3(i,j,1:3)
			U0(1:3,4) = U4(i,j,1:3)
			P0(3,1:4) = 0.5*g*(U0(1,1:4)**2.-zE(i,j,k)**2.)
			
			M0(1:3,1) = (/dvL(i,j,k),dvL(i,j,k)*duL(i,j,k)/dL(i,j,k),dvL(i,j,k)**2./dL(i,j,k)/)
			M0(1:3,4) = (/dvR(i,j,k),dvR(i,j,k)*duR(i,j,k)/dR(i,j,k),dvR(i,j,k)**2./dR(i,j,k)/)
			M0(1:3,2) = M0(1:3,1)+SI(1)*(U0(1:3,2)-U0(1:3,1))
			M0(1:3,3) = M0(1:3,4)-SI(3)*(U0(1:3,4)-U0(1:3,3))
			M0 = M0+P0
			
			ck = SI*dt/dy
			Ak(1) = WAF(ck(1),U1(i,j-1:j+1,1),U2(i,j-1:j+1,1))
			Ak(2) = WAF(ck(2),U2(i,j-1:j+1,2)/U2(i,j-1:j+1,1),U3(i,j-1:j+1,2)/U3(i,j-1:j+1,1))
			Ak(3) = WAF(ck(3),U3(i,j-1:j+1,1),U4(i,j-1:j+1,1))
			do n = 1,3
				M1(1:3,n) = sgn(ck(n))*Ak(n)*(M0(1:3,n+1)-M0(1:3,n))
			enddo
			
			G1(i,j) = 0.5*(M0(1,1)+M0(1,4))-0.5*sum(M1(1,1:3))
			G2(i,j) = 0.5*(M0(2,1)+M0(2,4))-0.5*sum(M1(2,1:3))
			G3(i,j) = 0.5*(M0(3,1)+M0(3,4))-0.5*sum(M1(3,1:3))
		enddo
	enddo
endif

END SUBROUTINE HLLC_FLUX
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE ADVECTION
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j
real	(kind=rk)	:: pw(3)
! advection term calculation

call HLLC_CONV(1)
call HLLC_FLUX(1)

call HLLC_CONV(2)
call HLLC_FLUX(2)

! W-E
i = feedxs-1
F1(i,j1:j9) = du(i,j1:j9)
F2(i,j1:j9) = u(i,j1:j9)*du(i,j1:j9)+0.5*g*(d(i,j1:j9)**2.-zE(i,j1:j9,1)**2.)
F3(i,j1:j9) = u(i,j1:j9)*dv(i,j1:j9)

i = i9-9
F1(i-1,j1:j9) = du(i,j1:j9)
F2(i-1,j1:j9) = u(i,j1:j9)*du(i,j1:j9)+0.5*g*(d(i,j1:j9)**2.-zE(i-1,j1:j9,1)**2.)
F3(i-1,j1:j9) = u(i,j1:j9)*dv(i,j1:j9)

! N-S
G1(i1:i9,j0) = 0.; G1(i1:i9,j9) = 0.
G2(i1:i9,j0) = 0.; G2(i1:i9,j9) = 0.
G3(i1:i9,j0) = 0.; G3(i1:i9,j9) = 0.

do i = i1,i9
	pw(1) = (v(i,j1)-2.*SQRT(g*d(i,j1)))**2./(4.*g)
	pw(2) = (v(i,j9)+2.*SQRT(g*d(i,j9)))**2./(4.*g)
	G3(i,j0) = 0.5*g*(pw(1)**2.-z(i,j1)**2.)
	G3(i,j9) = 0.5*g*(pw(2)**2.-z(i,j9)**2.)
enddo

END SUBROUTINE ADVECTION
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE HYDRO
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j
real	(kind=4)	:: Rk0,a3
! Hydrodynamics time evolution

call FRICTION
call SRC
call SPLIT
call RECON
call DIFFUSION
call SOURCE
call ADVECTION

FG1(i1:i9,j1:j9) = (F1(i1:i9,j1:j9)-F1(i1-1:i9-1,j1:j9))/dx+(G1(i1:i9,j1:j9)-G1(i1:i9,j1-1:j9-1))/dy
FG2(i1:i9,j1:j9) = (F2(i1:i9,j1:j9)-F2(i1-1:i9-1,j1:j9))/dx+(G2(i1:i9,j1:j9)-G2(i1:i9,j1-1:j9-1))/dy
FG3(i1:i9,j1:j9) = (F3(i1:i9,j1:j9)-F3(i1-1:i9-1,j1:j9))/dx+(G3(i1:i9,j1:j9)-G3(i1:i9,j1-1:j9-1))/dy

Sx(i1:i9,j1:j9) = Sbx(i1:i9,j1:j9)+Dfx(i1:i9,j1:j9)
Sy(i1:i9,j1:j9) = Sby(i1:i9,j1:j9)+Dfy(i1:i9,j1:j9)

d1(i1:i9,j1:j9) = d(i1:i9,j1:j9)-dt*FG1(i1:i9,j1:j9)
du1(i1:i9,j1:j9) = du(i1:i9,j1:j9)-dt*(FG2(i1:i9,j1:j9)-Sx(i1:i9,j1:j9))
dv1(i1:i9,j1:j9) = dv(i1:i9,j1:j9)-dt*(FG3(i1:i9,j1:j9)-Sy(i1:i9,j1:j9))

! friction term with new flow depth
do j = j1,j9
	do i = i1,i9
		Rk0 = MAX(0.1001,d1(i,j)/MIN(0.5*d1(i,j),ks(i,j)))
		a3 = (1.0-0.1/Rk0)*Ca(1)
		cf(i,j) = (a3*LOG(Ca(2)*Rk0))**(-2.)
	enddo
enddo
Sf(i1:i9,j1:j9) = 1.+dt*cf(i1:i9,j1:j9)/d1(i1:i9,j1:j9)**2.*HYPOT(du1(i1:i9,j1:j9),dv1(i1:i9,j1:j9))

! 2nd-order Runge-Kutta
if (rkn==1) then
	d0(i1:i9,j1:j9) = d(i1:i9,j1:j9)
	du0(i1:i9,j1:j9) = du(i1:i9,j1:j9)
	dv0(i1:i9,j1:j9) = dv(i1:i9,j1:j9)
	d(i1:i9,j1:j9) = d1(i1:i9,j1:j9)
	du(i1:i9,j1:j9) = du1(i1:i9,j1:j9)/Sf(i1:i9,j1:j9)
	dv(i1:i9,j1:j9) = dv1(i1:i9,j1:j9)/Sf(i1:i9,j1:j9)
elseif (rkn==2) then
	d(i1:i9,j1:j9) = 0.5*(d0(i1:i9,j1:j9)+d1(i1:i9,j1:j9))
	du(i1:i9,j1:j9) = 0.5*(du0(i1:i9,j1:j9)+du1(i1:i9,j1:j9)/Sf(i1:i9,j1:j9))
	dv(i1:i9,j1:j9) = 0.5*(dv0(i1:i9,j1:j9)+dv1(i1:i9,j1:j9)/Sf(i1:i9,j1:j9))
endif

END SUBROUTINE HYDRO
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SEDSPLIT
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j,k,n
! sediment flux splitting with TVD

do n = 1,2
	vc = 0.
	if (n==1) then
		vc = h
	elseif (n==2) then
		vc = d
	endif

	do k = 1,2
		call TVD(k)
		if (k==1) then
			fm(1:3,j1:j9) = vc(1:3,j1:j9)
			fp(1:3,j1:j9) = vc(2:4,j1:j9)
			fm(ix-3:ix-1,j1:j9) = vc(ix-3:ix-1,j1:j9)
			fp(ix-3:ix-1,j1:j9) = vc(ix-2:ix,j1:j9)
		elseif (k==2) then
			fm(i1:i9,1:2) = vc(i1:i9,1:2)
			fp(i1:i9,1:2) = vc(i1:i9,2:3)
			fm(i1:i9,iy-2:iy-1) = vc(i1:i9,iy-2:iy-1)
			fp(i1:i9,iy-2:iy-1) = vc(i1:i9,iy-1:iy)
		endif
		
		if (n==1) then
			hL(:,:,k) = fm; hR(:,:,k) = fp
		elseif (n==2) then
			dL(:,:,k) = fm; dR(:,:,k) = fp
		endif
	enddo
enddo

zL = hL-dL; zR = hR-dR

do j = j1,j9
	do i = 1,ix-1
		zE(i,j,1) = MAX(zL(i,j,1),zR(i,j,1))
	enddo
enddo

do j = 2,iy-2
	do i = i1,i9
		zE(i,j,2) = MAX(zL(i,j,2),zR(i,j,2))
	enddo
enddo

zE(i1:i9,2,2) = z(i1:i9,3); zE(i1:i9,j9,2) = z(i1:i9,j9)

! wetting/drying
do j = j1,j9
	do i = 3,ix-3
		if (dL(i,j,1).LE.sma .AND. dR(i,j,1).GT.sma) then
			if (zL(i,j,1) .GE. hR(i,j,1)) hL(i,j,1) = hR(i,j,1)
		endif
		if (dL(i,j,1).GT.sma .AND. dR(i,j,1).LE.sma) then
			if (zR(i,j,1) .GE. hL(i,j,1)) hR(i,j,1) = hL(i,j,1)
		endif
	enddo
enddo
do j = 3,iy-3
	do i = i1,i9
		if (dL(i,j,2).LE.sma .AND. dR(i,j,2).GT.sma) then
			if (zL(i,j,2) .GE. hR(i,j,2)) hL(i,j,2) = hR(i,j,2)
		endif
		if (dL(i,j,2).GT.sma .AND. dR(i,j,2).LE.sma) then
			if (zR(i,j,2) .GE. hL(i,j,2)) hR(i,j,2) = hL(i,j,2)
		endif
	enddo
enddo

do j = j1,j9
	do i = 1,ix-1
		dL(i,j,1) = MAX(0.,hL(i,j,1)-zE(i,j,1))
		dR(i,j,1) = MAX(0.,hR(i,j,1)-zE(i,j,1))
	enddo
enddo

do j = 2,iy-2
	do i = i1,i9
		dL(i,j,2) = MAX(0.,hL(i,j,2)-zE(i,j,2))
		dR(i,j,2) = MAX(0.,hR(i,j,2)-zE(i,j,2))
	enddo
enddo

hL = dL+zE; hR = dR+zE

END SUBROUTINE SEDSPLIT
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SEDMOV
USE variables
IMPLICIT none
integer		(kind=2)	:: i,j
real		(kind=rk)	:: m1,m2
! phit = angle between the flow direction and the bed slope direction
! alpha = angle between average particle path and x-axis
! dzds, dzdn = bed slope in sediment transport and transverse direction

call SEDSPLIT

cfs(i1:i9,j1:j9) = cf(i1:i9,j1:j9)
theta(i1:i9,j1:j9) = cf(i1:i9,j1:j9)*uv(i1:i9,j1:j9)**2./tszero

do j = j1,j9
	do i = i1,i9
		m2 = 2./vonk**2.*(1.-SQRT(cfs(i,j))/vonk)
		phit(i,j) = ATAN2(v(i,j),u(i,j))-ATAN(m2*d(i,j)*Cs(i,j))
	enddo
enddo

do j = j1,j9
	do i = i1,i9
		m1 = 9.*(dsm/d(i,j))**0.3*theta(i,j)**0.5
		alphax(i,j) = COS(phit(i,j))-dzdx(i,j)/m1
		alphay(i,j) = SIN(phit(i,j))-dzdy(i,j)/m1
		alpha(i,j) = ATAN2(alphay(i,j),alphax(i,j))
		dzds(i,j) = dzdx(i,j)*alphax(i,j)+dzdy(i,j)*alphay(i,j)
		dzdn(i,j) = dzdy(i,j)*alphax(i,j)-dzdx(i,j)*alphay(i,j)
		
		if (ISNAN(alpha(i,j))) then
			write(6,*) i,j,theta(i,j),cfs(i,j),cf(i,j),Cs(i,j),dzdx(i,j),dzdy(i,j),alpha(i,j)
			pause '"alpha" is a NaN'
		endif
	enddo
enddo

END SUBROUTINE SEDMOV
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SHIELDS
USE variables
IMPLICIT none
integer		(kind=2)	:: i,j
real		(kind=rk)   :: bs,bn,phi,Sc,k1(2)
! critical Shields correction for effect of local bed gradient
! (Soulsby, 1998; Wiberg & Smith, 1987)
! phi = angle of repose ~ 30 degrees
! beta = local bed slope

thetacg(i1:i9,j1:j9) = ksc(i1:i9,j1:j9)**thetac1

do j = j1,j9
	do i = i1,i9
		Sc = 0.4*ksg(i,j)/MAX(1.0,rd)
		phi = ACOS((dsm/Sc-0.02)/(dsm/Sc+1.0))
		bs = ATAN(dzds(i,j)); bn = ATAN(dzdn(i,j))
		
		if (bs.GE.0.) then
			phir(i,j) = MIN(MAX(phi,ABS(bn)),0.5*pi-bs)
		else
			phir(i,j) = MAX(MAX(phi,ABS(bn)),ABS(bs))
		endif
		
		k1(1) = SIN(bs+phir(i,j))/SIN(phir(i,j)); k1(2) = COS(bn)*SQRT(1.-(TAN(bn)/TAN(phir(i,j)))**2.)
		M01(i,j) = MAX(0.001,PRODUCT(k1))
		
		if (ISNAN(M01(i,j)) .OR. M01(i,j).LE.0.) then
			write(6,*) i,j,M01(i,j),Sc,phir(i,j),dzds(i,j),dzdn(i,j),bs,bn,k1
			pause '"theta_cr" is a NaN'
		endif
		
		thetacr(i,j) = thetacg(i,j)*M01(i,j)
	enddo
enddo

END SUBROUTINE SHIELDS
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SEDRPF
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j
real	(kind=rk)	:: qbcs,m1
! calculates ripple factor for near-bed shear stress correction

do j = j1,j9
	do i = i1,i9
		if (d(i,j).LE.ks(i,j)) then
			rpf(i,j) = 1.0
		else
			if (theta(i,j)*thetacg(i,j).GT.0. .AND. theta(i,j).GT.thetacg(i,j)) then
				qbcs = 17.*(SQRT(theta(i,j))-SQRT(thetacg(i,j)))*(theta(i,j)-thetacg(i,j))
				m1 = MAX(0.001,MIN(qbcs,1.0))
				rpf(i,j) = 0.5*(1.8+0.27*LOG10(m1))
			else
				rpf(i,j) = 1.0
			endif
		endif
		
		if (ISNAN(rpf(i,j))) then
			write(6,*) i,j,rpf(i,j),cfg(i,j),cfs(i,j),cf(i,j),theta(i,j),thetacg(i,j),thetacr(i,j),m1,qbcs
			pause '"rpf" is a NaN'
		endif
	enddo
enddo

END SUBROUTINE SEDRPF
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SEDCAP
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j
real	(kind=rk)	:: t0,sn(2)
! calculates sediment fluxes at cell centers
! saltation velocity (Sklar & Dietrich, 2004)
! bedload transport rate (Ashida & Michiue, 1972)

DO j = j1,j9
	DO i = i1,i9
		t0 = (cfg(i,j)/cf(i,j))**rpf(i,j)*theta(i,j)
		theta(i,j) = t0
		
		IF (t0*thetacr(i,j).GT.0. .AND. t0.GT.thetacr(i,j)) THEN
			qbc(i,j) = qszero*17.*(SQRT(t0)-SQRT(thetacr(i,j)))*(t0-thetacr(i,j))
			us(i,j) = SQRT(tszero)*1.56*(t0/thetacr(i,j)-1.)**0.56
			
			vbc(i,j) = qbc(i,j)/us(i,j)
			qb(i,j) = MAX(0.,MIN(1.,vb(i,j)/vbc(i,j)))*qbc(i,j)
			sn = (/alphax(i,j),alphay(i,j)/)
			qbe(i,j) = sn(1)*qb(i,j); qbn(i,j) = sn(2)*qb(i,j)
			usx(i,j) = sn(1)*us(i,j); usy(i,j) = sn(2)*us(i,j)
			qbx(i,j) = sn(1)*qbc(i,j); qby(i,j) = sn(2)*qbc(i,j)
		ELSE
			us(i,j) = 0.0; usx(i,j) = 0.0; usy(i,j) = 0.0; qbc(i,j) = 0.0; vbc(i,j) = 0.0
			qb(i,j) = 0.0; qbe(i,j) = 0.0; qbn(i,j) = 0.0; qbx(i,j) = 0.0; qby(i,j) = 0.0
		ENDIF
		
		if (qbc(i,j)==0.) then
			us(i,j) = 0.0; usx(i,j) = 0.0; usy(i,j) = 0.0; qbc(i,j) = 0.0; vbc(i,j) = 0.0
			qb(i,j) = 0.0; qbe(i,j) = 0.0; qbn(i,j) = 0.0; qbx(i,j) = 0.0; qby(i,j) = 0.0
		endif
		
		if (qb(i,j).LT.0. .OR. ISNAN(qbe(i,j)) .OR. ISNAN(qbn(i,j)) .OR. ISNAN(vbc(i,j))) then
			write(6,*) i,j,qbc(i,j),us(i,j),qbx(i,j),qby(i,j),alphax(i,j),alphay(i,j),vbc(i,j),vba(i,j),qbe(i,j),qbn(i,j), &
				theta(i,j),thetacr(i,j),cf(i,j),cfg(i,j),cft(i,j),cfs(i,j),qb(i,j)
			pause '"qb" is a NaN'
		endif
	ENDDO
ENDDO

i = feedxs
do j = j1,j9
	us(1:i-1,j) = us(i,j); us(i9:ix,j) = us(i9-1,j)
	qbc(1:i-1,j) = qbc(i,j); qbc(i9:ix,j) = qbc(i9-1,j)
	vbc(1:i-1,j) = vbc(i,j); vbc(i9:ix,j) = vbc(i9-1,j)
	qbn(1:i-1,j) = qbn(i,j); qbn(i9:ix,j) = qbn(i9-1,j)
	qbx(1:i-1,j) = qbx(i,j); qbx(i9:ix,j) = qbx(i9-1,j)
	qby(1:i-1,j) = qby(i,j); qby(i9:ix,j) = qby(i9-1,j)
	usx(1:i-1,j) = usx(i,j); usx(i9:ix,j) = usx(i9-1,j)
	usy(1:i-1,j) = usy(i,j); usy(i9:ix,j) = usy(i9-1,j)
enddo

do j = j1,j9
	do i = feedxs,2,-1
		qbe(i-1,j) = qbe(i,j)+(qbe(i,j)-qbe(i+1,j))/15.
	enddo
	do i = i9-10,ix-1
		qbe(i+1,j) = qbe(i,j)+(qbe(i,j)-qbe(i-1,j))/15.
	enddo
enddo

! N-S
qbn(i1:i9,j1-2:j1-1) = -qbn(i1:i9,j1+1:j1:-1); qbn(i1:i9,j9+1:j9+2) = -qbn(i1:i9,j9:j9-1:-1)
qby(i1:i9,j1-2:j1-1) = -qby(i1:i9,j1+1:j1:-1); qby(i1:i9,j9+1:j9+2) = -qby(i1:i9,j9:j9-1:-1)

END SUBROUTINE SEDCAP
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SEDLIM
USE variables
IMPLICIT none  
integer	(kind=2)	:: i,j
real	(kind=rk)	:: a2(2),a3
! calculate sediment flux at the cell faces using cell centered qbe and qbn

! x-direction sediment flux
vc = qbe
call TVD(1)
fm(1:3,j1:j9) = vc(1:3,j1:j9); fm(i9-1:i9+1,j1:j9) = vc(i9-1:i9+1,j1:j9)
fp(1:3,j1:j9) = vc(2:4,j1:j9); fp(i9-1:i9+1,j1:j9) = vc(i9:i9+2,j1:j9)

do j = j1,j9
	do i = i1,i9-1
		a2 = (/zd(i+1,j)-zd(i,j),qbe(i+1,j)-qbe(i,j)/)
		a3 = fdir(a2)
		fbx(i,j) = 0.5*((1.+a3)*fm(i,j)+(1.-a3)*fp(i,j))
		
		if (ISNAN(fbx(i,j))) then
			write(6,*) i,j,fbx(i,j),vba(i-3:i+3,j),qbe(i-3:i+3,j),fm(i,j),fp(i,j),a3
			pause '"fbx" is a NaN'
		endif
	enddo
enddo

! y-direction sediment flux
vc = qbn
call TVD(2)
fm(i1:i9,1:2) = vc(i1:i9,1:2); fm(i1:i9,iy-2:iy-1) = vc(i1:i9,iy-2:iy-1)
fp(i1:i9,1:2) = vc(i1:i9,2:3); fp(i1:i9,iy-2:iy-1) = vc(i1:i9,iy-1:iy)

do j = j1,j9-1
	do i = i1+1,i9-1
		a2 = (/zd(i,j+1)-zd(i,j),qbn(i,j+1)-qbn(i,j)/)
		a3 = fdir(a2)
		fby(i,j) = 0.5*((1.+a3)*fm(i,j)+(1.-a3)*fp(i,j))
		
		if (ISNAN(fby(i,j))) then
			write(6,*) i,j,fby(i,j),vba(i-3:i+3,j),qbn(i-3:i+3,j),fm(i,j),fp(i,j),a3
			pause '"fbx" is a NaN'
		endif
	enddo
enddo

do j = j1,j9
	do i = i1+1,i9-1
		if (vba(i,j) .LE. smd) then
			if (fbx(i-1,j).LT.0.0)	fbx(i-1,j) = 0.0
			if (fbx(i,j).GT.0.0)	fbx(i,j) = 0.0
			if (fby(i,j-1).LT.0.0)	fby(i,j-1) = 0.0
			if (fby(i,j).GT.0.0)	fby(i,j) = 0.0
		endif
	enddo
enddo

! set boundary values <<<<<<<<<<<<<<<<<<<<<<<<<<<
! W-E
fbx(1:feedxs-1,j1:j9) = 0.0
fbx(i9-1:i9+1,j1:j9) = SPREAD(fm(i9-1,j1:j9),1,3)
fby(1:i1,j1:j9) = SPREAD(fby(i1+1,j1:j9),1,3)
fby(i9:ix,j1:j9) = SPREAD(fby(i9-1,j1:j9),1,3)

! N-S
fby(i1:i9,1:j1-1) = 0.0; fby(i1:i9,j9:j9+1)= 0.0

END SUBROUTINE SEDLIM
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SEDFEED
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j,n
real	(kind=4)	:: vbain0,c0(ny,2),sigma,mu,lambda,amp,Vb1,Vb2,Vc1,Vc2
! calculates channel average sediment transport capacity
! adjusts initial critical Shields parameter coefficient
! determines initial sediment supply rate
! determines sediment feed distirbution pattern

if (kcont.GT.INT(50/dt) .AND. kcont.LE.INT(100/dt)) then
	if (zbini==0.) then
		Qc0 = 0.4*SUM(MINVAL(qbx(ixn(1):ixn(2),j1:j9),DIM=2))/ixn(3)*BC
	else
		Qc0 = SUM(MINVAL(qbx(ixn(1):ixn(2),j1:j9),DIM=2))/ixn(3)*BC
		vban(i1:i9,j1:j9) = vbc(i1:i9,j1:j9)+zba(i1:i9,j1:j9)
	endif
	
	if (Qc0 .GT. 1.001*qstot) then
		rc(1) = rc(1)*1.005
	elseif (Qc0 .LT. 0.999*qstot) then
		rc(1) = rc(1)*0.995
	endif
	
	if (kcont == INT(100/dt)) then
		Qs0 = qstot*scrate
		if (SUM(thetacg(ixn(1):ixn(2),j1:j9))/(ixn(3)*ny) .LT. 1.e-3) then
			write(6,*) Qc0*2.65e06,thetac0, &
				SUM(thetacg(ixn(1):ixn(2),j1:j9))/(ixn(3)*ny),SUM(thetacr(ixn(1):ixn(2),j1:j9))/(ixn(3)*ny), &
				MINVAL(ks(i1:i9,j1:j9)),MAXVAL(ks(i1:i9,j1:j9)),SUM(ks(i1:i9,j1:j9))/(nx*ny), &
				MINVAL(kst(i1:i9,j1:j9)),MAXVAL(kst(i1:i9,j1:j9)),SUM(kst(i1:i9,j1:j9))/(nx*ny), &
				MINVAL(ksf(i1:i9,j1:j9)),MAXVAL(ksf(i1:i9,j1:j9)),SUM(ksf(i1:i9,j1:j9))/(nx*ny)
			pause '"theta_cr" is too small'
		endif
	endif
endif

! sediment feed
vbain(i1:i9,j1:j9) = 0.0

if (kcont .GT. INT(100/dt)) then
	vbain0 = Qs0/(dx*dy)*dt
else
	if (zbini==0.) then
		vbain0 = 0.0
	else
		vbain0 = SUM(fbx(feedxs,j1:j9))/dx*dt
	endif
endif

if (feed_s==0 .OR. kcont.LE.INT(100/dt)) then
	vbain(1,j1:j9) = vbain0/ny
elseif (feed_s==1) then
	c0(:,1) = (/1:ny/)
	mu = sum(c0(:,1))/ny
	! normal distribution
	sigma = fsigma; lambda = flambda
	c0(:,2) = 0.5*lambda*exp(0.5*lambda*(2.*mu+lambda*sigma**2-2.*c0(:,1))) &
				*erfc((mu+lambda*sigma**2.-c0(:,1))/(sqrt(2.)*sigma))
	vbain(feedxs,j1:j9) = c0(:,2)/sum(c0(:,2))*vbain0
	do j = j1,j9
		if (ISNAN(vbain(feedxs,j))) then
			write(6,*) i,j,Qs0,vbain0,mu
			pause '"feed" is a NaN'
		endif
		if (vbain(feedxs,j).lt.0.0) pause '"feed" is a negative'
	enddo
elseif (feed_s==2) then
	vbain(feedxs,j1:j9) = qbc(ixn(2),j1:j9)*dt/dx
elseif (feed_s==3) then
	vbain(feedxs,j1:j9) = vbain0*d(feedxs,j1:j9)/SUM(d(feedxs,j1:j9))
elseif (feed_s==4 .AND.  kcont.GT.INT(100/dt)) then
	if (kcont .GT. INT(100/dt)) then
		Vb1 = SUM(hd(feedxs+40,j1:j9)); Vb2 = SUM(hd(i9-40,j1:j9))
		Vc1 = MINVAL(Pc(feedxs,j1:j9))
		if (Vc1 == 1.) then
			famp = MAX(0.01,famp*0.95)
		else
			if (Vb2 .LT. Vb1) then
				famp = MIN(1.0,famp*1.05)
			else
				famp = MAX(0.01,famp*0.95)
			endif
		endif
	endif
	
	c0(:,1) = (/1:ny/)
	mu = sum(c0(:,1))/ny; sigma = 5.0; amp = famp
	c0(:,2) = amp+0.5*(1.+erf((c0(:,1)-mu)/(sigma*sqrt(2.0))))
	if (SUM(ZD(nw,j1:20)).GT.SUM(ZD(nw,21:j9))) then
		vbain(1,j1:j9) = c0(:,2)/SUM(c0(:,2))*vbain0
	else
		vbain(1,j9:j1:-1) = c0(:,2)/SUM(c0(:,2))*vbain0
	endif
endif

fbx(feedxs-1,j1:j9) = vbain(1,j1:j9)*dx/dt; vbain(feedxs,j1:j9) = vbain(1,j1:j9)

END SUBROUTINE SEDFEED
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SEDTRANS
USE variables
IMPLICIT none
! sediment flux changes calculation

! calculate friction coefficient, shear stress
call FRICTION
! radius of curvature
call SRC
! dimensionless shear stress, transport angle
call SEDMOV
! dimensionless critical shear stress
call SHIELDS
! Ripple factor
call SEDRPF
! calculate sediment transport capacity
call SEDCAP
! calculate sediment transport rate at the cell faces
call SEDLIM
! sediment feed
call SEDFEED

! divergence of sediment flux
FG0(i1:i9,j1:j9) = (fbx(i1:i9,j1:j9)-fbx(i1-1:i9-1,j1:j9))/dx+(fby(i1:i9,j1:j9)-fby(i1:i9,j1-1:j9-1))/dy
FG0(1:feedxs-1,j1:j9) = 0.0

END SUBROUTINE SEDTRANS
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE MORPH
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j
real	(kind=4)	:: vba2,zba2,vbc2,vb2
! sediment transport time evolution
! bedload layer and alluvial layer update

call SedTrans

! vba = vban-dt*FG0+vbain
if (kcont.GE.INT(50/dt)) then
	vba(i1:i9,j1:j9) = vban(i1:i9,j1:j9)-dt*FG0(i1:i9,j1:j9)
else
	vba(i1:i9,j1:j9) = vban(i1:i9,j1:j9)
endif

call SEDCON

do j = j1,j9
	do i = i1,i9
		vba2 = vba(i,j)
		vbc2 = vbc(i,j)
		if (vba2.GT.vbc2) then
			vb2 = vbc2
			zba2 = vba2-vbc2
		else
			vb2 = MAX(0.0,vba2)
			zba2 = 0.0
		endif
		
		if (vb2.LT.smd) then
			vb(i,j) = 0.0
		else
			vb(i,j) = vb2
		endif
		zba(i,j) = zba2
	enddo
enddo

END SUBROUTINE MORPH
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SEDCON
USE variables
IMPLICIT none
integer	(kind=2)	:: i,j,k1(1),k2
real	(kind=4)	:: vba2,vbneg(4)
! negative sediment volume treatment in bedrock cells

vzb = vba
vzb(i1-2:i1-1,:) = spread(vzb(i1,:),1,2)
vzb(i9+1:i9+2,:) = spread(vzb(i9,:),1,2)
vzb(:,j1-2:j1-1) = spread(vzb(:,j1),2,2)
vzb(:,j9+1:j9+2) = spread(vzb(:,j9),2,2)
do j = j1,j9
	do i = i1,i9
		vba2 = vzb(i,j)
		if (vba2 .LT. smd) then
			vbneg = (/vzb(i-1,j),vzb(i,j-1),vzb(i+1,j),vzb(i,j+1)/)-smd
			k1 = MAXLOC(vbneg)
			k2 = k1(1)
			if ((vbneg(k2)+vba2) .GE. smd) then
				if (mod(k2,2) == 1) then
					vzb(i-2+k2,j) = vzb(i-2+k2,j)+vba2
				elseif (mod(k2,2) == 0) then
					vzb(i,j-3+k2) = vzb(i,j-3+k2)+vba2
				endif
				vzb(i,j) = smd
			endif
		endif
	enddo
enddo
vba(i1:i9,j1:j9) = vzb(i1:i9,j1:j9)

END SUBROUTINE SEDCON
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE SW2D
USE variables
IMPLICIT none
! HYDRO: use 2nd-order Runge-Kutta method for flow calculation
! MORPH: sediment transport and bed evolution calculation

do rkn = 1,2
	call HYDRO
	call BoundCond(rkn)
enddo
call err_test
rkn = 3
call MORPH
call BoundCond(rkn)

END SUBROUTINE SW2D
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

SUBROUTINE err_test
USE variables
IMPLICIT none    
integer	(kind=2)	:: i,j,nan_flag,err_negs
! check to make sure we do not have NaNs & negative depth

nan_flag = 0
err_negs = 0
do j = j1,j9
	do i = i1,i9
		if (isnan(d		(i,j))) nan_flag = nan_flag+1
		if (isnan(du	(i,j))) nan_flag = nan_flag+1
		if (isnan(dv	(i,j))) nan_flag = nan_flag+1
		if (isnan(u		(i,j))) nan_flag = nan_flag+1
		if (isnan(v		(i,j))) nan_flag = nan_flag+1
		if (isnan(zb	(i,j))) nan_flag = nan_flag+1
		if (isnan(vb	(i,j))) nan_flag = nan_flag+1
		if (isnan(vbc	(i,j))) nan_flag = nan_flag+1
		if (isnan(vba	(i,j))) nan_flag = nan_flag+1
		if (isnan(pc	(i,j))) nan_flag = nan_flag+1
		
		if (d(i,j) .LT. smd) err_negs = err_negs+1
		if (zb(i,j) .LT. 0.) err_negs = err_negs+1
	enddo
enddo

186 format(2X,2(i5,2X),36(ES15.8E2,2X))
if (nan_flag.GT.0 .OR. err_negs.GT.0) then
	num_errs = num_errs+1
	write(test1,'(i8.8)') num_errs
	if (nan_flag.GT.0) then
		nome = "va_NaNs"//test//".var"//char(0)
	elseif (err_negs.GT.0) then
		nome = "va_NEGs"//test//".var"//char(0)
	endif
	
	open(unit=unit_var,file=nome)
	write(unit_var,'(2x,e15.8)')time
	do i = i1,i9
		do j = j1,j9
			write(unit_var,186)i,j,d(i,j),h(i,j),z(i,j), &
				u(i,j),du(i,j),v(i,j),dv(i,j),ksg(i,j), &
				kst(i,j),ksf(i,j),cfg(i,j),cf(i,j),cff(i,j), &
				rpf(i,j),dvL(i,j,1),dvR(i,j,1),Sx(i,j),Sy(i,j), &
				theta(i,j),thetacr(i,j),qbc(i,j),qbx(i,j),qby(i,j), &
				us(i,j),qbe(i,j),qbn(i,j),fx(i,j),fy(i,j), &
				fbx(i,j),fby(i,j),FG0(i,j),vbc(i,j),vb(i,j), &
				vba(i,j),zb(i,j),pc(i,j)
		enddo
	enddo
	close(unit_var)
endif  
if (err_negs.GT.0) write(6,'(A,f15.6)') 'error prof. NEGs at time ', time
if (nan_flag.GT.0) write(6,'(A,f15.6)') 'At least 1 value is NaN at time ', time

if (err_negs.gt.0) then
	do j = j1,j9
		do i = i1,i9
			if (d(i,j).lt.0.0 .or. zb(i,j).lt.0.0) then
				write(6,*) i,j
				PAUSE 'negative flow or sediment depth'
			endif
		enddo
	enddo
endif

if (nan_flag.GT.0) PAUSE 'NaN'
if (num_errs.GT.100) PAUSE 'Error 100'

END SUBROUTINE err_Test
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________

PROGRAM CHANNEL2D

! *******************************************************************
! 2-Dimensional Morphodynamic Model
! 
! Written by Jongseok Cho and Peter A. Nelson
! Department of Civil and Environmental Engineering
! Colorado State University, Fort Collins, Colorado
! *******************************************************************

USE variables
IMPLICIT none
integer		(kind=2)	:: i,j
character	(len=5)		:: header(12)

call allocation
call initialization

171 format(A,f16.2,A,i10,A,i10,A,f10.8)
172 format(2X,12(A6,8X))

181 format(2X,2(i5,2X),12(ES15.8E2,2X))
182 format(2X,2(i5,2X),28(ES15.8E2,2X))

188 format(2X,12(ES12.5E2,2X))
189 format(2X,7(ES12.5E2,2X))


idout = 'chnl_'
header = (/'x','y','z','u','v','d','Pc','Vb','Vbc','zb','qbc','vbain'/)

call input
call setbed
time = 0.0

write(test,'(i8.8)') time*10
nome = trim(idout)//test//".out"//char(0)
open(unit=unit_conc,file=nome)
write(unit_conc,172) header
do i = i1,i9
	do j = j1,j9
		write(unit_conc,188) XI(i,j),YI(i,j),z(i,j),u(i,j),v(i,j), &
			d(i,j),pc(i,j),vb(i,j),vbc(i,j),zb(i,j), &
			qbc(i,j),vbain(i,j)
	enddo
enddo
close(unit_conc)

nome = "abnd_"//test//".out"//char(0)
open(unit=unit_conc,file=nome)
write(unit_conc,'(2X,E15.8)') time
do i = i1,i9
	do j = j1,j9
		write(unit_conc,188) XI(i,j),YI(i,j),zbr0(i,j),zra(i,j),zbr(i,j), &
			z(i,j),ksb(i,j)
	enddo
enddo
close(unit_conc)

! start time stepping
do kcont = 1,ntot
	time = kcont*dt+nc_ini*0.1
	! ******************************************************************
	call SW2D
	if (MOD(kcont,NINT(1/dt)) == 0.) then
		write(*,171,advance='YES') 'Time: ',time, ', kcont: ',kcont, &
			', ntot: ',ntot,' kcont/ntot: ',REAL(kcont)/ntot
	endif
	! ******************************************************************
	
	! OUTPUT
	if (kcont==np(npc)) then
		write(test,'(i8.8)') int(time*10)
		nome = trim(idout)//test//".out"//char(0)
		open(unit=unit_chnl,file=nome)
		write(unit_chnl,172) header
		do i = i1,i9
			do j = j1,j9
				write(unit_chnl,188) XI(i,j),YI(i,j),z(i,j),u(i,j),v(i,j), &
					d(i,j),pc(i,j),vb(i,j),vbc(i,j),zb(i,j), &
					qbc(i,j),vbain(i,j)
			enddo
		enddo
		close(unit_chnl)
		
		if (kcont==NINT(100/dt) .OR. MOD(npc,60)==0) then
			nome = "variab_"//test//".var"//char(0)
			open(unit=unit_var,file=nome)
			write(unit_var,'(2X,E15.8)') time
			do i = 2,ix-1
				do j = 2,iy-1
					write(unit_var,182) i,j,d(i,j),h(i,j),z(i,j), &
						du(i,j),dv(i,j),cf(i,j),Sbx(i,j),Sby(i,j), &
						Dfx(i,j),Dfy(i,j),theta(i,j),thetacr(i,j),alpha(i,j), &
						rpf(i,j),qbc(i,j),us(i,j),qbe(i,j),qbn(i,j), &
						fbx(i,j),fby(i,j),vb(i,j),vbc(i,j),pc(i,j), &
						zb(i,j),ksg(i,j),kst(i,j),ksf(i,j),cfg(i,j)
				enddo
			enddo
			close(unit_var)
		endif
		
		npc = npc+1
	endif
enddo

call deallocation

END PROGRAM CHANNEL2D
! _________________________________________________________________________________________________
! _________________________________________________________________________________________________