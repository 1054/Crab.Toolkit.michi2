c       ===========================================================================
 	PROGRAM GET_INFRARED_COLORS
c       ---------------------------------------------------------------------------
c               Authors :   E. da Cunha & S. Charlot
c       Latest revision :   Sep. 15th, 2010
c       ---------------------------------------------------------------------------
c       Computes AB magnitudes of a library of near-infrared-to-submm model
c       spectra for a given set of photometric bands, at a given redshift.
c
c       INPUTS: - model library [default: InfraredLIB_newmir.bin]
c               - cosmological parameters [default: 70,30,30]
c               - redshift
c               - filter file - define USER_FILTERS in

c       OUTPUT: .lbr file containing: physical parameters + magnitudes
c       
c       REVISION HISTORY:
c       - Aug. 2009 - L. Dunne:
c       program has been modified to fix bugs with the MIPS filters. For SHADES/SWIRE
c       MIPS fluxes are calibrated as for IRAC and so the BB part has been commented out.
c       Some integrations have been taken out of the main loop to speed things up.
c
c       Version: Feb 16, 2011
c       changed to include new IR SEDs extending out to the radio
c       number of points in SED: 6450 -> 6750
c       ===========================================================================

	implicit none
        character infile*80,outfile*80,user_filt*80
        integer niw,io,nfilt_use,nmax,iz,index
	parameter(nmax=50)
	character*10 filt_name(nmax)
	character filter_header*250
        integer i,nc,imod,k
	integer k_use(nmax),filt_id_use(nmax)
        real z,xaux(nmax),mags(nmax),irlums(3)
	real wl(6750),irsed(6750),irprop(9),redshift(161)
	real lambda_eff(nmax),lambda_rest(nmax)
	real xi_pah,xi_mir,xi_warm,xi_cold,tdust
	character*6 numz
	integer nfilt,filt_id(100),fit(100),ifilt
        real h,omega,omega_lambda,clambda,q,cosmol_c,t,tu,dm,dismod
        data h/70./,omega/0.30/,omega_lambda/0.70/

c       INPUT file: spectral library (dust emission)
	call getenv('IRLIB',infile)

c       INPUT redshift
        write (6,'(x,a,$)') 'Enter redshift = '
        read (5,*,end=1) z

c       OUTPUT file: catalogue of parameters & magnitudes @ redshift z
	write(numz,'(f6.4)') z
	outfile='infrared_dce08_z'//numz//'.lbr'
	close (30)
        open (30,file=outfile,status='unknown')

c       Cosmological parameters
        write (6,'(/x,a,f4.0,a,f5.3,a,$)') 'Enter Cosmology [default: 70.,0.3,0.7]'
10      write (6,'(x,a,f4.0,a,f5.2,a,f5.2,a,$)')  'Enter Ho [',h,'], Omega [',
     +          omega,'], Omega_lambda [',omega_lambda,'] = '
        call nread(xaux,nc,*10,*1)
        if (nc.gt.0.) then
           h=xaux(1)
           omega=xaux(2)
           omega_lambda=xaux(3)
        endif

c       Obtain cosmological constant and q
        clambda=cosmol_c(h,omega,omega_lambda,q)

c       Age of the universe at z
        tu=t(h,q,z,clambda)*1.e9
        write (6,'(x,a,f4.2,a,f6.2,a  )') 'Age of this universe at z=',z,
     +                                    ':',tu/1.e9,' Gyr'

c       Distance modulus
        dm=dismod(h,q,z)
        write (6,'(x,a,f4.2,a,f6.2,a  )') 'Distance modulus at z=',z,
     +                                    ':',dm

c       Read filter file
	call getenv('USER_FILTERS',user_filt)
	close(22)
	open(22,file=user_filt,status='old')
	do i=1,1
	   read(22,*)
	enddo
	io=0
	ifilt=0
	do while(io.eq.0)
	   ifilt=ifilt+1
	   read(22,*,iostat=io) filt_name(ifilt),lambda_eff(ifilt),filt_id(ifilt),fit(ifilt)
	enddo
	nfilt=ifilt-1

c       WHAT FILTERS TO USE HERE?
c       Compute optical colours for filters sampling the rest-frame
c       emission at wavelengths longer tahn 2.5 microns
c       (if lambda < 2.5 microns, contribution from dust emission is negligible)

c       count nr of filters to use: nfilt_use
	nfilt_use=0

	do i=1,nfilt
	   lambda_rest(i)=lambda_eff(i)/(1.+z)
	   if (lambda_rest(i).gt.2.5) then
	      nfilt_use=nfilt_use+1
	      k_use(nfilt_use)=i
	   endif
	enddo

	write(filter_header,*) (filt_name(k_use(k)),k=1,nfilt_use)

	do k=1,nfilt_use
	   filt_id_use(k)=filt_id(k_use(k))
	enddo
	  
c       Output File Header
c       Store cosmological parameters
c       List of model parameters and predicted fluxes
	   write (30,'(8(a10),9(a12),a4,a250)')
     +  '#    index',   ! 0: index in the InfraredLIB library
     +  '       fmu',   ! 1: fraction of Ld(tot) contributed by the diffuse ISM
     +  '   xiC_ISM',   ! 2: contribution of cold dust to total ISM luminosity
     +  '     TW_BC',   ! 3: temperature of warm dust in birth clouds
     +  '    TC_ISM',   ! 4: temperature of cold dust in ISM
     +  '  xiPAH_BC',   ! 5: contribution of PAH to total BC luminosity
     +  '  xiMIR_BC',   ! 6: contribution of MIR continuum to total BC luminosity
     +  '    xiW_BC',   ! 7: contribution of warm dust to total BC luminosity
     +  '     Mdust',   ! 8: total dust mass (in Msun)
     +  '      q_IR',   ! 9: radio/far-IR
     +  '      L_IR',   !10: IR luminosity (8 - 1000 microns)
     +  '     L_FIR',   !11: far-IR luminosity (42.5 - 122.5 microns)
     +  '  Tdust_av',   !12: luminosity-weighted dust temperature
     +  ' xiPAH_tot',   !13: Ld(PAH)/Ld(tot)
     +  ' xiMIR_tot',   !14: Ld(MIR)/Ld(tot)
     +  '   xiW_tot',   !15: Ld(warm)/Ld(tot)
     +  '   xiC_tot',   !16: Ld(cold)/Ld(tot)
     +  '    ',filter_header  !17--nfilt_use: model AB magnitudes at redshift z

c       ---------------------------------------------------------------------------       
c       Read .bin file (300 models)
c       Compute magnitude for each model
c       Store parameters + magnitudes of each model in output file

	close (29)
	open (29,file=infile,status='old',form='unformatted')

        read (29) niw,(wl(i),i=1,niw)
        write (*,*)  'Number of wavelength points = ',niw

	index=0
        do imod=1,25000
	   index=index+1
              read (29,end=1) (irprop(i),i=1,9)
	      read (29) (irsed(i),i=1,niw)
	      read (29) (irlums(i),i=1,3)
c     re-define xi parameters as contributions to the total IR luminosity
c     Re-define IR parameters: xi^tot
	 xi_pah=irprop(5)*(1.-irprop(1))+0.550*(1-irprop(2))*irprop(1)  ! xi_PAH^tot Ld(PAH)/Ld(tot)
         xi_mir=irprop(6)*(1.-irprop(1))+0.275*(1-irprop(2))*irprop(1)  ! xi_MIR^tot Ld(MIR)/Ld(tot)
         xi_warm=irprop(7)*(1.-irprop(1))+0.175*(1-irprop(2))*irprop(1) ! xi_W^tot Ld(warm)/Ld(tot)
         xi_cold=irprop(1)*irprop(2)                                    ! xi_C^tot Ld(cold)/Ld(tot)
c     compute average (luminosity-weighted) dust temperature
c     Tdust_average = ( xi_W^tot * T_W^BC + xi_C^tot * T_C^ISM + 0.07 * 45 * fmu) / (xi_W^tot + xi_C^tot + 0.07*fmu)
         tdust=xi_warm*irprop(3) + xi_cold*irprop(4) + 0.07*45.*irprop(1)
         tdust=tdust/(xi_warm+xi_cold+0.07*irprop(1))
                 call model_ab_color(z,wl,irsed,niw,nfilt_use,filt_id_use,mags)
c                write output file:
                 write (30,200) index,(irprop(i),i=1,9),irlums(2),irlums(3),
     +                          tdust,xi_pah,xi_mir,xi_warm,xi_cold,(mags(i),i=1,nfilt_use)
200              format(i10,0p7f10.3,1pe12.3,0pf10.3,1p2e12.3,0p5f12.4,0p25f10.4)
        enddo
1       stop

	end
c       ===========================================================================


c       ===========================================================================
	SUBROUTINE MODEL_AB_COLOR(z,x,yd,inw,nf,ifilt,mag)
c       ===========================================================================
c	Computes AB magnitude of each model at given z in each band
c       ---------------------------------------------------------------------------
c       z    : redshift
c       x    : wavelength vector (in A)
c       yd   : total dust SED L_lambda (in Lsun/A)
c       inw  : number of wavelength points
c       nf   : number of photometric bands (filters)
c       ifilt: array containing the indexes of the filters in response file
c       mags : magnitudes
c       ---------------------------------------------------------------------------
        implicit none
        integer nf
        integer i,inw,icall
	integer ifilt(25)
	real x(inw),yd(inw),fx(nf)
	real z,f_mean,dl,mag(nf)
	data icall/0/

	if (icall.eq.0) then
		icall=1
	endif

c	Compute flux through each of nb filters
        do i=1,nf
	   fx(i)=f_mean(ifilt(i),x,yd,inw,z)
        enddo

c       Compute absolute (k-shifted) AB magnitude
c       10 pc in units of Mpc
        dl=1.e-5
	do i=1,nf
	   mag(i)=-2.5*alog10(fx(i))-48.6
	   mag(i)=mag(i)+(5.*alog10(1.7684e+08*dl)) ![this factor is SQRT(4*pi*(3.0856E24)^2/Lsun)]
	   !deal with INFs because SED doesnt have optical
	   if (mag(i).gt.10000) then
	      mag(i)=-99.
	   endif
	enddo

        return
	end

c       ===========================================================================
	SUBROUTINE FILTER0
c       ===========================================================================
c	Reads filter response functions from (binary) file FILTERBIN.RES
c       ---------------------------------------------------------------------------
	INCLUDE 'filter.dec'
	character filtfile*80
	close (81)

c	Get file name from environment variable FILTERS
	call getenv('FILTERS',filtfile)
	open (81,file=filtfile,form='unformatted',status='old')

	l=largo(filtfile)
	write (6,'(/x,2a)') 'Reading Filter File: ',filtfile(1:l)
	read (81,err=1) nf,(ni(i),i=1,nf),(nl(i),i=1,nf),(np(i),i=1,nf),
     &		(fid(i),i=0,nf),ltot,(r(i,1),i=1,ltot),(r(i,2),i=1,ltot)
	close (81)
	iread=1
	write (6,'(i4,a,i3,a,$)') nf,' filters defined, out of ',mxft,' maximum'
	write (6,'(x,a)') '    ...done'
	return
1	stop 'Program exits because of error reading file FILTERBIN.RES'
	end


c       ===========================================================================
	FUNCTION F_MEAN(I,X,Y,N,Z)	
c       ===========================================================================
c       Returns mean flux >>> (Fnu) <<< in ith filter (/Hz)
c       As: f_nu=int(dnu Fnu Rnu/h*nu)/int(dnu Rnu/h*nu)
c       ie: f_nu=int(dlm Flm Rlm lm / c)/int(dlm Rlm/lm)

c	The s.e.d. is Y(X) with N data points. Z is the redshift
c	to be applied to the s.e.d.

c	Improved version that samples correctly the s.e.d. and the
c	filter response function at all points available in these arrays.

c       The new version of the F_MEAN routine requires the filter
c	response function to be interpolated at each point in the
c	sed which is not a point in the filter response. All points in
c	the filter response are also used. The format of the binary
c	file was changed to reduce its size. Filters are stored
c	sequentially in r(i,k), k=1 => wavelength, k=2 => response
c	function. Information about starting point and number of
c	points per filter is kept in file, as well as filter id label.

c       i    : index of the filter
c       x    : wavelength vector (in A)
c       y    : SED L_lambda (in Lsun/A)
c       n    : number of wavelength points
c       z    : redshift
c       ---------------------------------------------------------------------------
	INCLUDE 'filter.dec'
	integer pos_i(mxft),pos_f(mxft),ifilt(mxft)
	real x(n),y(n),LINEAR
	real xlam(40000),rlam(40000),xf(8000),rf(8000),rfa(8000),rfb(8000)
	real rfc(8000),rfd(8000),xeff(mxft)
	real*8 h,c,kappa,t0
	parameter(h=6.6261e-27) !erg s
	parameter(c=2.99792478e+10) !speed of light cgs
	parameter(kappa=1.3807e-16) !erg/grad
	parameter(t0=1.e+4) !K
	common /save_filter_2/ nfilt,ntot,pos_i,pos_f,ifilt,xlam,rlam
	data zzlast/-100./

	f_mean=0.
c	Check filter number
	if (i.gt.mxft) then
		write (6,*) 'Filter No.',i,' not available.'
		return
	elseif (i.le.0) then
c       [Introduced to allow computation of K-correction with routine ~/is/k_correct.f]
		f_mean=1.
		return
	endif

c	Read filter file
	if (iread.eq.0) Call FILTER0

c	Check z value
	if (z.ne.zzlast.or.ireset.gt.0) then
		ireset=0
		zzlast=z
		nfilt=0
		ntot=0
		pos_i(1)=1
		do k=1,mxft
		ifilt(k)=-10
		enddo
	endif


c	Extract ith-filter. Shift by (1+z)
	z1=1.+z
	m=0

c	Extract filter
c	Check if filter has been stored already
	do k=1,nfilt
	if (ifilt(k).eq.i) goto 2
	enddo

	do k=ni(i),nl(i)
	m=m+1
	xf(m)=r(k,1)/z1
	rf(m)=r(k,2)

	rfc(m)=rf(m)*xf(m)
	rfd(m)=rf(m)

	enddo

c       Effective wavelength of the filter:
	xeff(i)=z1*TRAPZ1(xf,rfc,m)/TRAPZ1(xf,rfd,m)

c	Add wavelength points in the sed.
	l=0
	do k=1,n
c	xz=x(k)
	xz=x(k)/z1
	if (xz.ge.xf(np(i))) then
		goto 1
	elseif (xz.gt.xf(1)) then
		m=m+1
		xf(m)=xz
c		Interpolate shifted filter at shifted wavelength of sed
		rf(m)=LINEAR(xz,xf,rf,np(i),l)
	endif
	enddo

c	Sort (xf,rf) arrays according to xf
1	call sort2(m,xf,rf)

c	Store sorted arrays
	do k=1,m
	ntot=ntot+1
	xlam(ntot)=xf(k)
	rlam(ntot)=rf(k)
	enddo

c	Store filter No. and ending point in sorted arrays
	nfilt=nfilt+1
	ifilt(nfilt)=i
	pos_f(nfilt)=ntot
	pos_i(nfilt+1)=ntot+1
	k=nfilt
	kt=pos_f(k)-pos_i(k)+1
	write(6,'(a,i2,a,i4,a,a)')'#',nfilt,', ID:',i,' -- ',(fid(i)(1:60))
	write(*,'(a,f15.4,a)') 'lambda_eff= ', xeff(i)*1.e-4, ' microns'

c	Use kth filter in extracted filter array
c	Interpolate sed at shifted wavelength of filter
2	m=0
	l=0	   
	do j=pos_i(k),pos_f(k)
	   m=m+1
c	   take wavelength in detector frame
	   xf(m)=xlam(j)*z1
	   rfa(m)=rlam(j)*xf(m)*LINEAR(xlam(j),x,y,n,l)
	   rfb(m)=rlam(j)
	   enddo

c       Compute flux below filter
	f_mean=TRAPZ1(xf,rfa,m)/TRAPZ1(xf,rfb,m)
	f_mean=xeff(i)*f_mean/2.997925e+18
c	F(lambda)*dlambda = F[lambda/(1+z)]*dlambda/(1+z)
	f_mean=f_mean/z1

	return
	end

c       ===========================================================================
	SUBROUTINE NREAD(X,NA,*,*)
c       ===========================================================================
c	This routine is useful when reading several parameters separated by
c	commas (not supported by current version of f77 in Sun 4/110).

c	Returns NA arguments in array X (read from the screen).
c	RETURN 1: Statement to execute if error while reading.
c	RETURN 2: Statement to execute if EOF detected.
c       ---------------------------------------------------------------------------
	parameter (npar=24)
	character b*132
	real x(npar)

c	Clear buffer
	do i=1,npar
	x(i)=0.
	enddo

c	Read string b
	read (5,'(a)',end=2) b

c	Number of characters read into b
c	na=index(b,' ')-1
	na=largo(b)
	if (na.eq.0) return

c	Adds "," at the end of string b to guarantee extraction of last value
	b(na+1:na+1)=','

c	Extract numbers from b
	n=0
	l1=0
	do i=1,100
	l1=l1+n+1
	n=index(b(l1:),',')-1
	if (n.eq.0) then
		x(i)=0.
	elseif (n.gt.0) then
      		read(b(l1:l1+n-1),'(G132.0)',err=1) x(i)
	else
		na=i-1
		return
	endif
	enddo
1	return 1
2	return 2

	ENTRY QREAD(V,K,*,*)

c	Emulates q format of VMS/Fortran

c	V = required value
c	K = number of characters typed
c	RETURN 1: Statement to execute if error while reading.
c	RETURN 2: Statement to execute if EOF detected.

c	Read string b and assign value to v
	read (5,'(a)',end=4) b
c	k = index(b,' ') - 1
	k=largo(b)
	v=0.
	if (k.gt.0) then
      		read(b(:k),'(G132.0)',err=3) v
	endif
	return
3	return 1
4	return 2
	end

c       ===========================================================================
	REAL*4 FUNCTION COSMOL_C(h,omega,omega_lambda,q)
c       ===========================================================================
c	Returns cosmological constant = cosmol_c and parameter q

c       Omega is entered by the user
c       omega=1.-omega_lambda
c       ---------------------------------------------------------------------------

c       Cosmological constant
        cosmol_c=omega_lambda/(3.*h**2)

c       Compute q=q0
        if (omega_lambda.eq.0.) then
                q=omega/2.
        else
                q=(3.*omega/2.) - 1.
        endif
	return
	end

c       ===========================================================================
        REAL*4 FUNCTION T(h, q, z, lamb)
c       ===========================================================================
c	Returns age of the Universe at redshift z
c       for the specified cosmological parameters
c       ---------------------------------------------------------------------------
	implicit none
	real*4 h, q, z, a, b, c, hh0, cosh, x, lamb
	real*4 aa, bb, epsr, s0, s, funq, omega0
	integer npts
	logical success
	external funq
	common /cosm/ omega0

c	H = Ho in km/sec/Mpc
c	Q = qo  (if problems with qo = 0, try 0.0001)

	a(q,z) = (sqrt(1. + ((2. * q) * z)) / (1. - (2. * q))) / (1. + z)
	b(q) = q / (abs((2. * q) - 1.) ** 1.5)
	c(q,z) = ((1. - (q * (1. - z))) / q) / (1. + z)

	cosh(x) = alog(x + sqrt((x ** 2) - 1.))
	hh0 = h * 0.001022
c       in (billion years)**(-1)

	if (lamb .ne. 0.0) then
         	omega0 = (2. * (q + 1.)) / 3.
         	aa = 0.
         	bb = 1. / (1. + z)
	        success=.false.
		s0=1.e-10
		npts=0
		do while (.not.success)
		  npts=npts+1
		  callmidpnt(funq,aa,bb,s,npts)
		  epsr=abs(s-s0)/s0
		  if (epsr.lt.1.e-4) then
			success=.true.
		  else
			s0=s
		  endif
		enddo
		t=s
	     else if (q .eq. 0.0) then
		t = 1. / (1. + z)
	else if (q .lt. 0.5) then
	   t = a(q,z) - (b(q) * cosh(c(q,z)))
	else if (q .eq. 0.5) then
	   t = (2. / 3.) / ((1. + z) ** 1.5)
	else
	   t = a(q,z) + (b(q) * cos(c(q,z)))
	end if

	t = t / hh0

	return
	end


c       ===========================================================================
        REAL*4 FUNCTION FUNQ(x)
c       ===========================================================================
c       for non-zero cosmological constant
c       ---------------------------------------------------------------------------
	real*4 x, omega0, omegainv
	common /cosm/ omega0
	omegainv = 1. / omega0
	funq = sqrt(omegainv) / (x*sqrt((omegainv-1.)+(1./(x**3.))))
	return
	end


c       ===========================================================================
	FUNCTION DISMOD(H0,Q0,Z)
c       ===========================================================================
c	Returns cosmological distance modulus
c       ---------------------------------------------------------------------------
c	H0 = Ho in km/sec/Mpc
c	Q0 = qo
c	DL = Luminosity distance in Mpc

	dismod = 5*alog10(dl(h0,q0,z)*1.E6/10.)

	return
	end


c       ===========================================================================
	REAL*4 FUNCTION DL(h,q,z)
c       ===========================================================================
c	Computes luminosity distance corresponding to a redshift z.
c	Uses Mattig formulae for qo both 0 and non 0
c	Revised January 1991 to implement cosmolgical constant
c	Ho in km/sec/Mpc, DL is in Mpc
c       ---------------------------------------------------------------------------
	implicit none
	real*4 h, q, z, d1, d2
	real*4 aa, bb, epsr, s, s0, funl
	real*4 dd1, dd2, omega0
	logical success
	integer npts
	external funl
	common /cosm/ omega0

	if (z.le.0.) then
	   dl=1.e-5		!10pc
	   return
	endif

	if (q .eq. 0) then
	   dl = ((3.e5 * z) * (1 + (z / 2.))) / h
	else if (q .gt. 0.) then
	   d1 = (q * z) + ((q - 1.) * (sqrt(1. + ((2. * q) * z)) - 1.))
	   d2 = ((h * q) * q) / 3.e5
	   dl = d1 / d2
	else if (q .lt. 0.) then
	   omega0 = (2. * (q + 1.)) / 3.
	   aa = 1.
	   bb = 1. + z
	   success=.false.
	   s0=1.e-10
	   npts=0
	   do while (.not.success)
	      npts=npts+1
	      call midpnt(funl,aa,bb,s,npts)
	      epsr=abs(s-s0)/s0
	      if (epsr.lt.1.e-4) then
		 success=.true.
	      else
		 s0=s
	      endif
	   enddo
	   dd1=s
	   dd2 = (3.e5 * (1. + z)) / (h * sqrt(omega0))
	   dl = dd1 * dd2
	end if

      return
      end


c       ===========================================================================
	REAL*4 FUNCTION FUNL(x)
c       ===========================================================================
c       For non-zero cosmological constant
c       ---------------------------------------------------------------------------
	real*4 x, omega0, omegainv
	common /cosm/ omega0
	omegainv = 1. / omega0
	funl = 1. / sqrt(((x ** 3.) + omegainv) - 1.)
	return
	end

c       ===========================================================================
	FUNCTION LARGO(A)
c       ===========================================================================
c	Returns significant length of string a
c       ---------------------------------------------------------------------------
	character*(*) a
	largo=0
	do i=1,len(a)
	if (a(i:i).ne.' ') largo=i
	enddo
	return
	end


c       ===========================================================================
	REAL FUNCTION LINEAR (X0,X,Y,N,I0)
c       ===========================================================================
c	Interpolates linearly the function y(x) at x=x0
c       ---------------------------------------------------------------------------
	REAL X(N),Y(N)
	data ierr/0/
	IF (I0) 21,20,20
 20	I=IPLACE(X0,X,MAX0(1,I0),N)
	GOTO 22
 21	I=-I0
 22	I0=I
	IF (I.EQ.0) GOTO 2
	IF (X0.EQ.X(I).OR.X(I).EQ.X(I+1)) GOTO 1
	LINEAR=Y(I) + (Y(I+1)-Y(I))*(X0-X(I))/(X(I+1)-X(I))
	RETURN
 1	LINEAR=Y(I)
	RETURN
 2	if (ierr.eq.0) then
	   ierr=1
	   write (6,3) X0
	   write (6,*) 'Error reported only once. It may occur more than once.'
	endif
 3	FORMAT (' --- LINEAR:  X0 = ',1PE10.3,' is outside X range ---')
	LINEAR=0.
	RETURN
	END


c       ===========================================================================
	SUBROUTINE SORT2(N,RA,RB)
c       ===========================================================================
c	Sorts array ra(n)
c       ---------------------------------------------------------------------------      
	DIMENSION RA(N),RB(N)
	L=N/2+1
	IR=N
 10	CONTINUE
        IF(L.GT.1)THEN
	   L=L-1
	   RRA=RA(L)
	   RRB=RB(L)
        ELSE
	   RRA=RA(IR)
	   RRB=RB(IR)
	   RA(IR)=RA(1)
	   RB(IR)=RB(1)
	   IR=IR-1
	   IF(IR.EQ.1)THEN
	      RA(1)=RRA
	      RB(1)=RRB
	      RETURN
	   ENDIF
        ENDIF
        I=L
        J=L+L
 20	IF(J.LE.IR)THEN
	   IF(J.LT.IR)THEN
	      IF(RA(J).LT.RA(J+1))J=J+1
	   ENDIF
	   IF(RRA.LT.RA(J))THEN
	      RA(I)=RA(J)
	      RB(I)=RB(J)
	      I=J
	      J=J+J
	   ELSE
	      J=IR+1
	   ENDIF
	   GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
	GO TO 10
	END

c       ===========================================================================
	SUBROUTINE MIDPNT(func,a,b,s,n)
c       ===========================================================================
	INTEGER n
	REAL a,b,s,func
	EXTERNAL func
	INTEGER it,j
	REAL ddel,del,sum,tnm,x
	if (n.eq.1) then
	   s=(b-a)*func(0.5*(a+b))
	else
	   it=3**(n-2)
	   tnm=it
	   del=(b-a)/(3.*tnm)
	   ddel=del+del
	   x=a+0.5*del
	   sum=0.
	   do 11 j=1,it
	      sum=sum+func(x)
	      x=x+ddel
	      sum=sum+func(x)
	      x=x+del
 11	   continue
	   s=(s+(b-a)*sum/tnm)/3.
	endif
	return
	END


c       ===========================================================================
	FUNCTION TRAPZ1 (X,Y,N)
c       ===========================================================================
c	Integrates function y(x)
c       n : number of points
c       ---------------------------------------------------------------------------  
	REAL X(N),Y(N)
	TRAPZ1=0.
	IF (N.LE.1) RETURN
	DO 1 J=2,N
 1	   TRAPZ1= TRAPZ1 + ABS(X(J)-X(J-1))*(Y(J)+Y(J-1))/2.
	   RETURN
 2	   TRAPZ1=Y(1)*X(1)/2.
	   RETURN
	   END


c       ===========================================================================
	INTEGER FUNCTION IPLACE(X0,XX,N1,N2)
c       ===========================================================================
c	Finds L such that XX(L).lt.X0.lt.XX(L+1)
c       n1 : start point
c       n2 : end point
c       ---------------------------------------------------------------------------  
	DIMENSION XX(N2)
	X(I)=XX(N1+I-1)
	N=N2-N1+1
	I=1
	IF (X(1).GT.X(2)) I=2
	L=0
 1	L=L+1
	IF (L.EQ.N) GOTO 4
	GOTO (2,3),I
 2	IF (X(L).LE.X0.AND.X0.LT.X(L+1)) GOTO 5
	GOTO 1
 3	IF (X(L).GE.X0.AND.X0.GT.X(L+1)) GOTO 5
	GOTO 1
 4	IF (X0.NE.X(N)) GOTO 6
 5	IPLACE=L+N1-1
	RETURN
 6	IPLACE=0
	RETURN
	END
