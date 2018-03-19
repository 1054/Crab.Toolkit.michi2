c       ===========================================================================
	PROGRAM MAKE_ZGRID
c       ---------------------------------------------------------------------------
c       Authors :   E. da Cunha & S. Charlot
c       Latest revision :   Sep. 16th, 2010
c       ---------------------------------------------------------------------------
c       Makes a grid of redshift adapted to current data set
c       OUTPUT: 'zlibs.dat' with the redhifts at which the model libraries 
c       will be computed using the codes 'get_optic_colors' and 'get_infrared_colors'
c
c       version Feb 16, 2011
c       galmax: 5 000 -> 20 000
c       read input file changed (added 'aux' for RA & Dec)
c       ===========================================================================

	implicit none
	character ans
	integer i,io,j,k,igrid,ngrid
	integer nmax,galmax,nzmax,n_obs !galmax: maximum number of galaxies in one input file
	parameter(nmax=200,galmax=20000)  !nmax: maxium number of photometric points/filters
	integer nfilt,filt_id(nmax),fit(nmax),ifilt
	character*29 gal_name(galmax)
	character*29 filt_name(nmax)
	real*8 redshift(galmax),lambda_eff(nmax)
	real*8 flux_obs(galmax,nmax),sigma(galmax,nmax)	
	parameter(nzmax=100000)
	integer fbin(nzmax)
	real*8 zmin,zmax,zgrid(nzmax)
	real*8 zi(nzmax),za(nzmax),dz,aux
c	character filters*80,obs*80 ! dzliu modified 20160720
	character filters*255,obs*255 


c       READ FILTER FILE: e.g. "filters.dat"
	call getenv('USER_FILTERS',filters)
	close(22)
	open(22,file=filters,status='old')
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
	close(22)

c       READ FILE WITH OBSERVATIONS:
	call getenv('USER_OBS',obs)
	close(20)
	open(20,file=obs,status='old')
	do i=1,1
	   read(20,*)
	enddo
	io=0
	n_obs=0
	do while(io.eq.0)
	   n_obs=n_obs+1
	   read(20,*,iostat=io) gal_name(n_obs),redshift(n_obs),
     +	   (flux_obs(n_obs,k),sigma(n_obs,k),k=1,nfilt)
	enddo
	n_obs=n_obs-1
	close(20)

c       OUTPUT FILE:
	   close(21)
	   open(21,file='zlibs.dat',status='unknown')


c       OPTION: redshift grid or use exact z values
c       write (6,'(x,a,$)') 'Build redshift grid [Y/N]?'
c       read (5,'(a)',end=1) ans
        
                ! dzliu TODO
        ans='N' ! dzliu modified, so no need to input things while running batch jobs.
                ! dzliu TODO

	if (ans.eq.'y'.or.ans.eq.'Y') then
c	   write (6,'(x,a,$)') 'Building grid with dz=0.1...'   ! dzliu commented

c       Build grid
c       zmin: minimum grid redshift
c       zmax: maximum grid redshift
c       dz: redshift interval
	   zmin=100.
	   zmax=0.
	   dz=0.1
	   do i=1,n_obs
c	      zmin=min(zmin,redshift(i))
c	      zmax=max(zmax,redshift(i))
	      zmin=redshift(i)-0.1 ! dzliu modified
	      zmax=redshift(i)+0.1 ! dzliu modified
	   enddo
	   write (*,*) 'dzliu zmin  ', zmin ! dzliu modified
	   write (*,*) 'dzliu zmax  ', zmax ! dzliu modified
	   write (*,*) 'dzliu dz    ', dz   ! dzliu modified
	   zmin=0.01*dnint(zmin*100)-dz
	   zmax=0.01*dnint(zmax*100)+dz
	
	   igrid=1
	   zgrid(igrid)=zmin
	   do while (zgrid(igrid).le.zmax)
	      igrid=igrid+1
	      zgrid(igrid)=zgrid(igrid-1)+dz
	      zi(igrid)=zgrid(igrid)-(dz/2)
	      za(igrid)=zgrid(igrid)+(dz/2)
	   enddo
	   ngrid=igrid
	   write (*,*) 'dzliu ngrid ', ngrid   ! dzliu modified

c	No use in building library in a bin with no galaxies
c       Make fbin array
c       if fbin(i)=0 -> no galaxies with redshift at zgrid(i)
c       --> don't make library
c       if fbin(i)=0 -> at least one galaxy with redshift at zgrid(i)
c       --> make library

	   do i=1,ngrid
	      fbin(i)=0 
	   enddo
	   
	   do j=1,n_obs
	      do i=1,ngrid
c		 if (redshift(j).ge.zi(i).and.redshift(j).lt.za(i)) then ! dzliu modified -- to loop all redshifts
		    fbin(i)=fbin(i)+1
c		    write (*,*) 'dzliu fbin(i) ', fbin(i)                 dzliu modified
c		 endif                                                   ! dzliu modified -- to loop all redshifts
	      enddo
	   enddo
	   	   
	   j=0
	   do i=1,ngrid
	      if (fbin(i).gt.0) then
		 j=j+1
		 write(21,211) j,zgrid(i)
		 write(*,*) 'dzliu zgrid(i) ', zgrid(i) ! dzliu modified
	      endif
	   enddo
 211	   format(i6,f8.4)
	   
	   else if (ans.eq.'n'.or.ans.eq.'N') then
	      write (6,'(x,a,$)') 'Building formatted file with galaxy redshifts...'

	      do i=1,n_obs
		 write(21,211) i,redshift(i)
	      enddo
		 
	   endif
	      
 1	   stop
	end
