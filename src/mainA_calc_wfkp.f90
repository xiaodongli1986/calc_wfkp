!This code is used for DR12v1 CMASS North

program ap_main

use mpi
use ap_cosmo_funs

	implicit none

! Settings
	!Input file...
        logical :: apply_radecrange = .true. ! If true then the data will be trimed by the above ranges 

! Variables
	real(dl) :: x,y,z,vx,vy,vz,mass,zcosmo,zobs,rat, ra,dec,r, ra2,dec2,r2, RMatrix(3,3,100),Point(3),NewPoint(3), VPoint(3),NewVPoint(3),  tmpx,tmpy,tmpz
        real(dl) :: thetax,thetay,thetaz, Matrix(3,3), Matrix2(3,3)
	character(len=char_len) :: inifile, inputfilename, inputfile, tmpstr1, tmpstr2, tmpstr3, outputfiles(100), cmdstr
	integer ::  i,j,nowi,k, outputunits(100),nowpatch, numgaltot, numgals(100), i_mock
! Variables used to check distances between points        
        real(dl) :: XYZ(3,10000000), P1s(3,10), P2s(3,10), P3s(3,10), P4s(3,10), maxdiff = 0.0, d1, d2
        integer, parameter :: i1 = 10000 ! For each path, the galaxies i1 to i1+10 will be checked
        logical :: just_check = .false.

	i = iargc()
	if(i .ne.1) then
		print *, 'Usage: EXE inifile'
		stop
	endif
	call getarg(1, inifile)
	call get_info_from_ini(inifile)
	
	print *
	if(only_do_vlosadd) stop
! Test of rotation matrix
!        thetax = 0.0; thetay = -1.41_dl; thetaz = -1.0*const_pi;
!        thetax = 1.11035138067; thetay = -1.88934879318; thetaz = -1.40741456079
!        Point(1) = 0.2; Point(2) = 0.3; Point(3) = 0.4
!        thetax = -0.12; thetay = -3.09_dl; thetaz = 1.25;
!        thetax = 1.0; thetay = 2.0; thetaz = 3.0
!        thetax = 0.17041976792201244; thetay = 0.00038080272519722277; thetaz = 1.6985941760660432
!        Point = (/0.3791855654012633_dl, -0.7605936250258089d0, 1.2236474563792248d0 /)
!        thetax = -0.8464486782744852; thetay = -2.292039663901131; thetaz = 0.08450248586292837
!        Point = (/2.6094154405931684d0, -1.2716103325698767d0, 0.7803131378334187d0/)
!        call RotateMatrix(thetax,thetay,thetaz, Matrix)
!        print *, 'Matrix: '
!        do i = 1, 3
!                print *, Matrix(i,1:3)
!        enddo
!        print *, 'Point:   ', real(Point)
!        call RotateM(Point,Matrix,NewPoint)
!        print *, 'NewPoint:', real(NewPoint)
!        call nizhen(Matrix,Matrix2,3)
!        call RotateM(NewPoint,Matrix2,Point)
!        print *, 'Rotate back:', real(Point)
!        stop

! Names/Units of output files 	
    do i_mock = imock0, imock1
	inputfilename = get_inputfilename(i_mock) 
        inputfile = trim(adjustl(inputdir))//trim(adjustl(inputfilename))
        write(*,'(A)') ' inputfile: '//trim(adjustl(inputfile))
        cmdstr = 'mkdir -p '//trim(adjustl(outputdir))//trim(adjustl(patchname)); call system(cmdstr)
	do i = 1, num_patch	
		write(tmpstr1,*) i
		outputfiles(i) = trim(adjustl(outputdir))//trim(adjustl(patchname))//'/'//trim(adjustl(inputfilename))//'.patch'//trim(adjustl(tmpstr1))
		outputunits(i) = 4028202+i
		open(file=outputfiles(i),unit=outputunits(i),action='write')
                write(*,'(A)') ' outputfile: '//trim(adjustl(outputfiles(i)))
	enddo
        
        if(just_check) then
                goto 300
        endif

! Rotate matrix
	tmpx =  const_pi
	do i = 1, num_patch
		call RotateMatrix(thetaxs(i), thetays(i), thetazs(i), RMatrix(1:3,1:3,i))
		print *, 'Rotation Matrix of patch (comparing it with python code for check)', i
		do j = 1,3
			print *, '   ', real(RMatrix(j,1:3,i))
		enddo
	enddo
        do i = 1, num_patch
                Matrix = RMatrix(1:3,1:3,i)
                call nizhen(Matrix, RMatrix(1:3,1:3,i), 3)
        enddo
        ! Matrix one
!        RMatrix(1,1,1) = 0.97589745_dl;  RMatrix(3,3,1) = RMatrix(1,1,1)
!        RMatrix(2,2,1) = 1.0;
!        RMatrix(1,3,1) = 0.21822962_dl;  RMatrix(3,1,1) = -RMatrix(1,3,1)
        ! Matrix two
!        RMatrix(1,1,2) = -9.75897449e-01;  RMatrix(3,3,2) = -RMatrix(1,1,1)
!        RMatrix(1,2,2) = 1.19509022e-16;   RMatrix(2,1,2) = -1.22460635e-16;
!        RMatrix(1,3,2) = 0.21822962_dl;  RMatrix(3,1,2) = -RMatrix(1,3,1)
        print *
	do i = 1, num_patch
		print *, 'Inversed Rotation Matrix of patch (comparing it with python code for check)', i
		do j = 1,3
			print *, '   ', real(RMatrix(j,1:3,i))
		enddo
	enddo
	
! Read in the data, select the patches, rotate them, and output them
	numgaltot = 0; numgals(1:num_patch)=0
	open(file=inputfile,unit=1,action='read')
	do while(.true.)
		read(1,*,end=100) x,y,z,vx,vy,vz,mass,zcosmo,zobs,rat
		numgaltot = numgaltot + 1
!                call random_number(tmpx); if(tmpx<0.99) cycle !!!! ONLY SCAN 1% SAMPLE!!!
                call xyz_to_radecr(x,y,z,ra,dec,r)
		Point(1)=x;Point(2)=y;Point(3)=z
		
		do nowpatch = 1, num_patch
			call RotateM(Point, RMatrix(1:3,1:3,nowpatch), NewPoint)
                        call xyz_to_radecr(NewPoint(1),NewPoint(2),NewPoint(3),ra2,dec2,r2)
                        if(.not.apply_radecrange.or.(apply_radecrange.and.ra2>ramin.and.dec2>decmin.and.ra2<ramax.and.dec2<decmax)) then
                        	VPoint(1)=vx;VPoint(2)=vy;VPoint(3)=vz ! We add the rotating of velocity...
                        	call RotateM(VPoint, RMatrix(1:3,1:3,nowpatch), NewVPoint)
			        write(outputunits(nowpatch),'(6(1x,f13.6), e14.7, 3(1x,f10.7), i10 )') &
				        NewPoint(1:3),NewVPoint(1:3),mass,zcosmo,zobs,rat, numgaltot
                                numgals(nowpatch)=numgals(nowpatch)+1
                        endif
		enddo
		
		cycle !!! Will just rotate all galaxies, and select those within our ra/dec range. Will not apply ra/dec selection before rotation (too complicate; time consumed not too much different)
		
		!!! Very careful: Determine the patch, and output!
                !!! Settings: modify the ra/dec ranges of the patch!
		! Patch 1 
		if((85.0<ra.and.ra<285.0) .and. (-5.0<dec.and.dec<90.0)) then
			nowpatch = 1 
			call RotateM(Point, RMatrix(1:3,1:3,nowpatch), NewPoint)
                        call xyz_to_radecr(NewPoint(1),NewPoint(2),NewPoint(3),ra2,dec2,r2)
                        if(.not.apply_radecrange.or.(apply_radecrange.and.ra2>ramin.and.dec2>decmin.and.ra2<ramax.and.dec2<decmax)) then
			        write(outputunits(nowpatch),'(6(1x,f13.6), e14.7, 3(1x,f10.7), i10 )') &
				        NewPoint(1:3),vx,vy,vz,mass,zcosmo,zobs,rat, numgaltot
                                numgals(nowpatch)=numgals(nowpatch)+1
                        endif
		endif
		! Patch 2
		if((ra<100.0 .or. 265.0<ra) .and. (-5.0<dec.and.dec<90.0)) then
			nowpatch = 2
			call RotateM(Point, RMatrix(1:3,1:3,nowpatch), NewPoint)
                        call xyz_to_radecr(NewPoint(1),NewPoint(2),NewPoint(3),ra2,dec2,r2)
                        if(.not.apply_radecrange.or.(apply_radecrange.and.ra2>ramin.and.dec2>decmin.and.ra2<ramax.and.dec2<decmax)) then
        			write(outputunits(nowpatch),'(6(1x,f13.6), e14.7, 3(1x,f10.7), i10 )') &
	        			NewPoint(1:3),vx,vy,vz,mass,zcosmo,zobs,rat, numgaltot
                                numgals(nowpatch)=numgals(nowpatch)+1
                        endif
		endif
		! Patch 3
		if((85.0<ra.and.ra<285.0) .and. (-90.0<dec.and.dec<5.0)) then
			nowpatch = 3
			call RotateM(Point, RMatrix(1:3,1:3,nowpatch), NewPoint)
                        call xyz_to_radecr(NewPoint(1),NewPoint(2),NewPoint(3),ra2,dec2,r2)
                        if(.not.apply_radecrange.or.(apply_radecrange.and.ra2>ramin.and.dec2>decmin.and.ra2<ramax.and.dec2<decmax)) then
        			write(outputunits(nowpatch),'(6(1x,f13.6), e14.7, 3(1x,f10.7), i10 )') &
	        			NewPoint(1:3),vx,vy,vz,mass,zcosmo,zobs,rat, numgaltot
                                numgals(nowpatch)=numgals(nowpatch)+1
                        endif
		endif
		! Patch 4
		if((ra<100.0 .or. 265.0<ra) .and. (-90.0<dec.and.dec<5.0)) then
			nowpatch = 4
			call RotateM(Point, RMatrix(1:3,1:3,nowpatch), NewPoint)
                        call xyz_to_radecr(NewPoint(1),NewPoint(2),NewPoint(3),ra2,dec2,r2)
                        if(.not.apply_radecrange.or.(apply_radecrange.and.ra2>ramin.and.dec2>decmin.and.ra2<ramax.and.dec2<decmax)) then
        			write(outputunits(nowpatch),'(6(1x,f13.6), e14.7, 3(1x,f10.7), i10 )') &
	        			NewPoint(1:3),vx,vy,vz,mass,zcosmo,zobs,rat, numgaltot
                                numgals(nowpatch)=numgals(nowpatch)+1
                        endif
		endif
		cycle
100		exit		

	enddo
	
	print *, 'Tot gal-#: ', numgaltot
	do i = 1, num_patch
		print *, '  gal-# in patch ', i,':', numgals(i)
	enddo
!        * Blue region: ra: (85,285); dec: (-5,90)
!        * Red  region: ra: (0,100)&(265,360); dec: (-5,90)
!        * Yellow region: ra: (85,285); dec: (-90,5)
!        * Green region:  ra: (0,100)&(265,360); dec: (-90,5)

!(2). Cut them down from the data;

!(3). Rotate them to the original ra, dec ranges;;;
!        Angles used to rotate away from the original region are:
!            [[0.0, 0.22, 0.0], [0.0, 0.22, np.pi*1.0], [0.0, -1.41, 0.0], [0.0, -1.41, np.pi]]
!        So, to rotate back, we shall use:
!            [[0.0, -0.22, 0.0], [0.0, -0.22, -np.pi*1.0], [0.0, 1.41, 0.0], [0.0, 1.41, -np.pi]]
		
		
! Closing output files	

300     continue
        close(1)
	do i = 1, num_patch
		close(outputunits(i))
	enddo

	cycle

        print *, 'Consistency Check (distances between galaxies shall not change)...'
	numgaltot = 0
	open(file=inputfile,unit=1,action='read')
	do while(.true.)
		read(1,*,end=101) x,y,z
		numgaltot = numgaltot + 1
                XYZ(1:3,numgaltot) = (/x,y,z/)
                cycle
101             exit                
        enddo
        close(1)
	do nowpatch = 1, num_patch	
                print *, 'Checking patch ', nowpatch, '...'
                j = 0; k = 0
                open(file=outputfiles(nowpatch),unit=outputunits(nowpatch),action='read')
                do while (.true.)
                        read(outputunits(nowpatch),*,end=102) x,y,z,vx,vy,vz,mass,zcosmo,zobs,rat, nowi
                        j = j+1
                        if(j>i1) then
                                k = j-i1;
                                if(k>10) exit
                                P1s(1,k) = x; P1s(2,k)=y; P1s(3,k)=z
                                P2s(1:3,k) = XYZ(1:3,nowi)
                        endif
                        cycle
102                     exit                        
                enddo
                maxdiff = 0.0;
                do i = 1, 9
                        do j = i+1, 10
                                d1 = distance(P1s(1:3,i),P1s(1:3,j))
                                d2 = distance(P2s(1:3,i),P2s(1:3,j))
!                                print *, '      Distances:', i,j, real(d1), real(d2)
                                maxdiff = max(maxdiff,abs(d1-d2))
                                d1 = P1s(1,i)**2.0 + P1s(2,i)**2.0 + P1s(3,i)**2.0
                                d2 = P2s(1,i)**2.0 + P2s(2,i)**2.0 + P2s(3,i)**2.0
!                                maxdiff = max(maxdiff,abs(sqrt(d1)-sqrt(d2)))
                        enddo
                enddo
                print *, 'Max-diff in distances before/after rotation: ', maxdiff
                close(outputunits(nowpatch))
	enddo
    enddo ! do loop of i_mock
end program ap_main
