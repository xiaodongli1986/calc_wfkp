
program main

use ap_cosmo_funs
implicit none

	character(len=char_len) :: inputfilename, tmpstr1, tmpstr2
	integer :: numarg, i, numz, ndat, ibin, numbin, col_weight, col_weight2, i1,i2, skiprow, maxcol
	logical :: hasweight, hasweight2, appendrlt, outputrlt, appendreplace
	real(dl) :: omegam, w, Area, zmin,zmax,zstep, binned_sumweight(10000), nbar(10000), zedges(10001), &
		redges(10001), Vol(10000), wfkp(10000), tmp(10000)
	real(dl) :: x, y, z, Red, Weight, z1,z2,wfkp1,wfkp2,nbar1,nbar2,nowwfkp,nownbar
	
	numarg = iargc()
	if(numarg.le.4) then
		print *, 'Usgae: ./calc_wfkp -om omegam -w w -inputfilename inputfilename -zmin zmin -zmax zmax -zstep zstep -Area Area -hasweight hasweight -hasweight2 hasweight2 -col_weight col_weight -col_weight col_weight2 -appendrlt appendrlt -outputrlt outputrlt -appendreplace appendreplace -skiprow skiprow'
		print *, 'if hasweight=True, use weight read in from col_weight'
		print *, 'if hasweight2=True, use multiplication of weights readin from col_weight, col_weight2...'
		print *, 'if outputrlt, output a file with nbar, wfkp'
		print *, 'if appendrlt, output a file with nbar, wfkp appended to the original file (will replace if appendreplace=True)'
		print *, 'skip first skiprow lines '
		print *, 'fmt of inputfile must be: x/y/z OR x/y/z/weight'
		stop
	endif

	omegam = 0.26
	w = -1.0
	Area = 41252.96124941928d0
	skiprow = 0
	appendrlt = .false.
	hasweight = .false.
	col_weight = 4
	hasweight2 = .false.
	col_weight2 = 4
	outputrlt = .false.
	appendreplace = .false.
	
	do i = 1, numarg
		if(mod(i,2).eq.0) cycle
		call getarg(i,tmpstr1)
		call getarg(i+1,tmpstr2)
		if(trim(adjustl(tmpstr1)).eq.'-om') then
			read(tmpstr2,*) omegam
		elseif(trim(adjustl(tmpstr1)).eq.'-w') then
			read(tmpstr2,*) w
		elseif(trim(adjustl(tmpstr1)).eq.'-inputfilename') then
			read(tmpstr2,'(A)') inputfilename
		elseif(trim(adjustl(tmpstr1)).eq.'-zmin') then
			read(tmpstr2,*) zmin
		elseif(trim(adjustl(tmpstr1)).eq.'-zmax') then
			read(tmpstr2,*) zmax
		elseif(trim(adjustl(tmpstr1)).eq.'-zstep') then
			read(tmpstr2,*) zstep
		elseif(trim(adjustl(tmpstr1)).eq.'-Area') then
			read(tmpstr2,*) Area
		elseif(trim(adjustl(tmpstr1)).eq.'-hasweight') then
			read(tmpstr2,*) hasweight
		elseif(trim(adjustl(tmpstr1)).eq.'-col_weight') then
			read(tmpstr2,*) col_weight
		elseif(trim(adjustl(tmpstr1)).eq.'-hasweight2') then
			read(tmpstr2,*) hasweight2
		elseif(trim(adjustl(tmpstr1)).eq.'-col_weight2') then
			read(tmpstr2,*) col_weight2
		elseif(trim(adjustl(tmpstr1)).eq.'-appendrlt') then
			read(tmpstr2,*) appendrlt
		elseif(trim(adjustl(tmpstr1)).eq.'-outputrlt') then
			read(tmpstr2,*) outputrlt
		elseif(trim(adjustl(tmpstr1)).eq.'-appendreplace') then
			read(tmpstr2,*) appendreplace
		elseif(trim(adjustl(tmpstr1)).eq.'-skiprow') then
			read(tmpstr2,*) skiprow
		else
			print *, 'Unkown argument: ', trim(adjustl(tmpstr1))
			print *, 'Usgae: ./calc_wfkp -om omegam -w w -inputfilename inputfilename -zmin zmin -zmax zmax -zstep zstep -Area Area -hasweight hasweight -hasweight2 hasweight2 -col_weight col_weight -col_weight col_weight2 -appendrlt appendrlt -outputrlt outputrlt -appendreplace appendreplace -skiprow skiprow'
			stop
		endif
	enddo

	! print information
	if(.true.) then
		print *, '#####################################'
		print *, ' Settings of calc_nbar_wfkp'
		print *, '  inputfile:  ', trim(adjustl(inputfilename))
		print *, '  omega, w = ', real(omegam), real(w)
		print *, '  zmin, zmax, zstep = ', real(zmin), real(zmax), real(zstep)
		print *, '  Area = ', Area
		print *, '  hasweight = ', hasweight
		print *, '  col_weight =', col_weight
		print *, '  hasweight2 = ', hasweight2
		print *, '  col_weight2 =', col_weight2
		print *, '  outputrlt  =', outputrlt
		print *, '  appendrlt  =', appendrlt
		if(appendrlt) then
			print *, '  Replance previous file  =', appendreplace
		endif
		print *, '#####################################'
	endif
	
	if(col_weight.lt.4) then
		print *, 'ERROR! xyz 1-3; col_weight shall be at least 4;'
		stop
	endif

	! preparation for cosmological calculation
  	gb_omegam = omegam; gb_w=w; gb_h = 0.72_dl
  	call cosmo_funs_init(.true.)
	call de_calc_comovr()
	
	! edges of redshifts
	numbin = ceiling((zmax-zmin)/zstep)
	if(numbin>10000) then
		print *, 'Too many bins (>10000)! increase the size of arrays!: ', numbin
		stop
	endif
	do i = 1, numbin
		zedges(i) = zmin+zstep*(i-1)
	enddo
	zedges(numbin+1) = zmax
!	write(*,'(A,<numbin+1>(f8.4))') 'Edges of redshifts: ', zedges(1:numbin+1)

	! edges of r
	do i = 1, numbin+1
		redges(i) = comov_r(dble(zedges(i)))
	enddo
!	write(*,'(A,<numbin+1>(e15.7))') 'Edges of Distances: ', redges(1:numbin+1)

	! Volumn
	do i = 1, numbin
		Vol(i) = vol_fun(redges(i), redges(i+1), Area)
	enddo
!	write(*,'(A,<numbin>(e15.7))') 'Volumes: ', Vol(1:numbin)
	
	! processing files
	call count_line_number(inputfilename,ndat)
	ndat = ndat-skiprow
	print *, '#-data = ', ndat, '...'
	print *, 'calculating sum_weight...'
	binned_sumweight = 0.0
	open(unit=1,file=inputfilename)
	do i = 1, skiprow
		read(1,*) tmpstr1
	enddo
	
	if(.not.hasweight2) then
		maxcol = col_weight
	else
		maxcol = max(col_weight,col_weight2)
	endif
		
	maxcol = max(col_weight,col_weight2)
	do i = 1, ndat
		if(hasweight) then
			read(1,*) tmp(1:maxcol)
			x = tmp(1)
			y = tmp(2)
			z = tmp(3)
			Weight = tmp(col_weight)
			if(hasweight2) then
				Weight = Weight*tmp(col_weight2)
			endif
		else
			read(1,*) x,y,z
			Weight = 1.0_dl
		endif
		Red = de_zfromintpl(dble(sqrt(x**2.0+y**2.0+z**2.0)))
		if(Red>zmin .and. Red<zmax) then
			ibin = ceiling((Red - zmin) / zstep)
			binned_sumweight(ibin) = binned_sumweight(ibin)+Weight
		endif
	enddo
	close(1)
	
	print *, sum(binned_sumweight)
	
	! nbar, wfkp
	do i = 1, numbin
		nbar(i) = binned_sumweight(i) / Vol(i)
		wfkp(i) = 1.0/(1.0+nbar(i)*20000.0d0)
	enddo
	
	! print & output information
	tmpstr1 = trim(adjustl(inputfilename))//'.nbar-wfkp-info'
	write(*,'(A,A)'), ' *** information of nbar/wfkp wrote to: ', trim(adjustl(tmpstr1))
	open(unit=2,file=tmpstr1)
!	write(*,'(A)')  '#zcen     zlow      zhigh     nbar         wfkp      shell_vol     num-gal'
	write(2,'(A)')  '#zcen     zlow      zhigh     nbar         wfkp      shell_vol     num-gal'
	do i = 1, numbin
!		write(*,'(3(f9.6),4e14.6)'), (zedges(i)+zedges(i+1))/2.0,zedges(i),zedges(i+1),nbar(i),wfkp(i),vol(i),binned_sumweight(i)
		write(2,'(3(f9.6),4e13.6)'), (zedges(i)+zedges(i+1))/2.0,zedges(i),zedges(i+1),nbar(i),wfkp(i),vol(i),binned_sumweight(i)
	enddo
	close(2)

	! output information of nbar, wkfp for each point...	
	if(.not.outputrlt .and. .not.appendrlt) stop
	
	open(unit=1,file=inputfilename)
	if(outputrlt) then
		tmpstr1 = trim(adjustl(inputfilename))//'.nbar-wfkp'
		open(unit=2,file=tmpstr1)
		write(*,'(A,A)'), ' *** nbar/wfkp of every point wrote to: ', trim(adjustl(tmpstr1))
	endif
	if(appendrlt) then
		tmpstr2 = trim(adjustl(inputfilename))//'.nbar-wfkp-appended'
		open(unit=3,file=tmpstr2)
		write(*,'(A,A)'), ' *** nbar/wfkp of every point append to: ', trim(adjustl(tmpstr2))
	endif
	do i = 1, skiprow
		read(1,'(A)') tmpstr1
		if(appendrlt) write(3,'(A)') trim(adjustl(tmpstr1))
	enddo
	do i = 1, ndat
		read(1,'(A)') tmpstr1
		read(tmpstr1,*) x,y,z
		Red = de_zfromintpl(dble(sqrt(x**2.0+y**2.0+z**2.0)))
		if(Red>zmin .and. Red<zmax) then
			ibin = ceiling((Red - zmin) / zstep)
			if(Red>(zedges(ibin)+zedges(ibin+1))*0.5) then
				i1 = ibin
			else
				i1 = ibin-1
			endif
			i1 = max(i1,1); i1 = min(i1, numbin-1); i2=i1+1
			z1 = (zedges(i1)+zedges(i1+1))*0.5; 
			z2 = (zedges(i2)+zedges(i2+1))*0.5
			nbar1 = nbar(i1); nbar2 = nbar(i2)
			wfkp1 = wfkp(i1); wfkp2 = wfkp(i2)
			nownbar = nbar1 + (nbar2-nbar1)/(z2-z1)*(red-z1)
			nowwfkp = wfkp1 + (wfkp2-wfkp1)/(z2-z1)*(red-z1)
		else
			if(Red<zmin) then
				nownbar = nbar(1); nowwfkp = wfkp(1)
			elseif(Red>zmax) then
				nownbar = nbar(numbin); nowwfkp = wfkp(numbin)
			else
				print *, 'Error! strange redshift: ', zmin, Red, zmax
				nownbar = nbar(1); nowwfkp = 0.0
			endif
		endif
		if(outputrlt) write(2,*) nownbar, nowwfkp
		if(appendrlt) write(3,'(A,2e15.7)') trim(adjustl(tmpstr1)), nownbar, nowwfkp
	enddo
	close(1); close(2); close(3)

	if(appendreplace .and.appendrlt) then
		print *, 'Replacing previous file by appended one...'
		tmpstr1 = trim(adjustl(inputfilename))//'.nbar-wfkp-appended'	
		call system('mv '//trim(adjustl(tmpstr1))//' '//trim(adjustl(inputfilename)))
	endif
end program main

		
	
	
