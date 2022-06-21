!==========================================================================================
!===============================================================[VERSION]
!@FeifeiTian
!2021.08.03
!===============================================================[NOTES]
!===============================================================[USE]
!./xas_ave.x xas_tt.dat
!===============================================================[DATAFILE]
!xas_tt.dat
!===============================================================[OUTFILE]
!xas_ave.dat
!xas_ave.log
![my.err]
!==========================================================================================
program main
implicit none

integer, parameter               :: q = 8
real (q), parameter              :: pi = 3.14159265359d0, bohr2angstrom=0.529177210903d0, eang2debye=4.799d0
integer                          :: clock_start, clock_end, clock_rate
integer (q)                      :: i, j, k, l, m, ios, temp_int0, input_tot
real (q)                         :: temp_real0, temp_real (3)
character (len=500)              :: temp_char0, temp_char (2), leng=" (65 ('='), ", infile, datafile, outfile='xas_ave.dat', errfile='my.err', logfile='xas_ave.log'
integer, parameter               :: errfile_unit=11, infile_unit=12, datafile_unit=13, outfile_unit=14
character (len=500), allocatable :: input_list(:)

integer (q)                     :: npiece, natom
!integer (q), allocatable        ::
!character (len=100), allocatable ::
real (q), allocatable           :: energy(:,:), intensity(:,:,:), eminn(:), emaxn(:), deltaEn(:), xas(:,:)
real (q)                        :: emin0, emax0, deltaE0
character (len=500)             :: title

open (errfile_unit, file=errfile, status='replace')
call system_clock (clock_start)
call get_command_argument(1,datafile)
!===============================================================[READ DATAFILE]
write (errfile_unit, leng//"'[READ DATAFILE]: "//trim(datafile)//"')")
open(datafile_unit, file = datafile, status = 'old')
temp_int0 = 0
natom = 0
do while (.true.)
        read (datafile_unit, *, iostat = ios) temp_char0
        if (ios < 0) exit
        if (temp_char0(1:1)=='#') then
            natom = natom + 1
            cycle
        endif
        temp_int0 = temp_int0 + 1
enddo
write(errfile_unit, '(A25,I25)') 'natom',natom
write(errfile_unit, '(A25,I25)') 'nline',temp_int0
npiece = temp_int0/natom
write(errfile_unit, '(A25,I25)') 'npiece',npiece

rewind (datafile_unit)
allocate( energy(natom,npiece), intensity(natom,npiece,6))
do i=1, natom
    read (datafile_unit, *, iostat = ios) temp_char0
    if (temp_char0(1:1)/='#' .or. ios < 0) then
        write (*,*) "ERROR!! See my.err!"
        write (errfile_unit, *) "ERROR: datafile not formated!"
        stop
    endif
    do j=1, npiece
        read (datafile_unit, '(a)', iostat = ios) temp_char0
        temp_char0 = adjustl(temp_char0)
        if (temp_char0(1:1)=='#' .or. ios < 0) then
            write (*,*) "ERROR!! See my.err!"
            write (errfile_unit, *) "ERROR: datafile not formated!"
            stop
        endif
        read (temp_char0, *) energy(i,j), intensity(i,j,1:6)
    enddo
enddo

rewind (datafile_unit)
read (datafile_unit, '(a)') title

close(datafile_unit)
!===============================================================[ave]
allocate( deltaEn(natom), eminn(natom), emaxn(natom) )
do i=1, natom
    eminn(i) = minval(energy(i,:))
    emaxn(i) = maxval(energy(i,:))
    deltaEn(i) = (emaxn(i) -eminn(i))/(npiece-1)
enddo
emin0 = maxval(eminn(:))
emax0 = minval(emaxn(:))
deltaE0 = (emax0-emin0)/(npiece-1)
allocate( xas(npiece,7) )
xas(:,:) = 0.0_q
do j=1, npiece
    xas(j,1) = emin0 + (j-1)*deltaE0
    do i=1, natom
        temp_real0 = (xas(j,1) -eminn(i)) /deltaEn(i)+1.0_q
        temp_int0 = floor(temp_real0)
        temp_real0 = temp_real0 - temp_int0
        if (temp_int0<1 .OR. temp_int0>npiece) then
            write (*,*) "ERROR!! See my.err!"
            write (errfile_unit, *) "ERROR: emin0 or emax0 error!"
            stop
        elseif (temp_int0 == npiece) then
            xas(j,2:7) =xas(j,2:7) +intensity(i,temp_int0,1:6)
        else
            xas(j,2:7) =xas(j,2:7) +intensity(i,temp_int0,1:6)*(1.0_q-temp_real0) +intensity(i,temp_int0+1,1:6)*temp_real0
        endif
    enddo
enddo
xas(:,2:7) = xas(:,2:7)/natom
!===============================================================[OUTDATA]
write(errfile_unit, leng//"'[OUTDATA]: "//trim(outfile)//"')")
open(outfile_unit, file = outfile, status = 'replace')
write(outfile_unit, *) trim(title)
do i = 1, npiece
        write(outfile_unit, '(7F20.12)') xas(i,1:7)
enddo
close(outfile_unit)
!===============================================================[FINALIZE
deallocate(energy, intensity, eminn, emaxn, deltaEn, xas)
call system_clock (clock_end, clock_rate)
write (errfile_unit, ' (A25, F)') 'Elapsed time (s)', (clock_end-clock_start)/real(clock_rate)
close(errfile_unit)
call rename(errfile, logfile)
end program main
