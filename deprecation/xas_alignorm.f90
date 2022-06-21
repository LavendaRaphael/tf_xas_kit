!==========================================================================================
!---------------------------------------------------------------[VERSION]
! 2021.09.08
! @FeifeiTian
!---------------------------------------------------------------[NOTES]
!---------------------------------------------------------------[USE]
! ./xas_alignorm.x xas_alignorm.in
!---------------------------------------------------------------[xas_alignorm.in]
! datafile      "*.dat"
! datacolumn    2
! e_align       535.0d0             # eV
! area          10.d0
! e_begin       530.d0               # eV
! e_end         545.d0               # eV
! peak_tolera   0.1                 # %
! peak_index    1
!---------------------------------------------------------------[DATAFILE]
! *.dat
!---------------------------------------------------------------[OUTFILE]
! xas_alignorm.dat
! xas_alignorm.log
! [my.err]
!==========================================================================================
program main
implicit none
!----------------------------------------------------------------[CommonVariable]
integer, parameter              :: DP = 8
real (DP), parameter            :: pi = 3.14159265359d0, bohr2angstrom=0.529177210903d0, eang2debye=4.799d0
integer                         :: clock_start, clock_end, clock_rate
integer (DP)                    :: i, j, k, l, m, ios, temp_int0, input_tot
real (DP)                       :: temp_real0, temp_real (3)
character (len=500)             :: temp_char0, temp_char (2)
character (len=500)             :: leng=" (65 ('='), ", infile, datafile
character (len=500)             :: outfile="xas_alignorm.dat", errfile='my.err', logfile="xas_alignorm.log"
integer, parameter              :: errfile_unit=11, infile_unit=12, datafile_unit=13, outfile_unit=14
character (len=500), allocatable:: input_list(:)
!---------------------------------------------------------------[InputVariable]
real (DP)                       :: e_align, area, e_begin, e_end, peak_tolera
integer (DP)                    :: datacolumn, peak_index
!---------------------------------------------------------------[InternalVariable] 
integer (DP)                    :: ei_tot, i_begin, i_end, columnmax
real (DP), allocatable          :: xas(:,:)
real (DP)                       :: xas_sum, e_onset
character (len=500)             :: title

open (errfile_unit, file=errfile, status='replace')
call system_clock (clock_start)
!---------------------------------------------------------------[READ INFILE]
call get_command_argument(1,infile)
write (errfile_unit, leng//"'[READ INFILE]: "//trim(infile)//"')")
open(infile_unit, file = infile, status = 'old')
input_tot = 8
allocate(input_list(input_tot))
input_list = (/'datafile','datacolumn','e_align','area','e_begin','e_end','peak_tolera','peak_index'/)
do i = 1, input_tot
do while (.true.)
    read (infile_unit, *, iostat = ios) temp_char0
    if (ios < 0) then
        write (errfile_unit, *) "ERROR: parameter '"//trim(input_list(i))//"' not found in "//trim(infile)//"!"
        stop
    endif
    if (temp_char0 /= input_list(i)) cycle
    backspace(infile_unit)
    exit
enddo
select case (i)
case (1)
    read (infile_unit, *) temp_char(1:2)
    write (errfile_unit, "(2A25)") trim(temp_char(1)), trim(temp_char(2))
    read (temp_char(2), *) datafile
case (2)
    read (infile_unit, *) temp_char(1:2)
    write (errfile_unit, "(2A25)") trim(temp_char(1)), trim(temp_char(2))
    read (temp_char(2), *) datacolumn
case (3)
    read (infile_unit, *) temp_char(1:2)
    write (errfile_unit, "(2A25)") trim(temp_char(1)), trim(temp_char(2))
    read (temp_char(2), *) e_align
case (4)
    read (infile_unit, *) temp_char(1:2)
    write (errfile_unit, "(2A25)") trim(temp_char(1)), trim(temp_char(2))
    read (temp_char(2), *) area
case (5)
    read (infile_unit, *) temp_char(1:2)
    write (errfile_unit, "(2A25)") trim(temp_char(1)), trim(temp_char(2))
    read (temp_char(2), *) e_begin
case (6)
    read (infile_unit, *) temp_char(1:2)
    write (errfile_unit, "(2A25)") trim(temp_char(1)), trim(temp_char(2))
    read (temp_char(2), *) e_end
case (7)
    read (infile_unit, *) temp_char(1:2)
    write (errfile_unit, "(2A25)") trim(temp_char(1)), trim(temp_char(2))
    read (temp_char(2), *) peak_tolera
case (8)
    read (infile_unit, *) temp_char(1:2)
    write (errfile_unit, "(2A25)") trim(temp_char(1)), trim(temp_char(2))
    read (temp_char(2), *) peak_index
endselect
rewind (infile_unit)
enddo
close (infile_unit)
!---------------------------------------------------------------[READ DATAFILE]
write (errfile_unit, leng//"'[READ DATAFILE]: "//trim(datafile)//"')")
open(datafile_unit, file = datafile, status = 'old')
ei_tot = 0
call sub_nlines(datafile_unit,ei_tot)
write(errfile_unit, '(A25,I25)') 'ei_tot',ei_tot

rewind (datafile_unit)
title=''
do while (.true.)
    read (datafile_unit, '(a)', iostat = ios) temp_char0
    if (ios < 0) exit
    temp_char0 = adjustl(temp_char0)
    if (temp_char0(1:1) == '#') then
        title = temp_char0
    else
        call sub_splitline(temp_char0, columnmax)
        exit
    endif
enddo
write(errfile_unit, '(A25,I25)') 'columnmax',columnmax

if (columnmax < datacolumn) then
    write (*,*) "ERROR !!! See 'my.err' for details."
    write (errfile_unit, *) " columnmax < datacolumn !!!"
    stop
endif

rewind (datafile_unit)
allocate( xas(ei_tot,columnmax))
i = 0
do while (.true.)
    read (datafile_unit, '(a)', iostat = ios) temp_char0
    if (ios < 0) exit
    temp_char0 = adjustl(temp_char0)
    if (temp_char0(1:1) == '#') cycle
    i = i + 1
    read(temp_char0, *) xas(i, 1:columnmax)
enddo

close(datafile_unit)
!---------------------------------------------------------------[ALIGN TO PRE-EDGE]
temp_int0=0
write (errfile_unit,'(A25)') "FIND PEAK"
do i=2, ei_tot-1
    if (xas(i,datacolumn) > xas(i-1,datacolumn) .AND. xas(i,datacolumn) > xas(i+1,datacolumn)) then
        temp_real0=xas(i,datacolumn)/maxval(xas(:,datacolumn))
        if (temp_real0 >peak_tolera) then
            temp_int0 = temp_int0 + 1
            write (errfile_unit,'(A25,F25.10,F25.10)') '',xas(i,1),temp_real0
            if (temp_int0==peak_index) e_onset = xas(i,1)
        endif
    endif
enddo
if (temp_int0<peak_index) then
    write (*,*) "ERROR !!! See 'my.err' for details."
    write (errfile_unit, *) "Find peak_num ",temp_int0,"< peak_index ",peak_index,"!!!"
    stop
endif
write(errfile_unit, '(A25,F25.10)') 'e_onset', e_onset
write(errfile_unit, '(A25,F25.10)') 'align_delta', e_align - e_onset
xas(:,1) = xas(:,1) - e_onset + e_align
!---------------------------------------------------------------[NORMALIZATION]
! e_begin = e_align + e_begin
! e_end = e_align + e_end

do i = 1, ei_tot
    if (xas(i,1) >= e_begin) exit
enddo
i_begin = i
write(errfile_unit, '(A25,I25)') 'i_begin', i_begin
xas_sum = 0.d0
do i= i_begin, ei_tot-1
    if (xas(i+1,1) > e_end) exit
    xas_sum = xas_sum + (xas(i+1,1)-xas(i,1))*(xas(i+1,datacolumn)+xas(i,datacolumn))/2.d0
enddo
i_end = i
write(errfile_unit, '(A25,I25)') 'i_end', i_end
write(errfile_unit, '(A25,F25.10)') 'xas_sum',xas_sum
write(errfile_unit, '(A25,F25.10)') 'norm_scale', area/xas_sum
! do i = i_begin, i_end
!     xas(i,2) = xas(i,2) / xas_sum * area
! enddo
xas(:, 2:columnmax) = xas(:, 2:columnmax) / xas_sum * area
!---------------------------------------------------------------[OUTDATA]
write(errfile_unit, leng//"'[OUTDATA]: "//trim(outfile)//"')")
open(outfile_unit, file = outfile, status = 'replace')
! do i = i_begin, i_end
!     write(outfile_unit, '(2F25.10)') xas(i,1:2)
! enddo
write(outfile_unit, *) trim(title)
write (temp_char0, *) columnmax
do i = 1, ei_tot
    write(outfile_unit, '('//temp_char0//'F20.12)') xas(i,1:columnmax)
enddo

close(outfile_unit)
!---------------------------------------------------------------[FINALIZE]
call system_clock (clock_end, clock_rate)
write (errfile_unit, ' (A25, F)') 'Elapsed time (s)', (clock_end-clock_start)/real(clock_rate)
close(errfile_unit)
call rename(errfile, logfile)
end program main
!================================================================================[读取文件行数]
subroutine sub_nlines(file_unit,nlines)
    implicit none
    integer, parameter                  :: q = 8
    
    integer, intent(in)                 :: file_unit
    
    integer(q), intent(out)             :: nlines
    
    integer(q)                          :: ios
    character(100)                      :: temp

    nlines = 0
    do while (.true.)
        read (file_unit, *, iostat = ios) temp
        if (ios < 0) exit
        if (temp(1:1) == '#') cycle    ! 忽略注释行
        nlines = nlines + 1
    enddo
    rewind(file_unit)
    
end subroutine sub_nlines
!=================================================================================[读取行参数个数]
subroutine sub_splitline(line, n)
    implicit none
    integer, parameter                  :: q = 8
    
    character(*), intent(in)            :: line
    
    integer(q), intent(out)             :: n
    
    integer(q)                          :: i, ios
    character(100)                      :: buff(10000)
    
    do n=1, 10000
        read(line,*,iostat=ios) buff(1:n)
        if (ios<0) exit
    enddo
    n=n-1

end subroutine sub_splitline   
