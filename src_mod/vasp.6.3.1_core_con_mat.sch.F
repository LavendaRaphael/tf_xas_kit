#include "symbol.inc"
      MODULE CORE_CON_MAT
      USE prec

      ! variables need for core-conduction band matrix elements
      LOGICAL,SAVE :: CH_LSPEC             ! controls whether corehole spectra
                                           ! are calculated or not
      REAL(q), SAVE :: CH_SIGMA            ! controls broadening of corehole
      REAL(q), SAVE :: CH_AMPLIFICATION    ! factor for amplification of ch-spectrum
                                           ! only set this value greater
                                           ! than one as a last resort
      LOGICAL, SAVE :: CH_JDOSCOREHOLE     ! if true matrix elements are set to 1
      INTEGER, SAVE :: CH_NEDOS            ! Number of frequency gridpoints in
                                           ! dielectric function

      CONTAINS

!*******************************************************************
! CALC_CORE_CON_MAT:
! this subroutine reads in neccessary control variables for the
! calculation of matrix elements between the AE orbitals of core and
! conduction states on a selected atom type
!
!*******************************************************************

      SUBROUTINE CORE_CON_MAT_READER(IU5,IU0)
      USE reader_tags

      INTEGER IU5, IU0
      LOGICAL LOPEN

! check if INCAR file exists
      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)

      CH_LSPEC = .FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5,'CH_LSPEC', &
                         CH_LSPEC, IERR, WRITEXMLINCAR)
      CH_JDOSCOREHOLE = .FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5,'CH_JDOSCOREHOLE', &
                         CH_JDOSCOREHOLE, IERR, WRITEXMLINCAR)
      CH_NEDOS = 1000
      CALL PROCESS_INCAR(LOPEN, IU0, IU5,'CH_NEDOS', &
                         CH_NEDOS, IERR, WRITEXMLINCAR)
      CH_SIGMA = 0.3_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5,'CH_SIGMA', &
                         CH_SIGMA, IERR, WRITEXMLINCAR)
      CH_AMPLIFICATION = 1.0_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5,'CH_AMPLIFICATION', &
                         CH_AMPLIFICATION, IERR, WRITEXMLINCAR)

      END SUBROUTINE CORE_CON_MAT_READER

!*******************************************************************
! CALC_CORE_CON_MAT:
! this subroutine calculates the matrix elements between the AE
! orbitals of core and conduction states on a selected atom type
!
!*******************************************************************

      SUBROUTINE CALC_CORE_CON_MAT(T_INFO,IO,KPOINTS,LATT_CUR, &
                                   P,WDES,W,ISYMM)
      USE cl            !N_CL,L_CL,...
      USE base          !IO, symmetry
      USE constant      !PI, EDEPS, HSQDTM, AUTOA
      USE lattice       !LATT_CUR
      USE main_mpi      !COMM_INB
      USE mkpoints      !KPOINTS
      USE mlr_optic     !MAX_ENERGY_UNOCCUPIED
      USE mlrf_main     !OMEGAMAX_OPTIC
      USE poscar        !T_INFO
      USE pseudo        !P
      USE wave          !W, WDES, CPROJ, NI_LOCAL

! internal variables:
      REAL(q) ::                 EMIN
      REAL(q), ALLOCATABLE ::    A(:,:)
      INTEGER ::                 ALPHA
      REAL(q), ALLOCATABLE ::    B(:,:)
      INTEGER ::                 BETA
      INTEGER ::                 CHANNELS
      REAL(q) ::                 DELTAE
      REAL(q) ::                 DELTAFUN
      REAL(q), ALLOCATABLE ::    EIG(:)
      REAL(q) ::                 EMAX
      REAL(q) ::                 ENERGY_DIFF
      COMPLEX(q), ALLOCATABLE :: EPSDD(:,:,:) ! density-density response function
      REAL(q) ::                 EPSILON_HELP
      REAL(q) ::                 FREQUENCY
      INTEGER ::                 I,J
      INTEGER ::                 I_CONBANDS
      INTEGER ::                 I_PROJ
      INTEGER ::                 ISPIN
      INTEGER ::                 ISTEP
      INTEGER ::                 KPOINT
      INTEGER ::                 L_CON
      INTEGER ::                 L_CON_RUN
      INTEGER, ALLOCATABLE ::    LC(:)
      INTEGER ::                 LHELP
      INTEGER ::                 LM_CON
      INTEGER ::                 LM_CORE
      INTEGER ::                 M_CON
      INTEGER ::                 M_CORE
      INTEGER ::                 MAXNL
      COMPLEX(q) ::              MECC1
      COMPLEX(q) ::              MECC2
      INTEGER, ALLOCATABLE ::    N(:)
      REAL(q),ALLOCATABLE ::     NABLA(:,:,:)
      INTEGER ::                 NI
      INTEGER ::                 NIS
      INTEGER ::                 NLCRUN
      INTEGER ::                 NL_CORE
      LOGICAL ::                 NL_NOTSET
      INTEGER ::                 NPRO_BASE
      INTEGER ::                 NT
      INTEGER ::                 OMEGA_FREQU
      REAL(q) ::                 PREFACTOR
      INTEGER ::                 RNMAX
      REAL(q) ::                 STEPFUN
      REAL(q), ALLOCATABLE ::    W_CORE(:,:)

!@FeifeiTian
!---------------------------------------------------->>
      COMPLEX(q), ALLOCATABLE :: EPSDD_TM(:,:,:,:,:)
      REAL(q), ALLOCATABLE ::    FREQUENCY_TM(:,:,:)
      REAL(q), ALLOCATABLE ::    OCC_TM(:,:,:)
      INTEGER ::                 NB_GLOBAL
      INTEGER, ALLOCATABLE ::    NCOUNT_TM(:,:,:)
      INTEGER ::                 IOME_IU1=-1
!---------------------------------------------------->>

! input variables
      TYPE (type_info)           T_INFO           !type_type_info from poscar.F
      TYPE (in_struct)           IO               !in_struct from base.F
      TYPE (kpoints_struct)      KPOINTS          !kpoints_struct from mkpoints.F
      TYPE (latt)                LATT_CUR         !latt from lattice.F
      TYPE (potcar)              P(1:T_INFO%NTYP)   !module pseudo from pseudo.F
!!!      TYPE (symmetry)            SYMM             !module base from base.F
      INTEGER ::                 ISYMM
      TYPE (wavespin)            W                !wavespin from wave.F
      TYPE (wavedes)             WDES             !wavedes from wave.F

! common blocks needed for symmetrization of 3x3 tensor
      COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
         GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL


! other variables
!    COMM_INB - variable from module main_mpi (main_mpi.F)
!    EDEPS, PI, HSQDTM, AUTOA - from module constant (constant.inc)
!    L_CL - variable from module cl (cl_shift.F)
!    N_CL - variable from module cl (cl_shift.F)
!    NT_CL - variable from module cl (cl_shift.F)
!    OMEGAMAX_OPTIC - from module mlrf_main (elinear_response.F)
!    Z_CL - variable from module cl (cl_shift.F)

      EXTERNAL                   DELSTP
      EXTERNAL                   M_sum_z
      EXTERNAL                   NABIJ_RADIAL
      EXTERNAL                   TSYM_CMPLX
      EXTERNAL                   XML_EPSILON_W

! start of profiling routines
      PROFILING_START('calc_core_con_mat')

! prefactor 4*pi^2*e^2*hbar^4/(m_electron^2) where HSQDTM=hbar^2/(2*m_electron)
! and EDEPS=4*pi*e^2, WDES%RSPIN is the spin multiplicity
      PREFACTOR = EDEPS*PI*(HSQDTM*2.0_q)**2/(LATT_CUR%OMEGA)* WDES%RSPIN

! closing handler to be sure before calling it again
      CALL CL_CLEAR_CORE_CONF()

! initialize core wavefunctions
      CALL CL_INIT_CORE_CONF(P(NT_CL),MAXNL)

! set number of grid points
      RNMAX = P(NT_CL)%R%NMAX

! set number of n,l channels for core WF (in principle we need one, but to keep structure we
! calculate all)
      CHANNELS = MAXNL + P(NT_CL)%LMAX

! count the size of LM_CON needed for the allocation of NABLA
! and determine LM_CORE
      LM_CON = 0
      DO L_CON = 1, P(NT_CL)%LMAX
         L_CON_RUN = P(NT_CL)%LPS(L_CON)
         DO M_CON = -L_CON_RUN, L_CON_RUN
            LM_CON = LM_CON + 1
         ENDDO
      ENDDO
      LM_CORE = 2*L_CL+1

! allocate arrays needed in SET_CORE_WF
      ALLOCATE(W_CORE(RNMAX,CHANNELS),A(RNMAX,CHANNELS), B(RNMAX,CHANNELS), &
               N(CHANNELS), LC(CHANNELS),EIG(CHANNELS),NABLA(1:3,1:LM_CORE,1:LM_CON))

! Calculate core wavefunctions
      CALL SET_CORE_WF(P(NT_CL)%RHOAE, P(NT_CL)%POTAE , P(NT_CL)%R, P(NT_CL)%ZCORE, P(NT_CL)%ZVALF_ORIG , &
                       W_CORE, N, LC, EIG, A_=A, B_=B)

      NL_CORE = 0
      NL_NOTSET = .TRUE.
! determine index of desired core state in W_CORE
      chan: DO I=1,CHANNELS
         IF ((N(I)==N_CL).AND.(LC(I)==L_CL)) THEN
! set index for main and angular quantum number of core electron
            NL_CORE = I
            NL_NOTSET  = .FALSE.
            EXIT chan
         ENDIF
      ENDDO chan

! if NL_CORE not found then problem
      IF (NL_NOTSET) THEN
         CALL vtutor%error("NL_CORE could not be set, ERROR, exiting...")
      ENDIF

! if CH_JDOSCOREHOLE then all matrix elements are considered equal 1
! and only joint DOS is calculated
      IF (CH_JDOSCOREHOLE) THEN
         NABLA = 1.0_q
      ELSE
         NABLA = 0.0_q
! calculation of matrix elements between AE core states and conduction states
! start loop over l and m of conduction band wavefunction
! atoms, principal and angular quantum numbers from which to excite are specified in
! NT_CL, N_CL, L_CL
         LM_CORE = 0
         coremloop: DO M_CORE = -L_Cl, L_CL
            LM_CORE = LM_CORE + 1
            LM_CON = 0
            lloop: DO L_CON = 1, P(NT_CL)%LMAX
               L_CON_RUN = P(NT_CL)%LPS(L_CON)
               mloop: DO M_CON = -L_CON_RUN, L_CON_RUN
                  LM_CON = LM_CON + 1
                  CALL NABIJ_RADIAL(P(NT_CL)%R,A(1,NL_CORE),P(NT_CL)%WAE(1,L_CON), &
                                    NABLA(1,LM_CORE,LM_CON),L_CL,L_CON_RUN,M_CORE,M_CON,1._q)
               ENDDO mloop
            ENDDO lloop
         ENDDO coremloop
      ENDIF

! check which core contains atom with core hole
! only this (these) cores should calculate
      NIS=1+SUM(T_INFO%NITYP(1:NT_CL-1))
! determine starting position of current atom (NT_CL) in first index of CPROJ
! should be zero if first atom is chosen
      NPRO_BASE = BASE_INDEX_IN_CPROJ(NIS)

! determine highest conduction state
      EMAX=MAX_ENERGY_UNOCCUPIED(WDES,W)*1.2_q
      IF (OMEGAMAX_OPTIC/=-1) THEN
         EMAX=OMEGAMAX_OPTIC
      ENDIF

! setting EMAX to that of the difference of EMAX and the core state
      EMAX = EMAX+ABS(EIG(NL_CORE))

! determining frequency grid
      EMIN=ABS(EIG(NL_CORE))*0.9_q
      DELTAE=(EMAX-EMIN)/(CH_NEDOS-1)

! allocate neccessary array for dielectric function
      ALLOCATE(EPSDD(CH_NEDOS,3,3))
      EPSDD = (0.0_q,0.0_q)

!@FeifeiTian
!--------------------------------------------------------------->>
      ALLOCATE( EPSDD_TM    ( WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN, 3, 3 ) ) 
      ALLOCATE( FREQUENCY_TM( WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN       ) )
      ALLOCATE( OCC_TM      ( WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN       ) )
      ALLOCATE( NCOUNT_TM   ( WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN       ) )
      EPSDD_TM = (0.0_q,0.0_q)
      FREQUENCY_TM = 0.0_q
      OCC_TM=0.0_q
      NCOUNT_TM=0
!---------------------------------------------------------------<<

! loop over atoms
      atomloop: DO NI=NIS,(NIS+T_INFO%NITYP(NT_CL)-1)
! do the job only on cores that have the atom with the core hole (NT_CL, should
! be only one core!!!)
         IF (NI_LOCAL(NI,WDES%COMM_INB).eq.0) cycle
! spin loop
         spinloop: DO ISPIN = 1, WDES%ISPIN
! frequency loop
            omegaloop: DO OMEGA_FREQU = 1, CH_NEDOS
               FREQUENCY = EMIN+DELTAE*(OMEGA_FREQU-1)
! if Frequency below tolerance set epsilon_imag to zero
               IF (FREQUENCY>=1.E-8_q) THEN
! k-point loop
                  kloop: DO KPOINT = 1, WDES%NKPTS
#ifdef MPI
                     IF (MOD(KPOINT-1,WDES%COMM_KINTER%NCPU).NE.(WDES%COMM_KINTER%NODE_ME-1)) CYCLE
#endif
! loop over conduction bands
                     condloop: DO I_CONBANDS = 1, WDES%NBANDS

!@FeifeiTian
!---------------------------------------------------------------------------------------------------------------------->>
                        NB_GLOBAL=(I_CONBANDS-1)*WDES%COMM_INTER%NCPU+WDES%COMM_INTER%NODE_ME
                        IF (NB_GLOBAL<1 .or. NB_GLOBAL>WDES%NB_TOT) THEN
                            CALL vtutor%error("STOP at 'core_con_mat.F', DEBUG id: TFF02")
                        ENDIF
                        FREQUENCY_TM( NB_GLOBAL,KPOINT,ISPIN) = REAL(W%CELEN(I_CONBANDS,KPOINT,ISPIN),q)-EIG(NL_CORE)
                        OCC_TM      ( NB_GLOBAL,KPOINT,ISPIN) =      W%FERWE(I_CONBANDS,KPOINT,ISPIN)
                        NCOUNT_TM   ( NB_GLOBAL,KPOINT,ISPIN) =    NCOUNT_TM( NB_GLOBAL,KPOINT,ISPIN) + 1
!----------------------------------------------------------------------------------------------------------------------<<

! if band is occupied then cycle
                        ENERGY_DIFF = REAL(W%CELEN(I_CONBANDS,KPOINT,ISPIN),q)-EIG(NL_CORE)-FREQUENCY
! calculate delta function
                        CALL DELSTP(KPOINTS%ISMEAR,ENERGY_DIFF/CH_SIGMA,DELTAFUN,STEPFUN)
                        DELTAFUN = DELTAFUN/CH_SIGMA 
! alternatively use Lorentzian
!                       DELTAFUN =AIMAG(1/(ENERGY_DIFF-(0.0_q,1.0_q)*CH_SIGMA))/PI
                        LM_CORE = 0
! loop over magnetic quantum numbers of core electron
                        lmcoreloop: DO M_CORE = -L_CL,L_CL
                           LM_CORE = LM_CORE + 1
! loop over directions of the dielectric tensor
                           alphaloop: DO ALPHA = 1, 3
                              betaloop: DO BETA = 1, 3
                                 MECC1=0
                                 MECC2=0
! loop over projector quantum numbers (all valence states)
! I_PROJ goes over all energy (or n), l and m channels
                                 projloop: DO I_PROJ = 1, P(NT_CL)%LMMAX
! MECC = wavefunction character times core-con-matrixelements
                                    MECC1 = MECC1 + W%CPROJ(NPRO_BASE+I_PROJ,I_CONBANDS,KPOINT,ISPIN) * &
                                            NABLA(ALPHA,LM_CORE,I_PROJ)
                                    MECC2 = MECC2 + W%CPROJ(NPRO_BASE+I_PROJ,I_CONBANDS,KPOINT,ISPIN) * &
                                            NABLA(BETA,LM_CORE,I_PROJ)
                                 ENDDO projloop
                                  ! we need to multiply by (1-Fermi_weight) to account for metals too
                                 EPSILON_HELP = (1.0_q-W%FERWE(I_CONBANDS,KPOINT,ISPIN)) * &
                                                   PREFACTOR/FREQUENCY**2.0_q * &
                                                   WDES%WTKPT(KPOINT)*DELTAFUN*&
                                                   MECC1*CONJG(MECC2)*CH_AMPLIFICATION
                                 EPSDD(OMEGA_FREQU,ALPHA,BETA)=EPSDD(OMEGA_FREQU,ALPHA,BETA) + &
                                                               CMPLX(0.0_q,EPSILON_HELP,q)

!@FeifeiTian
!------------------------------------------------------------------------------------------------------------------------------->>
                                 EPSILON_HELP = (1.0_q-W%FERWE(I_CONBANDS,KPOINT,ISPIN)) * &
                                                   PREFACTOR/FREQUENCY_TM(NB_GLOBAL,KPOINT,ISPIN) * &
                                                   WDES%WTKPT(KPOINT)*         &
                                                   MECC1*CONJG(MECC2)*CH_AMPLIFICATION
                                 EPSDD_TM(NB_GLOBAL,KPOINT,ISPIN,ALPHA,BETA)= &
                                                                EPSDD_TM(NB_GLOBAL,KPOINT,ISPIN,ALPHA,BETA) + &
                                                                CMPLX(0.0_q,EPSILON_HELP,q)
!-------------------------------------------------------------------------------------------------------------------------------<<

                              ENDDO betaloop
                           ENDDO alphaloop
                        ENDDO lmcoreloop
                     ENDDO condloop
                  ENDDO kloop
               ENDIF
            ENDDO omegaloop
         ENDDO spinloop
      ENDDO atomloop

! allreduce of parallel part (summing up of different contributions)
! reduction is over all nodes, this should be fine
      CALLMPI(  M_sum_z(WDES%COMM,EPSDD(1,1,1),CH_NEDOS*3*3) )

!@FeifeiTian
!----------------------------------------------------------------------------------------------->>
      EPSDD_TM=EPSDD_TM/CH_NEDOS
      CALLMPI( M_sum_z( WDES%COMM,    EPSDD_TM (1,1,1,1,1),(WDES%NB_TOT *WDES%NKPTS *WDES%ISPIN *3*3)))
      CALLMPI( M_sum_d( WDES%COMM,FREQUENCY_TM (1,1,1    ),(WDES%NB_TOT *WDES%NKPTS *WDES%ISPIN     )))
      CALLMPI( M_sum_d( WDES%COMM,      OCC_TM (1,1,1    ),(WDES%NB_TOT *WDES%NKPTS *WDES%ISPIN     )))
      CALLMPI( M_sum_i( WDES%COMM,   NCOUNT_TM (1,1,1    ),(WDES%NB_TOT *WDES%NKPTS *WDES%ISPIN     )))
      IF (ANY(NCOUNT_TM/=CH_NEDOS)) THEN
        CALL vtutor%error("STOP at 'core_con_mat.F', DEBUG id: TFF03. Check POSCAR whether or not only 1 atom was excited!")
      ENDIF
!-----------------------------------------------------------------------------------------------<<

! if symmetry present then symmetrize 3x3 tensor using TSYM
!!      IF (SYMM%ISYM>0) THEN
      IF (ISYMM>0) THEN
         DO OMEGA_FREQU = 1, CH_NEDOS
            CALL TSYM_CMPLX(EPSDD(OMEGA_FREQU,:,:),ISYMOP,NROTK,LATT_CUR%A)
         ENDDO

!@FeifeiTian
!------------------------------------------------------------------------------------------->>
         DO ISPIN = 1, WDES%ISPIN 
         DO KPOINT = 1, WDES%NKPTS                                                         
         DO NB_GLOBAL= 1, WDES%NB_TOT
            CALL TSYM_CMPLX(EPSDD_TM(NB_GLOBAL,KPOINT,ISPIN,:,:),ISYMOP,NROTK,LATT_CUR%A) 
         ENDDO                                                                              
         ENDDO                                                                              
         ENDDO
!-------------------------------------------------------------------------------------------<<

      ENDIF


! deallocate variables needed for SET_CORE_WF
      CALL CL_CLEAR_CORE_CONF()

! for write statement
1000 FORMAT(/ &
     &        '  frequency dependent IMAGINARY DIELECTRIC FUNCTION' &
     &       ,' (independent particle, no local field effects)',A,/&
     &        '     E(ev)  ', &
     &        4X,'X', 9X,'Y', 9X,'Z', 8X,'XY', 8X,'YZ', 8X,'ZX'/ &
     &        '  ----------------------------------------------------', &
     &        '----------------------------------------------------------'/ &
     &        (7F20.12))
1100 FORMAT(/ &
     &        '  frequency dependent      REAL DIELECTRIC FUNCTION' &
     &       ,' (independent particle, no local field effects)',A,/&
     &        '     E(ev)  ', &
     &        6X,'X', 11X,'Y', 11X,'Z', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
     &        '  ----------------------------------------------------', &
     &        '----------------------------------------------------------'/ &
     &        (7F20.12))

! writing to standard out
      IF (IO%IU6>=0) THEN
!         ISTEP=MIN(10,CH_NEDOS/40)
!         IF (KPOINTS%ISMEAR<0) ISTEP = 1
         ISTEP = 1
         WRITE(IO%IU6,1000) " density-density", &
            ((EMIN+DELTAE*(I-1)),AIMAG(EPSDD(I,1,1)),AIMAG(EPSDD(I,2,2)),AIMAG(EPSDD(I,3,3)), &
            AIMAG(EPSDD(I,1,2)),AIMAG(EPSDD(I,2,3)),AIMAG(EPSDD(I,3,1)), I=1,CH_NEDOS,ISTEP)
         WRITE(IO%IU6,1100) " density-density", &
            ((EMIN+DELTAE*(I-1)),REAL(EPSDD(I,1,1),q),REAL(EPSDD(I,2,2),q),REAL(EPSDD(I,3,3),q), &
            REAL(EPSDD(I,1,2),q),REAL(EPSDD(I,2,3),q),REAL(EPSDD(I,3,1),q),I=1,CH_NEDOS,ISTEP)
      ENDIF

! writing to XML file?
      CALL XML_EPSILON_W(DELTAE, REAL(EPSDD,q), AIMAG(EPSDD), CH_NEDOS)
!!!      CALL XML_EPSILON_W(DELTAE, REAL(EPSJJ,q), AIMAG(EPSJJ), CH_NEDOS)

! deallocate neccessary arrays
      DEALLOCATE(EPSDD,W_CORE,A,B,N,LC,EIG,NABLA)

!@FeifeiTian
!------------------------------------------------------------------------------------------->>
20402 FORMAT(/ &
     &        '#  Transition Matrix: KPOINT=',I5,' ISPIN=',I1,/ &
     &        '#     E(ev)',8X,'X', 9X,'Y', 9X,'Z', 8X,'XY', 8X,'YZ', 8X,'ZX', &
     &        8X,'NB_GLOBAL'/ &
     &        (7F20.12,I8))

    IF (IO%IU6>0) THEN
        IOME_IU1=20401
        OPEN(UNIT=IOME_IU1,FILE='MYCARXAS',STATUS='REPLACE')
        DO KPOINT=1, WDES%NKPTS
        DO ISPIN =1, WDES%ISPIN
            WRITE(IOME_IU1,20402) KPOINT,ISPIN, &
                ((FREQUENCY_TM(NB_GLOBAL,KPOINT,ISPIN    ) ,&
                AIMAG(EPSDD_TM(NB_GLOBAL,KPOINT,ISPIN,1,1)),AIMAG(EPSDD_TM(NB_GLOBAL,KPOINT,ISPIN,2,2)),&
                AIMAG(EPSDD_TM(NB_GLOBAL,KPOINT,ISPIN,3,3)),AIMAG(EPSDD_TM(NB_GLOBAL,KPOINT,ISPIN,1,2)),&
                AIMAG(EPSDD_TM(NB_GLOBAL,KPOINT,ISPIN,2,3)),AIMAG(EPSDD_TM(NB_GLOBAL,KPOINT,ISPIN,3,1)),&
                NB_GLOBAL),&
                NB_GLOBAL=1,WDES%NB_TOT)
        ENDDO
        ENDDO
        CLOSE(IOME_IU1)
    ENDIF
    DEALLOCATE(EPSDD_TM,FREQUENCY_TM, OCC_TM, NCOUNT_TM)
!--------------------------------------------------------------------------------------------<<

! End of profiling routines
      PROFILING_STOP('calc_core_con_mat')

      CONTAINS

! local function to get starting index of projectors
         FUNCTION BASE_INDEX_IN_CPROJ(ION_TO_FIND)
         IMPLICIT NONE
         INTEGER BASE_INDEX_IN_CPROJ
         INTEGER ION_TO_FIND     ! index of atom for which we seek the base index in CPROJ
         INTEGER LASTTYP, NPRO, NT, LMMAXC, NI, NIS_LOCAL

         LASTTYP=0
         NIS_LOCAL = 1
         BASE_INDEX_IN_CPROJ=0

         DO NT=1,NT_CL
            LMMAXC=P(NT)%LMMAX
            IF (LMMAXC/=0) THEN
               DO NI=NIS_LOCAL,(T_INFO%NITYP(NT)+NIS_LOCAL-1)
                  IF (NI==ION_TO_FIND) EXIT
! does this element reside on local node
                  IF (NI_LOCAL(NI,WDES%COMM_INB) /=0 ) THEN
                     BASE_INDEX_IN_CPROJ = BASE_INDEX_IN_CPROJ + LMMAXC
                  ENDIF
               ENDDO
            ENDIF
            NIS_LOCAL = NIS_LOCAL+T_INFO%NITYP(NT)
         ENDDO
         END FUNCTION BASE_INDEX_IN_CPROJ

      END SUBROUTINE CALC_CORE_CON_MAT

      END MODULE CORE_CON_MAT

