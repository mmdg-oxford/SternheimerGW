SUBROUTINE extract (filplot,plot_num)
  !-----------------------------------------------------------------------
  !
  !    This subroutine reads the data for the output file produced by pw.x
  !    extracts and calculates the desired quantity (rho, V, ...)
  !    writes it to a file for further processing or plotting
  !
  !    DESCRIPTION of the INPUT: see file Doc/INPUT_PP
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : bg
  USE ener,      ONLY : ef
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau
  USE gvect
  USE klist,     ONLY : two_fermi_energies
  USE vlocal,    ONLY : strf
  USE io_files,  ONLY : tmp_dir, prefix, trimcheck
  USE io_global, ONLY : ionode, ionode_id
  USE mp_global,     ONLY : nproc, nproc_pool, nproc_file, nproc_pool_file
  USE control_flags, ONLY : twfcollect
  USE noncollin_module, ONLY : i_cons
  USE paw_variables, ONLY : okpaw
  USE mp,        ONLY : mp_bcast

  IMPLICIT NONE
  CHARACTER(len=256), INTENT(out) :: filplot
  INTEGER, INTENT(out) :: plot_num

  INTEGER :: kpoint, kband, spin_component, ios
  LOGICAL :: lsign, needwf

  REAL(DP) :: emin, emax, sample_bias, z, dz, epsilon
  ! directory for temporary files
  CHARACTER(len=256) :: outdir

  NAMELIST / inputpp / outdir, prefix, plot_num, sample_bias, &
       spin_component, z, dz, emin, emax, kpoint, kband, &
       filplot, lsign, epsilon

  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  filplot = 'tmp.pp' 
  plot_num = -1
  spin_component = 0
  sample_bias = 0.01d0
  z = 1.d0
  dz = 0.05d0
  lsign=.FALSE.
  emin = -999.0d0
  emax = +999.0d0
  epsilon=1.d0
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck ( outdir )
     !
  END IF
  !
  call mp_bcast (ios, ionode_id)
  !
  IF ( ios /= 0) CALL errore ('postproc', 'reading inputpp namelist', ABS(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( plot_num, ionode_id )
  CALL mp_bcast( sample_bias, ionode_id )
  CALL mp_bcast( spin_component, ionode_id )
  CALL mp_bcast( z, ionode_id )
  CALL mp_bcast( dz, ionode_id )
  CALL mp_bcast( emin, ionode_id )
  CALL mp_bcast( emax, ionode_id )
  CALL mp_bcast( kband, ionode_id )
  CALL mp_bcast( kpoint, ionode_id )
  CALL mp_bcast( filplot, ionode_id )
  CALL mp_bcast( lsign, ionode_id )
  CALL mp_bcast( epsilon, ionode_id )
  !
  ! no task specified: do nothing and return
  !
  IF (plot_num == -1) return
  !
  IF (plot_num < 0 .OR. plot_num > 18) CALL errore ('postproc', &
          'Wrong plot_num', ABS (plot_num) )

  IF (plot_num == 7 .OR. plot_num == 13 .OR. plot_num==18) THEN
     IF  (spin_component < 0 .OR. spin_component > 3) CALL errore &
          ('postproc', 'wrong spin_component', 1)
  ELSE IF (plot_num == 10) THEN
     IF  (spin_component < 0 .OR. spin_component > 2) CALL errore &
          ('postproc', 'wrong spin_component', 2)
  ELSE
     IF (spin_component < 0 ) CALL errore &
         ('postproc', 'wrong spin_component', 3)
  END IF
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )

  needwf=(plot_num==3).or.(plot_num==4).or.(plot_num==5).or.(plot_num==7).or. &
         (plot_num==8).or.(plot_num==10)
  IF (nproc /= nproc_file .and. .not. twfcollect .and. needwf)  &
     CALL errore('postproc',&
     'pw.x run with a different number of processors. Use wf_collect=.true.',1)

  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect .and. needwf)  &
     CALL errore('postproc',&
     'pw.x run with a different number of pools. Use wf_collect=.true.',1)

  IF ( ( two_fermi_energies .or. i_cons /= 0) .AND. &
       ( plot_num==3 .or. plot_num==4 .or. plot_num==5 ) ) &
     CALL errore('postproc',&
     'Post-processing with constrained magnetization is not available yet',1)
!   IF (okpaw) CALL errore('postproc', &
!              'post-processing paw routines not yet tested',1)

  CALL openfil_pp ( )
  CALL struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, nr3, &
       strf, eigts1, eigts2, eigts3)
!  CALL init_us_1 ( )
  !
  ! The following line sets emax to its default value if not set
  ! It is done here because Ef must be read from file
  !
  IF (emax == +999.0d0) emax = ef
  IF (plot_num == 10) THEN
     emin = emin / 13.6058d0
     emax = emax / 13.6058d0
  END IF
  !
  !
  !   Now do whatever you want
  !
  !HL
  !CALL punch_plot (filplot, plot_num, sample_bias, z, dz, &
  !     emin, emax, kpoint, kband, spin_component, lsign, epsilon)
  !
END SUBROUTINE extract

