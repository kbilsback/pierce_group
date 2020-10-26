! Offline RRTMG
!0;95;c JKodros

program rrtmg_main

  USE MCICA_SUBCOL_GEN_SW, ONLY : MCICA_SUBCOL_SW
  USE PARKIND,        ONLY : IM=>KIND_IM, RB=>KIND_RB
  use parrrsw, only : nbndsw, ngptsw, naerec
  USE RRTMG_SW_RAD, ONLY: RRTMG_SW

  implicit none
  include '/usr/local/include/netcdf.inc'

  ! Met file
  character (len = *), parameter :: IN_FILE_NAME = '/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects_025x03125/AIE/setup_output/china_boff_w1415.nc'
  character (len=*), parameter :: OUT_FILE_NAME = '/pierce-scratch/kbilsback/geos_chem/projects/china_rescoal/radiative_effects_025x03125/AIE/RRTMG_output/china_boff_w1415.nc'

  ! CLOUD INPUTS
!  integer, parameter :: MNS=12, NCOL=3312, Z = 47, WVS=nbndsw, HRS=8
  integer, parameter :: MNS=12, NCOL=3312, Z = 30, WVS=nbndsw, HRS=8
  REAL(KIND=RB) :: PCENTER(NCOL, Z, MNS)
  REAL(KIND=RB) :: CLDFR(NCOL, Z, MNS)
  REAL(KIND=RB) :: CICEWP(NCOL, Z, MNS)
  REAL(KIND=RB) :: CLIQWP(NCOL, Z, MNS) 
  REAL(KIND=RB) :: REICE(NCOL, Z, MNS)
  REAL(KIND=RB) :: RELIQ(NCOL, Z, MNS)
  REAL(KIND=RB) :: MONTHS(MNS)

  ! RRTMG SW INPUTS
  REAL(KIND=RB) :: PLEV(NCOL, Z+1, MNS)
  REAL(KIND=RB) :: TSFC(NCOL, MNS)
  REAL(KIND=RB) :: TLAY(NCOL, Z, MNS)
  REAL(KIND=RB) :: TLEV(NCOL, Z+1, MNS)
  REAL(KIND=RB) :: SUNCOS(NCOL, HRS, MNS)
  !REAL(KIND=RB) :: SUNCOS(NCOL, MNS)
  !REAL(KIND=RB) :: ALBEDO(NCOL, MNS)

  ! netcdf input ids
  integer :: ncid, varid1, varid2, varid3, varid4, varid5, varid6
  integer :: varid7, varid8, varid9, varid10, varid11
  integer :: varid12, varid13, varid14, varid15, varid16, varid17, varid18
  integer :: varid19
  integer :: varid17a, varid17b, varid17c, varid17d
  integer :: retval, rec

  ! rad input options
  INTEGER(KIND=IM) :: ONECOL, LLPAR, ICLDMCL, SEEDSW, IRNG
  INTEGER(KIND=IM) :: INFLGSW, ICEFLGSW, LIQFLGSW
  ! Loops
  integer :: T, H, IL, w, zl

  ! UNUSED, YET NECESSARY, INPUTS
  REAL(KIND=RB) :: TAUCLD_SW(WVS, NCOL, Z)
  REAL(KIND=RB) :: SSACLD(WVS, NCOL, Z)
  REAL(KIND=RB) :: ASMCLD(WVS, NCOL, Z)
  REAL(KIND=RB) :: FSFCLD(WVS, NCOL, Z)

  ! MCICA SW OUTPUTS
  REAL(KIND=RB)         :: CLDFMCL_SW0(NGPTSW,NCOL,Z)
  REAL(KIND=RB)         :: CIWPMCL_SW0(NGPTSW,NCOL,Z)
  REAL(KIND=RB)         :: CLWPMCL_SW0(NGPTSW,NCOL,Z)
  REAL(KIND=RB)         :: TAUCMCL_SW0(NGPTSW,NCOL,Z)
  REAL(KIND=RB)         :: SSACMCL0(NGPTSW,NCOL,Z)
  REAL(KIND=RB)         :: ASMCMCL0(NGPTSW,NCOL,Z)
  REAL(KIND=RB)         :: FSFCMCL0(NGPTSW,NCOL,Z)
  REAL(KIND=RB)         :: RELQMCL0(NCOL,Z)
  REAL(KIND=RB)         :: REICMCL0(NCOL,Z)
     
  !REAL*8, ALLOCATABLE  :: CLDFMCL_SW_OUT(:,:,:,:)
  !REAL*8, ALLOCATABLE  :: CIWPMCL_SW_OUT(:,:,:,:)
  !REAL*8, ALLOCATABLE  :: CLWPMCL_SW_OUT(:,:,:,:)
  !REAL*8, ALLOCATABLE  :: TAUCMCL_SW_OUT(:,:,:,:)
  !REAL*8, ALLOCATABLE  :: SSACMCL_OUT(:,:,:,:)
  !REAL*8, ALLOCATABLE  :: ASMCMCL_OUT(:,:,:,:)
  !REAL*8, ALLOCATABLE  :: FSFCMCL_OUT(:,:,:,:)
  !REAL*8, ALLOCATABLE  :: REICMCL_OUT(:,:,:) 
  !REAL*8, ALLOCATABLE  :: RELQMCL_OUT(:,:,:) 

  ! RRTMG INPUT OPTIONS
  ! SW SOLAR VARIABLES
  REAL(KIND=RB)  :: ADJES=1.0              ! FLUX ADJUSTMENT FOR EARTH/SUN DISTANCE
  REAL(KIND=RB)  :: SCON=1368.22           ! SOLAR CONSTANT (W/M2)
  
  integer :: DOY    ! DOY

  ! CHEMICAL SPECIES
  REAL(KIND=RB) :: H2OVMR(NCOL, Z, MNS)
  REAL(KIND=RB) :: O3VMR(NCOL, Z)
  REAL(KIND=RB) :: CO2VMR(NCOL, Z)
  REAL(KIND=RB) :: CH4VMR(NCOL, Z)
  REAL(KIND=RB) :: N2OVMR(NCOL, Z)
  REAL(KIND=RB) :: O2VMR(NCOL, Z)

  !ALBEDO INPUTS
  REAL(KIND=RB) :: ALBDIRVIS(NCOL, MNS)
  REAL(KIND=RB) :: ALBDIFVIS(NCOL, MNS)
  REAL(KIND=RB) :: ALBDIRNIR(NCOL, MNS)
  REAL(KIND=RB) :: ALBDIFNIR(NCOL, MNS)

  ! AEROSOL INPUTS
  REAL(KIND=RB) :: TAUAER_SW(NCOL, Z, NBNDSW)
  REAL(KIND=RB) :: SSAAER(NCOL, Z, NBNDSW)
  REAL(KIND=RB) :: ASMAER(NCOL, Z, NBNDSW)
  REAL(KIND=RB) :: ECAER(1,Z,NAEREC)

  ! RADSW OUTPUTS
  REAL(KIND=RB) :: SWUFLX(NCOL,Z+1)  ! TOTAL SKY SHORTWAVE UPWARD FLUX (W/M2)
  REAL(KIND=RB) :: SWDFLX(NCOL,Z+1)  ! TOTAL SKY SHORTWAVE DOWNWARD FLUX (W/M2)
  REAL(KIND=RB) :: SWHR(NCOL,Z)      ! TOTAL SKY SHORTWAVE RADIATIVE HEATING RATE (K/D)
  REAL(KIND=RB) :: SWUFLXC(NCOL,Z+1) ! CLEAR SKY SHORTWAVE UPWARD FLUX (W/M2)
  REAL(KIND=RB) :: SWDFLXC(NCOL,Z+1) ! CLEAR SKY SHORTWAVE DOWNWARD FLUX (W/M2)
  REAL(KIND=RB) :: SWHRC(NCOL,Z)     ! CLEAR SKY SHORTWAVE RADIATIVE HEATING RATE (K/D)

  ! RADSW OUTPUT STORE
  REAL(KIND=RB) :: SW_UFLHRS(HRS, NCOL, Z+1)
  REAL(KIND=RB) :: SW_DFLHRS(HRS, NCOL, Z+1)
  REAL(KIND=RB) :: SW_UFLCHRS(HRS, NCOL, Z+1)
  REAL(KIND=RB) :: SW_DFLCHRS(HRS, NCOL, Z+1)

  REAL(KIND=RB) :: SW_UFLUX(MNS, NCOL, Z+1)
  REAL(KIND=RB) :: SW_DFLUX(MNS, NCOL, Z+1)
  REAL(KIND=RB) :: SW_UFLUXC(MNS, NCOL, Z+1)
  REAL(KIND=RB) :: SW_DFLUXC(MNS, NCOL, Z+1)
  
  
  
  ! OUT netCDF variables


  integer, parameter :: NDIMS=3
  integer, parameter :: NT=MNS, NC=NCOL, NZ=Z+1

  integer :: ncid_out, varid1_out, varid2_out, varid3_out, varid4_out, dimids(NDIMS)
  integer :: NT_dimid, NC_dimid, NZ_dimid, retvalo

  ONECOL = 1
  !LLPAR = 47
  LLPAR = Z
  ICLDMCL = 2
  SEEDSW = NGPTSW + 1
  !SEEDSW = 1
  IRNG = 1

  INFLGSW = 2
  ICEFLGSW = 2
  LIQFLGSW = 1
  
  !open file
  retval = nf_open(IN_FILE_NAME, NF_NOWRITE, ncid)

  retval = nf_inq_varid(ncid, "pcenter", varid1)
  retval = nf_inq_varid(ncid, "cldfr", varid2)
  retval = nf_inq_varid(ncid, "cicewp", varid3)
  retval = nf_inq_varid(ncid, "cliqwp", varid4)
  retval = nf_inq_varid(ncid, "reice", varid5)
  retval = nf_inq_varid(ncid, "reliq", varid6)
  retval = nf_inq_varid(ncid, "month", varid11)
  retval = nf_inq_varid(ncid, "plev", varid12)
  retval = nf_inq_varid(ncid, "tsfc", varid13)
  retval = nf_inq_varid(ncid, "tlay", varid14)
  retval = nf_inq_varid(ncid, "tlev", varid15)
  retval = nf_inq_varid(ncid, "suncos", varid16)
  retval = nf_inq_varid(ncid, "albdrvis", varid17a)
  retval = nf_inq_varid(ncid, "albdfvis", varid17b)
  retval = nf_inq_varid(ncid, "albdrnir", varid17c)
  retval = nf_inq_varid(ncid, "albdfnir", varid17d)
  retval = nf_inq_varid(ncid, "o3vmr", varid18)
  retval = nf_inq_varid(ncid, "h2ovmr", varid19)

  ! Read data
  retval = nf_get_var(ncid, varid1, PCENTER)  
  retval = nf_get_var(ncid, varid2, CLDFR)
  retval = nf_get_var(ncid, varid3, CICEWP)
  retval = nf_get_var(ncid, varid4, CLIQWP)
  retval = nf_get_var(ncid, varid5, REICE)
  retval = nf_get_var(ncid, varid6, RELIQ)

  retval = nf_get_var(ncid, varid11, MONTHS)
  
  retval = nf_get_var(ncid, varid12, PLEV)
  retval = nf_get_var(ncid, varid13, TSFC)
  retval = nf_get_var(ncid, varid14, TLAY)
  retval = nf_get_var(ncid, varid15, TLEV)
  retval = nf_get_var(ncid, varid16, SUNCOS)

  retval = nf_get_var(ncid, varid17a, ALBDIRVIS)
  retval = nf_get_var(ncid, varid17b, ALBDIFVIS)
  retval = nf_get_var(ncid, varid17c, ALBDIRNIR)
  retval = nf_get_var(ncid, varid17d, ALBDIFNIR)
  
  retval = nf_get_var(ncid, varid18, O3VMR)
  retval = nf_get_var(ncid, varid19, H2OVMR)

  retval = nf_close(ncid)
  
  print *, "WEEEELLLPPP I READ IN SOME NETCDF FILES"
  
  !! CHEMCIAL SPECIES
  !H2OVMR(:,:) = 100E-4
  !O3VMR(:,:) = 2.0E-6
  CO2VMR(:,:) = 3.90E-4
  CH4VMR(:,:) = 1800.0E-9
  N2OVMR(:,:) = 315.0E-9
  O2VMR(:,:) = 0.209

  !H2OVMR(:,:,:) = 0.0
  !O3VMR(:,:) = 0.0
  !CO2VMR(:,:) = 0.0
  !CH4VMR(:,:) = 0.0
  !N2OVMR(:,:) = 0.0
  !O2VMR(:,:) = 0.0

  !! ALBEDO
  !ALBDIRVIS(:, :) = ALBEDO(:,:)
  !ALBDIFVIS(:,:) = ALBEDO(:,:)
  !ALBDIFVIS(:,:) = 0.01
  !ALBDIRNIR(:, :) = ALBEDO(:,:)
  !ALBDIFNIR(:, :) = ALBEDO(:,:)
  !ALBDIFVIS(:,:) = 0.01
  

  !! Zero aerosol inputs for indirect effect
  TAUAER_SW(:,:,:) = 0.0
  SSAAER(:,:,:) = 1.0
  ASMAER(:,:,:) = 0.0
  ECAER(:,:,:) = 0.0

  !! ZERO UNUSED CLOUD INPUTS
  TAUCLD_SW(:,:,:) = 0.
  SSACLD(:,:,:) = 0.
  ASMCLD(:,:,:) = 0.
  FSFCLD(:,:,:) = 0.

  !ALLOCATE( CLDFMCL_SW_OUT( MNS, NGPTSW, NCOL, Z ))
  !ALLOCATE( CIWPMCL_SW_OUT( MNS, NGPTSW, NCOL, Z ))
  !ALLOCATE( CLWPMCL_SW_OUT( MNS, NGPTSW, NCOL, Z ))
  !ALLOCATE( TAUCMCL_SW_OUT( MNS, NGPTSW, NCOL, Z ))
  !ALLOCATE( SSACMCL_OUT( MNS, NGPTSW, NCOL, Z ))
  !ALLOCATE( ASMCMCL_OUT( MNS, NGPTSW, NCOL, Z ))
  !ALLOCATE( FSFCMCL_OUT( MNS, NGPTSW, NCOL, Z ))
  !ALLOCATE( REICMCL_OUT( MNS, NCOL, Z ))
  !ALLOCATE( RELQMCL_OUT( MNS, NCOL, Z ))

 do T = 1, MNS
 !do T = 12, 12
     print*, 'MONTH: ', T
     !do IL = 1, NCOL
     !do IL = 1, 1
        
     !DOY = 3*30 + T*15
     !DOY = 3*30 + T*30 -15  
     DOY = T*30 -15  
        
     CALL MCICA_SUBCOL_SW(NCOL, LLPAR, ICLDMCL&
             ,SEEDSW, IRNG,PCENTER(:,:,T)&
             ,CLDFR(:,:,T), CICEWP(:,:,T), CLIQWP(:,:,T)&

             ,REICE(:,:,T), RELIQ(:,:,T)&
             ,TAUCLD_SW(:,:,:),SSACLD(:,:,:)&
             ,ASMCLD(:,:,:),FSFCLD(:,:,:)&
             ,CLDFMCL_SW0, CIWPMCL_SW0, CLWPMCL_SW0&
             ,REICMCL0, RELQMCL0, TAUCMCL_SW0, SSACMCL0&
             ,ASMCMCL0, FSFCMCL0)

     !print*, "I made it here"
        
     !CLDFMCL_SW_OUT(T,:,:,:)   = CLDFMCL_SW0(:,:,:)
     !CIWPMCL_SW_OUT(T,:,:,:)   = CIWPMCL_SW0(:,:,:)
     !CLWPMCL_SW_OUT(T,:,:,:)   = CLWPMCL_SW0(:,:,:)
     !TAUCMCL_SW_OUT(T,:,:,:)   = TAUCMCL_SW0(:,:,:)
     !SSACMCL_OUT(T,:,:,:)    = SSACMCL0(:,:,:)
     !ASMCMCL_OUT(T,:,:,:)    = ASMCMCL0(:,:,:)
     !FSFCMCL_OUT(T,:,:,:)    = FSFCMCL0(:,:,:)
     !REICMCL_OUT(T,:,:)  = REICMCL0(:,:)
     !RELQMCL_OUT(T,:,:)  = RELQMCL0(:,:)
      
        !print*, "CLDFMCL_SW0: ", CLDFMCL_SW0(:,1,5)
        !print*, "CLWPMCL_SW0: ",CLWPMCL_SW0(:,1,5)
        !print*, "RELQMCL0: ", RELQMCL0(1,:)   
        !print*, 'NCOL: ', NCOL
        !print*, 'LLPAR: ',LLPAR
        !print*, 'ICLDMCL: ',ICLDMCL
        !print*, 'PCENTER: ',PCENTER(1,1,T)
        !print*, 'PLEV: ',PLEV(1,1,T)
        !print*, 'TLAY: ',TLAY(1,1,T)
        !print*, 'TLEV: ',TLEV(1,1,T)
        !print*, 'TSFC: ',TSFC(1,T)
        !print*, 'H2OVMR: ',H2OVMR(1,1,T)
        !print*, 'O3VMR: ',O3VMR(1,1)
        !print*, 'CO2VMR: ',CO2VMR(IL,:)
        !print*, 'CH4VMR: ',CH4VMR(IL,:)
        !print*, 'N2OVMR: ',N2OVMR(IL,:)
        !print*, 'O2VMR: ',O2VMR(IL,:)
        !print*, 'ALBDIRVIS: ',ALBDIRVIS(IL,T)
        !print*, 'ALBDIFVIS: ',ALBDIFVIS(IL,T)
        !print*, 'ALBDIRNIR: ',ALBDIRNIR(IL,T)
        !print*, 'ALBDIFNIR: ',ALBDIFNIR(IL,T)
        !print*, 'SUNCOS: ',SUNCOS(246,7,T)
        !print*, 'ADJES: ',ADJES
        !print*, 'DOY: ',DOY
        !print*, 'SCON: ',SCON
        !print*, 'INFLGSW: ',INFLGSW
        !print*, 'ICEFLGSW: ',ICEFLGSW
        !print*, 'LIQFLGSW: ',LIQFLGSW
     do H=1,HRS  
      print*, 'HOUR: ', H
   
        CALL RRTMG_SW(NCOL, LLPAR, ICLDMCL&
             ,PCENTER(:,:,T), PLEV(:,:,T), TLAY(:,:,T), TLEV(:,:,T), TSFC(:,T)&
             ,H2OVMR(:,:,T), O3VMR(:,:), CO2VMR(:,:), CH4VMR(:,:), N2OVMR(:,:), O2VMR(:,:)&
             ,ALBDIRVIS(:,T), ALBDIFVIS(:,T), ALBDIRNIR(:,T), ALBDIFNIR(:,T)&
             ,SUNCOS(:,H,T), ADJES, DOY, SCON&
             ,INFLGSW, ICEFLGSW, LIQFLGSW&
             ,CLDFMCL_SW0(:,:,:), TAUCMCL_SW0(:,:,:)&
             ,SSACMCL0(:,:,:), ASMCMCL0(:,:,:), FSFCMCL0(:,:,:)&
             ,CIWPMCL_SW0(:,:,:), CLWPMCL_SW0(:,:,:)&
             ,REICMCL0(:,:), RELQMCL0(:,:)&
             ,TAUAER_SW(:,:,:), SSAAER(:,:,:), ASMAER(:,:,:)&
             ,ECAER, SWUFLX, SWDFLX, SWHR, SWUFLXC, SWDFLXC, SWHRC)

        print*, "After RRTMG"
        !print*, 'SWUFLUX', SWUFLX(1,:)
        
        SW_UFLHRS(H,:,:) = SWUFLX(:,:)
        SW_DFLHRS(H,:,:) = SWDFLX(:,:)

        SW_UFLCHRS(H,:,:) = SWUFLXC(:,:)
        SW_DFLCHRS(H,:,:) = SWDFLXC(:,:)
        
        !print*, "SW_UFLUX", SW_UFLUX(T,IL,:)
     end do
     SW_UFLUX(T,:,:) = sum(SW_UFLHRS,1)/HRS
     SW_DFLUX(T,:,:) = sum(SW_DFLHRS,1)/HRS
     SW_UFLUXC(T,:,:) = sum(SW_UFLCHRS,1)/HRS
     SW_DFLUXC(T,:,:) = sum(SW_DFLCHRS,1)/HRS
  end do
print *, 'WELLLLPPP CALLED RRTMG'


!! WIRTE OUTPUT NC FILE
retvalo = nf_create(OUT_FILE_NAME, NF_CLOBBER, ncid_out)

retvalo = nf_def_dim(ncid_out, "month", NT, NT_dimid)
retvalo = nf_def_dim(ncid_out, "column", NC, NC_dimid)
retvalo = nf_def_dim(ncid_out, "layer", NZ, NZ_dimid)


dimids(1) = NT_dimid
dimids(2) = NC_dimid
dimids(3) = NZ_dimid


retvalo = nf_def_var(ncid_out, "SWUFLUX", NF_DOUBLE, NDIMS, dimids, varid1_out)
retvalo = nf_def_var(ncid_out, "SWDFLUX", NF_DOUBLE, NDIMS, dimids, varid2_out)
retvalo = nf_def_var(ncid_out, "SWUFLUC", NF_DOUBLE, NDIMS, dimids, varid3_out)
retvalo = nf_def_var(ncid_out, "SWDFLUC", NF_DOUBLE, NDIMS, dimids, varid4_out)

retvalo = nf_enddef(ncid_out)

retvalo = nf_put_var_double(ncid_out, varid1_out, SW_UFLUX) 
retvalo = nf_put_var_double(ncid_out, varid2_out, SW_DFLUX) 
retvalo = nf_put_var_double(ncid_out, varid3_out, SW_UFLUXC) 
retvalo = nf_put_var_double(ncid_out, varid4_out, SW_DFLUXC) 

retvalo = nf_close(ncid_out)

print*, "WEEEELLLPP WROTE TO NETCDF"


end program rrtmg_main
