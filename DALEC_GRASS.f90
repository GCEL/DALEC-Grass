! -----------------------------------------------------------------------------------------------------------
! POOLS:   1.labile 2.foliar 3.root                        ! PARAMETERS:
!         *4.wood   5.litter 6.som                         !
! ------------------------------------------               ! 1.Decomposition rate
! FLUXES:  1.GPP                                           ! 2.Fraction of GPP respired 
! (daily)  2.temprate                                      ! 3.GSI sens for leaf growth
!          3.respiration_auto                              ! 4.NPP belowground allocation parameter
!          4.leaf production                               ! 5.GSI max leaf turnover 
!          5.labile production                             ! 6.TOR roots
!          6.root production                               ! 7.TOR litter
!          7.aboveground production                        ! 8.TOR SOM
!          8.labile consumption -> leaves                  ! 9.Temp factor Q10 (1.2-1.6)
!          9.leaffall factor                               ! 10.Photosynthetic N use efficiency
!          10.leaf litter production                       ! 11.GSI max labile turnover
!         *11.woodlitter production                        ! 12.GSI min temperature threshold (K)
!          12.rootlitter production                        ! 13.GSI max temperature threshold (K)
!          13.respiration het litter                       ! 14.GSI min photoperiod threshold (sec)
!          14.respiration het som                          ! 15.LCA - g.C.leaf_m-2
!          15.litter2som                                   ! 16.C labile (initialization)
!          16.labrelease factor(leaf growth)               ! 17.C foliar (initialization)
!         *17.carbon flux due to fire                      ! 18.C roots  (initialization)
!          18.growing season index                         ! 19.C litter (initialization)
!          19.animal manure C soil input (per time step)   ! 20.GSI max photoperiod threshold (sec)
!          20.animal resp co2 (per time step)              ! 21.GSI min VPD threshold (Pa) 
!          21.animal ch4 (per time step)                   ! 22.GSI max VPD threshold (Pa)
! ------------------------------------------               ! 23.critical GPP for LAI increase (gC.m-2.day-1)
! MET:     1.run day                                       ! 24.GSI senstivity for leaf senescence 
!          2.min T (C)                                     ! 25.GSI - have I just left a growing state (>1)
!          3.max T (C)                                     ! 26.GSI - initial GSI value 
!          4.Radiation (MJ.m-2)                            ! 27.DM min lim for grazing (kg.DM.ha-1)
!          5.CO2 (ppm)                                     ! 28.DM min lim for cutting (kg.DM.ha-1)
!          6.DOY                                           ! 29.leaf-vs-stem allocation factor 
!         *7.lagged precip                                 ! 30.C SOM (initialization)
!          8.cutting/grazing :                             ! 31.DM demand % of animal weight
!            - spatial mode = lai removed (m2.m-2)         ! 32.Post-grazing labile loss
!            - field mode = LSU.ha-1                       ! 33.Post-cut labile loss 
!         *9.burnt area fraction                           
!          10.21-day avg min T (K)        
!          11.21-day avg photoperiod (sec) 
!          12.21-day avg VPD (Pa)         
!         *13.Forest mgmt after clearing
!         *14.Mean T
!
! * not used / empty 
! * this code cannot be used on spatial mode (met driver #8)
! -----------------------------------------------------------------------------------------------------------

module CARBON_MODEL_MOD

implicit none

public :: CARBON_MODEL           &
         ,acm                    &
         ,linear_model_gradient  

! ACM related parameters
double precision, parameter :: pi = 3.1415927
double precision, parameter :: deg_to_rad = pi/180d0

! local variables for GSI phenology model
double precision :: Tfac,Photofac,VPDfac        & ! oC, seconds, Pa
                   ,tmp,gradient                & 
                   ,fol_turn_crit,lab_turn_crit &
                   ,gsi_history(22),just_grown

integer :: gsi_lag_remembered 

double precision, allocatable, dimension(:) :: tmp_x, tmp_m


contains

!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start,finish,deltat,lat,met,pars &
                         ,nodays,nopars,nomet,nopools,nofluxes   &
                         ,LAI,GPP,NEE,POOLS,FLUXES,REMOVED_C,version_code)

    ! The Data Assimilation Linked Ecosystem Carbon - Growing Season
    ! Index - Forest Rotation (DALEC_GSI_FR) model. 
    ! The subroutine calls the Aggregated Canopy Model to simulate GPP and 
    ! partitions between various ecosystem carbon pools. These pools are
    ! subject to turnovers / decompostion resulting in ecosystem phenology and fluxes of CO2

    implicit none

    ! declare input variables
    integer, intent(in) :: start    &
                          ,finish   & 
                          ,nodays   & ! number of days in simulation
                          ,nopars   & ! number of paremeters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nopools  & ! number of model pools
                          ,nofluxes & ! number of model fluxes
                          ,version_code

    double precision, intent(in) :: deltat(nodays)    & ! time step in decimal days
                                   ,lat               & ! site latitude (degrees)
                                   ,met(nomet,nodays) & ! met drivers
                                   ,pars(nopars)        ! number of parameters


    double precision, intent(out) :: LAI(nodays) & ! leaf area index
                                    ,GPP(nodays) & ! Gross primary productivity
                                    ,NEE(nodays)   ! net ecosystem exchange of CO2

    double precision, intent(out) :: POOLS((nodays+1),nopools) ! vector of ecosystem pools
 
    double precision, intent(out) :: FLUXES(nodays,nofluxes) ! vector of ecosystem fluxes

    double precision, intent(out) :: REMOVED_C(2,nodays) ! vector of removed C (grazed,cut)
    
    ! declare general local variables
    double precision :: gpppars(12)        & ! ACM inputs (LAI+met)
                       ,constants(10)        ! parameters for ACM

    integer :: f,n,test,m

    double precision :: foliage_frac_res    & 
                       ,roots_frac_death    &
                       ,labile_loss,foliar_loss       &
                       ,roots_loss                    &
                       ,labile_residue,foliar_residue &
                       ,roots_residue                 & 
                       ,labile_frac_res               &
                       ,tot_abg_exp,fol_frac,lab_frac &
                       ,f_root,NPP                


    integer :: gsi_lag

    gpppars(4) = 2.0  ! g N leaf_m-2
    gpppars(7) = lat
    gpppars(9) = -2.0 ! leafWP-soilWP
    gpppars(10) = 1.0 ! hydraulic resistance
    gpppars(11) = pi

    ! assign acm parameters
    constants(1)=pars(10) 
    constants(2)=0.0156935
    constants(3)=4.22273
    constants(4)=208.868
    constants(5)=0.0453194
    constants(6)=0.37836
    constants(7)=7.19298
    constants(8)=0.011136
    constants(9)=2.1001
    constants(10)=0.789798

    ! post-removal residues and root death | 0:none 1:all
    foliage_frac_res  = 0.05  ! fraction of removed foliage that goes to litter
    labile_frac_res   = 0.05  ! fraction of removed labile that goes to litter
    roots_frac_death  = 0.01  ! fraction of roots that dies and goes to litter

    if (start == 1) then

        ! assigning initial conditions
        POOLS(1,1) = pars(16)
        POOLS(1,2) = pars(17)
        POOLS(1,3) = pars(18)
        POOLS(1,4) = 0 ! no wood in grasslands
        POOLS(1,5) = pars(19)
        POOLS(1,6) = pars(30)

        if (.not.allocated(tmp_x)) then
            allocate(tmp_x(22),tmp_m(nodays))
            do f = 1, 22
               tmp_x(f) = f
            end do
            do n = 1, nodays
              if (sum(deltat(1:n)) < 21) then
                  tmp_m(n) = n-1
              else
                 m = 0 ; test = 0
                 do while (test < 21)
                    m=m+1 ; test = sum(deltat((n-m):n))
                    if (m > (n-1)) then 
                      test = 21 
                    endif
                 end do
                 tmp_m(n) = m
               endif 
            end do
            gsi_lag_remembered = max(2,maxval(nint(tmp_m)))
        end if
        gsi_history = pars(24)-1d0
        just_grown = pars(25)

    endif 

    gsi_lag = gsi_lag_remembered 
    fol_turn_crit=pars(24)-1d0
    lab_turn_crit=pars(3)-1d0

    ! 
    ! Begin looping through each time step
    ! 

    do n = start, finish
      
      ! calculate LAI value
      LAI(n) = POOLS(n,2) / pars(15)

      ! load next met / lai values for ACM
      gpppars(1)=LAI(n)   ! LAI
      gpppars(2)=met(3,n) ! max temp
      gpppars(3)=met(2,n) ! min temp
      gpppars(5)=met(5,n) ! co2
      gpppars(6)=ceiling(met(6,n)-(deltat(n)*0.5)) ! doy
      gpppars(8)=met(4,n) ! radiation

      ! GPP (gC.m-2.day-1)
      if (LAI(n) > 0.) then
         FLUXES(n,1) = acm(gpppars,constants)
      else
         FLUXES(n,1) = 0.
      endif 

      ! temprate (i.e. T modified rate of metabolic activity))
      FLUXES(n,2) = exp(pars(9)*0.5*(met(3,n)+met(2,n)))

      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(n,3) = FLUXES(n,1) * pars(2)
      
      ! NPP 
      NPP = FLUXES(n,1) - FLUXES(n,3)
      
      ! fraction for allocation to roots
      f_root = 1 - exp(-1*pars(4)*LAI(n))
      if (f_root < 0.2) then
         f_root = 0.2
      endif
      if (f_root > 0.8) then 
         f_root = 0.8
      endif
      
      ! allocation to roots 
      FLUXES(n,6) = NPP * f_root
      
      ! C left for aboveground allocation                        
      FLUXES(n,7) = NPP - FLUXES(n,6)               
      
      ! allocation of ABG C to leaf   
      FLUXES(n,4) = FLUXES(n,7) * (1 - (pars(29)*(LAI(n)/6)))
      
      ! allocation of ABG C to labile/stem              
      FLUXES(n,5) = FLUXES(n,7) * (pars(29)*(LAI(n)/6))
      
      ! labile consumption
      FLUXES(n,8) = 0.0

      ! Calculate the Growing Season Index based on Jolly et al. 
      ! doi: 10.1111/j.1365-2486.2005.00930.x doi:10.1029/2010JG001545.
      ! It is the product of 3 limiting factors for temperature, photoperiod and
      ! vapour pressure deficit that grow linearly from 0 to 1 between a calibrated 
      ! min and max value. Photoperiod, VPD and avgTmin are direct input

      ! temperature limitation
      Tfac = ( met(10,n)-pars(12)) / (pars(13)-pars(12) )
      Tfac = min(1d0,max(0d0,Tfac))
      ! photoperiod limitation
      Photofac = ( met(11,n)-pars(14)) / (pars(20)-pars(14) )
      Photofac = min(1d0,max(0d0,Photofac))
      ! VPD limitation
      VPDfac = 1.0 - ( (met(12,n)-pars(21)) / (pars(22)-pars(21)) )
      VPDfac = min(1d0,max(0d0,VPDfac))

      ! calculate and store the GSI index
      FLUXES(n,18) = Tfac * Photofac * VPDfac

      m = tmp_m(n)
      if (n == 1) then
          gsi_history(gsi_lag) = FLUXES(n,18)
      else
          gsi_history((gsi_lag-m):gsi_lag) = FLUXES((n-m):n,18)
      endif

      ! calculate gradient
      gradient = linear_model_gradient(tmp_x(1:(gsi_lag)),gsi_history(1:gsi_lag),gsi_lag)
      ! adjust gradient to daily rate
      if (deltat(n) > 1) then 
        if (nint((sum(deltat((n-m+1):n))) / (gsi_lag-1)) == 0) then 
          gradient =0 
        else 
          gradient = gradient / nint((sum(deltat((n-m+1):n))) / (gsi_lag-1))
        endif 
      endif
       
      gsi_lag_remembered = gsi_lag

      FLUXES(n,9) = 0d0  ! leaf turnover
      FLUXES(n,16) = 0d0 ! leaf growth

      ! update foliage and labile conditions based on gradient calculations
      if (gradient < fol_turn_crit .or. FLUXES(n,18) == 0) then
         ! we are in a decending condition so foliar turnover
         FLUXES(n,9) = pars(5)*(1.0-FLUXES(n,18))
         just_grown = 0.5
      else if (gradient > lab_turn_crit) then
         ! we are in a assending condition so labile turnover
         FLUXES(n,16) = pars(11)*FLUXES(n,18)
         just_grown = 1.5
         ! check carbon return
         tmp = POOLS(n,1)*(1d0-(1d0-FLUXES(n,16))**deltat(n))/deltat(n)
         tmp = (POOLS(n,2)+tmp)/pars(15)
         gpppars(1)=tmp
         tmp = acm(gpppars,constants)
         ! determine if increase in LAI leads to an improvement in GPP greater
         ! than critical value, if not then no labile turnover allowed
         if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(25) ) then
             FLUXES(n,16) = 0d0
         endif
      else

         if (just_grown >= 1.0) then
            FLUXES(n,16) = pars(11)*FLUXES(n,18)
            tmp = POOLS(n,1)*(1d0-(1d0-FLUXES(n,16))**deltat(n))/deltat(n)
            tmp = (POOLS(n,2)+tmp)/pars(15)
            gpppars(1)=tmp
            tmp = acm(gpppars,constants)
            if ( ((tmp - FLUXES(n,1))/FLUXES(n,1)) < pars(23) ) then
                FLUXES(n,16) = 0d0
            endif
         end if 
      endif 

      ! FLUXES WITH TIME DEPENDENCIES

      ! labile release = P_labile * (1-(1-leafgrowth)**deltat)/deltat
      FLUXES(n,8)  = POOLS(n,1)*(1.-(1.-FLUXES(n,16))**deltat(n))/deltat(n)
      ! leaf litter production = P_foliar * (1-(1-leaffall)**deltat)/deltat  
      FLUXES(n,10) = POOLS(n,2)*(1.-(1.-FLUXES(n,9))**deltat(n))/deltat(n)
      !  wood litter production
      FLUXES(n,11) = 0
      !  root litter production = P_root * (1-(1-rootTOR)**deltat)/deltat  
      FLUXES(n,12) = POOLS(n,3)*(1.-(1.-pars(6))**deltat(n))/deltat(n)

      ! FLUXES WITH TEMPERATURE AND TIME DEPENDENCIES

      ! resp het litter = P_litter * (1-(1-GPP_respired*litterTOR)**deltat)/deltat  
      FLUXES(n,13) = POOLS(n,5)*(1.-(1.-FLUXES(n,2)*pars(7))**deltat(n))/deltat(n)
      ! resp het som = P_som * (1-(1-GPP_respired*somTOR)**deltat)/deltat
      FLUXES(n,14) = POOLS(n,6)*(1.-(1.-FLUXES(n,2)*pars(8))**deltat(n))/deltat(n)
      ! litter to som = P_litter * (1-(1-dec_rate*temprate)**deltat)/deltat
      FLUXES(n,15) = POOLS(n,5)*(1.-(1.-pars(1)*FLUXES(n,2))**deltat(n))/deltat(n)

      ! NEE = resp_auto + resp_het_litter + resp_het_som - GPP
      NEE(n) = (FLUXES(n,3) + FLUXES(n,13) + FLUXES(n,14)) - FLUXES(n,1)
      ! GPP 
      GPP(n) = FLUXES(n,1)

      ! update pools for next timestep

      ! labile pool = labile_pool[†-1] + (lab_prod - lab_cons)*deltat
      POOLS(n+1,1) = POOLS(n,1) + (FLUXES(n,5)-FLUXES(n,8))*deltat(n)
      !  foliar pool = foliar_pool[†-1] + (leaf_prod - leaf_litter_prod + lab_prod2)*deltat
      POOLS(n+1,2) = POOLS(n,2) + (FLUXES(n,4)-FLUXES(n,10) + FLUXES(n,8))*deltat(n)
      ! wood pool
      POOLS(n+1,4) = 0.0
      ! root pool = root_pool[†-1] + (root_prod - root_litter_prod)*deltat
      POOLS(n+1,3) = POOLS(n,3) + (FLUXES(n,6)-FLUXES(n,12))*deltat(n)
      ! litter pool = litter_pool[†-1] + (leaf_litter_prod + root_litter_prod - resp_het_litter - litter2som)*deltat
      POOLS(n+1,5) = POOLS(n,5) + (FLUXES(n,10)+FLUXES(n,12)-FLUXES(n,13)-FLUXES(n,15))*deltat(n)
      ! som pool = som_pool[†-1] + (litter2som - resp_het_som + wood_litter_prod)
      POOLS(n+1,6) = POOLS(n,6) + (FLUXES(n,15)-FLUXES(n,14)+FLUXES(n,11))*deltat(n)


      if (version_code .EQ. 2) then 

        ! CUTTING (if : AGB > cutting limit & met(8,n) = 100)
        if ( ((POOLS(n+1,2)+POOLS(n+1,1)) .GE. (pars(28)*0.475*0.1)) .AND. (met(8,n) .EQ. 100) ) then
                                                                   
            ! direct C losses
            labile_loss  = POOLS(n+1,1) * pars(33)
            foliar_loss  = max(0.,POOLS(n+1,2) - (pars(27)*0.475*0.1 + labile_loss))
            roots_loss   = 0 

            ! fraction of harvest wasted 
            labile_residue = labile_loss * labile_frac_res
            foliar_residue = foliar_loss * foliage_frac_res

            ! extracted carbon via cutting
            REMOVED_C(2,n) = (labile_loss-labile_residue) + (foliar_loss-foliar_residue)

            ! update pools 
            POOLS(n+1,1) = max(0.,POOLS(n+1,1)-labile_loss)
            POOLS(n+1,2) = max(0.,POOLS(n+1,2)-foliar_loss)
            POOLS(n+1,3) = max(0.,POOLS(n+1,3)-roots_loss)
            POOLS(n+1,4) = 0.0
            POOLS(n+1,5) = max(0., POOLS(n+1,5) + (labile_residue+foliar_residue+roots_loss))
            POOLS(n+1,6) = max(0., POOLS(n+1,6))

        endif ! end cutting 

        ! GRAZING (if LSU.ha-1 > 0 & AGB > limit )
        if ( (met(8,n) > 0.) .AND. (met(8,n) .NE. 100 ) .AND. ((POOLS(n+1,2)+POOLS(n+1,1)) .GE. (pars(27)*0.475*0.1)) ) then
            
            ! direct C losses
            labile_loss  = POOLS(n+1,1) * pars(32)
            ! Remove demand (g.C.m-2) from foliage : LSU.per.ha * 650 kg_weight * 2.5% * convert_kg.DM.ha-1_to_g.C.m-2 
            foliar_loss  = max(0.,met(8,n) * 650 * pars(31) * 0.047619 - labile_loss)
            roots_loss   = 0 ! POOLS(n+1,3) * roots_frac_death

            ! fraction of harvest wasted 
            labile_residue = labile_loss * labile_frac_res
            foliar_residue = foliar_loss * foliage_frac_res

            ! extracted carbon via grazing (if grass remains  > grazing limit DM )
            if ( ((POOLS(n+1,2)+POOLS(n+1,1))-foliar_loss-labile_loss) .GE. (pars(27)*0.475*0.1) ) then
                
                ! extracted carbon via grazing
                REMOVED_C(1,n) = (labile_loss-labile_residue) + (foliar_loss-foliar_residue)

                ! animal manure-C production
                FLUXES(n,19) = REMOVED_C(1,n) * 0.32
                
                ! animal respiration CO2-C
                FLUXES(n,20) = REMOVED_C(1,n) * 0.54 

                ! animal CH4-C 
                FLUXES(n,21) = REMOVED_C(1,n) * 0.04 

                ! update pools 
                POOLS(n+1,1) = max(0.,POOLS(n+1,1)-labile_loss)
                POOLS(n+1,2) = max(0.,POOLS(n+1,2)-foliar_loss)
                POOLS(n+1,3) = max(0.,POOLS(n+1,3)-roots_loss)
                POOLS(n+1,4) = 0.0
                POOLS(n+1,5) = max(0., POOLS(n+1,5) + (labile_residue+foliar_residue+roots_loss) + FLUXES(n,19) )
                POOLS(n+1,6) = max(0., POOLS(n+1,6))
            endif 

            ! extracted carbon via grazing (if grass remains  < grazing limit DM )
            if ( (REMOVED_C(1,n) .EQ. 0.0) .AND. ((POOLS(n+1,2)+POOLS(n+1,1))-foliar_loss-labile_loss) < (pars(27)*0.475*0.1) ) then

                ! direct C losses
                labile_loss  = POOLS(n+1,1) * pars(32)
                foliar_loss  = POOLS(n+1,2) - (pars(27)*0.475*0.1 + labile_loss)
                roots_loss   = 0 ! POOLS(n+1,3) * roots_frac_death

                if (foliar_loss > 0.0 ) then

                  ! fraction of harvest wasted 
                  labile_residue = labile_loss * labile_frac_res
                  foliar_residue = foliar_loss * foliage_frac_res

                  ! extracted carbon via grazing
                  REMOVED_C(1,n) = (labile_loss-labile_residue) + (foliar_loss-foliar_residue)

                  ! animal manure-C production
                  FLUXES(n,19) = REMOVED_C(1,n) * 0.32
                  
                  ! animal respiration CO2-C
                  FLUXES(n,20) = REMOVED_C(1,n) * 0.54 

                  ! animal CH4-C 
                  FLUXES(n,21) = REMOVED_C(1,n) * 0.04 

                  ! update pools 
                  POOLS(n+1,1) = max(0.,POOLS(n+1,1)-labile_loss)
                  POOLS(n+1,2) = max(0.,POOLS(n+1,2)-foliar_loss)
                  POOLS(n+1,3) = max(0.,POOLS(n+1,3)-roots_loss)
                  POOLS(n+1,4) = 0.0
                  POOLS(n+1,5) = max(0., POOLS(n+1,5) + (labile_residue+foliar_residue+roots_loss) + FLUXES(n,19) )
                  POOLS(n+1,6) = max(0., POOLS(n+1,6))
                endif 

            endif 

        endif 
      
      endif

    end do


  end subroutine CARBON_MODEL
  
  !
  !------------------------------------------------------------------
  !
  
  double precision function acm(drivers,constants)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    double precision, intent(in) :: drivers(12) & ! ACM inputs
                                   ,constants(10) ! ACM parameters

    double precision :: gc, pn, pd, pp, qq, ci, e0, dayl, cps, dec, nit &
                       ,trange, sinld, cosld,aob  &
                       ,mint,maxt,radiation,co2,lai,doy,lat &
                       ,deltaWP,Rtot,NUE,temp_exponent,dayl_coef &
                       ,dayl_const,hydraulic_exponent,hydraulic_temp_coef &
                       ,co2_comp_point,co2_half_sat,lai_coef,lai_const

    gc=0.0 ; pp=0.0 ; qq=0.0 ; ci=0.0 ; e0=0.0 ; dayl=0.0 ; cps=0.0 ; dec=0.0 ; nit=1.0

    lai  = drivers(1)
    maxt = drivers(2)
    mint = drivers(3)
    nit  = drivers(4)   
    co2  = drivers(5)
    doy  = drivers(6)
    lat = drivers(7)
    radiation = drivers(8)
    deltaWP = drivers(9)
    Rtot = drivers(10)

    NUE = constants(1)
    dayl_coef = constants(2)
    co2_comp_point = constants(3) 
    co2_half_sat = constants(4)
    dayl_const = constants(5)
    hydraulic_temp_coef = constants(6)
    lai_coef = constants(7)
    temp_exponent = constants(8)
    lai_const = constants(9)
    hydraulic_exponent = constants(10)

    ! determine temperature range 
    trange = 0.5*(maxt-mint)
    ! daily canopy conductance, of CO2 or H2O? 
    gc = abs(deltaWP)**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot+trange))
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn = lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp = pn/gc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    ci = 0.5*(co2+qq-pp+sqrt(((co2+qq-pp)*(co2+qq-pp))-4.0*(co2*qq-pp*co2_comp_point)))
    ! limit maximum quantium efficiency by leaf area, hyperbola
    e0 = lai_coef*(lai*lai)/((lai*lai)+lai_const)
    ! calculate day length (hours)
    dec = - asin( sin( 23.45 * deg_to_rad ) * cos( 2.0 * pi * ( doy + 10.0 ) / 365.0 ) )
    sinld = sin( lat*deg_to_rad ) * sin( dec )
    cosld = cos( lat*deg_to_rad ) * cos( dec )
    aob = max(-1.0,min(1.0,sinld / cosld))
    dayl = 12.0 * ( 1.0 + 2.0 * asin( aob ) / pi )

    ! calculate CO2 limited rate of photosynthesis
    pd=gc*(co2-ci)
    ! calculate combined light and CO2 limited photosynthesis
    cps=e0*radiation*pd/(e0*radiation+pd)
    ! correct for day length variation
    acm=cps*(dayl_coef*dayl+dayl_const)

    return

  end function acm
  
  !
  !------------------------------------------------------------------
  !
  
  double precision function linear_model_gradient(x,y,interval)

    ! Function to calculate the gradient of a linear model for a given depentent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    implicit none

    ! declare input variables
    integer :: interval
    double precision, dimension(interval) :: x,y

    ! declare local variables
    double precision :: sum_x, sum_y, sumsq_x,sum_product_xy

    ! calculate the sum of x
    sum_x = sum(x)
    ! calculate the sum of y
    sum_y = sum(y)
    ! calculate the sum of squares of x
    sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    sum_product_xy = sum(x*y)
    ! calculate the gradient
    linear_model_gradient = ( (interval*sum_product_xy) - (sum_x*sum_y) ) &
                          / ( (interval*sumsq_x) - (sum_x*sum_x) )

    return

  end function linear_model_gradient

!
!--------------------------------------------------------------------
!
end module CARBON_MODEl_MOD
