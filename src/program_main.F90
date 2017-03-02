!-----------------------------!
! EXAMPLE 3 IN ZHAO JSC 2014  !
!-----------------------------!
#if example == 3

#if   time == 1
   INCLUDE "program_adid1_disc.f90"
#elif time == 2
   INCLUDE "program_adipr_mib2d.f90"
#elif time == 3
   INCLUDE "program_cn_disc.f90"
#elif time == 4
   INCLUDE "program_ie_disc.f90"
#endif

!-----------------------------!
! EXAMPLE 5 IN ZHAO JSC 2014  !
!-----------------------------!
#elif example == 5

#if   time == 1
   INCLUDE "program_adid1_leaves.f90"
#elif time == 2
   INCLUDE "program_adipr_mib2d.f90"
#elif time == 3
   INCLUDE "program_cn_leaves.f90"
#elif time == 4
   INCLUDE "program_ie_leaves.f90"
#endif

!-----------------------------!
! NEW EXAMPLE 6               !
!-----------------------------!
#elif example == 6

#if   time == 2
   INCLUDE "program_adipr_mib2d.f90"
#endif

!-----------------------------!
! NEW EXAMPLE 7               !
!-----------------------------!
#elif example == 7

#if   time == 2
   INCLUDE "program_adipr_mib2d.f90"
#endif

#endif
