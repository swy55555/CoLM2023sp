MODULE YOS_CMF_DIAG
!==========================================================
!* PURPOSE: Shared diagnostic variables for CaMa-Flood
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
    
! Modified by Zhongwang Wei @ SYSU 2022.11.20: add water re-infiltration calculation 

! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
USE PARKIND1, ONLY: JPIM, JPRB, JPRM, JPRD
IMPLICIT NONE
SAVE
!================================================
!*** Inst. diagnostics 
REAL(KIND=JPRB),ALLOCATABLE,TARGET :: D2DIAG(:,:,:)        !! Storage array for average diagnostics

REAL(KIND=JPRB),POINTER         :: D2RIVINF(:,:)           !! river      inflow   [m3/s] (from upstream)
REAL(KIND=JPRB),POINTER         :: D2RIVDPH(:,:)           !! river      depth    [m]
REAL(KIND=JPRB),POINTER         :: D2RIVVEL(:,:)           !! flow velocity       [m/s]

REAL(KIND=JPRB),POINTER         :: D2FLDINF(:,:)           !! floodplain inflow   [m3/s]
REAL(KIND=JPRB),POINTER         :: D2FLDDPH(:,:)           !! floodplain depth    [m]
REAL(KIND=JPRB),POINTER         :: D2FLDFRC(:,:)           !! flooded    fractipn [m2/m2]
REAL(KIND=JPRB),POINTER         :: D2FLDARE(:,:)           !! flooded    area     [m2]

REAL(KIND=JPRB),POINTER         :: D2PTHOUT(:,:)           !! flood path outflow   [m3/s]
REAL(KIND=JPRB),POINTER         :: D2PTHINF(:,:)           !! flood path inflow   [m3/s]

REAL(KIND=JPRB),POINTER         :: D2SFCELV(:,:)           !! water surface elev  [m]    (elevtn - rivhgt + rivdph)
REAL(KIND=JPRB),POINTER         :: D2OUTFLW(:,:)           !! total outflow       [m3/s] (rivout + fldout)
REAL(KIND=JPRB),POINTER         :: D2STORGE(:,:)           !! total storage       [m3]   (rivsto + fldsto)

REAL(KIND=JPRB),POINTER         :: D2OUTINS(:,:)           !! instantaneous discharge [m3/s] (unrouted runoff)
REAL(KIND=JPRB),POINTER         :: D2WEVAPEX(:,:)          !! Evaporation water extracted
REAL(KIND=JPRB),POINTER         :: D2WINFILTEX(:,:)          !! Infiltration water extracted

INTEGER(KIND=JPIM)              :: N2DIAG                  !! number of 2D diagnostics

!================================================
!*** Average diagnostics 
REAL(KIND=JPRB),ALLOCATABLE,TARGET :: D2DIAG_AVG(:,:,:)    !! Storage array for average diagnostics (nseqmax,1,variable)

REAL(KIND=JPRB),POINTER         :: D2RIVOUT_AVG(:,:)       !! average river       discharge
REAL(KIND=JPRB),POINTER         :: D2OUTFLW_AVG(:,:)       !! average total outflow       [m3/s] (rivout + fldout)  !! bugfix v362
REAL(KIND=JPRB),POINTER         :: D2FLDOUT_AVG(:,:)       !! average floodplain  discharge
REAL(KIND=JPRB),POINTER         :: D2RIVVEL_AVG(:,:)       !! average flow velocity
REAL(KIND=JPRB),POINTER         :: D2PTHOUT_AVG(:,:)       !! flood pathway net outflow (2D)

REAL(KIND=JPRB),POINTER         :: D2GDWRTN_AVG(:,:)       !! average ground water return flow
REAL(KIND=JPRB),POINTER         :: D2RUNOFF_AVG(:,:)       !! average input runoff
REAL(KIND=JPRB),POINTER         :: D2ROFSUB_AVG(:,:)       !! average input sub-surface runoff
REAL(KIND=JPRB),POINTER         :: D2WEVAPEX_AVG(:,:)      !! average extracted water evaporation
REAL(KIND=JPRB),POINTER         :: D2WINFILTEX_AVG(:,:)    !! average extracted Infiltration water

INTEGER(KIND=JPIM)              :: N2DIAG_AVG              !! Number of 2D diagnostics averages
REAL(KIND=JPRB)                 :: NADD                    !! sum DT to calculate average
!*** Average diagnostics (1D)
REAL(KIND=JPRB),ALLOCATABLE  :: D1PTHFLW_AVG(:,:)          !! bifurcation channel flow (1D, not 2D variable)

!================================================
!*** Daily max diagnostics 
REAL(KIND=JPRB),ALLOCATABLE,TARGET :: D2DIAG_MAX(:,:,:)    !! Storage array for maximum diagnostics (nseqmax,1,variable)

REAL(KIND=JPRB),POINTER         :: D2OUTFLW_MAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)
REAL(KIND=JPRB),POINTER         :: D2STORGE_MAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)
REAL(KIND=JPRB),POINTER         :: D2RIVDPH_MAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)

INTEGER(KIND=JPRB)              :: N2DIAG_MAX              !! Number of 2D diagnostics maximum

!================================================
!*** Global total
! discharge calculation budget
REAL(KIND=JPRD)                 :: P0GLBSTOPRE              !! global water storage      [m3] (befre flow calculation)
REAL(KIND=JPRD)                 :: P0GLBSTONXT              !! global water storage      [m3] (after flow calculation)
REAL(KIND=JPRD)                 :: P0GLBSTONEW              !! global water storage      [m3] (after runoff input)
REAL(KIND=JPRD)                 :: P0GLBRIVINF              !! global inflow             [m3] (rivinf + fldinf)
REAL(KIND=JPRD)                 :: P0GLBRIVOUT              !! global outflow            [m3] (rivout + fldout)

! stage calculation budget
REAL(KIND=JPRD)                 :: P0GLBSTOPRE2             !! global water storage      [m3] (befre stage calculation)
REAL(KIND=JPRD)                 :: P0GLBSTONEW2             !! global water storage      [m3] (after stage calculation)
REAL(KIND=JPRD)                 :: P0GLBRIVSTO              !! global river storage      [m3]
REAL(KIND=JPRD)                 :: P0GLBFLDSTO              !! global floodplain storage [m3]
REAL(KIND=JPRD)                 :: P0GLBLEVSTO              !! global protected-side storage [m3] (levee scheme)
REAL(KIND=JPRD)                 :: P0GLBFLDARE              !! global flooded area       [m2]

!================================================
!*** dam variable
REAL(KIND=JPRB),POINTER         :: D2DAMINF_AVG(:,:)       !! average reservoir inflow [m3/s]  !!!added

!================================================
!!!*** levee variables
REAL(KIND=JPRB),POINTER         :: D2LEVDPH(:,:)           !! flood depth in protected side (water depth betwen river & levee)



END MODULE YOS_CMF_DIAG
