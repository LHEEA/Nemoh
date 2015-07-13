!--------------------------------------------------------------------------------------
!
!   Copyright 2014 Ecole Centrale de Nantes, 1 rue de la Noë, 44300 Nantes, France
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License. 
!
!   Contributors list:
!   - G. Delhommeau
!   - J. Singh
!   - A. Babarit
!
!--------------------------------------------------------------------------------------
  MODULE FIC_COM

    IMPLICIT NONE
    INTEGER:: IIR,JJZ
    REAL::XR,XZ,APD1X,APD1Z,APD2X,APD2Z
    COMMON/FIC/IIR,JJZ,XR(676),XZ(124),APD1X(676,124),APD1Z(676,124),APD2X(676,124),APD2Z(676,124)

    END MODULE FIC_COM