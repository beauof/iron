!> \file
!> $Id: opencmiss_c.f90 542 2009-06-03 17:16:22Z chrispbradley $
!> \author Chris Bradley
!> \brief The top level OpenCMISS module for C bindings.
!>
!> \mainpage OpenCMISS Documentation
!>
!> An open source interactive computer program for Continuum Mechanics, Image analysis, Signal processing and System
!> Identification. Target usage: Bioengineering application of finite element analysis, boundary element and collocation
!> techniques.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>
!>
!> The top level OpenCMISS module for C. This module is the buffer module between the OpenCMISS library and user C code.
MODULE OPENCMISS_C

  USE ISO_C_BINDING
  USE ISO_VARYING_STRING
  USE OPENCMISS
   
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  INTEGER(C_INT), PARAMETER :: CMISSNoError = 0
  INTEGER(C_INT), PARAMETER :: CMISSPointerIsNULL = -1
  INTEGER(C_INT), PARAMETER :: CMISSPointerNotNULL = -2
  INTEGER(C_INT), PARAMETER :: CMISSCouldNotAllocatePointer = -3
  INTEGER(C_INT), PARAMETER :: CMISSErrorConvertingPointer = -4
  
  !Module types

  !Module variables

  !Interfaces

  PUBLIC CMISSBasisTypeFinaliseC, CMISSBasisTypeInitialiseC

  PUBLIC CMISSBoundaryConditionsTypeFinaliseC, CMISSBoundaryConditionsTypeInitialiseC

  PUBLIC CMISSControlLoopTypeFinaliseC, CMISSControlLoopTypeInitialiseC

  PUBLIC CMISSCoordinateSystemTypeFinaliseC,CMISSCoordinateSystemTypeInitialiseC

  PUBLIC CMISSDecompositionTypeFinaliseC, CMISSDecompositionTypeInitialiseC

  PUBLIC CMISSEquationsTypeFinaliseC, CMISSEquationsTypeInitialiseC

  PUBLIC CMISSEquationsSetTypeFinaliseC, CMISSEquationsSetTypeInitialiseC

  PUBLIC CMISSRegionTypeFinaliseC, CMISSRegionTypeInitialiseC
  
  PUBLIC CMISSFinaliseC,CMISSInitialiseCNum,CMISSInitialiseCPtr

CONTAINS

  !
  !================================================================================================================================
  !

  !>Copys/converts a C string (array of characters) to a Fortran String (length of characters)
  SUBROUTINE CMISSC2FString(Cstring,Fstring)
    !Argument variables
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Cstring(:)
    CHARACTER(LEN=*), INTENT(OUT) :: Fstring
    !Local variables
    INTEGER(C_INT) :: i,LENGTH

    IF(LEN(Fstring)>=SIZE(Cstring,1)-1) THEN
      LENGTH=SIZE(Cstring,1)-1
    ELSE
      LENGTH=LEN(Fstring)
    ENDIF
    Fstring=""
    DO i=1,LENGTH
      IF(Cstring(i)==C_NULL_CHAR) THEN
        EXIT
      ELSE
        Fstring(i:i)=Cstring(i)
      ENDIF
    ENDDO !i
    
    RETURN
    
  END SUBROUTINE CMISSC2FSTRING
   
  !
  !================================================================================================================================
  !

  !>Copys/converts a  Fortran String (length of characters) to a C string (array of characters)
  SUBROUTINE CMISSF2CString(Fstring,Cstring)
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: Fstring
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Cstring(:)
    !Local variables
    INTEGER(C_INT) :: i,LENGTH

    IF(SIZE(Cstring,1)>LEN_TRIM(Fstring)) THEN
      LENGTH=LEN_TRIM(Fstring)
    ELSE
      LENGTH=SIZE(Cstring,1)-1
    ENDIF
    DO i=1,LENGTH     
      Cstring(i)=Fstring(i:i)
    ENDDO !i
    !Null terminate the string
    Cstring(LENGTH+1)=C_NULL_CHAR
    
    RETURN
    
  END SUBROUTINE CMISSF2CSTRING
   
  !
  !================================================================================================================================
  !

  !>Finalises a CMISSBasisType object for C.
  FUNCTION CMISSBasisTypeFinaliseC(BasisTypePtr) BIND(C, NAME= "CMISSBasisTypeFinalise")

    !Argument Variables
    TYPE(C_PTR), INTENT(INOUT) :: BasisTypePtr !<C pointer to CMISSBasisType object to finalise.
    !Function Variable
    INTEGER(C_INT) :: CMISSBasisTypeFinaliseC !<Error code.
    !Local Variables
    TYPE(CMISSBasisType), POINTER :: BasisType

    CMISSBasisTypeFinaliseC = CMISSNoError

    IF(C_ASSOCIATED(BasisTypePtr)) THEN
      CALL C_F_POINTER(BasisTypePtr,BasisType)
      IF(ASSOCIATED(BasisType)) THEN
        CALL CMISSBasisTypeFinalise(BasisType,CMISSBasisTypeFinaliseC)
        DEALLOCATE(BasisType)
        BasisTypePtr = C_NULL_PTR
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSBasisTypeFinaliseC

  !
  !============================================================================
  !

  !>Initialises a CMISSBasisType object for C.

  FUNCTION CMISSBasisTypeInitialiseC (BasisTypePtr) BIND(C, NAME = "CMISSBasisTypeInitialise")

	!Argument variables
  TYPE(C_PTR), INTENT(INOUT) :: BasisTypePtr !<C pointer to CMISSBasisType object to initialise.
	!Function variable
  INTEGER(C_INT) :: CMISSBasisTypeInitialiseC !<Error code.
	!Local Variables
  INTEGER(C_INT) :: Err
    TYPE (CMISSBasisType), POINTER :: BasisType

    IF(C_ASSOCIATED(BasisTypePtr)) THEN
      CMISSBasisTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY (BasisType)
      ALLOCATE (BasisType, STAT = Err)
      IF (Err /= 0) THEN
        CMISSBasisTypeInitialiseC =CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSBasisTypeInitialise (BasisType,CMISSBasisTypeInitialiseC)
        BasisTypePtr = C_LOC (BasisType)
        CMISSBasisTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF
    RETURN

  END FUNCTION CMISSBasisTypeInitialiseC

  !
  !============================================================================
  !

  !>Finalises a CMISSBoundaryConditionsType object for C.

  FUNCTION CMISSBoundaryConditionsTypeFinaliseC(BoundaryConditionsTypePtr) BIND (C, NAME = "CMISSBoundaryConditionsTypeFinalise")

    !Argument Variables
    TYPE (C_PTR), INTENT(INOUT) :: BoundaryConditionsTypePtr !<C pointer to CMISSBoundaryConditionsType object to finalise.
    !Function Variable
    INTEGER(C_INT) :: CMISSBoundaryConditionsTypeFinaliseC !<Error Code.
    !Local Variables
    TYPE(CMISSBoundaryConditionsType), POINTER :: BoundaryConditionsType

    CMISSBoundaryConditionsTypeFinaliseC = CMISSNoError
    IF (C_ASSOCIATED(BoundaryConditionsTypePtr)) THEN
      CALL C_F_POINTER (BoundaryConditionsTypePtr, BoundaryConditionsType)
      IF(ASSOCIATED(BoundaryConditionsType)) THEN
        CALL CMISSBoundaryConditionsTypeFinalise (BoundaryConditionsType, CMISSBoundaryConditionsTypeFinaliseC)
        DEALLOCATE (BoundaryConditionsType)
        BoundaryConditionsTypePtr = C_NULL_PTR
      ELSE
        CMISSBoundaryConditionsTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSBoundaryConditionsTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSBoundaryConditionsTypeFinaliseC

  !
  !============================================================================
  !

  !>Initialises a CMISSBoundaryConditionsType object for C.

  FUNCTION CMISSBoundaryConditionsTypeInitialiseC (BoundaryConditionsTypePtr) BIND (C, NAME = &
  & "CMISSBoundaryConditionsTypeInitialise")

    !Argument variables
    TYPE(C_PTR), INTENT (INOUT) :: BoundaryConditionsTypePtr !<C pointer to the CMISSBoundaryConditionsType object to be initialised.
    !Function variables
    INTEGER(C_INT) :: CMISSBoundaryConditionsTypeInitialiseC !<Error Code.
    !Local Variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSBoundaryConditionsType), POINTER :: BoundaryConditionsType

    IF (C_ASSOCIATED(BoundaryConditionsTypePtr)) THEN
      CMISSBoundaryConditionsTypeInitialiseC = CMISSPointerNotNull
    ELSE
      NULLIFY (BoundaryConditionsType)
      ALLOCATE(BoundaryConditionsType, STAT = Err)
      IF (Err /= 0) THEN
        CMISSBoundaryConditionsTypeInitialiseC =CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsType, CMISSBoundaryConditionsTypeInitialiseC)
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSBoundaryConditionsTypeInitialiseC

  !
  !===========================================================================
  !

  !>Finalises a CMISSControlLoopType object for C.

  FUNCTION CMISSControlLoopTypeFinaliseC (ControlLoopTypePtr) BIND (C, NAME = "CMISSControlLoopTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT (INOUT) :: ControlLoopTypePtr !<C pointer to the CMISSControlLoopType object to be finalised.
    !Function variables
    INTEGER(C_INT) :: CMISSControlLoopTypeFinaliseC !<Error Code.
    !Local variables
    TYPE(CMISSControlLoopType), POINTER :: ControlLoopType

    CMISSControlLoopTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(ControlLoopTypePtr)) THEN
      CALL C_F_POINTER(ControlLoopTypePtr, ControlLoopType)
      IF(ASSOCIATED(ControlLoopType)) THEN
        CALL CMISSControlLoopTypeFinalise(ControlLoopType,CMISSControlLoopTypeFinaliseC)
        DEALLOCATE(ControlLoopType)
        ControlLoopTypePtr = C_NULL_PTR
      ELSE
        CMISSControlLoopTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSControlLoopTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSControlLoopTypeFinaliseC

  !
  !===========================================================================
  !

  !>Initialises a CMISSControlLoopType object for C.

  FUNCTION CMISSControlLoopTypeInitialiseC (ControlLoopTypePtr) BIND(C, NAME = "CMISSControlLoopTypeInitialise")

    !Argument variables
    TYPE(C_PTR), INTENT (INOUT) :: ControlLoopTypePtr !<C pointer to the CMISSControlLoopType object to be intialised.
    !Function variables
    INTEGER(C_INT) :: CMISSControlLoopTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE (CMISSControlLoopType), POINTER :: ControlLoopType

    IF(C_ASSOCIATED(ControlLoopTypePtr)) THEN
      CMISSControlLoopTypeInitialiseC=CMISSPointerNotNULL
    ELSE
      NULLIFY(ControlLoopType)
      ALLOCATE(ControlLoopType,STAT=Err)
      IF(Err/=0) THEN
        CMISSControlLoopTypeInitialiseC=CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSControlLoopTypeInitialise(ControlLoopType,CMISSControlLoopTypeInitialiseC)
        ControlLoopTypePtr=C_LOC(ControlLoopType)
        CMISSControlLoopTypeInitialiseC=CMISSNoError
      ENDIF
    ENDIF


    RETURN

  END FUNCTION CMISSControlLoopTypeInitialiseC

  !
  !===========================================================================
  !

  !>Finalises a CMISSCoordinateSystemType object for C.
  FUNCTION CMISSCoordinateSystemTypeFinaliseC(CoordinateSystemTypePtr)  BIND(C,NAME="CMISSCoordinateSystemTypeFinalise")
    
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: CoordinateSystemTypePtr !<C pointer to the CMISSCoordinateSystemType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSCoordinateSystemTypeFinaliseC !<Error Code.
    !Local variables
    TYPE(CMISSCoordinateSystemType), POINTER :: CoordinateSystemType
    
    CMISSCoordinateSystemTypeFinaliseC=CMISSNoError
    IF(C_ASSOCIATED(CoordinateSystemTypePtr)) THEN
      CALL C_F_POINTER(CoordinateSystemTypePtr,CoordinateSystemType)
      IF(ASSOCIATED(CoordinateSystemType)) THEN
        CALL CMISSCoordinateSystemTypeFinalise(CoordinateSystemType,CMISSCoordinateSystemTypeFinaliseC)
        DEALLOCATE(CoordinateSystemType)
        CoordinateSystemTypePtr=C_NULL_PTR
      ELSE
        CMISSCoordinateSystemTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSCoordinateSystemTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN
    
  END FUNCTION CMISSCoordinateSystemTypeFinaliseC
 
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSCoordinateSystemType object for C.
  FUNCTION CMISSCoordinateSystemTypeInitialiseC(CoordinateSystemTypePtr)  BIND(C,NAME="CMISSCoordinateSystemTypeInitialise")
    
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: CoordinateSystemTypePtr !<C pointer to the CMISSCoordinateSystemType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSCoordinateSystemTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSCoordinateSystemType), POINTER :: CoordinateSystemType
    
    IF(C_ASSOCIATED(CoordinateSystemTypePtr)) THEN
      CMISSCoordinateSystemTypeInitialiseC=CMISSPointerNotNULL
    ELSE
      NULLIFY(CoordinateSystemType)
      ALLOCATE(CoordinateSystemType,STAT=Err)
      IF(Err/=0) THEN
        CMISSCoordinateSystemTypeInitialiseC=CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystemType,CMISSCoordinateSystemTypeInitialiseC)
        CoordinateSystemTypePtr=C_LOC(CoordinateSystemType)
        CMISSCoordinateSystemTypeInitialiseC=CMISSNoError
      ENDIF
    ENDIF

    RETURN
    
  END FUNCTION CMISSCoordinateSystemTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSDecompositionType object.
  FUNCTION CMISSDecompositionTypeFinaliseC(DecompositionTypePtr) BIND(C, NAME = "CMISSDecompositionTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: DecompositionTypePtr !<C pointer to the CMISSDecompositionType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionTypeFinaliseC !Error Code.
    !Local Variables
    TYPE(CMISSDecompositionType), POINTER :: DecompositionType

    CMISSDecompositionTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(DecompositionTypePtr)) THEN
      CALL C_F_POINTER (DecompositionTypePtr, DecompositionType)
      IF(ASSOCIATED(DecompositionType)) THEN
        CALL CMISSDecompositionTypeFinalise (DecompositionType, CMISSDecompositionTypeFinaliseC)
        DEALLOCATE(DecompositionType)
        DecompositionTypePtr = C_NULL_PTR
      ELSE
        CMISSDecompositionTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSDecompositionType object.
  FUNCTION CMISSDecompositionTypeInitialiseC(DecompositionTypePtr)  BIND(C, NAME = "CMISSDecompositionTypeInitialise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: DecompositionTypePtr !<C pointer to the CMISSDecompositionType object to be initialised.
    !Function variables
    INTEGER(C_INT) ::  CMISSDecompositionTypeInitialiseC !<Error Code.
    !Local Variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSDecompositionType), POINTER :: DecompositionType

    IF(C_ASSOCIATED(DecompositionTypePtr)) THEN
      CMISSDecompositionTypeInitialiseC = CMISSPointerNotNull
    ELSE
      NULLIFY(DecompositionType)
      ALLOCATE(DecompositionType, STAT=Err)
      IF (Err/=0) THEN
        CMISSDecompositionTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSDecompositionTypeInitialise(DecompositionType, CMISSDecompositionTypeInitialiseC)
        DecompositionTypePtr = C_LOC(DecompositionType)
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSEquationsType object for C.
  FUNCTION CMISSEquationsTypeFinaliseC (EquationsTypePtr) BIND(C, NAME = "CMISSEquationsTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: EquationsTypePtr !<C pointer to the CMISSEquationsType object to be finalised.
    !Function variables
    INTEGER(C_INT) :: CMISSEquationsTypeFinaliseC !<Error Code.
    !Local variables
    TYPE(CMISSEquationsType), POINTER :: EquationsType

    CMISSEquationsTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(EquationsTypePtr)) THEN
      CALL C_F_POINTER(EquationsTypePtr, EquationsType)
      IF (ASSOCIATED(EquationsType)) THEN
        CALL CMISSEquationsTypeFinalise (EquationsType, CMISSEquationsTypeFinaliseC)
        DEALLOCATE(EquationsType)
        EquationsTypePtr = C_NULL_PTR
      ELSE
        CMISSEquationsTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSEquationsTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

      RETURN

    END FUNCTION CMISSEquationsTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSEquationsType object for C.
  FUNCTION CMISSEquationsTypeInitialiseC (EquationsTypePtr) BIND(C, NAME = "CMISSEquationsTypeInitialise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: EquationsTypePtr  !<C pointer to the CMISSEquationsType object to be intialised.
    !Function variables
    INTEGER(C_INT) :: CMISSEquationsTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSEquationsType), POINTER :: EquationsType

    IF(C_ASSOCIATED(EquationsTypePTR)) THEN
      CMISSEquationsTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(EquationsType)
      ALLOCATE(EquationsType, STAT= Err)
      IF(Err/=0) THEN
        CMISSEquationsTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSEquationsTypeInitialise(EquationsType, CMISSEquationsTypeInitialiseC)
        EquationsTypePtr = C_LOC(EquationsType)
        CMISSEquationsTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN
  END FUNCTION CMISSEquationsTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSEquationsSetType object for C.
  FUNCTION CMISSEquationsSetTypeFinaliseC (EquationsSetTypePtr) BIND(C, NAME = "CMISSEquationsSetTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: EquationsSetTypePtr  !<C pointer to the CMISSEquationsSetType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSEquationsSetTypeFinaliseC !<Error Code.
    !Local variables
    TYPE(CMISSEquationsSetType), POINTER :: EquationsSetType

    CMISSEquationsSetTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(EquationsSetTypePtr)) THEN
      CALL C_F_POINTER(EquationsSetTypePtr, EquationsSetType)
      IF (ASSOCIATED(EquationsSetType)) THEN
        CALL CMISSEquationsSetTypeFinalise(EquationsSetType, CMISSEquationsSetTypeFinaliseC)
        DEALLOCATE(EquationsSetType)
        EquationsSetTypePtr = C_NULL_PTR
      ELSE
        CMISSEquationsSetTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSEquationsSetTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN
  END FUNCTION CMISSEquationsSetTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSEquationsSetType object for C.
  FUNCTION CMISSEquationsSetTypeInitialiseC(EquationsSetTypePtr) BIND(C, NAME = "CMISSEquationsSetTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: EquationsSetTypePtr  !<C pointer to the CMISSEquationsSetType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSEquationsSetTypeInitialiseC !<Error Code.
    !Local variable
    INTEGER(C_INT) :: Err
    TYPE(CMISSEquationsSetType), POINTER :: EquationsSetType

    IF(C_ASSOCIATED(EquationsSetTypePtr)) THEN
      CMISSEquationsSetTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY (EquationsSetType)
      ALLOCATE(EquationsSetType, STAT = Err)
      IF(Err/=0) THEN
        CMISSEquationsSetTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSEquationsSetTypeInitialise(EquationsSetType, CMISSEquationsSetTypeInitialiseC)
        EquationsSetTypePtr = C_LOC(EquationsSetType)
        CMISSEquationsSetTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSEquationsSetTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSFieldType object for C.
  FUNCTION CMISSFieldTypeFinaliseC(FieldTypePtr) BIND(C,NAME="CMISSFieldTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: FieldTypePtr !<C pointer to the CMISSFieldType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: FieldType

    CMISSFieldTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(FieldTypePtr)) THEN
      CALL C_F_POINTER (FieldTypePtr, FieldType)
      IF(ASSOCIATED(FieldType)) THEN
        CALL CMISSFieldTypeFinalise(FieldType, CMISSFieldTypeFinaliseC)
        DEALLOCATE(FieldType)
        FieldTypePtr = C_NULL_PTR
     ELSE
        CMISSFieldTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSFieldType object for C.
  FUNCTION CMISSFieldTypeInitialiseC(FieldTypePtr) BIND(C, NAME = "CMISSFieldTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: FieldTypePtr !<C pointer to the CMISSFieldType object to be finalised.
    !Function Variable
    INTEGER(C_INT) :: CMISSFieldTypeInitialiseC !<Error Code
    !Local variable
    INTEGER(C_INT) :: Err
    TYPE(CMISSFieldType), POINTER :: FieldType

    IF(C_ASSOCIATED(FieldTypePtr)) THEN
      CMISSFieldTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(FieldType)
      ALLOCATE(FieldType, STAT = Err)
      IF(Err/=0) THEN
        CMISSFieldTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSFieldTypeInitialise (FieldType, CMISSFieldTypeInitialiseC)
        FieldTypePtr = C_LOC(FieldType)
        CMISSFieldTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN
  END FUNCTION CMISSFieldTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Creates a pointer to a CMISSFieldsType object for an object reference for C.
  FUNCTION CMISSFieldsTypeCreateC(RegionPtr, FieldsPtr) BIND(C, NAME="CMISSFieldsTypeCreate")

  !Argument variables
  TYPE(C_PTR), INTENT(IN) :: RegionPtr !<C pointer to the region to get the fields from.
  TYPE(C_PTR), INTENT(INOUT) :: FieldsPtr !<C pointer to the fields attached to the specified region.
  !Function variable
  INTEGER(C_INT) :: CMISSFieldsTypeCreateC !<Error Code.
  !Local variables
  INTEGER(C_INT) :: Err
  TYPE(CMISSRegionType), POINTER :: Region
  TYPE(CMISSFieldsType), POINTER :: Fields

  IF(C_ASSOCIATED(RegionPtr)) THEN
    IF (C_ASSOCIATED(FieldsPtr)) THEN
      CMISSFieldsTypeCreateC = CMISSPointerNotNULL
    ELSE
      NULLIFY(Fields)
      ALLOCATE(Fields, STAT= Err)
      IF(Err/=0) THEN
        CMISSFieldsTypeCreateC = CMISSCouldNotAllocatePointer
      ELSE
        CALL C_F_POINTER (RegionPtr, Region)
        IF(ASSOCIATED(Region)) THEN
          CALL CMISSFieldsTypeCreate(Region, Fields, CMISSFieldsTypeCreateC)
          FieldsPtr = C_LOC(Fields)
        ELSE
          CMISSFieldsTypeCreateC = CMISSErrorConvertingPointer
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  RETURN

END FUNCTION CMISSFieldsTypeCreateC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSFieldsType object for C.

  FUNCTION CMISSFieldsTypeFinaliseC (FieldsTypePtr) BIND(C,NAME="CMISSFieldsTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: FieldsTypePtr !<C pointer to the CMISSFieldsType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldsTypeFinaliseC !<Error code.
    !Local variables
    TYPE(CMISSFieldsType), POINTER :: FieldsType

    CMISSFieldsTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(FieldsTypePtr)) THEN
      CALL C_F_POINTER(FieldsTypePtr, FieldsType)
      IF(ASSOCIATED(FieldsType)) THEN
        CALL CMISSFieldsTypeFinalise(FieldsType, CMISSFieldsTypeFinaliseC)
        DEALLOCATE(FieldsType)
        FieldsTypePtr = C_NULL_PTR
      ELSE
        CMISSFieldsTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldsTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldsTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSFieldsType object for C.

  FUNCTION CMISSFieldsTypeInitialiseC(FieldsTypePtr) BIND(C,NAME = "CMISSFieldsTypeInitialise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: FieldsTypePtr !<C pointer to the CMISSFieldsType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldsTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSFieldsType), POINTER :: FieldsType

    IF(C_ASSOCIATED(FieldsTypePtr)) THEN
      CMISSFieldsTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY (FieldsType)
      ALLOCATE(FieldsType, STAT = Err)
      IF(Err/=0) THEN
        CMISSFieldsTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSFieldsTypeInitialise(FieldsType, CMISSFieldsTypeInitialiseC)
        FieldsTypePtr = C_LOC(FieldsType)
        CMISSFieldsTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldsTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSGeneratedMeshType object for C.

  FUNCTION CMISSGeneratedMeshTypeFinaliseC(GeneratedMeshTypePtr) BIND(C, NAME = "CMISSGeneratedMeshTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: GeneratedMeshTypePtr !<C pointer to the CMISSGeneratedMeshType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMeshType

    CMISSGeneratedMeshTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshTypePtr)) THEN
      CALL C_F_POINTER(GeneratedMeshTypePtr, GeneratedMeshType)
      IF(ASSOCIATED(GeneratedMeshType)) THEN
        CALL CMISSGeneratedMeshTypeFinalise(GeneratedMeshType, CMISSGeneratedMeshTypeFinaliseC)
        DEALLOCATE(GeneratedMeshType)
        GeneratedMeshTypePtr = C_NULL_PTR
      ELSE
        CMISSGeneratedMeshTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSGeneratedMeshType object for C.

  FUNCTION CMISSGeneratedMeshTypeInitialiseC (GeneratedMeshTypePtr) BIND(C, NAME = "CMISSGeneratedMeshTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: GeneratedMeshTypePtr !<C pointer to the CMISSGeneratedMeshType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMeshType

    IF(C_ASSOCIATED(GeneratedMeshTypePtr)) THEN
      CMISSGeneratedMeshTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY (GeneratedMeshType)
      ALLOCATE(GeneratedMeshType, STAT = Err)
      IF(Err/=0) THEN
        CMISSGeneratedMeshTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSGeneratedMeshTypeInitialise(GeneratedMeshType, CMISSGeneratedMeshTypeInitialiseC)
        GeneratedMeshTypePtr = C_LOC(GeneratedMeshType)
        CMISSGeneratedMeshTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSHistoryType object for C.

  FUNCTION CMISSHistoryTypeFinaliseC(HistoryTypePtr) BIND(C, NAME = "CMISSHistoryTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: HistoryTypePtr !<C pointer to the CMISSHistoryType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSHistoryTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSHistoryType), POINTER :: HistoryType

    CMISSHistoryTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(HistoryTypePtr)) THEN
      CALL C_F_POINTER(HistoryTypePtr, HistoryType)
      IF(ASSOCIATED(HistoryType)) THEN
        CALL CMISSHistoryTypeFinalise(HistoryType, CMISSHistoryTypeFinaliseC)
        DEALLOCATE(HistoryType)
        HistoryTypePtr = C_NULL_PTR
      ELSE
        CMISSHistoryTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSHistoryTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSHistoryTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSHistoryType object for C.
  FUNCTION CMISSHistoryTypeInitialiseC(HistoryTypePtr) BIND(C,NAME="CMISSHistoryTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: HistoryTypePtr !<C pointer to the CMISSHistoryType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSHistoryTypeInitialiseC !<Error Code.
    !Local varaibles
    INTEGER(C_INT) :: Err
    TYPE(CMISSHistoryType), POINTER :: HistoryType

    IF(C_ASSOCIATED(HistoryTypePtr)) THEN
      CMISSHistoryTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(HistoryType)
      ALLOCATE(HistoryType, STAT=Err)
      IF(Err/=0) THEN
        CMISSHistoryTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSHistoryTypeInitialise(HistoryType, CMISSHistoryTypeInitialiseC)
        HistoryTypePtr = C_LOC(HistoryType)
        CMISSHistoryTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSHistoryTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSMeshType object for C.
  FUNCTION CMISSMeshTypeFinaliseC(MeshTypePtr) BIND(C,NAME="CMISSMeshTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: MeshTypePtr !<C pointer to the CMISSMeshType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSMeshType), POINTER :: MeshType

    CMISSMeshTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(MeshTypePtr))THEN
      CALL C_F_POINTER(MeshTypePtr, MeshType)
      IF(ASSOCIATED(MeshType)) THEN
        CALL CMISSMeshTypeFinalise(MeshType, CMISSMeshTypeFinaliseC)
        DEALLOCATE(MeshType)
        MeshTypePtr = C_NULL_PTR
      ELSE
        CMISSMeshTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN
  END FUNCTION CMISSMeshTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSMeshType object for C.
  FUNCTION CMISSMeshTypeInitialiseC (MeshTypePtr)BIND(C,NAME= "CMISSMeshTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: MeshTypePtr
    !Function variable
    INTEGER(C_INT) :: CMISSMeshTypeInitialiseC
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSMeshType), POINTER :: MeshType

    IF(C_ASSOCIATED(MeshTypePtr)) THEN
      CMISSMeshTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(MeshType)
      ALLOCATE(MeshType, STAT=Err)
      IF(Err/=0) THEN
        CMISSMeshTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSMeshTypeInitialise(MeshType, CMISSMeshTypeInitialiseC)
        MeshTypePtr = C_LOC(MeshType)
        CMISSMeshTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSMeshTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSMeshElementsType object for C.

  FUNCTION CMISSMeshElementsTypeFinaliseC(MeshElementsTypePtr) BIND(C, NAME="CMISSMeshElementsTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: MeshElementsTypePtr !<C pointer to the CMISSMeshElementsType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSMeshElementsType), POINTER :: MeshElementsType

    CMISSMeshElementsTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsTypePtr)) THEN
      CALL C_F_POINTER(MeshElementsTypePtr, MeshElementsType)
        IF(ASSOCIATED(MeshElementsType)) THEN
          CALL CMISSMeshElementsTypeFinalise(MeshElementsType, CMISSMeshElementsTypeFinaliseC)
          DEALLOCATE(MeshElementsType)
          MeshElementsTypePtr = C_NULL_PTR
      ELSE
        CMISSMeshElementsTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

      RETURN

    END FUNCTION CMISSMeshElementsTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSMeshElementsType object for C.
  FUNCTION CMISSMeshElementsTypeInitialiseC(MeshElementsTypePtr) BIND(C, NAME = "CMISSMeshElementsTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT (INOUT) :: MeshElementsTypePtr !<C pointer to the CMISSMeshElementsType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSMeshElementsType), POINTER :: MeshElementsType

    IF(C_ASSOCIATED(MeshElementsTypePtr)) THEN
      CMISSMeshElementsTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(MeshElementsType)
      ALLOCATE(MeshElementsType, STAT = Err)
      IF(Err/=0) THEN
        CMISSMeshElementsTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSMeshElementsTypeInitialise(MeshElementsType, CMISSMeshElementsTypeInitialiseC)
        MeshElementsTypePtr = C_LOC(MeshElementsType)
        CMISSMeshElementsTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSNodesType object for C.
  FUNCTION CMISSNodesTypeFinaliseC(NodesTypePtr) BIND(C, NAME="CMISSNodesTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: NodesTypePtr !<C pointer to the CMISSNodesType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSNodesType), POINTER :: NodesType

    CMISSNodesTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(NodesTypePtr)) THEN
      CALL C_F_POINTER(NodesTypePtr, NodesType)
      IF(ASSOCIATED(NodesType)) THEN
        CALL CMISSNodesTypeFinalise(NodesType, CMISSNodesTypeFinaliseC)
        DEALLOCATE(NodesType)
        NodesTypePtr = C_NULL_PTR
      ELSE
        CMISSNodesTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSNodesTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSNodesTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSNodesType object for C.
  FUNCTION CMISSNodesTypeInitialiseC(NodesTypePtr) BIND(C, NAME= "CMISSNodesTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: NodesTypePtr !<C pointer to the CMISSNodesType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSNodesType), POINTER :: NodesType

    IF(C_ASSOCIATED(NodesTypePtr)) THEN
      CMISSNodesTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(NodesType)
      ALLOCATE(NodesType, STAT = Err)
      IF(Err/=0) THEN
        CMISSNodesTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSNodesTypeInitialise(NodesType, CMISSNodesTypeInitialiseC)
        NodesTypePtr = C_NULL_PTR
        CMISSNodesTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSNodesTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSProblemType object for C.
  FUNCTION CMISSProblemTypeFinaliseC(ProblemTypePtr) BIND(C, NAME="CMISSProblemTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: ProblemTypePtr !<C pointer to the CMISSProblemType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSProblemType), POINTER :: ProblemType

    CMISSProblemTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(ProblemTypePtr)) THEN
      CALL C_F_POINTER(ProblemTypePtr, ProblemType)
      IF(ASSOCIATED(ProblemType)) THEN
        CALL CMISSProblemTypeFinalise(ProblemType, CMISSProblemTypeFinaliseC)
        DEALLOCATE(ProblemType)
        ProblemTypePtr = C_NULL_PTR
      ELSE
        CMISSProblemTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSProblemType object for C.
  FUNCTION CMISSProblemTypeInitialiseC(ProblemTypePtr) BIND(C, NAME = "CMISSProblemTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: ProblemTypePtr !<C pointer to the CMISSProblemType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSProblemType), POINTER :: ProblemType

    IF(C_ASSOCIATED(ProblemTypePtr)) THEN
      CMISSProblemTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(ProblemType)
      ALLOCATE(ProblemType, STAT=Err)
      IF(Err/=0) THEN
        CMISSProblemTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSProblemTypeInitialise(ProblemType, CMISSProblemTypeInitialiseC)
        ProblemTypePtr = C_LOC(ProblemType)
        CMISSProblemTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSProblemTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSQuadratureType object for C.
  FUNCTION CMISSQuadratureTypeFinaliseC(QuadratureTypePtr) BIND(C, NAME="CMISSQuadratureTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: QuadratureTypePtr !<C pointer to the CMISSQuadratureType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSQuadratureTypeFinaliseC !<Error Code.
    !Local variables
    TYPE(CMISSQuadratureType), POINTER :: QuadratureType

    CMISSQuadratureTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(QuadratureTypePtr)) THEN
      CALL C_F_POINTER(QuadratureTypePtr, QuadratureType)
      IF(ASSOCIATED(QuadratureType)) THEN
        CALL CMISSQuadratureTypeFinalise(QuadratureType, CMISSQuadratureTypeFinaliseC)
        DEALLOCATE(QuadratureType)
        QuadratureTypePtr = C_NULL_PTR
      ELSE
        CMISSQuadratureTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSQuadratureTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN
  END FUNCTION CMISSQuadratureTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSQuadratureType object for C.
  FUNCTION CMISSQuadratureTypeInitialiseC(QuadratureTypePtr) BIND(C,NAME="CMISSQuadratureTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: QuadratureTypePtr !<C pointer to the CMISSQuadratureType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSQuadratureTypeInitialiseC !<Error Code.
    !Local variable
    INTEGER(C_INT) :: Err
    TYPE(CMISSQuadratureType), POINTER :: QuadratureType

    IF(C_ASSOCIATED(QuadratureTypePtr)) THEN
      CMISSQuadratureTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(QuadratureType)
      ALLOCATE(QuadratureType, STAT=Err)
      IF(Err/=0) THEN
        CMISSQuadratureTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSQuadratureTypeInitialise(QuadratureType, CMISSQuadratureTypeInitialiseC)
        QuadratureTypePtr = C_LOC(QuadratureType)
        CMISSQuadratureTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

  END FUNCTION CMISSQuadratureTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSRegionType object.
  FUNCTION CMISSRegionTypeFinaliseC(RegionTypePtr) BIND(C, NAME="CMISSRegionTypeFinalise")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: RegionTypePtr !<C pointer to the CMISSRegionType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSRegionTypeFinaliseC !<Error code.
    !Local variable
    TYPE(CMISSRegionType), POINTER :: RegionType

    CMISSRegionTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(RegionTypePtr)) THEN
      CALL C_F_POINTER(RegionTypePtr, RegionType)
      IF(ASSOCIATED(RegionType)) THEN
        CALL CMISSRegionTypeFinalise(RegionType, CMISSRegionTypeFinaliseC)
        DEALLOCATE(RegionType)
        RegionTypePtr = C_NULL_PTR
      ELSE
        CMISSRegionTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSRegionTypeFinaliseC
 
  !
  !================================================================================================================================
  !

  !>Initialises a CMISSRegionType object for C.
  FUNCTION CMISSRegionTypeInitialiseC(RegionTypePtr) BIND(C,NAME="CMISSRegionTypeInitialise")
    
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: RegionTypePtr  !<C pointer to the CMISSRegionType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSRegionTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSRegionType), POINTER :: RegionType
    
    IF(C_ASSOCIATED(RegionTypePtr)) THEN
      CMISSRegionTypeInitialiseC=CMISSPointerNotNULL
    ELSE
      NULLIFY(RegionType)
      ALLOCATE(RegionType,STAT=Err)
      IF(Err/=0) THEN
        CMISSRegionTypeInitialiseC=CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSRegionTypeInitialise(RegionType,CMISSRegionTypeInitialiseC)
        RegionTypePtr=C_LOC(RegionType)
        CMISSRegionTypeInitialiseC=CMISSNoError
      ENDIF
    ENDIF

    RETURN
    
  END FUNCTION CMISSRegionTypeInitialiseC
  
  !
  !================================================================================================================================
  !

  !>Finalises a CMISSSolverType object for C.
  FUNCTION CMISSSolverTypeFinaliseC(SolverTypePtr) BIND(C, NAME="CMISSSolverTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: SolverTypePtr !<C pointer to the CMISSSolverType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSSolverTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSSolverType), POINTER :: SolverType

    CMISSSolverTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(SolverTypePtr)) THEN
      CALL C_F_POINTER(SolverTypePtr, SolverType)
      IF(ASSOCIATED(SolverType)) THEN
        CALL CMISSSolverTypeFinalise(SolverType, CMISSSolverTypeFinaliseC)
        DEALLOCATE(SolverType)
        SolverTypePtr = C_NULL_PTR
      ELSE
        CMISSSolverTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSSolverTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSSolverTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSSolverType object for C.
  FUNCTION CMISSSolverTypeInitialiseC(SolverTypePtr) BIND(C, NAME="CMISSSolverTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: SolverTypePtr !<C pointer to the CMISSSolverType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSSolverTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSSolverType), POINTER :: SolverType

    IF(C_ASSOCIATED(SolverTypePtr)) THEN
      CMISSSolverTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(SolverType)
      ALLOCATE(SolverType, STAT=Err)
      IF(Err/=0) THEN
        CMISSSolverTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSSolverTypeInitialise(SolverType, CMISSSolverTypeInitialiseC)
        SolverTypePtr = C_LOC(SolverType)
        CMISSSolverTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSSolverTypeInitialiseC

  !
  !================================================================================================================================
  !

  !>Finalises a CMISSSolverEquationsType object for C.
  FUNCTION CMISSSolverEquationsTypeFinaliseC (SolverEquationsTypePtr) BIND(C, NAME="CMISSSolverEquationsTypeFinalise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: SolverEquationsTypePtr !<C pointer to the CMISSSolverEquationsType object to be finalised.
    !Function variable
    INTEGER(C_INT) :: CMISSSolverEquationsTypeFinaliseC !<Error Code.
    !Local variable
    TYPE(CMISSSolverEquationsType), POINTER :: SolverEquationsType

    CMISSSolverEquationsTypeFinaliseC = CMISSNoError
    IF(C_ASSOCIATED(SolverEquationsTypePtr)) THEN
      CALL C_F_POINTER(SolverEquationsTypePtr, SolverEquationsType)
      IF(ASSOCIATED(SolverEquationsType)) THEN
        CALL CMISSSolverEquationsTypeFinalise(SolverEquationsType, CMISSSolverEquationsTypeFinaliseC)
        DEALLOCATE(SolverEquationsType)
        SolverEquationsTypePtr = C_NULL_PTR
      ELSE
        CMISSSolverEquationsTypeFinaliseC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSSolverEquationsTypeFinaliseC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSSolverEquationsTypeFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises a CMISSSolverEquationsType object for C.
  FUNCTION CMISSSolverEquationsTypeInitialiseC(SolverEquationsTypePtr) BIND(C, NAME="CMISSSolverEquationsTypeInitialise")

    !Argument variable
    TYPE(C_PTR), INTENT(INOUT) :: SolverEquationsTypePtr !<C pointer to the CMISSSolverEquationsType object to be initialised.
    !Function variable
    INTEGER(C_INT) :: CMISSSolverEquationsTypeInitialiseC !<Error Code.
    !Local variables
    INTEGER(C_INT) :: Err
    TYPE(CMISSSolverEquationsType), POINTER :: SolverEquationsType

    IF(C_ASSOCIATED(SolverEquationsTypePtr)) THEN
      CMISSSolverEquationsTypeInitialiseC = CMISSPointerNotNULL
    ELSE
      NULLIFY(SolverEquationsType)
      ALLOCATE(SolverEquationsType, STAT=Err)
      IF(Err/=0) THEN
        CMISSSolverEquationsTypeInitialiseC = CMISSCouldNotAllocatePointer
      ELSE
        CALL CMISSSolverEquationsTypeInitialise(SolverEquationsType, CMISSSolverEquationsTypeInitialiseC)
        SolverEquationsTypePtr = C_LOC(SolverEquationsType)
        CMISSSolverEquationsTypeInitialiseC = CMISSNoError
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSSolverEquationsTypeInitialiseC


!!==================================================================================================================================
!!
!! ANALYTIC_ANALYSIS_ROUTINES
!!
!!==================================================================================================================================

  !>Output the analytic error analysis for a field specified by a user number compared to the analytic values parameter set for C.
  FUNCTION CMISSAnalyticAnalysisOutputNumberC(RegionUserNumber, FieldUserNumber, FileNameSize, FileName) BIND(C, NAME = &
  & "CMISSAnalyticAnalysisOutputNumber")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to
    INTEGER(C_INT),VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to calculate the analytic error analysis for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FileNameSize !<The file name size, for C.
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: FileName(FileNameSize) !<If not empty, the filename to output the analytic analysis to, for C. If empty, the analysis will be output to the standard output.
    !Function variable
    INTEGER(C_INT) :: CMISSAnalyticAnalysisOutputNumberC !<Error Code.
    !Local variables
    CHARACTER(LEN=FileNameSize-1) :: FFileName

    CALL CMISSC2FString(FileName, FFileName)
    CALL CMISSAnalyticAnalysisOutputNumber(RegionUserNumber, FieldUserNumber, FFilename, CMISSAnalyticAnalysisOutputNumberC)

    RETURN

  END FUNCTION CMISSAnalyticAnalysisOutputNumberC

  !
  !================================================================================================================================
  !

  !>Output the analytic error analysis for a field identified by an object compared to the analytic values parameter set for C.
  FUNCTION CMISSAnalyticAnalysisOutputPtrC(FieldPtr,FileNameSize, FileName) BIND(C, NAME = "CMISSAnalyticAnalysisOutputObj")

    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: FieldPtr !<The dependent field to calculate the analytic error analysis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FileNameSize !<The file name size, for C.
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: FileName(FileNameSize) !<If not empty, the filename to output the analytic analysis to, for C. If empty, the analysis will be output to the standard output.
    !Function variable
    INTEGER(C_INT) :: CMISSAnalyticAnalysisOutputPtrC !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN = FileNameSize-1) :: FFileName

    CMISSAnalyticAnalysisOutputPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER (FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSC2FString(Filename, FFileName)
        CALL CMISSAnalyticAnalysisOutputObj(Field, FFileName, CMISSAnalyticAnalysisOutputPtrC)
      ELSE
        CMISSAnalyticAnalysisOutputPtrC=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSAnalyticAnalysisOutputPtrC=CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSAnalyticAnalysisOutputPtrC

  !
  !================================================================================================================================
  !

  !>Get absolute error value for the node in a field specified by a user number compared to the analytic value for C.
  FUNCTION CMISSAnalyticAnalysisNodeAbsoluteErrorGetCNum(RegionUserNumber,FieldUserNumber,DerivativeNumber,NodeNumber, &
    & ComponentNumber,VariableNumber,Value) BIND(C, NAME = "CMISSAnalyticAnalysisNodeAbsoluteErrorGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field for analytic error analysis.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to calculate the analytic error analysis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<Derivative number of the field to calculate the analytic error analysis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NodeNumber !<Node number of the field to calculate the analytic error analysis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<Component number of the field to calculate the analytic error analysis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableNumber !<Variable number of the field to calculate the analytic error analysis for, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The absolute error of the field to calculate the analytic error analysis for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSAnalyticAnalysisNodeAbsoluteErrorGetCNum
    !Local variable

    CALL CMISSAnalyticAnalysisNodeAbsoluteErrorGet(RegionUserNumber,FieldUserNumber,DerivativeNumber,NodeNumber, &
      & ComponentNumber,VariableNumber,Value, CMISSAnalyticAnalysisNodeAbsoluteErrorGetCNum)

    RETURN

  END FUNCTION CMISSAnalyticAnalysisNodeAbsoluteErrorGetCNum

  !
  !================================================================================================================================
  !

  !>Get absolute error value for the node in a field identified by an object compared to the analytic value for C.
  FUNCTION CMISSAnalyticAnalysisNodeAbsoluteErrorGetPtrC(FieldPtr,DerivativeNumber,NodeNumber,ComponentNumber,VariableNumber, &
    & Value) BIND(C, NAME = "CMISSAnalyticAnalysisNodeAbsoluteErrorGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the dependent field to calculate the analytic error analysis for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<Derivative number of the field to calculate the analytic error analysis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NodeNumber !<Node number of the field to calculate the analytic error analysis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<Component number of the field to calculate the analytic error analysis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableNumber !<Variable number of the field to calculate the analytic error analysis for, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The absolute error of the field to calculate the analytic error analysis for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSAnalyticAnalysisNodeAbsoluteErrorGetPtrC
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSAnalyticAnalysisNodeAbsoluteErrorGetPtrC = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSAnalyticAnalysisNodeAbsoluteErrorGet(Field,DerivativeNumber,NodeNumber,ComponentNumber,VariableNumber, &
          & Value, CMISSAnalyticAnalysisNodeAbsoluteErrorGetPtrC)
      ELSE
        CMISSAnalyticAnalysisNodeAbsoluteErrorGetPtrC = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSAnalyticAnalysisNodeAbsoluteErrorGetPtrC = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSAnalyticAnalysisNodeAbsoluteErrorGetPtrC


 !stop


  !
  !================================================================================================================================
  !

  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************
  !***********************************************************************************************************




!!==================================================================================================================================
!!
!! CMISS_ROUTINES
!!
!!==================================================================================================================================

  !>Finalises CMISS for C.
  FUNCTION CMISSFinaliseC() BIND(C,NAME="CMISSFinalise")
  
    !Argument variables
    !Function variable
    INTEGER(C_INT) :: CMISSFinaliseC
    !Local variables

    CALL CMISSFinalise(CMISSFinaliseC)

    RETURN
    
  END FUNCTION CMISSFinaliseC

  !
  !================================================================================================================================
  !

  !>Initialises CMISS for C returning a user number to the world coordinate system and region.
  FUNCTION CMISSInitialiseCNum(WorldCoordinateSystemUserNumber,WorldRegionUserNumber) BIND(C,NAME="CMISSInitialiseNum")
  
    !Argument variables
    INTEGER(C_INT), INTENT(OUT) :: WorldCoordinateSystemUserNumber
    INTEGER(C_INT), INTENT(OUT) :: WorldRegionUserNumber
    !Function variable
    INTEGER(C_INT) :: CMISSInitialiseCNum
    !Local variables
  
    CALL CMISSInitialise(WorldCoordinateSystemUserNumber,WorldRegionUserNumber,CMISSInitialiseCNum)

    RETURN
    
    END FUNCTION CMISSInitialiseCNum

  !
  !================================================================================================================================
  !

  !>Initialises CMISS for C returning pointers to the world coordinate system and region.
  FUNCTION CMISSInitialiseCPtr(WorldCoordinateSystemPtr,WorldRegionPtr) BIND(C,NAME="CMISSInitialise")
  
    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: WorldCoordinateSystemPtr
    TYPE(C_PTR), INTENT(INOUT) :: WorldRegionPtr
    !Function variable
    INTEGER(C_INT) :: CMISSInitialiseCPtr
    !Local variables
    TYPE(CMISSCoordinateSystemType), POINTER :: WorldCoordinateSystem
    TYPE(CMISSRegionType), POINTER :: WorldRegion

    CMISSInitialiseCPtr=CMISSCoordinateSystemTypeInitialiseC(WorldCoordinateSystemPtr)
    IF(CMISSInitialiseCPtr==CMISSNoError) THEN
      CMISSInitialiseCPtr=CMISSRegionTypeInitialiseC(WorldRegionPtr)
      IF(CMISSInitialiseCPtr==CMISSNoError) THEN
        CALL C_F_POINTER(WorldCoordinateSystemPtr,WorldCoordinateSystem)
        IF(ASSOCIATED(WorldCoordinateSystem)) THEN
          CALL C_F_POINTER(WorldRegionPtr,WorldRegion)
          IF(ASSOCIATED(WorldRegion)) THEN        
            CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,CMISSInitialiseCPtr)
          ELSE
            CMISSInitialiseCPtr=CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSInitialiseCPtr=CMISSErrorConvertingPointer
        ENDIF
      ENDIF
    ENDIF

    RETURN
    
  END FUNCTION CMISSInitialiseCPtr

!!==================================================================================================================================
!!
!! FIELD_ROUTINES
!!
!!==================================================================================================================================

  !>Returns the interpolation type for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentInterpolationGetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & InterpolationType) BIND(C,NAME = "CMISSFieldComponentInterpolationGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the interpolation type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the interpolation type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the interpolation type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the interpolation type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) ::  CMISSFieldComponentInterpolationGetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldComponentInterpolationGet(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & InterpolationType, CMISSFieldComponentInterpolationGetCNum)

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationGetCNum

  !!
  !!==================================================================================================================================
  !!

  !>Returns the interpolation type for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentInterpolationGetCPtr(FieldPtr,VariableType,ComponentNumber,InterpolationType) BIND(C, NAME = &
  & "CMISSFieldComponentInterpolationGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the interpolation type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the interpolation type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the interpolation type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentInterpolationGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentInterpolationGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentInterpolationGet(Field, VariableType,ComponentNumber,InterpolationType, &
        & CMISSFieldComponentInterpolationGetCPtr)
      ELSE
        CMISSFieldComponentInterpolationGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentInterpolationGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationGetCPtr

  !!
  !!==================================================================================================================================
  !!

  !>Sets/changes the interpolation type for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentInterpolationSetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & InterpolationType) BIND(C, NAME = "CMISSFieldComponentInterpolationSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the interpolation type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentInterpolationSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentInterpolationSet(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,InterpolationType,&
    & CMISSFieldComponentInterpolationSetCNum)

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationSetCNum


  !!
  !!==================================================================================================================================
  !!

  !>Sets/changes the interpolation type for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentInterpolationSetCPtr(FieldPtr,VariableType,ComponentNumber,InterpolationType) &
   & BIND(C, NAME = "CMISSFieldComponentInterpolationSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the interpolation type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the interpolation type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: InterpolationType !<The interpolation type for C. \see OPENCMISS_FieldInterpolationTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentInterpolationSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentInterpolationSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentInterpolationSet(Field, VariableType, ComponentNumber,InterpolationType, &
        & CMISSFieldComponentInterpolationSetCPtr)
      ELSE
        CMISSFieldComponentInterpolationSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentInterpolationSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentInterpolationSetCPtr


  !!
  !!==================================================================================================================================
  !!

  !>Returns the character string label for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentLabelGetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,LabelSize,Label) &
  & BIND(C, NAME = "CMISSFieldComponentLabelGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !< Label size
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable component character string label to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelGetCNum !<Error Code.
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel

    CALL CMISSFieldComponentLabelGet(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,FLabel, &
    & CMISSFieldComponentLabelGetCNum)
    CALL CMISSF2CString(Flabel,Label)

    RETURN

  END FUNCTION CMISSFieldComponentLabelGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentLabelGetCPtr(FieldPtr,VariableType,ComponentNumber,LabelSize,Label) BIND(C, NAME = &
  & "CMISSFieldComponentLabelGet")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the label for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !< Label size
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable component character string label to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelGetCPtr !<Error Code.
    !Local variable
    CHARACTER(LEN=LabelSize-1) :: FLabel
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentLabelGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentLabelGet(Field,VariableType,ComponentNumber,FLabel,CMISSFieldComponentLabelGetCPtr)
        CALL CMISSF2CString(FLabel, Label)
      ELSE
        CMISSFieldComponentLabelGetCPtr=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentLabelGetCPtr=CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentLabelGetCPtr

  !
  !================================================================================================================================
  !
  !>Sets/changes the character string label for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentLabelSetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,LabelSize,Label)&
  & BIND(C,NAME= "CMISSFieldComponentLabelSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the label to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !< Label size
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field variable component character string label to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelSetCNum !<Error Code.
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel

    CALL CMISSC2FString(Label,Flabel)
    CALL CMISSFieldComponentLabelSet(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,FLabel, &
    & CMISSFieldComponentLabelSetCNum)

    RETURN

  END FUNCTION CMISSFieldComponentLabelSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the character string label for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentLabelSetCPtr(FieldPtr,VariableType,ComponentNumber,LabelSize,Label) BIND(C, NAME= &
  & "CMISSFieldComponentLabelSetC")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the label to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the label to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Label Size
    CHARACTER(LEN = 1, KIND=C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field variable component character string label to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentLabelSetCPtr !<Error Code.
    !Local variable
    CHARACTER(LEN=LabelSize-1) :: FLabel
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentLabelSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSC2FString(Label, FLabel)
        CALL CMISSFieldComponentLabelSet(Field,VariableType,ComponentNumber,FLabel,CMISSFieldComponentLabelSetCPtr)
      ELSE
        CMISSFieldComponentLabelSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentLabelSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN
  END FUNCTION CMISSFieldComponentLabelSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number for a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentMeshComponentGetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & MeshComponent) BIND(C, NAME = "CMISSFieldComponentMeshComponentGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the mesh component number for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the mesh component number for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the mesh component number for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the mesh component number for, for C.
    INTEGER(C_INT), INTENT(OUT) :: MeshComponent !<The mesh component number to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentMeshComponentGet(RegionUserNumber, FieldUserNumber, VariableType, ComponentNumber, MeshComponent,&
    & CMISSFieldComponentMeshComponentGetCNum)

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentMeshComponentGetCPtr(FieldPtr,VariableType,ComponentNumber,MeshComponent) BIND(C, &
  & NAME = "CMISSFieldComponentMeshComponentGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the mesh component number for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the mesh component number for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the mesh component number for, for C.
    INTEGER(C_INT), INTENT(OUT) :: MeshComponent !<The mesh component number to get, for C.
    !Function Variables
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentMeshComponentGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentMeshComponentGet(Field, VariableType, ComponentNumber, MeshComponent, &
        & CMISSFieldComponentMeshComponentGetCPtr)
      ELSE
        CMISSFieldComponentMeshComponentGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentMeshComponentGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentGetCPtr

  !
  !================================================================================================================================
  !
  !>Sets/changes the mesh component number for a field variable component for a field identified by a user number.
  FUNCTION CMISSFieldComponentMeshComponentSetCNum(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber, &
    & MeshComponent) BIND(C, NAME = "CMISSFieldComponentMeshComponentSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the mesh component number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the mesh component number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the mesh component number to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the mesh component number to, for C.
    INTEGER(C_INT), INTENT(IN) :: MeshComponent !<The mesh component number to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentMeshComponentSet(RegionUserNumber,FieldUserNumber,VariableType,ComponentNumber,MeshComponent, &
    & CMISSFieldComponentMeshComponentSetCNum)

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number for a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldComponentMeshComponentSetCPtr(FieldPtr,VariableType,ComponentNumber,MeshComponent) BIND(C, &
  & NAME = "CMISSFieldComponentMeshComponentSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the mesh component number to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the mesh component number to. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to set the mesh component number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponent !<The mesh component number to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentMeshComponentSetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentMeshComponentSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentMeshComponentSet(Field,VariableType,ComponentNumber,MeshComponent, &
        & CMISSFieldComponentMeshComponentSetCPtr)
      ELSE
        CMISSFieldComponentMeshComponentSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentMeshComponentSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentMeshComponentSetCPtr
  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to an integer constant value for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentValuesInitialiseIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseIntgNum")

        !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseIntgCNum !<Error Code.
    !Local variable

    CALL CMISSFieldComponentValuesInitialiseIntg(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, &
    & ComponentNumber, Value, CMISSFieldComponentValuesInitialiseIntgCNum)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseIntgCNum

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to an integer constant value for a field identified by an object for C.
  FUNCTION CMISSFieldComponentValuesInitialiseIntgCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) &
  & BIND(C, NAME ="CMISSFieldComponentValuesInitialiseIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseIntg(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseIntgCPtr)
      ELSE
        CMISSFieldComponentValuesInitialiseIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseIntgCPtr

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a single precision constant value for a field identified by a user number for C.
  FUNCTION CMISSFieldComponentValuesInitialiseSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for , for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentValuesInitialiseSP(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, &
    & ComponentNumber, Value, CMISSFieldComponentValuesInitialiseSPCNum)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseSPCNum


  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a single precision constant value for a field identified by an object, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseSPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
  & "CMISSFieldComponentValuesInitialiseSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseSPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseSP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseSPCPtr)
      ELSE
        CMISSFieldComponentValuesInitialiseSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseSPCPtr

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a double precision constant value for a field identified by a user number, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentValuesInitialiseDP(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType,&
    & ComponentNumber, Value, CMISSFieldComponentValuesInitialiseDPCNum)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseDPCNum

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a double precision constant value for a field identified by an object, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseDPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldComponentValuesInitialiseDP")

      !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseDPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseDP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseDPCPtr)
      ELSE
        CMISSFieldComponentValuesInitialiseDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseDPCPtr

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a logical constant value for a field identified by a user number, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldComponentValuesInitialiseLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to initialise the field variable component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldComponentValuesInitialiseL(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, ComponentNumber,&
    & Value, CMISSFieldComponentValuesInitialiseLCNum)

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseLCNum

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a logical constant value for a field identified by an object, for C.
  FUNCTION CMISSFieldComponentValuesInitialiseLCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldComponentValuesInitialiseL")

    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to initialise the field variable component for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to initialise the field variable component for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The value to initialise the parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldComponentValuesInitialiseLCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldComponentValuesInitialiseLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldComponentValuesInitialiseL(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldComponentValuesInitialiseLCPtr)
      ELSE
        CMISSFieldComponentValuesInitialiseLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldComponentValuesInitialiseLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldComponentValuesInitialiseLCPtr

  !
  !================================================================================================================================
  !

  !>Returns the data type for a field variable for a field identified by a user number, for C.
  FUNCTION CMISSFieldDataTypeGetCNum(RegionUserNumber,FieldUserNumber,VariableType,DataType) BIND(C,NAME= &
  & "CMISSFieldDataTypeGetNum")
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the data type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the data type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the data type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DataType !<The field variable data type to get, for C. \see OPENCMISS_FieldDataTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDataTypeGetCNum !<Error Code.
    !Local variables

    CALL  CMISSFieldDataTypeGet(RegionUserNumber,FieldUserNumber,VariableType,DataType,CMISSFieldDataTypeGetCNum)

    RETURN

  END FUNCTION  CMISSFieldDataTypeGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the data type for a field variable for a field identified by an object, for C.
  FUNCTION CMISSFieldDataTypeGetCPtr(FieldPtr,VariableType,DataType) BIND(C, NAME = "CMISSFieldDataTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<The field to get the data type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the data type for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DataType !<The field variable data type to get, for C. \see OPENCMISS_FieldDataTypes
    !Function Variables
    INTEGER(C_INT) :: CMISSFieldDataTypeGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDataTypeGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDataTypeGet(Field, VariableType, DataType, CMISSFieldDataTypeGetCPtr)
      ELSE
        CMISSFieldDataTypeGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDataTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDataTypeGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a field variable for a field identified by a user number, for C.
  FUNCTION CMISSFieldDataTypeSetCNum(RegionUserNumber,FieldUserNumber,VariableType,DataType) BIND(C, NAME = &
  & "CMISSFieldDataTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the data type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the data type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the data type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DataType !<The field variable data type to set, for C. \see OPENCMISS_FieldDataTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDataTypeSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldDataTypeSet(RegionUserNumber, FieldUserNumber,VariableType, DataType,CMISSFieldDataTypeSetCNum)

    RETURN

  END FUNCTION CMISSFieldDataTypeSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a field variable for a field identified by an object, for C.
  FUNCTION CMISSFieldDataTypeSetCPtr(FieldPtr,VariableType,DataType) BIND(C, NAME = "CMISSFieldDataTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the data type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to set the data type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DataType !<The field variable data type to set, for C. \see OPENCMISS_FieldDataTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDataTypeSetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDataTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDataTypeSet(Field, VariableType, DataType,CMISSFieldDataTypeSetCPtr)
      ELSE
        CMISSFieldDataTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDataTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDataTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the DOF order type for a field variable for a field identified by a user number, for C.
  FUNCTION CMISSFieldDOFOrderTypeGetCNum(RegionUserNumber,FieldUserNumber,VariableType,DOFOrderType) BIND( &
  & C, NAME = "CMISSFieldDOFOrderTypeGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the DOF order type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the DOF order type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the DOF order type for, for C.\see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DOFOrderType !<The field variable DOF Order type to get, for C.  \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeGetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldDOFOrderTypeGet(RegionUserNumber, FieldUserNumber, VariableType, DOFOrderType, &
    & CMISSFieldDOFOrderTypeGetCNum)

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the DOF Order type for a field variable for a field identified by an object, for C.
  FUNCTION CMISSFieldDOFOrderTypeGetCPtr(FieldPtr,VariableType,DOFOrderType) BIND(C, NAME = "CMISSFieldDOFOrderTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer of the field to get the DOF Order type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to get the DOF Order type for, for C \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: DOFOrderType !<The field variable DOF Order type to get, for C. \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeGetCPtr
 !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDOFOrderTypeGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDOFOrderTypeGet(Field, VariableType, DOFOrderType,CMISSFieldDOFOrderTypeGetCPtr)
      ELSE
        CMISSFieldDOFOrderTypeGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDOFOrderTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the DOF order type for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldDOFOrderTypeSetCNum(RegionUserNumber,FieldUserNumber,VariableType,DOFOrderType) BIND(C, &
  & NAME = "CMISSFieldDOFOrderTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to set the DOF Order type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the DOF Order type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the DOF Order type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DOFOrderType !<The field variable DOF Order type to set, for C. \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeSetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldDOFOrderTypeSet(RegionUserNumber, FieldUserNumber, VariableType, DOFOrderType, &
    & CMISSFieldDOFOrderTypeSetCNum)

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the DOF Order type for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldDOFOrderTypeSetCPtr(FieldPtr,VariableType,DOFOrderType) BIND(C, NAME = "CMISSFieldDOFOrderTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the DOF Order type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the DOF Order type to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DOFOrderType !<The field variable DOF Order type to set, for C. \see OPENCMISS_FieldDOFOrderTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDOFOrderTypeSetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDOFOrderTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDOFOrderTypeSet(Field, VariableType, DOFOrderType, CMISSFieldDOFOrderTypeSetCPtr)
      ELSE
        CMISSFieldDOFOrderTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDOFOrderTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDOFOrderTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a field identified by a user number for C.
  FUNCTION CMISSFieldCreateFinishCNum(RegionUserNumber,FieldUserNumber) BIND(C, NAME = "CMISSFieldCreateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to finish the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to finish the creation of, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldCreateFinishCNum !<Error Code.

    CALL CMISSFieldCreateFinish(RegionUserNumber,FieldUserNumber,CMISSFieldCreateFinishCNum)

    RETURN

  END FUNCTION CMISSFieldCreateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a field identified by an object for C.
  FUNCTION CMISSFieldCreateFinishCPtr(FieldPtr) BIND(C, NAME = "CMISSFieldCreateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finish the creation of.
    !Function variable
    INTEGER(C_INT) ::CMISSFieldCreateFinishCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldCreateFinish(Field,CMISSFieldCreateFinishCPtr)
      ELSE
        CMISSFieldCreateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts the creation of a field identified by a user number for C.
  FUNCTION CMISSFieldCreateStartCNum(FieldUserNumber,RegionUserNumber) BIND(C, NAME = "CMISFieldCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber!<The user number of the region containing the field to start the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the field to start the creation of, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldCreateStartCNum !<Error Code.
    !Local variable

    CALL CMISSFieldCreateStart(FieldUserNumber, RegionUserNumber, CMISSFieldCreateStartCNum)

    RETURN

  END FUNCTION CMISSFieldCreateStartCNum


  !
  !================================================================================================================================
  !

  !>Starts the creation of a field identified by an object for C.
  FUNCTION CMISSFieldCreateStartCPtr(FieldUserNumber,RegionPtr,FieldPtr) BIND(C, NAME ="CMISSFieldCreateStart")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber  !<The user number of the field to start the creation of, for C.
    TYPE(C_PTR), INTENT(IN) :: RegionPtr !<C pointer to the region to create the field on.
    TYPE(C_PTR), INTENT(IN) :: FieldPtr !<C pointer to the created field.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldCreateStartCPtr !<Error Code.
    !Local variable
    TYPE(CMISSRegionType), POINTER :: Region
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr, Region)
      IF(ASSOCIATED(Region)) THEN
        IF(C_ASSOCIATED(FieldPtr)) THEN
          CALL C_F_POINTER(FieldPtr, Field)
          IF(ASSOCIATED(Field)) THEN
            CALL CMISSFieldCreateStart(FieldUserNumber, Region, Field, CMISSFieldCreateStartCPtr)
          ELSE
            CMISSFieldCreateStartCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldCreateStartCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldCreateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldCreateStartCPtr


  !
  !================================================================================================================================
  !

  !>Returns the dependent type for a field identified by a user number for C.
  FUNCTION CMISSFieldDependentTypeGetCNum(RegionUserNumber,FieldUserNumber,DependentType) BIND(C, NAME = &
  & "CMISSFieldDependentTypeGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to get the dependent type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to get the dependent type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: DependentType !<The field dependent type to get, for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldDependentTypeGet(RegionUserNumber, FieldUserNumber, DependentType, CMISSFieldDependentTypeGetCNum)

    RETURN

  END FUNCTION CMISSFieldDependentTypeGetCNum


  !
  !================================================================================================================================
  !

  !>Returns the dependent type for a field identified by an object for C.
  FUNCTION CMISSFieldDependentTypeGetCPtr(FieldPtr,DependentType) BIND(C, NAME = "CMISSFieldDependentTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the dependent type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: DependentType !<The field dependent type for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeGetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDependentTypeGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF (ASSOCIATED(Field)) THEN
        CALL CMISSFieldDependentTypeGet(Field, DependentType, CMISSFieldDependentTypeGetCPtr)
      ELSE
        CMISSFieldDependentTypeGetCPtr= CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDependentTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDependentTypeGetCPtr


  !
  !================================================================================================================================
  !

  !>Sets/changes the dependent type for a field identified by a user number for C.
  FUNCTION CMISSFieldDependentTypeSetCNum(RegionUserNumber,FieldUserNumber,DependentType) BIND(C, NAME  = &
  & "CMISSFieldDependentTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to set the dependent type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the dependent type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DependentType !<The field dependent type for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldDependentTypeSet(RegionUserNumber,FieldUserNumber,DependentType,CMISSFieldDependentTypeSetCNum)

    RETURN

  END FUNCTION CMISSFieldDependentTypeSetCNum


  !
  !================================================================================================================================
  !

  !>Sets/changes the dependent type for a field identified by an object for C.
  FUNCTION CMISSFieldDependentTypeSetCPtr(FieldPtr,DependentType) BIND(C, NAME = "CMISSFieldDependentTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the dependent type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DependentType !<The field dependent type, for C. \see OPENCMISS_FieldDependentTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDependentTypeSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDependentTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDependentTypeSet(Field, DependentType, CMISSFieldDependentTypeSetCPtr)
      ELSE
        CMISSFieldDependentTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDependentTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDependentTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Destroys a field identified by a user number for C.
  FUNCTION CMISSFieldDestroyCNum(RegionUserNumber,FieldUserNumber) BIND(C, NAME = "CMISSFieldDestroyNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to destroy for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to destroy for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDestroyCNum !<Error Code.
    !Local variable

    CALL CMISSFieldDestroy(RegionUserNumber, FieldUserNumber, CMISSFieldDestroyCNum)

    RETURN

  END FUNCTION CMISSFieldDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys a field identified by an object for C.
  FUNCTION CMISSFieldDestroyCPtr(FieldPtr) BIND(C,NAME = "CMISSFieldDestroy")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to destroy for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDestroyCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF (ASSOCIATED(Field)) THEN
        CALL CMISSFieldDestroy(Field, CMISSFieldDestroyCPtr)
      ELSE
        CMISSFieldDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Returns the dimension for a field identified by a user number for C.
  FUNCTION CMISSFieldDimensionGetCNum(RegionUserNumber,FieldUserNumber,VariableType,DIMENSION) BIND(C, NAME = &
  & "CMISSFieldDimensionGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the dimension for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the dimension for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the dimension for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: Dimension !<The field dimension for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionGetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldDimensionGet(RegionUserNumber, FieldUserNumber, VariableType, DIMENSION, CMISSFieldDimensionGetCNum)

    RETURN

  END FUNCTION CMISSFieldDimensionGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the dimension for a field identified by an object for C.
  FUNCTION CMISSFieldDimensionGetCPtr(FieldPtr,VariableType,Dimension) BIND(C, NAME = "CMISSFieldDimensionGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the dimension for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the dimension for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: Dimension !<The field dimension for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionGetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDimensionGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDimensionGet(Field, VariableType, Dimension, CMISSFieldDimensionGetCPtr)
      ELSE
        CMISSFieldDimensionGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDimensionGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDimensionGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the dimension for a field identified by a user number for C.
  FUNCTION CMISSFieldDimensionSetCNum(RegionUserNumber,FieldUserNumber,VariableType,Dimension) BIND(C, NAME = &
  & "CMISSFieldDimensionSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the dimension to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the dimension to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the dimension to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: Dimension !<The field dimension for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionSetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldDimensionSet(RegionUserNumber, FieldUserNumber, VariableType, Dimension, CMISSFieldDimensionSetCNum)

    RETURN

  END FUNCTION CMISSFieldDimensionSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the dimension for a field identified by an object for C.
  FUNCTION CMISSFieldDimensionSetCPtr(FieldPtr,VariableType,Dimension) BIND(C, NAME= "CMISSFieldDimensionSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the dimension to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the dimension to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: Dimension !<The field dimension, for C. \see OPENCMISS_FieldDimensionTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldDimensionSetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldDimensionSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldDimensionSet(Field, VariableType, Dimension, CMISSFieldDimensionSetCPtr)
      ELSE
        CMISSFieldDimensionSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldDimensionSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldDimensionSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the geometric field for a field identified by a user number for C.
  FUNCTION CMISSFieldGeometricFieldGetCNum(RegionUserNumber,FieldUserNumber,GeometricFieldUserNumber) BIND(C, NAME = &
  & "CMISSFieldGeometricFieldGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the geometric field for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the geometric field for, for C.
    INTEGER(C_INT), INTENT(OUT) :: GeometricFieldUserNumber !<The field geometric field user number, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldGeometricFieldGet(RegionUserNumber, FieldUserNumber, GeometricFieldUserNumber, &
    & CMISSFieldGeometricFieldGetCNum)

    RETURN

  END FUNCTION CMISSFieldGeometricFieldGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the geometric field for a field identified by an object for C.
  FUNCTION CMISSFieldGeometricFieldGetCPtr(FieldPtr,GeometricFieldPtr) BIND(C, NAME = "CMISSFieldGeometricFieldGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the geometric field to, for C.
    TYPE(C_PTR), INTENT(OUT) :: GeometricFieldPtr !<C pointer to the geometric field for the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSFieldType), POINTER :: GeometricField

    CMISSFieldGeometricFieldGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldGeometricFieldGet(Field, GeometricField, CMISSFieldGeometricFieldGetCPtr)
        IF(ASSOCIATED(GeometricField)) THEN
          GeometricFieldPtr = C_LOC(GeometricField)
        ELSE
          CMISSFieldGeometricFieldGetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldGeometricFieldGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldGeometricFieldGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldGeometricFieldGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the geometric field for a field identified by a user number for C.
  FUNCTION CMISSFieldGeometricFieldSetCNum(RegionUserNumber,FieldUserNumber,GeometricFieldUserNumber) BIND(C, &
  & NAME = "CMISSFieldGeometricFieldSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region corresponding to the field to set the geometric field to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !< The user number for the field to set the geometric field to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeometricFieldUserNumber !<The field geometric field user number to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldSetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldGeometricFieldSetNum(RegionUserNumber,FieldUserNumber,GeometricFieldUserNumber,&
    & CMISSFieldGeometricFieldSetCNum)

    RETURN

  END FUNCTION CMISSFieldGeometricFieldSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the geometric field for a field identified by an object for C.
  FUNCTION CMISSFieldGeometricFieldSetCPtr(FieldPtr,GeometricFieldPtr) BIND(C, NAME = "CMISSFieldGeometricFieldSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the geometric field to, for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeometricFieldPtr !<C pointer to the geometric field for the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldGeometricFieldSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSFieldType), POINTER :: GeometricField

    CMISSFieldGeometricFieldSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(GeometricFieldPtr)) THEN
          CALL C_F_POINTER(GeometricFieldPtr, GeometricField)
          IF(ASSOCIATED(GeometricField)) THEN
            CALL CMISSFieldGeometricFieldSet(Field, GeometricField, CMISSFieldGeometricFieldSetCPtr)
          ELSE
            CMISSFieldGeometricFieldSetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldGeometricFieldSetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldGeometricFieldSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
     CMISSFieldGeometricFieldSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldGeometricFieldSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field identified by a user number for C.
  FUNCTION CMISSFieldLabelGetCNum(RegionUserNumber,FieldUserNumber,LabelSize,Label) BIND(C, NAME = "CMISSFieldLabelGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Label Size
    CHARACTER(LEN = 1, KIND = C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelGetCNum !<Error Code.
    !Local variables
    CHARACTER(LEN = LabelSize - 1) :: FLabel

    CALL CMISSFieldLabelGet(RegionUserNumber,FieldUserNumber,FLabel,CMISSFieldLabelGetCNum)
    CALL CMISSF2CString(FLabel, Label)

    RETURN

  END FUNCTION CMISSFieldLabelGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field identified by an object for C.
  FUNCTION CMISSFieldLabelGetCPtr(FieldPtr,LabelSize,Label) BIND(C, NAME = "CMISSFieldLabelGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Label size
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN = LabelSize -1) :: FLabel

    CMISSFieldLabelGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldLabelGet(Field, FLabel, CMISSFieldLabelGetCPtr)
        CALL CMISSF2CString(FLabel, Label)
      ELSE
        CMISSFieldLabelGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldLabelGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldLabelGetCPtr

  !!
  !================================================================================================================================
  !!

  !>Sets/changes the character string label for a field identified by a user number for C.
  FUNCTION CMISSFieldLabelSetCNum(RegionUserNumber,FieldUserNumber,LabelSize, Label) BIND(C, NAME = "CMISSFieldLabelSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size.
    CHARACTER(LEN = 1, KIND= C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelSetCNum !<Error Code.
    !Local variables
    CHARACTER(LEN = LabelSize-1) :: FLabel

    CALL CMISSC2FString (Label, FLabel)
    CALL CMISSFieldLabelSet(RegionUserNumber, FieldUserNumber, FLabel,CMISSFieldLabelSetCNum)

    RETURN

  END FUNCTION CMISSFieldLabelSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the character string label for a field identified by an object for C.
  FUNCTION CMISSFieldLabelSetCPtr(FieldPtr,LabelSize, Label) BIND(C, NAME = "CMISSFieldLabelSet")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field character string label for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldLabelSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN=LabelSize-1) :: FLabel

    CMISSFieldLabelSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSC2FString(Label, FLabel)
        CALL CMISSFieldLabelSet(Field,FLabel,CMISSFieldLabelSetCPtr)
      ELSE
        CMISSFieldLabelSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldLabelSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldLabelSetCPtr


  !
  !================================================================================================================================
  !

  !>Returns the mesh decomposition for a field identified by a user number for C.

  FUNCTION CMISSFieldMeshDecompositionGetCNum(RegionUserNumber,FieldUserNumber,DecompositionUserNumber) &
    & BIND(C, NAME = "CMISSFieldMeshDecompositionGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), INTENT(OUT) :: DecompositionUserNumber !<The field mesh decomposition user number for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldMeshDecompositionGet(RegionUserNumber,FieldUserNumber,DecompositionUserNumber,&
    & CMISSFieldMeshDecompositionGetCNum)

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the mesh decomposition for a field identified by an object for C.
  FUNCTION CMISSFieldMeshDecompositionGetCPtr(FieldPtr,MeshDecompositionPtr) BIND(C, NAME = "CMISSFieldMeshDecompositionGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the mesh decomposition for, for C.
    TYPE(C_PTR), INTENT(OUT) :: MeshDecompositionPtr !<C pointer to the field mesh decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSDecompositionType), POINTER :: MeshDecomposition

    CMISSFieldMeshDecompositionGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldMeshDecompositionGet(Field, MeshDecomposition, CMISSFieldMeshDecompositionGetCPtr)
        IF(ASSOCIATED(MeshDecomposition)) THEN
          MeshDecompositionPtr = C_LOC(MeshDecomposition)
        ELSE
          CMISSFieldMeshDecompositionGetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldMeshDecompositionGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldMeshDecompositionGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh decomposition for a field identified by a user number for C.
  FUNCTION CMISSFieldMeshDecompositionSetCNum(RegionUserNumber,FieldUserNumber,MeshUserNumber,DecompositionUserNumber) &
  & BIND(C, NAME = "CMISSFieldMeshDecompositionSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the mesh decomposition to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number for the mesh to set the mesh decomposititon to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The field mesh decomposition user number for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldMeshDecompositionSet(RegionUserNumber,FieldUserNumber,MeshUserNumber,DecompositionUserNumber,&
    & CMISSFieldMeshDecompositionSetCNum)

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh decomposition for a field identified by an object for C.
  FUNCTION CMISSFieldMeshDecompositionSetCPtr(FieldPtr,MeshDecompositionPtr) BIND(C,NAME = "CMISSFieldMeshDecompositionSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the mesh decomposition to, for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshDecompositionPtr !<C pointer to the mesh decomposition for the field to set.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldMeshDecompositionSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSDecompositionType), POINTER :: MeshDecomposition

    CMISSFieldMeshDecompositionSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(MeshDecompositionPtr)) THEN
          CALL C_F_POINTER(MeshDecompositionPtr, MeshDecomposition)
          IF (ASSOCIATED(MeshDecomposition)) THEN
            CALL CMISSFieldMeshDecompositionSet(Field, MeshDecomposition, CMISSFieldMeshDecompositionSetCPtr)
          ELSE
            CMISSFieldMeshDecompositionSetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldMeshDecompositionSetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldMeshDecompositionSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldMeshDecompositionSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldMeshDecompositionSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the number of componenets for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfComponentsGetCNum(RegionUserNumber,FieldUserNumber,VariableType,NumberOfComponents) &
  & BIND(C, NAME = "CMISSFieldNumberOfComponentsGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType
    INTEGER(C_INT), INTENT(OUT) :: NumberOfComponents
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsGetCNum
    !Local variables

    CALL CMISSFieldnumberOfComponentsGet(RegionUserNumber, FieldUserNumber, VariableType, NumberOfComponents, &
    & CMISSFieldNumberOfComponentsGetCNum)

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the number of components for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfComponentsGetCPtr(FieldPtr,VariableType,NumberOfComponents) BIND(C, NAME = &
  & "CMISSFieldNumberOfComponentsGet")

    !Argument variables
    TYPE(C_PTR),VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the number of components for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the number of components for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), INTENT(OUT) :: NumberOfComponents !<The number of components in the field variable for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfComponentsGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfComponentsGet(Field, VariableType, NumberOfComponents, CMISSFieldNumberOfComponentsGetCPtr)
      ELSE
        CMISSFieldNumberOfComponentsGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfComponentsGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsGetCPtr


  !
  !================================================================================================================================
  !

  !>Sets/changes the number of componenets for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfComponentsSetCNum(RegionUserNumber,FieldUserNumber,VariableType,NumberOfComponents) &
  & BIND(C, NAME = "CMISSFieldNumberOfComponentsSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to set the number of components to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the number of components to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the number of components to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfComponents !<The number of components in the field variable for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldNumberOfComponentsSet(RegionUserNumber,FieldUserNumber,VariableType,NumberOfComponents, &
    & CMISSFieldNumberofComponentsSetCNum)

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of components for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfComponentsSetCPtr(FieldPtr,VariableType,NumberOfComponents) BIND(C, NAME = &
  & "CMISSFieldNumberOfComponentsSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the number of components for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to set the number of components for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfComponents !<The number of components in the field variables for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfComponentsSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfComponentsSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfComponentsSet(Field, VariableType, NumberOfComponents, CMISSFieldNumberOfComponentsSetCPtr)
      ELSE
        CMISSFieldNumberOfComponentsSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfComponentsSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfComponentsSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the number of variables for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfVariablesGetCNum(RegionUserNumber,FieldUserNumber,NumberOfVariables) BIND(C, NAME = &
  & "CMISSFieldNumberOfVariablesGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the number of variables for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to get the number of variables for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfVariables !<The number of variables in the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldNumberOfVariablesGet(RegionUserNumber, FieldUserNumber, NumberOfVariables, &
    & CMISSFieldNumberOfVariablesGetCNum)

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the number of variables for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfVariablesGetCPtr(FieldPtr,NumberOfVariables) BIND(C, NAME = "CMISSFieldNumberOfVariablesGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the number of variables for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfVariables !<The number of variables in the field for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfVariablesGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfVariablesGet(Field, NumberOfVariables, CMISSFieldNumberOfVariablesGetCPtr)
      ELSE
        CMISSFieldNumberOfVariablesGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfVariablesGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of variables for a field identified by a user number for C.
  FUNCTION CMISSFieldNumberOfVariablesSetCNum(RegionUserNumber,FieldUserNumber,NumberOfVariables) BIND(C, &
  & NAME = "CMISSFieldNumberOfVariablesSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the number of variables to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the number of variables to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfVariables !<The number of variables set to the field, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesSetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldNumberofVariablesSet(RegionUserNumber, FieldUserNumber, NumberOfVariables, &
    & CMISSFieldNumberOfVariablesSetCNum)

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of variables for a field identified by an object for C.
  FUNCTION CMISSFieldNumberOfVariablesSetCPtr(FieldPtr,NumberOfVariables) BIND(C, NAME = "CMISSFieldNumberOfVariablesSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfVariables
    !Function variables
    INTEGER(C_INT) :: CMISSFieldNumberOfVariablesSetCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldNumberOfVariablesSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldNumberOfVariablesSet(Field, NumberOfVariables, CMISSFieldNumberOfVariablesSetCPtr)
      ELSE
        CMISSFieldNumberOfVariablesSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldNumberOfVariablesSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldNumberOfVariablesSetCPtr


  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddConstantIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME= "CMISSFieldParameterSetAddConstantIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantIntgCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddConstantIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value, CMISSFieldParameterSetAddConstantIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantIntgCNum

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for the constant of the field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetAddConstantIntgCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME= "CMISSFieldParameterSetAddConstantIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value  !<The integer value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantIntgCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantIntgCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantIntg(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantIntgCPtr)
      ELSE
        CMISSFieldParameterSetAddConstantIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantIntgCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.

  FUNCTION CMISSFieldParameterSetAddConstantSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddConstantSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantSPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddConstantSP (RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber,&
    & Value,CMISSFieldParameterSetAddConstantSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantSPCNum

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddConstantSPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
  & "CMISSFieldParameterSetAddConstantSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to add the constant to the field parameter set for, for C
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantSP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantSPCPtr)
      ELSE
        CMISSFieldParameterSetAddConstantSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantSPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddConstantDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddConstantDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number for the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantDPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddConstantDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber,&
    & Value,CMISSFieldParameterSetAddConstantDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantDPCNum

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddConstantDPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME = "CMISSFieldParameterSetAddConstantDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to add the constant to the field parameter set for, for C
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantDP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantDPCPtr)
      ELSE
        CMISSFieldParameterSetAddConstantDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantDPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddConstantLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddConstantLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetAddConstantL(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, ComponentNumber, &
    & Value, CMISSFieldParameterSetAddConstantLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddConstantLCNum

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for the constant of the field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetAddConstantLCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
  & NAME = "CMISSFieldParameterSetAddConstantL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the constant to the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the constant to the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to add the constant to the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddConstantLCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddConstantLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddConstantL(Field, VariableType, FieldSetType, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddConstantLCPtr)
      ELSE
        CMISSFieldParameterSetAddConstantLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddConstantLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN


  END FUNCTION CMISSFieldParameterSetAddConstantLCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an element in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddElementIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementIntgCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value, CMISSFieldParameterSetAddElementIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementIntgCNum

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an element in the given parameter set for field variable component for a field identified by an object, for C.
  FUNCTION CMISSFieldParameterSetAddElementIntgCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementIntgCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementIntg(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddElementIntgCPtr)
      ELSE
        CMISSFieldParameterSetAddElementIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementIntgCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an element in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddElementSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementSPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementSP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value,CMISSFieldParameterSetAddElementSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementSPCNum

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an element in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddElementSPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_FLOAT), VALUE, INTENT(IN) :: Value  !<The single precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementSPCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementSP(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddElementSPCPtr)
      ELSE
        CMISSFieldParameterSetAddElementSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementSPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an element in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddElementDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementDPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
    & ComponentNumber,Value,CMISSFieldParameterSetAddElementDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementDPCNum


    !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an element in the given parameter set for field variable component for a field identified by an object, for C.
  FUNCTION CMISSFieldParameterSetAddElementDPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    REAL(C_DOUBLE), VALUE, INTENT(IN) :: Value  !<The double precision value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementDPCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementDP(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
        & CMISSFieldParameterSetAddElementDPCPtr)
      ELSE
        CMISSFieldParameterSetAddElementDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementDPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to an element in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddElementLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddElementLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the element in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementLCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddElementL(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
    & ComponentNumber,Value,CMISSFieldParameterSetAddElementLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementLCNum

    !
  !================================================================================================================================
  !

  !>Adds the given logical value to an element in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddElementLCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
  & BIND(C, NAME = "CMISSFieldParameterSetAddElementL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the element in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the element in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to add the value to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the element in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the element in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddElementLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddElementLCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddElementL(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value, &
         & CMISSFieldParameterSetAddElementLCPtr)
      ELSE
        CMISSFieldParameterSetAddElementLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddElementLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddElementLCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an node in the given parameter set for field variable component for a field identified by a user number, for C.
  FUNCTION CMISSFieldParameterSetAddNodeIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value  !<The integer value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeIntgCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber,&
    & UserNodeNumber,ComponentNumber,Value,CMISSFieldParameterSetAddNodeIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeIntgCNum

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to an node in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddNodeIntgCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value  !<The integer value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeIntgCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeIntg(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetAddNodeIntgCPtr)
      ELSE
        CMISSFieldParameterSetAddNodeIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeIntgCPtr



  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an node in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddNodeSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value  !<The single precision value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeSPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeSP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value,CMISSFieldParameterSetAddNodeSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeSPCNum

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to an node in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddNodeSPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value  !<The single precision value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeSPCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeSP(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetAddNodeSPCPtr)
      ELSE
        CMISSFieldParameterSetAddNodeSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeSPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an node in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddNodeDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value  !<The double precision value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeDPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value,CMISSFieldParameterSetAddNodeDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeDPCNum

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to an node in the given parameter set for field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetAddNodeDPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value  !<The double precision value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeDPCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeDP(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
        & ComponentNumber, Value, CMISSFieldParameterSetAddNodeDPCPtr)
      ELSE
        CMISSFieldParameterSetAddNodeDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeDPCPtr

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to an node in the given parameter set for field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetAddNodeLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber, &
    & UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to add the value to the node in the field parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeLCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetAddNodeL(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber,&
    & ComponentNumber,Value,CMISSFieldParameterSetAddNodeLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeLCNum

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to an node in the given parameter set for field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetAddNodeLCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber,ComponentNumber, &
    & Value) BIND(C, NAME = "CMISSFieldParameterSetAddNodeL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr  !<C pointer to the field to add the value to the node in the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to add the value to the node in the field parameter set for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The node derivative number of the node to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to add the value to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber  !<The component number of the field variable to add the value to the node in the field parameter set for, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value  !<The logical value to add to the node in the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetAddNodeLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetAddNodeLCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetAddNodeL(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, ComponentNumber, &
        & Value, CMISSFieldParameterSetAddNodeLCPtr)
      ELSE
        CMISSFieldParameterSetAddNodeLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetAddNodeLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetAddNodeLCPtr

  !
  !================================================================================================================================
  !

  !>Creates a new parameter set of type set type for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetCreateCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, &
  & NAME = "CMISSFieldParameterSetCreateNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to create the parameter set on for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to create the parameter set on for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to create the parameter set on for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to create, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetCreateCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetCreate(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & CMISSFieldParameterSetCreateCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetCreateCNum

  !
  !================================================================================================================================
  !

  !>Creates a new parameter set of type set type for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetCreateCPtr(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = "CMISSFieldParameterSetCreate")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to create the parameter set on for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to create the parameter set on for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to create, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetCreateCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetCreateCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetCreate(Field,VariableType,FieldSetType,CMISSFieldParameterSetCreateCPtr)
      ELSE
        CMISSFieldParameterSetCreateCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetCreateCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetCreateCPtr

  !
  !================================================================================================================================
  !

  !>Destroys the specified parameter set type for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetDestroyCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, NAME= &
  & "CMISSFieldParameterSetDestroyNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to destroy the parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to destroy the parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to destroy the parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to destroy, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDestroyCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetDestroy(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & CMISSFieldParameterSetDestroyCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys the specified parameter set type for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetDestroyCPtr(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = "CMISSFieldParameterSetDestroy")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to destroy the parameter set for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to destroy the parameter set for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to destroy, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDestroyCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetDestroy(Field,VariableType,FieldSetType,CMISSFieldParameterSetDestroyCPtr)
      ELSE
        CMISSFieldParameterSetDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local integer data array for a field identified by an user number, for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ParametersSize, &
  & ParametersPtr) BIND(C, NAME = "CMISSFieldParameterSetDataGetIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetIntgCNum
    !Local variables
    INTEGER(C_INT), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataGetIntgCNum = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetIntgCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSFieldParameterSetDataGetIntg(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType, Parameters, &
        & CMISSFieldParameterSetDataGetIntgCNum)
      ParametersSize = Size(Parameters,1)
      ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetIntgCNum

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local integer data array for a field identified by an object for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetIntgCPtr(FieldPtr,VariableType,FieldSetType,ParametersSize,ParametersPtr) BIND(C, NAME = &
  & "CMISSFieldParameterSetDataGetIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetIntgCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    INTEGER(C_INT), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataGetIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetIntgCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(FieldPtr)) THEN
        CALL C_F_POINTER(FieldPtr, Field)
        IF(ASSOCIATED(Field)) THEN
          CALL CMISSFieldParameterSetDataGetIntg(Field, VariableType, FieldSetType, Parameters, &
            & CMISSFieldParameterSetDataGetIntgCPtr)
          ParametersSize = Size(Parameters,1)
          ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
          ELSE
          CMISSFieldParameterSetDataGetIntgCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSFieldParameterSetDataGetIntgCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetIntgCPtr

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local single precision data array for a field identified by an user number for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
  &ParametersSize,ParametersPtr) BIND(C, NAME = "CMISSFieldParameterSetDataGetSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetSPCNum !<Error code.
    !Local variables
    REAL(C_FLOAT), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataGetSPCNum = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetSPCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSFieldParameterSetDataGetSP(RegionUserNumber,FieldUserNumber, VariableType, FieldSetType, Parameters, &
      & CMISSFieldParameterSetDataGetSPCNum)
      ParametersSize = Size(Parameters, 1)
      ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetSPCNum

  !
  !================================================================================================================================
  !
  !>Returns a pointer to the specified field parameter set local single precision data array for a field identified by an object. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.

  FUNCTION CMISSFieldParameterSetDataGetSPCPtr(FieldPtr,VariableType,FieldSetType,ParametersSize,ParametersPtr) BIND(C, NAME = &
  & "CMISSFieldParameterSetDataGetSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetSPCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    REAL(C_FLOAT), POINTER :: Parameters (:)

    CMISSFieldParameterSetDataGetSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetSPCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(FieldPtr)) THEN
        CALL C_F_POINTER(FieldPtr,Field)
        IF(ASSOCIATED(Field)) THEN
          CALL CMISSFieldParameterSetDataGetSP(Field, VariableType, FieldSetType, Parameters, CMISSFieldParameterSetDataGetSPCPtr)
          ParametersSize = Size(Parameters, 1)
          ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
        ELSE
          CMISSFieldParameterSetDataGetSPCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSFieldParameterSetDataGetSPCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetSPCPtr

 !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local double precision data array for a field identified by an user number, for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ParametersSize,&
    &ParametersPtr) BIND(C, NAME = "CMISSFieldParameterSetDataGetDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetDPCNum !<Error code.
    !Local variables
    REAL(C_DOUBLE), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataGetDPCNum = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetDPCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSFieldParameterSetDataGetDP(RegionUserNumber,FieldUserNumber, VariableType, FieldSetType, Parameters, &
      & CMISSFieldParameterSetDataGetDPCNum)
      ParametersSize = Size(Parameters, 1)
      ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetDPCNum

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local double precision data array for a field identified by an object, for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetDPCPtr(FieldPtr,VariableType,FieldSetType,ParametersSize,ParametersPtr) BIND(C, NAME = &
    & "CMISSFieldParameterSetDataGetDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetDPCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    REAL(C_DOUBLE), POINTER :: Parameters (:)

    CMISSFieldParameterSetDataGetDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetDPCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(FieldPtr)) THEN
        CALL C_F_POINTER(FieldPtr,Field)
        IF(ASSOCIATED(Field)) THEN
          CALL CMISSFieldParameterSetDataGetDP(Field, VariableType, FieldSetType, Parameters, CMISSFieldParameterSetDataGetDPCPtr)
          ParametersSize = Size(Parameters, 1)
          ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
        ELSE
          CMISSFieldParameterSetDataGetDPCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSFieldParameterSetDataGetDPCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetDPCPtr

!
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local logical data array for a field identified by an user number, for C. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ParametersSize, &
  & ParametersPtr) BIND(C, NAME = "CMISSFieldParameterSetDataGetLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to get the parameter set data for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetLCNum !<Error code.
    !Local variables
    LOGICAL(C_BOOL), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataGetLCNum = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetLCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSFieldParameterSetDataGetL(RegionUserNumber,FieldUserNumber, VariableType, FieldSetType, Parameters, &
      & CMISSFieldParameterSetDataGetLCNum)
      ParametersSize = Size(Parameters, 1)
      ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetLCNum

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set local logical data array for a field identified by an object. The pointer must be restored with a call to OPENCMISS::CMISSFieldParameterSetDataRestore call. Note: the values can be used for read operations but a field parameter set update or add calls must be used to change any values.
  FUNCTION CMISSFieldParameterSetDataGetLCPtr(FieldPtr,VariableType,FieldSetType,ParametersSize,ParametersPtr) BIND(C, NAME = &
    & "CMISSFieldParameterSetDataGetL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to get, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataGetLCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    LOGICAL(C_BOOL), POINTER :: Parameters (:)

    CMISSFieldParameterSetDataGetLCPtr = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataGetLCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(FieldPtr)) THEN
        CALL C_F_POINTER(FieldPtr,Field)
        IF(ASSOCIATED(Field)) THEN
          CALL CMISSFieldParameterSetDataGetL(Field, VariableType, FieldSetType, Parameters, CMISSFieldParameterSetDataGetLCPtr)
          ParametersSize = Size(Parameters, 1)
          ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
        ELSE
          CMISSFieldParameterSetDataGetLCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSFieldParameterSetDataGetLCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataGetLCPtr

  !
  !================================================================================================================================
  !

  !>Restores the specified field variable parameter set local integer array that was obtained with an OPENCMISS::CMISSFieldParameterSetDataGet call for a field that is specified with an user number, for C.
  FUNCTION CMISSFieldParameterSetDataRestoreIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ParametersSize, &
    & ParametersPtr) BIND(C, NAME = "CMISSFieldParameterSetDataRestoreIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to restore the parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to restore the parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to restore the parameter set data for. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to restore. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataRestoreIntgCNum
    !Local variables
    INTEGER(C_INT), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataRestoreIntgCNum = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataRestoreIntgCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSFieldParameterSetDataRestoreIntg(RegionUserNumber,FieldUserNumber, VariableType, FieldSetType, Parameters, &
        & CMISSFieldParameterSetDataRestoreIntgCNum)
      ParametersSize = Size(Parameters, 1)
      ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataRestoreIntgCNum

  !
  !================================================================================================================================
  !

  !>Restores the specified field variable parameter set local integer array that was obtained with an OPENCMISS::CMISSFieldParameterSetDataGet call for a field that is specified with an object.
  FUNCTION CMISSFieldParameterSetDataRestoreIntgCPtr(FieldPtr,VariableType,FieldSetType,ParametersSize,ParametersPtr) BIND(C, &
    & NAME = "CMISSFieldParameterSetDataRestoreIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to restore the field parameter restore data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to restore the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to restore, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data to restore.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataRestoreIntgCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    INTEGER(C_INT), POINTER :: Parameters (:)

    CMISSFieldParameterSetDataRestoreIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataRestoreIntgCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(FieldPtr)) THEN
        CALL C_F_POINTER(FieldPtr,Field)
        IF(ASSOCIATED(Field)) THEN
          CALL CMISSFieldParameterSetDataRestoreIntg(Field, VariableType, FieldSetType, Parameters, &
            & CMISSFieldParameterSetDataRestoreIntgCPtr)
          ParametersSize = Size(Parameters, 1)
          ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
        ELSE
          CMISSFieldParameterSetDataRestoreIntgCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSFieldParameterSetDataRestoreIntgCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataRestoreIntgCPtr

  !
  !================================================================================================================================
  !

  !>Restores the specified field variable parameter set local single precision array that was obtained with an OPENCMISS::CMISSFieldParameterSetDataGet call for a field that is specified with an user number.
  FUNCTION CMISSFieldParameterSetDataRestoreSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, ParametersSize, &
    & ParametersPtr) BIND(C, NAME = "CMISSFieldParameterSetDataRestoreSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to restore the parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to restore the parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to restore the parameter set data for. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to restore. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataRestoreSPCNum
    !Local variables
    REAL(C_FLOAT), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataRestoreSPCNum = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataRestoreSPCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSFieldParameterSetDataRestoreSP(RegionUserNumber,FieldUserNumber, VariableType, FieldSetType, Parameters, &
        & CMISSFieldParameterSetDataRestoreSPCNum)
      ParametersSize = Size(Parameters, 1)
      ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataRestoreSPCNum


  !
  !================================================================================================================================
  !

  !>>Restores the specified field variable parameter set local single precision array that was obtained with an OPENCMISS::CMISSFieldParameterSetDataGet call for a field that is specified with an object.
  FUNCTION CMISSFieldParameterSetDataRestoreSPCPtr(FieldPtr,VariableType,FieldSetType,ParametersSize,ParametersPtr) BIND(C, NAME &
    & = "CMISSFieldParameterSetDataRestoreSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to restore the field parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to restore the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to restore, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data to restore.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataRestoreSPCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    REAL(C_FLOAT), POINTER :: Parameters (:)

    CMISSFieldParameterSetDataRestoreSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataRestoreSPCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(FieldPtr)) THEN
        CALL C_F_POINTER(FieldPtr,Field)
        IF(ASSOCIATED(Field)) THEN
          CALL CMISSFieldParameterSetDataRestoreSP(Field, VariableType, FieldSetType, Parameters, &
            & CMISSFieldParameterSetDataRestoreSPCPtr)
          ParametersSize = Size(Parameters, 1)
          ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
        ELSE
          CMISSFieldParameterSetDataRestoreSPCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSFieldParameterSetDataRestoreSPCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataRestoreSPCPtr

  !
  !================================================================================================================================
  !

  !>Restores the specified field variable parameter set local double precision array that was obtained with an OPENCMISS::CMISSFieldParameterSetDataGet call for a field that is specified with an user number.
  FUNCTION CMISSFieldParameterSetDataRestoreDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ParametersSize, &
    & ParametersPtr) BIND(C, NAME = "CMISSFieldParameterSetDataRestoreDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to restore the parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to restore the parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to restore the parameter set data for. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to restore. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataRestoreDPCNum
    !Local variables
    REAL(C_DOUBLE), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataRestoreDPCNum = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataRestoreDPCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSFieldParameterSetDataRestoreDP(RegionUserNumber,FieldUserNumber, VariableType, FieldSetType, Parameters, &
        & CMISSFieldParameterSetDataRestoreDPCNum)
      ParametersSize = Size(Parameters, 1)
      ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataRestoreDPCNum

  !
  !================================================================================================================================
  !

  !>>Restores the specified field variable parameter set local double precision array that was obtained with an OPENCMISS::CMISSFieldParameterSetDataGet call for a field that is specified with an object.
  FUNCTION CMISSFieldParameterSetDataRestoreDPCPtr(FieldPtr,VariableType,FieldSetType,ParametersSize,ParametersPtr) BIND(C, &
    & NAME = "CMISSFieldParameterSetDataRestoreDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to restore the field parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to restore the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to restore, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data to restore.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataRestoreDPCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    REAL(C_DOUBLE), POINTER :: Parameters (:)

    CMISSFieldParameterSetDataRestoreDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataRestoreDPCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(FieldPtr)) THEN
        CALL C_F_POINTER(FieldPtr,Field)
        IF(ASSOCIATED(Field)) THEN
          CALL CMISSFieldParameterSetDataRestoreDP(Field, VariableType, FieldSetType, Parameters, &
            & CMISSFieldParameterSetDataRestoreDPCPtr)
          ParametersSize = Size(Parameters, 1)
          ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
        ELSE
          CMISSFieldParameterSetDataRestoreDPCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSFieldParameterSetDataRestoreDPCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataRestoreDPCPtr

  !
  !================================================================================================================================
  !

  !>Restores the specified field variable parameter set local logical array that was obtained with an OPENCMISS::CMISSFieldParameterSetDataGet call for a field that is specified with an user number.
  FUNCTION CMISSFieldParameterSetDataRestoreLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ParametersSize, &
    & ParametersPtr) BIND(C, NAME = "CMISSFieldParameterSetDataRestoreLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to restore the parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to restore the parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to restore the parameter set data for. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to restore. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataRestoreLCNum
    !Local variables
    LOGICAL(C_BOOL), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataRestoreLCNum = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataRestoreLCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSFieldParameterSetDataRestoreL(RegionUserNumber,FieldUserNumber, VariableType, FieldSetType, Parameters, &
        & CMISSFieldParameterSetDataRestoreLCNum)
      ParametersSize = Size(Parameters, 1)
      ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataRestoreLCNum

  !
  !================================================================================================================================
  !

  !>Restores the specified field variable parameter set local logical array that was obtained with an OPENCMISS::CMISSFieldParameterSetDataGet call for a field that is specified with an object.
  FUNCTION CMISSFieldParameterSetDataRestoreLCPtr(FieldPtr,VariableType,FieldSetType,ParametersSize,ParametersPtr) BIND(C, NAME = &
    & "CMISSFieldParameterSetDataRestoreL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to restore the field parameter set data for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to restore the parameter set data for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the parameter set data to restore, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), INTENT(OUT) :: ParametersSize !<Size of parameter set data, for C.
    TYPE(C_PTR), INTENT(OUT) :: ParametersPtr !<C pointer to the parameter set data to restore.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetDataRestoreLCPtr
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    REAL(C_DOUBLE), POINTER :: Parameters(:)

    CMISSFieldParameterSetDataRestoreLCPtr = CMISSNoError
    IF(C_ASSOCIATED(ParametersPtr)) THEN
      CMISSFieldParameterSetDataRestoreLCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(FieldPtr)) THEN
        CALL C_F_POINTER(FieldPtr,Field)
        IF(ASSOCIATED(Field)) THEN
          CALL CMISSFieldParameterSetDataRestoreL(Field, VariableType, FieldSetType, Parameters, &
            & CMISSFieldParameterSetDataRestoreLCPtr)
          ParametersSize = Size(Parameters, 1)
          ParametersPtr = C_LOC(Parameters(1)) !Point to first element as fortran pointers to arrays are not interoperable. This assumes that the parameters array is sequential in memory
        ELSE
          CMISSFieldParameterSetDataRestoreLCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSFieldParameterSetDataRestoreLCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetDataRestoreLCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantIntgCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber,Value, &
      & CMISSFieldParameterSetGetConstantIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantIntgCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantIntgCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetGetConstantIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantIntgCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantIntg(Field,VariableType,FieldSetType,ComponentNumber,Value, &
          & CMISSFieldParameterSetGetConstantIntgCPtr)
      ELSE
        CMISSFieldParameterSetGetConstantIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantIntgCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantSP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
      & Value,CMISSFieldParameterSetGetConstantSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantSPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantSPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetGetConstantSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantSPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantSP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetConstantSPCPtr)
      ELSE
        CMISSFieldParameterSetGetConstantSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantSPCPtr


  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
      & Value,CMISSFieldParameterSetGetConstantDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantDPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantDPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetGetConstantDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantDPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantDP(Field, VariableType, FieldSetType, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetConstantDPCPtr)
      ELSE
        CMISSFieldParameterSetGetConstantDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantDPCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified constant of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetConstantLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetConstantLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetConstantL(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
      & Value,CMISSFieldParameterSetGetConstantLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantLCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified constant of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetConstantLCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetGetConstantL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the constant value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the constant value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the constant value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the constant value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetConstantLCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetConstantLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetConstantL(Field, VariableType, FieldSetType, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetConstantLCPtr)
      ELSE
        CMISSFieldParameterSetGetConstantLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetConstantLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetConstantLCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetElementIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementIntgCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
      & ComponentNumber,Value,CMISSFieldParameterSetGetElementIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementIntgCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified element of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetElementIntgCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
    & BIND(C, NAME = "CMISSFieldParameterSetGetElementIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementIntg(Field, VariableType, FieldSetType,UserElementNumber, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetElementIntgCPtr)
      ELSE
        CMISSFieldParameterSetGetElementIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementIntgCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetElementSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementSP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
      & ComponentNumber,Value,CMISSFieldParameterSetGetElementSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementSPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified element of a field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetGetElementSPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
    & BIND(C, NAME = "CMISSFieldParameterSetGetElementSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementSP(Field, VariableType, FieldSetType,UserElementNumber, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetElementSPCPtr)
      ELSE
        CMISSFieldParameterSetGetElementSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementSPCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetElementDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
      & ComponentNumber,Value,CMISSFieldParameterSetGetElementDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementDPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified element of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetElementDPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
    & BIND(C, NAME = "CMISSFieldParameterSetGetElementDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementDP(Field, VariableType, FieldSetType,UserElementNumber, ComponentNumber, Value, &
          & CMISSFieldParameterSetGetElementDPCPtr)
      ELSE
        CMISSFieldParameterSetGetElementDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementDPCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified element of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetElementLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetElementLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetElementL(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,UserElementNumber, &
      & ComponentNumber,Value,CMISSFieldParameterSetGetElementLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementLCNum


  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified element of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetElementLCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value) &
    & BIND(C, NAME = "CMISSFieldParameterSetGetElementL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the element value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the element value from the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the element value from, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number to get the value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the element value from the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetElementLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER  :: Field

    CMISSFieldParameterSetGetElementLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetElementL(Field, VariableType, FieldSetType, UserElementNumber, ComponentNumber, Value,  &
          & CMISSFieldParameterSetGetElementLCPtr)
      ELSE
        CMISSFieldParameterSetGetElementLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetElementLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetElementLCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeIntgCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeIntg(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeIntgCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set an integer value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeIntgCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), INTENT(OUT) :: Value !<The integer value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeIntgCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeIntg(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value, CMISSFieldParameterSetGetNodeIntgCPtr)
      ELSE
        CMISSFieldParameterSetGetNodeIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeIntgCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    REAL(C_FLOAT),INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeSP(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
      &UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeSPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a single precision value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeSPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    REAL(C_FLOAT), INTENT(OUT) :: Value !<The single precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeSPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeSP(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value, CMISSFieldParameterSetGetNodeSPCPtr)
      ELSE
        CMISSFieldParameterSetGetNodeSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeSPCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    REAL(C_DOUBLE),INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeDP(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeDPCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a double precision value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeDPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    REAL(C_DOUBLE), INTENT(OUT) :: Value !<The double precision value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeDPCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeDP(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value, CMISSFieldParameterSetGetNodeDPCPtr)
      ELSE
        CMISSFieldParameterSetGetNodeDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeDPCPtr

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified node and derivative of a field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetGetNodeLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the nodal value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    LOGICAL(C_BOOL),INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetGetNodeL(RegionUserNumber, FieldUserNumber, VariableType, FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber,Value, CMISSFieldParameterSetGetNodeLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeLCNum

  !
  !================================================================================================================================
  !

  !>Returns from the given parameter set a logical value for the specified node and derivative of a field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetGetNodeLCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetGetNodeL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the nodal value from the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to get the nodal value from the field parameter set, for C.. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to get the nodal value from, for C.. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number to get the value from the field parameter set, for C..
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to get the nodal value from the field parameter set, for C..
    LOGICAL(C_BOOL), INTENT(OUT) :: Value !<The logical value from the field parameter set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetGetNodeLCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetGetNodeLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetGetNodeL(Field, VariableType, FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value, CMISSFieldParameterSetGetNodeLCPtr)
      ELSE
        CMISSFieldParameterSetGetNodeLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetGetNodeLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetGetNodeLCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantIntgCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantIntg(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,&
      & ComponentNumber,Value,CMISSFieldParameterSetUpdateConstantIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantIntgCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the constant of the field variable component for a field identified by an object.
  FUNCTION CMISSFieldParameterSetUpdateConstantIntgCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetUpdateConstantIntg")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantIntg(Field, VariableType,FieldSetType,ComponentNumber,Value,&
          & CMISSFieldParameterSetUpdateConstantIntgCPtr)
      ELSE
        CMISSFieldParameterSetUpdateConstantIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantIntgCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantSP(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,ComponentNumber,&
      & Value,CMISSFieldParameterSetUpdateConstantSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantSPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantSPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
    & NAME = "CMISSFieldParameterSetUpdateConstantSP")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantSP(Field, VariableType,FieldSetType,ComponentNumber,Value, &
          & CMISSFieldParameterSetUpdateConstantSPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateConstantSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantSPCPtr
  
  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantDP(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,&
      & ComponentNumber,Value,CMISSFieldParameterSetUpdateConstantDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantDPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantDPCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, &
    & NAME = "CMISSFieldParameterSetUpdateConstantDP")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantDP(Field, VariableType,FieldSetType,ComponentNumber,Value,&
          & CMISSFieldParameterSetUpdateConstantDPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateConstantDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantDPCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the constant of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateConstantLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the constant value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateConstantL(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,ComponentNumber, &
      & Value,CMISSFieldParameterSetUpdateConstantLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantLCNum


  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the constant of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateConstantLCPtr(FieldPtr,VariableType,FieldSetType,ComponentNumber,Value) BIND(C, NAME = &
    & "CMISSFieldParameterSetUpdateConstantL")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the constant value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the constant value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the constant value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the constant value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateConstantLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateConstantLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateConstantL(Field, VariableType,FieldSetType,ComponentNumber,Value,&
          & CMISSFieldParameterSetUpdateConstantLCPtr)
      ELSE
        CMISSFieldParameterSetUpdateConstantLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateConstantLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateConstantLCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    INTEGER(C_INT), INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementIntgCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementIntg(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,&
      & UserElementNumber,ComponentNumber,Value,CMISSFieldParameterSetUpdateElementIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementIntgCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementIntgCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateElementIntg")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    INTEGER(C_INT), INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementIntg(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value,&
          & CMISSFieldParameterSetUpdateElementIntgCPtr)
      ELSE
        CMISSFieldParameterSetUpdateElementIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementIntgCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementSPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementSP(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
      & ComponentNumber,Value,CMISSFieldParameterSetUpdateElementSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementSPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementSPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateElementSP")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementSP(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value, &
          & CMISSFieldParameterSetUpdateElementSPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateElementSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementSPCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementDPCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementDP(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
      & ComponentNumber,Value,CMISSFieldParameterSetUpdateElementDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementDPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementDPCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateElementDP")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementDP(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value,&
          & CMISSFieldParameterSetUpdateElementDPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateElementDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementDPCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the element of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateElementLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & UserElementNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the element value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementLCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateElementL(RegionUserNumber, FieldUserNumber,VariableType,FieldSetType,UserElementNumber,&
      & ComponentNumber,Value,CMISSFieldParameterSetUpdateElementLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementLCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the element of the field variable component for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateElementLCPtr(FieldPtr,VariableType,FieldSetType,UserElementNumber,ComponentNumber, &
    & Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateElementL")

    !Argument variables.
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to update the element value for the field parameter set.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the element value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the element value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserElementNumber !<The user element number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the element value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateElementLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateElementLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateElementL(Field, VariableType,FieldSetType,UserElementNumber,ComponentNumber,Value, &
          & CMISSFieldParameterSetUpdateElementLCPtr)
      ELSE
        CMISSFieldParameterSetUpdateElementLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateElementLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateElementLCPtr

  !
  !================================================================================================================================
  !

  !>Finishes the parameter set update for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateFinishCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, NAME = &
    & "CMISSFieldParameterSetUpdateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to finish the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to finish the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to finish the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to finish the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateFinishCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateFinish(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
      & CMISSFieldParameterSetUpdateFinishCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the parameter set update for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateFinishCPtr(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = &
    & "CMISSFieldParameterSetUpdateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to finish the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to finish the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateFinishCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateFinish(Field,VariableType,FieldSetType,CMISSFieldParameterSetUpdateFinishCPtr)
      ELSE
        CMISSFieldParameterSetUpdateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeIntgCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeIntgNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeIntgCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeIntg(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeIntgCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeIntgCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the node and derivative of the field variable component for a field identified by an object for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeIntgCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeIntg")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: Value !<The integer value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeIntgCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeIntgCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeIntg(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber,Value,CMISSFieldParameterSetUpdateNodeIntgCPtr)
      ELSE
        CMISSFieldParameterSetUpdateNodeIntgCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeIntgCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeIntgCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeSPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeSPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeSPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeSP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeSPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeSPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the node and derivative of the field variable component for a field identified by an object, for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeSPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeSP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    REAL(C_FLOAT), INTENT(IN) :: Value !<The single precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeSPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeSPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeSP(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value,CMISSFieldParameterSetUpdateNodeSPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateNodeSPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeSPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeSPCPtr

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeDPCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value)  BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeDPNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeDPCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeDP(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber, &
      & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeDPCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeDPCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the node and derivative of the field variable component for a field identified by an object for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeDPCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeDP")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    REAL(C_DOUBLE), INTENT(IN) :: Value !<The double precision value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeDPCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeDPCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeDP(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value,CMISSFieldParameterSetUpdateNodeDPCPtr)
      ELSE
        CMISSFieldParameterSetUpdateNodeDPCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeDPCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeDPCPtr

  !
  !================================================================================================================================
  !,

  !>Updates the given parameter set with the given logical value for the node and derivative of the field variable component for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateNodeLCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, &
    & DerivativeNumber,UserNodeNumber,ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeLNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to update the nodal value for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeLCNum !<Error Code.
    !Local variable

    CALL CMISSFieldParameterSetUpdateNodeL(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType, DerivativeNumber,&
      & UserNodeNumber, ComponentNumber, Value, CMISSFieldParameterSetUpdateNodeLCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeLCNum

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the node and derivative of the field variable component for a field identified by an object for C.

  FUNCTION CMISSFieldParameterSetUpdateNodeLCPtr(FieldPtr,VariableType,FieldSetType,DerivativeNumber,UserNodeNumber, &
    & ComponentNumber,Value) BIND(C, NAME = "CMISSFieldParameterSetUpdateNodeL")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to finishe the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to update the nodal value for the field parameter set, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type of the field to update the nodal value for, for C. \see OPENCMISS_FieldParameterSetTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: DerivativeNumber !<The derivative number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: UserNodeNumber !<The user node number of the field variable component to update for the field parameter set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ComponentNumber !<The component number of the field variable to update the nodal value for the field parameter set, for C.
    LOGICAL(C_BOOL), INTENT(IN) :: Value !<The logical value to update the field parameter set to, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateNodeLCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateNodeLCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateNodeL(Field,VariableType,FieldSetType, DerivativeNumber, UserNodeNumber, &
          & ComponentNumber, Value,CMISSFieldParameterSetUpdateNodeLCPtr)
      ELSE
        CMISSFieldParameterSetUpdateNodeLCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateNodeLCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateNodeLCPtr

  !
  !================================================================================================================================
  !

  !>Starts the parameter set update for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldParameterSetUpdateStartCNum(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType) BIND(C, &
    & NAME = "CMISSFieldParameterSetUpdateStartNum")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber!<The user number of the region containing the field to start the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to start the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to start the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to start the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateStartCNum !<Error Code.
    !Local variables

    CALL CMISSFieldParameterSetUpdateStart(RegionUserNumber,FieldUserNumber,VariableType,FieldSetType,&
      & CMISSFieldParameterSetUpdateStartCNum)

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts the parameter set update for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldParameterSetUpdateStartCPtr(FieldPtr,VariableType,FieldSetType) BIND(C, NAME = &
    & "CMISSFieldParameterSetUpdateStart")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to start the parameter set update for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type of the field to start the parameter set update for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldSetType !<The parameter set type to start the update for, for C. \see OPENCMISS_FieldParameterSetTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldParameterSetUpdateStartCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldParameterSetUpdateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldParameterSetUpdateStart(Field,VariableType,FieldSetType,CMISSFieldParameterSetUpdateStartCPtr)
      ELSE
        CMISSFieldParameterSetUpdateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldParameterSetUpdateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldParameterSetUpdateStartCPtr

  !
  !================================================================================================================================
  !

  !>Returns the scaling type for a field identified by a user number for C.
  FUNCTION CMISSFieldScalingTypeGetCNum(RegionUserNumber,FieldUserNumber,ScalingType) BIND(C, NAME = &
    & "CMISSFieldScalingTypeGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the scaling type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the scaling type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ScalingType !<The field scaling type to get, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeGetCNum !<Error Code.
    !Local variables

    CALL CMISSFieldScalingTypeGet(RegionUserNumber,FieldUserNumber,ScalingType,CMISSFieldScalingTypeGetCNum)

    RETURN

  END FUNCTION CMISSFieldScalingTypeGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the scaling type for a field identified by an object for C.
  FUNCTION CMISSFieldScalingTypeGetCPtr(FieldPtr,ScalingType) BIND(C, NAME = "CMISSFieldScalingTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the scaling type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ScalingType !<The field scaling type to get, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeGetCPtr
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldScalingTypeGetCPtr =CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr,Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldScalingTypeGet(Field, ScalingType, CMISSFieldScalingTypeGetCPtr)
      ELSE
        CMISSFieldScalingTypeGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldScalingTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldScalingTypeGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the scaling type for a field identified by a user number for C.
  FUNCTION CMISSFieldScalingTypeSetCNum(RegionUserNumber,FieldUserNumber,ScalingType) BIND(C, NAME = &
    & "CMISSFieldScalingTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the scaling type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the scaling type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ScalingType !<The field scaling type to set, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeSetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldScalingTypeSet(RegionUserNumber,FieldUserNumber,ScalingType,CMISSFieldScalingTypeSetCNum)

    RETURN

  END FUNCTION CMISSFieldScalingTypeSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the scaling type for a field identified by an object for C.
  FUNCTION CMISSFieldScalingTypeSetCPtr(FieldPtr,ScalingType) BIND(C, NAME = "CMISSFieldScalingTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the scaling type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ScalingType !<The field scaling type to set, for C. \see OPENCMISS_FieldScalingTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldScalingTypeSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldScalingTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldScalingTypeSet(Field, ScalingType, CMISSFieldScalingTypeSetCPtr)
      ELSE
        CMISSFieldScalingTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldScalingTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldScalingTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the field type for a field identified by a user number for C.
  FUNCTION CMISSFieldTypeGetCNum(RegionUserNumber,FieldUserNumber,FieldType) BIND(C, NAME = "CMISSFieldTypeGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the field type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the field type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: FieldType !<The field type to get, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeGetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldTypeGet(RegionUserNumber,FieldUserNumber,FieldType,CMISSFieldTypeGetCNum)

    RETURN

  END FUNCTION CMISSFieldTypeGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the type for a field identified by an object for C.
  FUNCTION CMISSFieldTypeGetCPtr(FieldPtr,FieldType) BIND(C, NAME = "CMISSFieldTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: FieldType !<The field type to get, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldTypeGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldTypeGet(Field, FieldType, CMISSFieldTypeGetCPtr)
      ELSE
        CMISSFieldTypeGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldTypeGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the field type for a field identified by a user number for C.
  FUNCTION CMISSFieldTypeSetCNum(RegionUserNumber,FieldUserNumber,FieldType) BIND(C, NAME = "CMISSFieldTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the field type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the field type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldType !<The field type to set, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeSetCNum !<Error Code.
    !Local variable

    CALL CMISSFieldTypeSet(RegionUserNumber,FieldUserNumber,FieldType,CMISSFieldTypeSetCNum)

    RETURN

  END FUNCTION CMISSFieldTypeSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the type for a field identified by an object for C.
  FUNCTION CMISSFieldTypeSetCPtr(FieldPtr,FieldType) BIND(C, NAME = "CMISSFieldTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the field type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldType !<The field type to set, for C. \see OPENCMISS_FieldTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldTypeSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field

    CMISSFieldTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldTypeSet(Field, FieldType, CMISSFieldTypeSetCPtr)
      ELSE
        CMISSFieldTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field variable for a field identified by a user number.
  FUNCTION CMISSFieldVariableLabelGetCCNum(RegionUserNumber,FieldUserNumber,VariableType,LabelSize,Label) BIND(C, &
    & NAME = "CMISSFieldVariableLabelGetCNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable character string label, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableLabelGetCCNum !<Error Code.
    !Local variable
    CHARACTER(LEN = LabelSize-1) :: FLabel

    CALL CMISSFieldVariableLabelGetCNum(RegionUserNumber,FieldUserNumber,VariableType,FLabel,CMISSFieldVariableLabelGetCCNum)
    CALL CMISSF2CString(FLabel,Label)

    RETURN

  END FUNCTION CMISSFieldVariableLabelGetCCNum

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a field variable for a field identified by an object.
  FUNCTION CMISSFieldVariableLabelGetCCPtr(FieldPtr,VariableType,LabelSize,Label) BIND(C, NAME = "CMISSFieldVariableLabelGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to get the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The field variable character string label, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableLabelGetCCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN = LabelSize - 1) :: FLabel


    CMISSFieldVariableLabelGetCCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldVariableLabelGet(Field, VariableType,FLabel, CMISSFieldVariableLabelGetCCPtr)
        CALL CMISSF2CString(FLabel, Label)
      ELSE
        CMISSFieldVariableLabelGetCCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldVariableLabelGetCCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldVariableLabelGetCCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the character string label for a field variable for a field identified by a user number for C.
  FUNCTION CMISSFieldVariableLabelSetCCNum(RegionUserNumber,FieldUserNumber,VariableType,LabelSize,Label) BIND(C, &
    & NAME ="CMISSFieldVariableLabelSetCNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number for the field to set the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to set the label for, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN = 1, KIND = C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field variable character string label, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableLabelSetCCNum !<Error Code.
    !Local variable
    CHARACTER(LEN = LabelSize -1) :: FLabel

    CALL CMISSC2FLabel(Label, FLabel)
    CALL CMISSFieldVariableLabelSet(RegionUserNumber,FieldUserNumber,VariableType,FLabel,CMISSFieldVariableLabelSetCCNum)

    RETURN

  END FUNCTION CMISSFieldVariableLabelSetCCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the character string label for a field variable for a field identified by an object for C.
  FUNCTION CMISSFieldVariableLabelSetCCPtr(FieldPtr,VariableType,LabelSize,Label) BIND(C, NAME = "CMISSFieldVariableLabelSet")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the label to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: VariableType !<The variable type for the field to set the label to, for C. \see OPENCMISS_FieldVariableTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The label size
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN) :: Label(LabelSize) !<The field variable character string label, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableLabelSetCCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    CHARACTER(LEN = LabelSize - 1) :: FLabel


    CMISSFieldVariableLabelSetCCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSC2FString(Label, FLabel)
        CALL CMISSFieldVariableLabelSet(Field, VariableType,FLabel, CMISSFieldVariableLabelSetCCPtr)
      ELSE
        CMISSFieldVariableLabelSetCCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldVariableLabelSetCCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldVariableLabelSetCCPtr

  !
  !================================================================================================================================
  !

  !>Returns the field variable types for a field identified by a user number for C.
  FUNCTION CMISSFieldVariableTypesGetCNum(RegionUserNumber,FieldUserNumber,VariableTypesSize,VariableTypesPtr) BIND(C, NAME = &
    & "CMISSFieldVariableTypesGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to get the field variable types for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to get the field variable types for, for C.
    INTEGER(C_INT), INTENT(OUT) :: VariableTypesSize !<Size of the variable types array.
    TYPE(C_PTR),  INTENT(OUT) :: VariableTypesPtr !<VariableTypes(variable_idx). C pointer to the field variable types for the variable_idx'th field variable. \see OPENCMISS_FieldVariableTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableTypesGetCNum !<Error Code.
    !Local variables
    INTEGER(C_INT), POINTER :: VariableTypes(:)

    CALL CMISSFieldVariableTypesGet(RegionUserNumber,FieldUserNumber,VariableTypes, CMISSFieldVariableTypesGetCNum)
    VariableTypesSize = Size(VariableTypes)
    VariableTypesPtr = C_LOC(VariableTypes(1))

    RETURN

  END FUNCTION CMISSFieldVariableTypesGetCNum

    !
  !================================================================================================================================
  !

  !>Returns the variable types for a field identified by an object for C.
  FUNCTION CMISSFieldVariableTypesGetCPtr(FieldPtr,VariableTypesSize,VariableTypesPtr) BIND(C, NAME = "CMISSFieldVariableTypesGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to get the field variable types for, for C.
    INTEGER(C_INT), INTENT(OUT) :: VariableTypesSize !<Size of the field variable types array.
    TYPE(C_PTR), INTENT(OUT) :: VariableTypesPtr !<VariableTypes(variable_idx). The field variable types for the variable_idx'th field variable, for C. \see OPENCMISS_FieldVariableTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableTypesGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    INTEGER(C_INT), POINTER :: VariableTypes(:)

    CMISSFieldVariableTypesGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        CALL CMISSFieldVariableTypesGet(Field, VariableTypes, CMISSFieldVariableTypesGetCPtr)
        VariableTypesSize = Size(VariableTypes)
        VariableTypesPtr = C_LOC(VariableTypes(1))
      ELSE
        CMISSFieldVariableTypesGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldVariableTypesGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldVariableTypesGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the field variable types for a field identified by a user number for C.
  FUNCTION CMISSFieldVariableTypesSetCNum(RegionUserNumber,FieldUserNumber,VariableTypesPtr) BIND(C, NAME = &
    & "CMISSFieldVariableTypesSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to set the field variable types to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to set the field variable types to, for C.
    TYPE(C_PTR), INTENT(IN) :: VariableTypesPtr !<VariableTypes(variable_idx). The field variable types for the variable_idx'th field variable, for C. \see OPENCMISS_FieldVariableTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableTypesSetCNum !<Error Code.
    !Local variables
    INTEGER(C_INT), POINTER :: VariableTypes

    IF(C_ASSOCIATED(VariableTypesPtr)) THEN
      CALL C_F_POINTER(VariableTypesPtr, VariableTypes)
      IF(ASSOCIATED(VariableTypes)) THEN
        CALL CMISSFieldVariableTypesSetNum(RegionUserNumber,FieldUserNumber,VariableTypes,CMISSFieldVariableTypesSetCNum)
      ELSE
        CMISSFieldVariableTypesSetCNum = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldVariableTypesSetCNum = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldVariableTypesSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the variable types for a field identified by an object for C.
  FUNCTION CMISSFieldVariableTypesSetCPtr(FieldPtr,VariableTypesPtr) BIND(C, NAME = "CMISSFieldVariableTypesSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldPtr !<C pointer to the field to set the field variable types to, for C.
    TYPE(C_PTR), INTENT(IN) :: VariableTypesPtr !<VariableTypes(variable_idx). The field variable types for the variable_idx'th field variable, for C. \see OPENCMISS_FieldVariableTypes
    !Function variable
    INTEGER(C_INT) :: CMISSFieldVariableTypesSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldType), POINTER :: Field
    INTEGER(C_INT), POINTER :: VariableTypes

    CMISSFieldVariableTypesSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(VariableTypesPtr)) THEN
          CALL C_F_POINTER(VariableTypesPtr, VariableTypes)
          IF(ASSOCIATED(VariableTypes)) THEN
            CALL CMISSFieldVariableTypesSetObj(Field, VariableTypes, CMISSFieldVariableTypesSetCPtr)
          ELSE
            CMISSFieldVariableTypesSetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSFieldVariableTypesSetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSFieldVariableTypesSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldVariableTypesSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldVariableTypesSetCPtr

!!==================================================================================================================================
!!
!! FIELD_IO_ROUTINES
!!
!!==================================================================================================================================

  !>Export element information for fields set identified by an object for C. \todo number method
  FUNCTION CMISSFieldIOElementsExportCCPtr(FieldsPtr,FileNameSize,FileName,MethodSize,Method) BIND(C, NAME = &
    & "CMISSFieldIOElementsExportC")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldsPtr !<The fields to export the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FileNameSize !< Size of the file name to export the elements to for C.
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN) :: FileName(FileNameSize) !<The file name to export the elements to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MethodSize !<Size of the export method name for C.
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN):: Method(MethodSize) !<The export method to use for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldIOElementsExportCCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldsType), POINTER :: Fields
    CHARACTER(LEN = FileNameSize -1 ) :: FFileName
    CHARACTER(LEN = MethodSize - 1) :: FMethod

    CMISSFieldIOElementsExportCCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldsPtr)) THEN
      CALL C_F_POINTER(FieldsPtr, Fields)
      IF(ASSOCIATED(Fields)) THEN
        CALL CMISSC2FString(FileName, FFileName)
        CALL CMISSC2FString(Method, FMethod)
        CALL CMISSFieldIOElementsExportC(Fields,FileName,Method, CMISSFieldIOElementsExportCCPtr)
      ELSE
        CMISSFieldIOElementsExportCCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldIOElementsExportCCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldIOElementsExportCCPtr

  !
  !================================================================================================================================
  !

  !>Export nodal information for fields set identified by an object for C. \todo number method
  FUNCTION CMISSFieldIONodesExportCCPtr(FieldsPtr,FileNameSize,FileName,MethodSize,Method) BIND(C, NAME = &
    &  "CMISSFieldIONodesExportC")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: FieldsPtr !<The fields to export the nodes for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FileNameSize !< Size of the file name to export the nodes to for C.
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN) :: FileName(FileNameSize) !<The file name to export the nodes to for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MethodSize !<Size of the export method name for C.
    CHARACTER(LEN=1, KIND = C_CHAR), INTENT(IN):: Method(MethodSize) !<The export method to use for C.
    !Function variable
    INTEGER(C_INT) :: CMISSFieldIONodesExportCCPtr !<Error Code.
    !Local variables
    TYPE(CMISSFieldsType), POINTER :: Fields
    CHARACTER(LEN = FileNameSize -1 ) :: FFileName
    CHARACTER(LEN = MethodSize - 1) :: FMethod

    CMISSFieldIONodesExportCCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldsPtr)) THEN
      CALL C_F_POINTER(FieldsPtr, Fields)
      IF(ASSOCIATED(Fields)) THEN
        CALL CMISSC2FString(FileName, FFileName)
        CALL CMISSC2FString(Method, FMethod)
        CALL CMISSFieldIONodesExportC(Fields,FileName,Method, CMISSFieldIONodesExportCCPtr)
      ELSE
        CMISSFieldIONodesExportCCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSFieldIONodesExportCCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSFieldIONodesExportCCPtr

!!==================================================================================================================================
!!
!! GENERATED_MESH_ROUTINES
!!
!!==================================================================================================================================

  !>Returns the basis for a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshBasisGetCNum(GeneratedMeshUserNumber,BasisUserNumber) BIND(C, NAME = &
    & "CMISSGeneratedMeshBasisGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to get the basis for, for C.
    INTEGER(C_INT), INTENT(OUT) :: BasisUserNumber !<The user number of the basis to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshBasisGetCNum !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshBasisGet(GeneratedMeshUserNumber,BasisUserNumber,CMISSGeneratedMeshBasisGetCNum)

    RETURN

  END FUNCTION CMISSGeneratedMeshBasisGetCNum

    !
  !================================================================================================================================
  !

  !>Returns the basis for a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshBasisGetCPtr(GeneratedMeshPtr,BasisPtr) BIND(C, NAME = "CMISSGeneratedMeshBasisGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr!<C pointer to the generated mesh to get the basis for.
    TYPE(C_PTR), INTENT(INOUT) :: BasisPtr !<C pointer to the basis to get.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshBasisGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh
    TYPE(CMISSBasisType), POINTER :: Basis

    CMISSGeneratedMeshBasisGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        IF(C_ASSOCIATED(BasisPtr)) THEN
          CALL C_F_POINTER(BasisPtr, Basis)
          IF(ASSOCIATED(Basis)) THEN
            CALL CMISSGeneratedMeshBasisGet(GeneratedMesh, Basis, CMISSGeneratedMeshBasisGetCPtr)
          ELSE
            CMISSGeneratedMeshBasisGetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSGeneratedMeshBasisGetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSGeneratedMeshBasisGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshBasisGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshBasisGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the basis for a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshBasisSetCNum(GeneratedMeshUserNumber,BasisUserNumber) BIND(C, NAME = &
    & "CMISSGeneratedMeshBasisSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to set the basis to, for C.
    INTEGER(C_INT), INTENT(IN) :: BasisUserNumber !<The user number of the basis to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshBasisSetCNum !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshBasisSet(GeneratedMeshUserNumber,BasisUserNumber,CMISSGeneratedMeshBasisSetCNum)

    RETURN

  END FUNCTION CMISSGeneratedMeshBasisSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the basis for a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshBasisSetCPtr(GeneratedMeshPtr,BasisPtr) BIND(C, NAME = "CMISSGeneratedMeshBasisSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr!<C pointer to the generated mesh to set the basis to.
    TYPE(C_PTR), INTENT(INOUT) :: BasisPtr !<C pointer to the basis to set.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshBasisSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh
    TYPE(CMISSBasisType), POINTER :: Basis

    CMISSGeneratedMeshBasisSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        IF(C_ASSOCIATED(BasisPtr)) THEN
          CALL C_F_POINTER(BasisPtr, Basis)
          IF(ASSOCIATED(Basis)) THEN
            CALL CMISSGeneratedMeshBasisSet(GeneratedMesh, Basis, CMISSGeneratedMeshBasisSetCPtr)
          ELSE
            CMISSGeneratedMeshBasisSetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSGeneratedMeshBasisSetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSGeneratedMeshBasisSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshBasisSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshBasisSetCPtr

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshCreateFinishCNum(GeneratedMeshUserNumber,MeshUserNumber) BIND(C, NAME = &
    & "CMISSGeneratedMeshCreateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to set the basis to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to generate, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshCreateFinishCNum !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshCreateFinish(GeneratedMeshUserNumber,MeshUserNumber,CMISSGeneratedMeshCreateFinishCNum)

    RETURN

  END FUNCTION CMISSGeneratedMeshCreateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a generated mesh identified by an object.
  FUNCTION CMISSGeneratedMeshCreateFinishCPtr(GeneratedMeshPtr,MeshUserNumber,MeshPtr) BIND(C, NAME = &
    & "CMISSGeneratedMeshCreateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr!<C pointer to the generated mesh to finish the creation of.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to generate, for C.
    TYPE(C_PTR), INTENT(INOUT) :: MeshPtr !<C pointer to the generated mesh.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshCreateFinishCPtr !<Error Code.
    !Local variables
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSGeneratedMeshCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        IF(C_ASSOCIATED(MeshPtr)) THEN
          CALL C_F_POINTER(MeshPtr, Mesh)
          IF(ASSOCIATED(Mesh)) THEN
            CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh, MeshUserNumber, Mesh, CMISSGeneratedMeshCreateFinishCPtr)
          ELSE
            CMISSGeneratedMeshCreateFinishCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSGeneratedMeshCreateFinishCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSGeneratedMeshCreateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshCreateStartCNum(GeneratedMeshUserNumber,RegionUserNumber) BIND(C, NAME = &
    & "CMISSGeneratedMeshCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to create, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region to create the generated mesh in, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshCreateStartCNum !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,RegionUserNumber,CMISSGeneratedMeshCreateStartCNum)

    RETURN

  END FUNCTION CMISSGeneratedMeshCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts the creation of a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshCreateStartCPtr(GeneratedMeshUserNumber,RegionPtr,GeneratedMeshPtr) BIND(C, NAME = &
    & "CMISSGeneratedMeshCreateStart")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to create, for C.
    TYPE(C_PTR), INTENT(INOUT) :: GeneratedMeshPtr !<C pointer to the generated mesh to finish the creation of.
    TYPE(C_PTR), INTENT(INOUT) :: RegionPtr !<C pointer to the region to created generated mesh in.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshCreateStartCPtr !<Error Code.
    !Local variables
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh
    TYPE(CMISSregionType), POINTER :: Region

    CMISSGeneratedMeshCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        IF(C_ASSOCIATED(RegionPtr)) THEN
          CALL C_F_POINTER(RegionPtr, Region)
          IF(ASSOCIATED(Region)) THEN
            CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber, Region, GeneratedMesh, CMISSGeneratedMeshCreateStartCPtr)
          ELSE
            CMISSGeneratedMeshCreateStartCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSGeneratedMeshCreateStartCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSGeneratedMeshCreateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshCreateStartCPtr

  !
  !================================================================================================================================
  !

  !>Destroys a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshDestroyCNum(GeneratedMeshUserNumber) BIND(C, NAME = "CMISSGeneratedMeshDestroyNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to destroy, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshDestroyCNum !<Error Code.
    !Local variableC

    CALL CMISSGeneratedMeshDestroy(GeneratedMeshUserNumber,CMISSGeneratedMeshDestroyCNum)

    RETURN

  END FUNCTION CMISSGeneratedMeshDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshDestroyCPtr(GeneratedMeshPtr) BIND(C, NAME = "CMISSGeneratedMeshDestroy")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: GeneratedMeshPtr !<C pointer to the generated mesh to destroy.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshDestroyCPtr !<Error Code.
    !Local variables
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshDestroy(GeneratedMesh, CMISSGeneratedMeshDestroyCPtr)
      ELSE
        CMISSGeneratedMeshDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshDestroyCPtr


  !
  !================================================================================================================================
  !

  !>Returns the type of a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshTypeGetCNum(GeneratedMeshUserNumber,GeneratedMeshType) BIND(C, NAME = &
    & "CMISSGeneratedMeshTypeGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to get the type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: GeneratedMeshType !<On return, the type of the generated mesh to get, for C. \see OPENCMISS_GeneratedMeshTypes
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeGetCNum !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshTypeGet(GeneratedMeshUserNumber, GeneratedMeshType, CMISSGeneratedMeshTypeGetCNum)

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the type of a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshTypeGetCPtr(GeneratedMeshPtr,GeneratedMeshType) BIND(C, NAME = "CMISSGeneratedMeshTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr !<C pointer to the generated mesh to get the type for.
    INTEGER(C_INT), INTENT(OUT) :: GeneratedMeshType !<On return, the type of the generated mesh to get, for C. \see OPENCMISS_GeneratedMeshTypes
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeGetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshTypeGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshTypeGet(GeneratedMesh, GeneratedMeshType, CMISSGeneratedMeshTypeGetCPtr)
      ELSE
        CMISSGeneratedMeshTypeGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshTypeSetCNum(GeneratedMeshUserNumber,GeneratedMeshType) BIND(C, NAME = &
    & "CMISSGeneratedMeshTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to set the type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshType !<On return, the type of the generated mesh to set, for C. \see OPENCMISS_GeneratedMeshTypes
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeSetCNum !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshTypeSet(GeneratedMeshUserNumber, GeneratedMeshType, CMISSGeneratedMeshTypeSetCNum)

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a generated mesh identified by an object.
  FUNCTION CMISSGeneratedMeshTypeSetCPtr(GeneratedMeshPtr,GeneratedMeshType) BIND(C, NAME = "CMISSGeneratedMeshTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: GeneratedMeshPtr !<C pointer to the generated mesh to set the type to.
    INTEGER(C_INT), INTENT(OUT) :: GeneratedMeshType !<On return, the type of the generated mesh to set, for C. \see OPENCMISS_GeneratedMeshTypes
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshTypeSetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
      CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
      IF(ASSOCIATED(GeneratedMesh)) THEN
        CALL CMISSGeneratedMeshTypeSet(GeneratedMesh, GeneratedMeshType, CMISSGeneratedMeshTypeSetCPtr)
      ELSE
        CMISSGeneratedMeshTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Calculates and sets the geometric field parameters for a generated mesh identified by a user number for C.
  FUNCTION CMISSGeneratedMeshGeometricParametersCalculateCNum(RegionUserNumber,FieldUserNumber,GeneratedMeshUserNumber) &
    &  BIND(C, NAME = "CMISSGeneratedMeshGeometricParametersCalculateNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the field to calculate the geometric parameters for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: FieldUserNumber !<The user number of the field to calculate the geometric parameters for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GeneratedMeshUserNumber !<The user number of the generated mesh to calculate the geometric parameters for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshGeometricParametersCalculateCNum !<Error Code.
    !Local variable

    CALL CMISSGeneratedMeshGeometricParametersCalculate(RegionUserNumber,FieldUserNumber,GeneratedMeshUserNumber, &
      & CMISSGeneratedMeshGeometricParametersCalculateCNum)

    RETURN

  END FUNCTION CMISSGeneratedMeshGeometricParametersCalculateCNum

  !
  !================================================================================================================================
  !

  !>Calculates and sets the geometric field parameters for a generated mesh identified by an object for C.
  FUNCTION CMISSGeneratedMeshGeometricParametersCalculateCPtr(FieldPtr,GeneratedMeshPtr)  BIND(C, NAME = &
    & "CMISSGeneratedMeshGeometricParametersCalculate")

    !Argument variables
    TYPE(C_PTR), INTENT(INOUT) :: FieldPtr !<The field to calculate the geometric parameters for, for C.
    TYPE(C_PTR), INTENT(IN) :: GeneratedMeshPtr !<The generated mesh to calculate the geometric parameters for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSGeneratedMeshGeometricParametersCalculateCPtr !<Error Code.
    !Local variable
    TYPE(CMISSFieldType), POINTER :: Field
    TYPE(CMISSGeneratedMeshType), POINTER :: GeneratedMesh

    CMISSGeneratedMeshGeometricParametersCalculateCPtr = CMISSNoError
    IF(C_ASSOCIATED(FieldPtr)) THEN
      CALL C_F_POINTER(FieldPtr, Field)
      IF(ASSOCIATED(Field)) THEN
        IF(C_ASSOCIATED(GeneratedMeshPtr)) THEN
          CALL C_F_POINTER(GeneratedMeshPtr, GeneratedMesh)
          IF(ASSOCIATED(GeneratedMesh)) THEN
            CALL CMISSGeneratedMeshGeometricParametersCalculate(Field,GeneratedMesh, &
              & CMISSGeneratedMeshGeometricParametersCalculateCPtr)
          ELSE
            CMISSGeneratedMeshGeometricParametersCalculateCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSGeneratedMeshGeometricParametersCalculateCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSGeneratedMeshGeometricParametersCalculateCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSGeneratedMeshGeometricParametersCalculateCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSGeneratedMeshGeometricParametersCalculateCPtr

!!==================================================================================================================================
!!
!! MESH_ROUTINES
!!
!!==================================================================================================================================

  !>Finishes the creation of a domain decomposition for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionCreateFinishCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber) BIND(C, NAME = &
    & "CMISSDecompositionCreateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to finish the decomposition for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to finish the decomposition for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to finish for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionCreateFinishCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionCreateFinish(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
      & CMISSDecompositionCreateFinishCNum)

    RETURN

  END FUNCTION CMISSDecompositionCreateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a domain decomposition for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionCreateFinishCPtr(DecompositionPtr) BIND(C, NAME = "CMISSDecompositionCreateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to finish creating.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionCreateFinishCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionCreateFinish(Decomposition, CMISSDecompositionCreateFinishCPtr)
      ELSE
        CMISSDecompositionCreateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts the creation of a domain decomposition for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionCreateStartCNum(DecompositionUserNumber,RegionUserNumber,MeshUserNumber)  BIND(C, NAME = &
    & "CMISSDecompositionCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to create for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to create the decomposition for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to create the decomposition for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionCreateStartCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionCreateStart(DecompositionUserNumber,RegionUserNumber,MeshUserNumber, &
      & CMISSDecompositionCreateStartCNum)

    RETURN

  END FUNCTION CMISSDecompositionCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts the creation of a domain decomposition for a decomposition identified by an object.
  FUNCTION CMISSDecompositionCreateStartCPtr(DecompositionUserNumber,MeshPtr,DecompositionPtr)  BIND(C, NAME = &
    & "CMISSDecompositionCreateStart")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to finish for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to create the decomposition for.
    TYPE(C_PTR), INTENT(OUT) :: DecompositionPtr !<C pointer to the decomposition to create.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionCreateStartCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSDecompositionCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSDecompositionCreateStart(DecompositionUserNumber, Mesh, Decomposition, CMISSDecompositionCreateStartCPtr)
        IF(ASSOCIATED(Decomposition)) THEN
          DecompositionPtr = C_LOC(Decomposition)
        ELSE
          CMISSDecompositionCreateStartCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSDecompositionCreateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionCreateStartCPtr

  !
  !================================================================================================================================
  !

  !>Destroys a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionDestroyCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber)  BIND(C, NAME = &
    & "CMISSDecompositionDestroyNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to destroy the decomposition for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to destroy the decomposition for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to destroy for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionDestroyCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionDestroy(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, CMISSDecompositionDestroyCNum)

    RETURN

  END FUNCTION CMISSDecompositionDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionDestroyCPtr(DecompositionPtr) BIND(C, NAME = "CMISSDecompositionDestroy")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to destroy.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionDestroyCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionDestroy(Decomposition, CMISSDecompositionDestroyCPtr)
      ELSE
        CMISSDecompositionDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Calculates the element domains for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionElementDomainCalculateCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber) BIND(C, NAME =&
    & "CMISSDecompositionElementDomainCalculateNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to calculate the element domains for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to calculate the element domains, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to calculate the element domains for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainCalculateCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionElementDomainCalculate(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
      & CMISSDecompositionElementDomainCalculateCNum)

    RETURN

  END FUNCTION CMISSDecompositionElementDomainCalculateCNum

  !
  !================================================================================================================================
  !

  !>Calculates the element domains for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionElementDomainCalculateCPtr(DecompositionPtr) BIND(C, NAME = &
    & "CMISSDecompositionElementDomainCalculate")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to calculate the element domains.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainCalculateCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionElementDomainCalculateCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionElementDomainCalculate(Decomposition, CMISSDecompositionElementDomainCalculateCPtr)
      ELSE
        CMISSDecompositionElementDomainCalculateCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionElementDomainCalculateCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionElementDomainCalculateCPtr

  !
  !================================================================================================================================
  !

  !>Returns the domain for a given element in a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionElementDomainGetCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & ElementUserNumber,Domain) BIND(C, NAME = "CMISSDecompositionElementDomainGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the element domain for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the element domain for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to get the element domain for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementUserNumber !<The user number of the element to get the domain for, for C.
    INTEGER(C_INT), INTENT(OUT) :: Domain !<The computational domain of the element, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainGetCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionElementDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,ElementUserNumber,Domain,&
      & CMISSDecompositionElementDomainGetCNum)

    RETURN

  END FUNCTION CMISSDecompositionElementDomainGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the domain for a given element in a decomposition identified by an object.
  FUNCTION CMISSDecompositionElementDomainGetCPtr(DecompositionPtr,ElementUserNumber,Domain)  BIND(C, NAME = &
    & "CMISSDecompositionElementDomainGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to get the domain for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementUserNumber !<The user number of the element to get the domain for, for C.
    INTEGER(C_INT), INTENT(OUT) :: Domain !<The computational domain of the element, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionElementDomainGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionElementDomainGet(Decomposition,ElementUserNumber,Domain,CMISSDecompositionElementDomainGetCPtr)
      ELSE
        CMISSDecompositionElementDomainGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionElementDomainGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionElementDomainGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the domain for a given element in a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionElementDomainSetCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & ElementUserNumber,Domain) BIND(C, NAME = "CMISSDecompositionElementDomainSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the element domain to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the element domain to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to set the element domain to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementUserNumber !<The user number of the element to set the domain to, for C.
    INTEGER(C_INT), INTENT(IN) :: Domain !<The computational domain of the element to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainSetCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionElementDomainSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,ElementUserNumber,Domain,&
      & CMISSDecompositionElementDomainSetCNum)

    RETURN

  END FUNCTION CMISSDecompositionElementDomainSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the domain for a given element in a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionElementDomainSetCPtr(DecompositionPtr,ElementUserNumber,Domain)  BIND(C, NAME = &
    & "CMISSDecompositionElementDomainSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to set the domain to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementUserNumber !<The user number of the element to set the domain to, for C.
    INTEGER(C_INT), INTENT(IN) :: Domain !<The computational domain of the element to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionElementDomainSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionElementDomainSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionElementDomainSet(Decomposition,ElementUserNumber,Domain,CMISSDecompositionElementDomainSetCPtr)
      ELSE
        CMISSDecompositionElementDomainSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionElementDomainSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionElementDomainSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number used for the decomposition of a mesh for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionMeshComponentGetCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & MeshComponentNumber) BIND(C, NAME = "CMISSDecompositionMeshComponentGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the mesh component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the mesh component for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to get the mesh component for, for C.
    INTEGER(C_INT), INTENT(OUT) :: MeshComponentNumber !<The mesh component number for the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionMeshComponentGetCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionMeshComponentGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,MeshComponentNumber,&
      & CMISSDecompositionMeshComponentGetCNum)

    RETURN

  END FUNCTION CMISSDecompositionMeshComponentGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the mesh component number used for the decomposition of a mesh for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionMeshComponentGetCPtr(DecompositionPtr,MeshComponentNumber) BIND(C, NAME = &
    & "CMISSDecompositionMeshComponentGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to get the mesh component for.
    INTEGER(C_INT), INTENT(OUT) :: MeshComponentNumber !<The mesh component number for the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionMeshComponentGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionMeshComponentGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionMeshComponentGet(Decomposition,MeshComponentNumber,CMISSDecompositionMeshComponentGetCPtr)
      ELSE
        CMISSDecompositionMeshComponentGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionMeshComponentGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionMeshComponentGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number used for the decomposition of a mesh for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionMeshComponentSetCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & MeshComponentNumber) BIND(C, NAME = "CMISSDecompositionMeshComponentSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the mesh component to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the mesh component to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to set the mesh component to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number for the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionMeshComponentSetCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionMeshComponentSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,MeshComponentNumber,&
      & CMISSDecompositionMeshComponentSetCNum)

    RETURN

  END FUNCTION CMISSDecompositionMeshComponentSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number used for the decomposition of a mesh for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionMeshComponentSetCPtr(DecompositionPtr,MeshComponentNumber) BIND(C, NAME = &
    & "CMISSDecompositionMeshComponentSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to set the mesh component to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number for the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionMeshComponentSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionMeshComponentSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionMeshComponentSet(Decomposition,MeshComponentNumber,CMISSDecompositionMeshComponentSetCPtr)
      ELSE
        CMISSDecompositionMeshComponentSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionMeshComponentSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionMeshComponentSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the number of domains for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionNumberOfDomainsGetCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & NumberOfDomains) BIND(C, NAME = "CMISSDecompositionNumberOfDomainsGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the number of domains for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the number of domains for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to get the number of domains for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfDomains !<The number of domains in the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionNumberOfDomainsGetCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionNumberOfDomainsGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,NumberOfDomains,&
      & CMISSDecompositionNumberOfDomainsGetCNum)

    RETURN

  END FUNCTION CMISSDecompositionNumberOfDomainsGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the number of domains for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionNumberOfDomainsGetCPtr(DecompositionPtr,NumberOfDomains) BIND(C, NAME = &
    & "CMISSDecompositionNumberOfDomainsGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to get the number of domains for.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfDomains !<The number of domains in the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionNumberOfDomainsGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionNumberOfDomainsGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionNumberOfDomainsGet(Decomposition,NumberOfDomains,CMISSDecompositionNumberOfDomainsGetCPtr)
      ELSE
        CMISSDecompositionNumberOfDomainsGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionNumberOfDomainsGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionNumberOfDomainsGetCPtr

    !
  !================================================================================================================================
  !

  !>Sets/changes the number of domains for a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionNumberOfDomainsSetCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
    & NumberOfDomains) BIND(C, NAME = "CMISSDecompositionNumberOfDomainsSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the number of domains to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the number of domains to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to set the number of domains to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfDomains !<The number of domains in the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionNumberOfDomainsSetCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionNumberOfDomainsSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,NumberOfDomains,&
      & CMISSDecompositionNumberOfDomainsSetCNum)

    RETURN

  END FUNCTION CMISSDecompositionNumberOfDomainsSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of domains for a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionNumberOfDomainsSetCPtr(DecompositionPtr,NumberOfDomains) BIND(C, NAME = &
    & "CMISSDecompositionNumberOfDomainsSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to set the number of domains to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfDomains !<The number of domains in the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionNumberOfDomainsSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionNumberOfDomainsSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,CMISSDecompositionNumberOfDomainsSetCPtr)
      ELSE
        CMISSDecompositionNumberOfDomainsSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionNumberOfDomainsSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionNumberOfDomainsSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the type of a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionTypeGetCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,DecompositionType) BIND(C, NAME&
    & = "CMISSDecompositionTypeGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the decomposition type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the decomposition type for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to get the decomposition type for, for C.
    INTEGER(C_INT), INTENT(OUT) :: DecompositionType !<The type of the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionTypeGetCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionTypeGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,DecompositionType,&
      & CMISSDecompositionTypeGetCNum)

    RETURN

  END FUNCTION CMISSDecompositionTypeGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the type of a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionTypeGetCPtr(DecompositionPtr,DecompositionType) BIND(C, NAME = "CMISSDecompositionTypeGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to get the decomposition type for.
    INTEGER(C_INT), INTENT(OUT) :: DecompositionType !<The type of the decomposition for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionTypeGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionTypeGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionTypeGet(Decomposition,DecompositionType,CMISSDecompositionTypeGetCPtr)
      ELSE
        CMISSDecompositionTypeGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionTypeGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionTypeGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a decomposition identified by a user number for C.
  FUNCTION CMISSDecompositionTypeSetCNum(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,DecompositionType) BIND(C, NAME&
    & = "CMISSDecompositionTypeSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the decomposition type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the decomposition type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionUserNumber !<The user number of the decomposition to set the decomposition type to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: DecompositionType !<The type of the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionTypeSetCNum !<Error Code.
    !Local variable

    CALL CMISSDecompositionTypeSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,DecompositionType,&
      &CMISSDecompositionTypeSetCNum)

    RETURN

  END FUNCTION CMISSDecompositionTypeSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of a decomposition identified by an object for C.
  FUNCTION CMISSDecompositionTypeSetCPtr(DecompositionPtr,DecompositionType) BIND(C, NAME = "CMISSDecompositionTypeSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: DecompositionPtr !<C pointer to the decomposition to set the decomposition type to.
    INTEGER(C_INT), INTENT(OUT) :: DecompositionType !<The type of the decomposition to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSDecompositionTypeSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSDecompositionType), POINTER :: Decomposition

    CMISSDecompositionTypeSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(DecompositionPtr)) THEN
      CALL C_F_POINTER(DecompositionPtr, Decomposition)
      IF(ASSOCIATED(Decomposition)) THEN
        CALL CMISSDecompositionTypeSet(Decomposition,DecompositionType,CMISSDecompositionTypeSetCPtr)
      ELSE
        CMISSDecompositionTypeSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSDecompositionTypeSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSDecompositionTypeSetCPtr

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a mesh for a mesh identified by a user number for C.
  FUNCTION CMISSMeshCreateFinishCNum(RegionUserNumber,MeshUserNumber) BIND(C, NAME = "CMISSMeshCreateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to finish the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to finish the creation of, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshCreateFinishCNum !<Error Code.
    !Local variable

    CALL CMISSMeshCreateFinish(RegionUserNumber,MeshUserNumber,CMISSMeshCreateFinishCNum)

    RETURN

  END FUNCTION CMISSMeshCreateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a mesh for a mesh identified by an object for C.
  FUNCTION CMISSMeshCreateFinishCPtr(MeshPtr) BIND(C, NAME = "CMISSMeshCreateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to finish creating.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshCreateFinishCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshCreateFinish(Mesh,CMISSMeshCreateFinishCPtr)
      ELSE
        CMISSMeshCreateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts the creation of a mesh for a mesh identified by a user number for C.
  FUNCTION CMISSMeshCreateStartCNum(MeshUserNumber,RegionUserNumber,NumberOfDimensions) BIND(C, NAME = &
    &"CMISSMeshCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to start the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to start the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfDimensions !<The number of dimensions for the mesh, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshCreateStartCNum !<Error Code.
    !Local variable

    CALL CMISSMeshCreateStart(MeshUserNumber,RegionUserNumber,NumberOfDimensions,CMISSMeshCreateStartCNum)

    RETURN

  END FUNCTION CMISSMeshCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts the creation of a mesh for a mesh identified by an object for C.
  FUNCTION CMISSMeshCreateStartCPtr(MeshUserNumber,RegionPtr,NumberOfDimensions,MeshPtr) BIND(C, NAME = "CMISSMeshCreateStart")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to start the creation of, for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr !<C pointer to the region containing the mesh to start the creation of.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfDimensions !<The number of dimensions for the mesh, for C.
    TYPE(C_PTR), INTENT(OUT) :: MeshPtr !<C pointer to the created mesh.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshCreateStartCPtr !<Error Code.
    !Local variables
    TYPE(CMISSRegionType), POINTER :: Region
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr, Region)
      IF(ASSOCIATED(Region)) THEN
        IF(C_ASSOCIATED(MeshPtr)) THEN
          CALL C_F_POINTER(MeshPtr, Mesh)
          IF(ASSOCIATED(Mesh)) THEN
            CALL CMISSMeshCreateStart(MeshUserNumber,Region,NumberOfDimensions,Mesh,CMISSMeshCreateStartCPtr)
          ELSE
            CMISSMeshCreateStartCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSMeshCreateStartCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSMeshCreateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshCreateStartCPtr

  !
  !================================================================================================================================
  !

  !>Destroys a mesh identified by a user number for C.
  FUNCTION CMISSMeshDestroyCNum(RegionUserNumber,MeshUserNumber) BIND(C, NAME = "CMISSMeshDestroyNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to destroy, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to destroy, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshDestroyCNum !<Error Code.
    !Local variable

    CALL CMISSMeshDestroy(RegionUserNumber,MeshUserNumber,CMISSMeshDestroyCNum)

    RETURN

  END FUNCTION CMISSMeshDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys a mesh identified by an object for C.
  FUNCTION CMISSMeshDestroyCPtr(MeshPtr) BIND(C, NAME = "CMISSMeshDestroy")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to destroy.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshDestroyCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshDestroy(Mesh,CMISSMeshDestroyCPtr)
      ELSE
        CMISSMeshDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Returns the number of components in a mesh identified by a user number for C.
  FUNCTION CMISSMeshNumberOfComponentsGetCNum(RegionUserNumber,MeshUserNumber,NumberOfComponents) BIND(C, NAME = &
    &"CMISSMeshNumberOfComponentsGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the number of components for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the number of components for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfComponents !<The number of components in the mesh for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfComponentsGetCNum !<Error Code.
    !Local variable

    CALL CMISSMeshNumberOfComponentsGet(RegionUserNumber,MeshUserNumber,NumberOfComponents,&
      &CMISSMeshNumberOfComponentsGetCNum)

    RETURN

  END FUNCTION CMISSMeshNumberOfComponentsGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the number of components in a mesh identified by an object for C.
  FUNCTION CMISSMeshNumberOfComponentsGetCPtr(MeshPtr,NumberOfComponents) BIND(C, NAME = "CMISSMeshNumberOfComponentsGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to get the number of components for.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfComponents !<The number of components in the mesh for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfComponentsGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshNumberOfComponentsGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshNumberOfComponentsGet(Mesh,NumberOfComponents,CMISSMeshNumberOfComponentsGetCPtr)
      ELSE
        CMISSMeshNumberOfComponentsGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshNumberOfComponentsGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshNumberOfComponentsGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of components in a mesh identified by a user number for C.
  FUNCTION CMISSMeshNumberOfComponentsSetCNum(RegionUserNumber,MeshUserNumber,NumberOfComponents) BIND(C, NAME =&
    & "CMISSMeshNumberOfComponentsSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the number of components to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the number of components to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfComponents !<The number of components in the mesh to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfComponentsSetCNum !<Error Code.
    !Local variable

    CALL CMISSMeshNumberOfComponentsSet(RegionUserNumber,MeshUserNumber,NumberOfComponents,&
      &CMISSMeshNumberOfComponentsSetCNum)

    RETURN

  END FUNCTION CMISSMeshNumberOfComponentsSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of components in a mesh identified by an object for C.
  FUNCTION CMISSMeshNumberOfComponentsSetCPtr(MeshPtr,NumberOfComponents) BIND(C, NAME = "CMISSMeshNumberOfComponentsSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to set the number of components to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfComponents !<The number of components in the mesh to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfComponentsSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshNumberOfComponentsSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshNumberOfComponentsSet(Mesh,NumberOfComponents,CMISSMeshNumberOfComponentsSetCPtr)
      ELSE
        CMISSMeshNumberOfComponentsSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshNumberOfComponentsSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshNumberOfComponentsSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the number of elements in a mesh identified by a user number for C.
  FUNCTION CMISSMeshNumberOfElementsGetCNum(RegionUserNumber,MeshUserNumber,NumberOfElements) BIND(C, NAME = &
    &"CMISSMeshNumberOfElementsGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the number of elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the number of elements for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfElements !<The number of elements in the mesh for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfElementsGetCNum !<Error Code.
    !Local variable

    CALL CMISSMeshNumberOfElementsGet(RegionUserNumber,MeshUserNumber,NumberOfElements,CMISSMeshNumberOfElementsGetCNum)

    RETURN

  END FUNCTION CMISSMeshNumberOfElementsGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the number of elements in a mesh identified by an object for C.
  FUNCTION CMISSMeshNumberOfElementsGetCPtr(MeshPtr,NumberOfElements) BIND(C, NAME = "CMISSMeshNumberOfElementsGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to get the number of elements for.
    INTEGER(C_INT), INTENT(OUT) :: NumberOfElements !<The number of elements in the mesh for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfElementsGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshNumberOfElementsGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshNumberOfElementsGet(Mesh,NumberOfElements,CMISSMeshNumberOfElementsGetCPtr)
      ELSE
        CMISSMeshNumberOfElementsGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshNumberOfElementsGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshNumberOfElementsGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of elements in a mesh identified by a user number for C.
  FUNCTION CMISSMeshNumberOfElementsSetCNum(RegionUserNumber,MeshUserNumber,NumberOfElements) BIND(C, NAME = &
    &"CMISSMeshNumberOfElementsSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the number of elements to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the number of elements to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfElements !<The number of elements in the mesh to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfElementsSetCNum !<Error Code.
    !Local variable

    CALL CMISSMeshNumberOfElementsSet(RegionUserNumber,MeshUserNumber,NumberOfElements,CMISSMeshNumberOfElementsSetCNum)

    RETURN

  END FUNCTION CMISSMeshNumberOfElementsSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of elements in a mesh identified by an object for C.
  FUNCTION CMISSMeshNumberOfElementsSetCPtr(MeshPtr,NumberOfElements) BIND(C, NAME = "CMISSMeshNumberOfElementsSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to set the number of elements to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfElements !<The number of elements in the mesh to set for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshNumberOfElementsSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh

    CMISSMeshNumberOfElementsSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        CALL CMISSMeshNumberOfElementsSet(Mesh,NumberOfElements,CMISSMeshNumberOfElementsSetCPtr)
      ELSE
        CMISSMeshNumberOfElementsSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshNumberOfElementsSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshNumberOfElementsSetCPtr

  !
  !================================================================================================================================
  !

  !>Finishes creating elements for a mesh component of a mesh identified by a user number for C.
  FUNCTION CMISSMeshElementsCreateFinishCNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber) BIND(C, NAME = &
    &"CMISSMeshElementsCreateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to finish creating the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to finish creating the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to finish creating the elements for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsCreateFinishCNum !<Error Code.
    !Local variable

    CALL CMISSMeshElementsCreateFinish(RegionUserNumber,MeshUserNumber,MeshComponentNumber,&
      &CMISSMeshElementsCreateFinishCNum)

    RETURN

  END FUNCTION CMISSMeshElementsCreateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes creating elements for a mesh component of a mesh identified by an object for C.
  FUNCTION CMISSMeshElementsCreateFinishCPtr(MeshElementsPtr) BIND(C, NAME = "CMISSMeshElementsCreateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshElementsPtr !<C pointer the mesh elements to finish creating.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsCreateFinishCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements

    CMISSMeshElementsCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsPtr)) THEN
      CALL C_F_POINTER(MeshElementsPtr, MeshElements)
      IF(ASSOCIATED(MeshElements)) THEN
        CALL CMISSMeshElementsCreateFinish(MeshElements,CMISSMeshElementsCreateFinishCPtr)
      ELSE
        CMISSMeshElementsCreateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts creating elements for a mesh component of a mesh identified by a user number for C.
  FUNCTION CMISSMeshElementsCreateStartCNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,BasisUserNumber) BIND(C, NAME = &
    &"CMISSMeshElementsCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to start creating the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to start creating the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to start creating the elements for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: BasisUserNumber !<The user number of the default basis to use for the elements, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsCreateStartCNum !<Error Code.
    !Local variable

    CALL CMISSMeshElementsCreateStart(RegionUserNumber,MeshUserNumber,MeshComponentNumber,BasisUserNumber,&
      &CMISSMeshElementsCreateStartCNum)

    RETURN

  END FUNCTION CMISSMeshElementsCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts creating elements for a mesh component of a mesh identified by an object, for C.
  FUNCTION CMISSMeshElementsCreateStartCPtr(MeshPtr,MeshComponentNumber,BasisPtr,MeshElementsPtr) BIND(C, NAME = &
    &"CMISSMeshElementsCreateStart")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshPtr !<C pointer to the mesh to start the creation of elements for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to start creating the elements for, for C.
    TYPE(C_PTR), INTENT(IN) :: BasisPtr !<C pointer to the default basis to use for the elements.
    TYPE(C_PTR), INTENT(OUT) :: MeshElementsPtr !<C pointer to the created mesh elements.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsCreateStartCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshType), POINTER :: Mesh
    TYPE(CMISSBasisType), POINTER :: Basis
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements

    CMISSMeshElementsCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshPtr)) THEN
      CALL C_F_POINTER(MeshPtr, Mesh)
      IF(ASSOCIATED(Mesh)) THEN
        IF(C_ASSOCIATED(BasisPtr)) THEN
          CALL C_F_POINTER(BasisPtr, Basis)
          IF(ASSOCIATED(Basis)) THEN
            CALL CMISSMeshElementsCreateStart(Mesh, MeshComponentNumber, Basis, MeshElements, &
                  & CMISSMeshElementsCreateStartCPtr)
            IF(ASSOCIATED(MeshElements)) THEN
              MeshElementsPtr = C_LOC(MeshElements)
            ELSE
              CMISSMeshElementsCreateStartCPtr = CMISSPointerIsNULL
            ENDIF
          ELSE
            CMISSMeshElementsCreateStartCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSMeshElementsCreateStartCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSMeshElementsCreateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsCreateStartCPtr

  !
  !================================================================================================================================
  !

  !>Returns the basis for an element in a mesh identified by an user number, for C. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsBasisGetCNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber, &
      & BasisUserNumber) BIND(C, NAME = "CMISSMeshElementsBasisGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the basis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the basis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to get the basis for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to get the basis for, for C.
    INTEGER(C_INT), INTENT(OUT) :: BasisUserNumber !<The user number of the basis for the element, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsBasisGetCNum !<Error Code.
    !Local variable

    CALL CMISSMeshElementsBasisGet(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber,BasisUserNumber,&
      &CMISSMeshElementsBasisGetCNum)

    RETURN

  END FUNCTION CMISSMeshElementsBasisGetCNum

 !
  !================================================================================================================================
  !

  !>Returns the basis for an element in a mesh identified by an object, for C. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsBasisGetCPtr(MeshElementsPtr,GlobalElementNumber,BasisPtr) BIND(C, NAME = " CMISSMeshElementsBasisGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshElementsPtr !<C pointer to the mesh elements to get, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to get the basis for, for C.
    TYPE(C_PTR), INTENT(IN) :: BasisPtr !<C pointer to the basis for the element to get, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsBasisGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements
    TYPE(CMISSBasisType), POINTER :: Basis

    CMISSMeshElementsBasisGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsPtr)) THEN
      CALL C_F_POINTER(MeshElementsPtr, MeshElements)
      IF(ASSOCIATED(MeshElements)) THEN
        IF(C_ASSOCIATED(BasisPtr)) THEN
          CALL C_F_POINTER(BasisPtr, Basis)
          IF(ASSOCIATED(Basis)) THEN
            CALL CMISSMeshElementsBasisGet(MeshElements, GlobalElementNumber, Basis, &
              & CMISSMeshElementsBasisGetCPtr)
          ELSE
            CMISSMeshElementsBasisGetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSMeshElementsBasisGetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSMeshElementsBasisGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsBasisGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsBasisGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the basis for an element in a mesh identified by an user number, for C. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsBasisSetCNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber, &
    & BasisUserNumber) BIND(C, NAME = "CMISSMeshElementsBasisSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the basis to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the basis to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to set the basis to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to set the basis to, for C.
    INTEGER(C_INT), INTENT(IN) :: BasisUserNumber !<The user number of the basis for the element to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsBasisSetCNum !<Error Code.
    !Local variable

    CALL CMISSMeshElementsBasisSet(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber,BasisUserNumber,&
      &CMISSMeshElementsBasisSetCNum)

    RETURN

  END FUNCTION CMISSMeshElementsBasisSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the basis for an element in a mesh identified by an object, for C. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsBasisSetCPtr(MeshElementsPtr,GlobalElementNumber,BasisPtr) BIND(C, NAME = " CMISSMeshElementsBasisSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshElementsPtr !<C pointer to the mesh elements to set, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to set the basis for, for C.
    TYPE(C_PTR), INTENT(IN) :: BasisPtr !<C pointer to the basis for the element to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsBasisSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements
    TYPE(CMISSBasisType), POINTER :: Basis

    CMISSMeshElementsBasisSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsPtr)) THEN
      CALL C_F_POINTER(MeshElementsPtr, MeshElements)
      IF(ASSOCIATED(MeshElements)) THEN
        IF(C_ASSOCIATED(BasisPtr)) THEN
          CALL C_F_POINTER(BasisPtr, Basis)
          IF(ASSOCIATED(Basis)) THEN
            CALL CMISSMeshElementsBasisSet(MeshElements, GlobalElementNumber, Basis, &
              & CMISSMeshElementsBasisSetCPtr)
          ELSE
            CMISSMeshElementsBasisSetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSMeshElementsBasisSetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSMeshElementsBasisSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsBasisSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsBasisSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the element nodes for an element in a mesh identified by an user number for C. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsNodesGetCNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber, &
    & ElementUserNodesSize,ElementUserNodesPtr) BIND(C, NAME = "CMISSMeshElementsNodesGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the element nodes for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the element nodes for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to get the element nodes for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to get the element nodes for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ElementUserNodesSize !<Size of the element user node number array.
    TYPE(C_PTR), INTENT(OUT) :: ElementUserNodesPtr !<C pointer to location of the user node number.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsNodesGetCNum !<Error Code.
    !Local variable
    INTEGER(C_INT), POINTER :: ElementUserNodes(:)

    CMISSMeshElementsNodesGetCNum = CMISSNoError
    IF(C_ASSOCIATED(ElementUserNodesPtr)) THEN
      CMISSMeshElementsNodesGetCNum = CMISSPointerNotNULL
    ELSE
      CALL CMISSMeshElementsNodesGet(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber,ElementUserNodes, &
      &CMISSMeshElementsNodesGetCNum)
      ElementUserNodesSize = Size(ElementUserNodes,1)
      ElementUserNodesPtr = C_LOC(ElementUserNodes(1))
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsNodesGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the element nodes for an element in a mesh identified by an object, for C. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsNodesGetCPtr(MeshElementsPtr,GlobalElementNumber,ElementUserNodesSize,ElementUserNodesPtr) BIND(C, &
    &NAME = "CMISSMeshElementsNodesGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshElementsPtr !<The mesh elements to get the element nodes for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to get the element nodes for.
    INTEGER(C_INT), INTENT(OUT) :: ElementUserNodesSize !<Size of the element user node number array.
    TYPE(C_PTR), INTENT(OUT) :: ElementUserNodesPtr !<C pointer to location of the user node number.
    !Function variables
    INTEGER(C_INT) :: CMISSMeshElementsNodesGetCPtr !<Error code.
    !Local variables
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements
    INTEGER(C_INT), POINTER :: ElementUserNodes(:)

    CMISSMeshElementsNodesGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsPtr)) THEN
      CALL C_F_POINTER(MeshElementsPtr, MeshElements)
      IF(ASSOCIATED(MeshElements)) THEN
        IF(C_ASSOCIATED(ElementUserNodesPtr)) THEN
          CMISSMeshElementsNodesGetCPtr = CMISSPointerNotNULL
        ELSE
          CALL CMISSMeshElementsNodesGet(MeshElements, GlobalElementNumber, ElementUserNodes, &
            & CMISSMeshElementsNodesGetCPtr)
          ElementUserNodesSize = Size(ElementUserNodes,1)
          ElementUserNodesPtr = C_LOC(ElementUserNodes(1))
        ENDIF
      ELSE
        CMISSMeshElementsNodesGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsNodesGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsNodesGetCPtr



  !
  !================================================================================================================================
  !

  !>Sets/changes the element nodes for an element in a mesh identified by an user number, for C. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsNodesSetCNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber, &
    & ElementUserNodesSize,ElementUserNodesPtr) BIND(C, NAME = "CMISSMeshElementsNodesSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the element nodes to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the element nodes to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number of the mesh to set the element nodes to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to set the element nodes to, for C.
    INTEGER(C_INT), INTENT(IN) :: ElementUserNodesSize(1)
    TYPE(C_PTR), INTENT(IN) :: ElementUserNodesPtr !<C pointer to location of the user node number.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsNodesSetCNum !<Error Code.
    !Local variable
    INTEGER(C_INT), POINTER :: ElementUserNodes(:)

    CMISSMeshElementsNodesSetCNum = CMISSNoError
    IF(C_ASSOCIATED(ElementUserNodesPtr)) THEN
      CMISSMeshElementsNodesSetCNum = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(ElementUserNodesPtr)) THEN
        CALL C_F_POINTER(ElementUserNodesPtr, ElementUserNodes, ElementUserNodesSize)
        IF(ASSOCIATED(ElementUserNodes)) THEN
          CALL CMISSMeshElementsNodesSet(RegionUserNumber,MeshUserNumber,MeshComponentNumber,GlobalElementNumber,ElementUserNodes, &
            &CMISSMeshElementsNodesSetCNum)
        ELSE
          CMISSMeshElementsNodesSetCNum = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSMeshElementsNodesSetCNum = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsNodesSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the element nodes for an element in a mesh identified by an object for C. \todo should the global element number be a user number?
  FUNCTION CMISSMeshElementsNodesSetCPtr(MeshElementsPtr,GlobalElementNumber,ElementUserNodesSize,ElementUserNodesPtr) BIND(C, &
    & NAME="CMISSMeshElementsNodesSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshElementsPtr !<The mesh elements to get the element nodes for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: GlobalElementNumber !<The global element number to get the element nodes for.
    INTEGER(C_INT), INTENT(IN) :: ElementUserNodesSize(1)
    TYPE(C_PTR), INTENT(IN) :: ElementUserNodesPtr !<C pointer to location of the user node number.
    !Function variables
    INTEGER(C_INT) :: CMISSMeshElementsNodesSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements
    INTEGER(C_INT), POINTER :: ElementUserNodes(:)

    CMISSMeshElementsNodesSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsPtr)) THEN
      CALL C_F_POINTER(MeshElementsPtr, MeshElements)
      IF(ASSOCIATED(MeshElements)) THEN
        IF(C_ASSOCIATED(ElementUserNodesPtr)) THEN
        CALL C_F_POINTER(ElementUserNodesPtr, ElementUserNodes, ElementUserNodesSize)
          IF(ASSOCIATED(ElementUserNodes)) THEN
            CALL CMISSMeshElementsNodesSet(MeshElements, GlobalElementNumber, ElementUserNodes, &
              & CMISSMeshElementsNodesSetCPtr)
          ELSE
            CMISSMeshElementsNodesSetCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSMeshElementsNodesSetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSMeshElementsNodesSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsNodesSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsNodesSetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the user number for an element in a mesh identified by an user number, for C.
  FUNCTION CMISSMeshElementsUserNumberGetCNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,ElementGlobalNumber, &
    & ElementUserNumber) BIND(C, NAME = "CMISSMeshElementsUserNumberGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to get the element user number for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to get the element user number for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number to get the element user number for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementGlobalNumber !<The global element number to get the element user number for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ElementUserNumber !<The element user number, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSMeshElementsUserNumberGetCNum
    !Local variables

    CALL CMISSMeshElementsUserNumberGet(RegionUserNumber,MeshUserNumber,MeshComponentNumber,ElementGlobalNumber,ElementUserNumber, &
      & CMISSMeshElementsUserNumberGetCNum)

    RETURN

  END FUNCTION CMISSMeshElementsUserNumberGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the element user number for an element in a mesh identified by an object, for C.
  FUNCTION CMISSMeshElementsUserNumberGetCPtr(MeshElementsPtr,ElementGlobalNumber,ElementUserNumber) BIND(C, NAME = &
    & "CMISSMeshElementsUserNumberGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshElementsPtr !<C pointer to the mesh elements to get the element nodes for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementGlobalNumber !<The global element number to get the element user number for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ElementUserNumber !<The element user number, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsUserNumberGetCPtr
    !Local variables
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements

    CMISSMeshElementsUserNumberGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsPtr)) THEN
      CALL C_F_POINTER(MeshElementsPtr, MeshElements)
      IF(ASSOCIATED(MeshElements)) THEN
        CALL CMISSMeshElementsUserNumberGet(MeshElements,ElementGlobalNumber,ElementUserNumber,CMISSMeshElementsUserNumberGetCPtr)
      ELSE
        CMISSMeshElementsUserNumberGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsUserNumberGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsUserNumberGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the user number for an element in a mesh identified by an user number, for C.
  FUNCTION CMISSMeshElementsUserNumberSetCNum(RegionUserNumber,MeshUserNumber,MeshComponentNumber,ElementGlobalNumber, &
    & ElementUserNumber)  BIND(C, NAME = "CMISSMeshElementsUserNumberSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the mesh to set the element user number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshUserNumber !<The user number of the mesh to set the element user number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: MeshComponentNumber !<The mesh component number to set the element user number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementGlobalNumber !<The global element number to set the element user number to, for C.
    INTEGER(C_INT), INTENT(IN) :: ElementUserNumber !<The element user number to set, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSMeshElementsUserNumberSetCNum
    !Local variables

    CALL CMISSMeshElementsUserNumberSet(RegionUserNumber,MeshUserNumber,MeshComponentNumber,ElementGlobalNumber,ElementUserNumber, &
      & CMISSMeshElementsUserNumberSetCNum)

    RETURN

  END FUNCTION CMISSMeshElementsUserNumberSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the element user number for an element in a mesh identified by an object, for C.
  FUNCTION CMISSMeshElementsUserNumberSetCPtr(MeshElementsPtr,ElementGlobalNumber,ElementUserNumber) BIND(C, NAME = &
    & "CMISSMeshElementsUserNumberSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: MeshElementsPtr !<C pointer to the mesh elements to get the element nodes for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ElementGlobalNumber !<The global element number to set the element user number to, for C.
    INTEGER(C_INT), INTENT(IN) :: ElementUserNumber !<The element user number to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSMeshElementsUserNumberSetCPtr
    !Local variables
    TYPE(CMISSMeshElementsType), POINTER :: MeshElements

    CMISSMeshElementsUserNumberSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(MeshElementsPtr)) THEN
      CALL C_F_POINTER(MeshElementsPtr, MeshElements)
      IF(ASSOCIATED(MeshElements)) THEN
        CALL CMISSMeshElementsUserNumberSet(MeshElements,ElementGlobalNumber,ElementUserNumber,CMISSMeshElementsUserNumberSetCPtr)
      ELSE
        CMISSMeshElementsUserNumberSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSMeshElementsUserNumberSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSMeshElementsUserNumberSetCPtr

!!==================================================================================================================================
!!
!! NODE_ROUTINES
!!
!!==================================================================================================================================

  !>Finishes the process of creating nodes in a region for nodes identified by user number, for C.
  FUNCTION CMISSNodesCreateFinishCNum(RegionUserNumber) BIND(C, NAME= "CMISSNodesCreateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the nodes to finish the creation of, for C.
    !Function variables
    INTEGER(C_INT) :: CMISSNodesCreateFinishCNum !<Error Code.
    !Local variables

    CALL CMISSNodesCreateFinish(RegionUserNumber, CMISSNodesCreateFinishCNum)

    RETURN

  END FUNCTION CMISSNodesCreateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a nodes in a region for nodes identified by an object for C.
  FUNCTION CMISSNodesCreateFinishCPtr(NodesPtr) BIND(C, NAME = "CMISSNodesCreateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: NodesPtr !<C pointer to the nodes to finish creating.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesCreateFinishCPtr !<Error Code.
    !Local variables
    TYPE(CMISSNodesType), POINTER :: Nodes

    CMISSNodesCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(NodesPtr)) THEN
      CALL C_F_POINTER(NodesPtr, Nodes)
      IF(ASSOCIATED(Nodes)) THEN
        CALL CMISSNodesCreateFinish(Nodes, CMISSNodesCreateFinishCPtr)
      ELSE
        CMISSNodesCreateFinishCPtr= CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSNodesCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSNodesCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts the process of creating nodes in a region for nodes identified by user number.
  FUNCTION CMISSNodesCreateStartCNum(RegionUserNumber,NumberOfNodes) BIND(C, NAME = "CMISSNodesCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the nodes to start the creation of, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfNodes !<The number of nodes to create, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesCreateStartCNum
    !Local variable

    CALL CMISSNodesCreateStart(RegionUserNumber, NumberOfNodes, CMISSNodesCreateStartCNum)

    RETURN

  END FUNCTION CMISSNodesCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts the creation of a nodes in a region for nodes identified by an object for C.
  FUNCTION CMISSNodesCreateStartCPtr(RegionPtr,NumberOfNodes,NodesPtr) BIND(C, NAME = "CMISSNodesCreateStart")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr !<C pointer to the region to start the creation of nodes on.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NumberOfNodes !<The number of nodes to create for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: NodesPtr !<C pointer to the nodes created.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesCreateStartCPtr !<Error Code.
    !Local variables
    TYPE(CMISSRegionType), POINTER :: Region
    TYPE(CMISSNodesType), POINTER :: Nodes

    CMISSNodesCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(NodesPtr)) THEN
      CALL C_F_POINTER(NodesPtr, Nodes)
      IF(ASSOCIATED(Nodes)) THEN
        IF(C_ASSOCIATED(RegionPtr)) THEN
          CALL C_F_POINTER(RegionPtr, Region)
          IF(ASSOCIATED(Region)) THEN
            CALL CMISSNodesCreateStart(Region, NumberOfNodes, Nodes, CMISSNodesCreateStartCPtr)
          ELSE
            CMISSNodesCreateStartCPtr = CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSNodesCreateStartCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSNodesCreateStartCPtr= CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSNodesCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSNodesCreateStartCPtr

  !
  !================================================================================================================================
  !

  !>Destroys the nodes in a region for nodes identified by user number for C.
  FUNCTION CMISSNodesDestroyCNum(RegionUserNumber) BIND(C, NAME = "CMISSNodesDestroyNum")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the nodes to destroy, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesDestroyCNum !<Error Code.
    !Local variables

    CALL CMISSNodesDestroy(RegionUserNumber,CMISSNodesDestroyCNum)

    RETURN

  END FUNCTION CMISSNodesDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys the nodes in a region for nodes identified by an object for C.
  FUNCTION CMISSNodesDestroyCPtr(NodesPtr) BIND(C, NAME = "CMISSNodesDestroy")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: NodesPtr !<C pointer to the nodes to destroy.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesDestroyCPtr !<Error Code.
    !Local variable
    TYPE(CMISSNodesType), POINTER :: Nodes

    CMISSNodesDestroyCPtr=CMISSNoError
    IF(C_ASSOCIATED(NodesPtr)) THEN
      CALL C_F_POINTER(NodesPtr,Nodes)
      IF(ASSOCIATED(Nodes)) THEN
        CALL CMISSNodesDestroy(Nodes, CMISSNodesDestroyCPtr)
      ELSE
        CMISSNodesDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSNodesDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSNodesDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Returns the character label for a node in a set of nodes identified by user number for C. \todo should this be user number??
  FUNCTION CMISSNodesLabelGetCCNum(RegionUserNumber,NodeGlobalNumber,LabelSize, Label) BIND(C, NAME = "CMISSNodesLabelGetCNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the nodes to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NodeGlobalNumber !<The global number of the nodes to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Size of the label to get, for C.
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The label for the node, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesLabelGetCCNum !<Error Code.
    !Local variables
    CHARACTER(LEN = LabelSize-1) :: FLabel

    CALL CMISSNodesLabelGetC(RegionUserNumber,NodeGlobalNumber, FLabel, CMISSNodesLabelGetCCNum)
    CALL CMISSF2CString(FLabel, Label)

    RETURN

  END FUNCTION CMISSNodesLabelGetCCNum

    !
  !================================================================================================================================
  !

  !>Returns the character label for a node in a set of nodes identified by an object for C. \todo should this be user number??
  FUNCTION CMISSNodesLabelGetCCPtr(NodesPtr,NodeGlobalNumber,LabelSize, Label) BIND(C, NAME = "CMISSNodesLabelGetC")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: NodesPtr !<C pointer to the nodes to get the label for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NodeGlobalNumber !<The global number of the nodes to get the label for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The size of the label to get, for C.
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize) !<The label of the node, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesLabelGetCCPtr !<Error Code.
    !Local variables
    TYPE(CMISSNodesType), POINTER :: Nodes
    CHARACTER(LEN= LabelSize-1) :: FLabel

    CMISSNodesLabelGetCCPtr=CMISSNoError
    IF(C_ASSOCIATED(NodesPtr)) THEN
      CALL C_F_POINTER(NodesPtr,Nodes)
      IF(ASSOCIATED(Nodes)) THEN
        CALL CMISSNodesLabelGetC(Nodes, NodeGlobalNumber, FLabel, CMISSNodesLabelGetCCPtr)
        CALL CMISSF2CString(FLabel, Label)
      ELSE
        CMISSNodesLabelGetCCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSNodesLabelGetCCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSNodesLabelGetCCPtr


  !
  !================================================================================================================================
  !

  !>Sets/changes the character label for a node in a set of nodes identified by user number for C. \todo should this be user number??
  FUNCTION CMISSNodesLabelSetCCNum(RegionUserNumber,NodeGlobalNumber,LabelSize, Label) BIND(C, NAME = "CMISSNodesLabelSetCNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the nodes to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NodeGlobalNumber !<The global number of the nodes to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<Size of the label to set, for C.
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: Label(LabelSize) !<The label for the node to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesLabelSetCCNum !<Error Code.
    !Local variables
    CHARACTER(LEN = LabelSize-1) :: FLabel

    CALL CMISSC2FString(Label, FLabel)
    CALL CMISSNodesLabelSetC(RegionUserNumber,NodeGlobalNumber, FLabel, CMISSNodesLabelSetCCNum)


    RETURN

  END FUNCTION CMISSNodesLabelSetCCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the character label for a node in a set of nodes identified by an object for C. \todo should this be user number??
  FUNCTION CMISSNodesLabelSetCCPtr(NodesPtr,NodeGlobalNumber,LabelSize, Label) BIND(C, NAME = "CMISSNodesLabelSetC")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: NodesPtr !<C pointer to the nodes to set the label to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NodeGlobalNumber !<The global number of the nodes to set the label to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize !<The size of the label to set, for C.
    CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: Label(LabelSize) !<The label of the node to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesLabelSetCCPtr !<Error Code.
    !Local variables
    TYPE(CMISSNodesType), POINTER :: Nodes
    CHARACTER(LEN= LabelSize-1) :: FLabel

    CMISSNodesLabelSetCCPtr=CMISSNoError
    IF(C_ASSOCIATED(NodesPtr)) THEN
      CALL C_F_POINTER(NodesPtr,Nodes)
      IF(ASSOCIATED(Nodes)) THEN
        CALL CMISSC2FString(Label, FLabel)
        CALL CMISSNodesLabelSetC(Nodes, NodeGlobalNumber, FLabel, CMISSNodesLabelSetCCPtr)
      ELSE
        CMISSNodesLabelSetCCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSNodesLabelSetCCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSNodesLabelSetCCPtr

  !
  !================================================================================================================================
  !

  !>Returns the user number for a node in a set of nodes identified by user number for C.
  FUNCTION CMISSNodesUserNumberGetCNum(RegionUserNumber,NodeGlobalNumber,NodeUserNumber) BIND(C, NAME = &
    & "CMISSNodesUserNumberGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the nodes to get the node user number for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NodeGlobalNumber !<The global number of the nodes to get the node user number for, for C.
    INTEGER(C_INT), INTENT(OUT) :: NodeUserNumber !<The user number for the node, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesUserNumberGetCNum !<Error Code.
    !Local variables

    CALL CMISSNodesUserNumberGet(RegionUserNumber, NodeGlobalNumber,NodeUserNumber, CMISSNodesUserNumberGetCNum)

    RETURN

  END FUNCTION CMISSNodesUserNumberGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the user number for a node in a set of nodes identified by an object for C. \todo should this be user number??
  FUNCTION CMISSNodesUserNumberGetCPtr(NodesPtr,NodeGlobalNumber,NodeUserNumber) BIND(C, NAME = "CMISSNodesUserNumberGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: NodesPtr !<C pointer to the nodes to get the user number for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NodeGlobalNumber !<The global number of the nodes to set the node user number to, for C.
    INTEGER(C_INT), INTENT(OUT) :: NodeUserNumber !<The user number for the node to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesUserNumberGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSNodesType), POINTER :: Nodes

    CMISSNodesUserNumberGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(NodesPtr)) THEN
      CALL C_F_POINTER(NodesPtr, Nodes)
      IF(ASSOCIATED(Nodes)) THEN
        CALL CMISSNodesUserNumberGet(Nodes, NodeGlobalNumber, NodeUserNumber, CMISSNodesUserNumberGetCPtr)
      ELSE
        CMISSNodesUserNumberGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSNodesUserNumberGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSNodesUserNumberGetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the user number for a node in a set of nodes identified by user number, for C.
  FUNCTION CMISSNodesUserNumberSetCNum(RegionUserNumber,NodeGlobalNumber,NodeUserNumber) BIND(C, NAME = &
    & "CMISSNodesUserNumberSetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber !<The user number of the region containing the nodes to set the node user number to, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NodeGlobalNumber !<The global number of the nodes to set the node user number to, for C.
    INTEGER(C_INT), INTENT(OUT) :: NodeUserNumber !<The user number for the node to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesUserNumberSetCNum !<Error Code.
    !Local variables

    CALL CMISSNodesUserNumberSet(RegionUserNumber, NodeGlobalNumber,NodeUserNumber, CMISSNodesUserNumberSetCNum)

    RETURN

  END FUNCTION CMISSNodesUserNumberSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the user number for a node in a set of nodes identified by an object for C. \todo should this be user number??
  FUNCTION CMISSNodesUserNumberSetCPtr(NodesPtr,NodeGlobalNumber,NodeUserNumber) BIND(C, NAME = "CMISSNodesUserNumberSet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: NodesPtr !<C pointer to the nodes to set the user number to.
    INTEGER(C_INT), VALUE, INTENT(IN) :: NodeGlobalNumber !<The global number of the nodes to set the node user number to, for C.
    INTEGER(C_INT), INTENT(OUT) :: NodeUserNumber  !<The user number for the node to set, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSNodesUserNumberSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSNodesType), POINTER :: Nodes

    CMISSNodesUserNumberSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(NodesPtr)) THEN
      CALL C_F_POINTER(NodesPtr, Nodes)
      IF(ASSOCIATED(Nodes)) THEN
        CALL CMISSNodesUserNumberSet(Nodes, NodeGlobalNumber, NodeUserNumber, CMISSNodesUserNumberSetCPtr)
      ELSE
        CMISSNodesUserNumberSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSNodesUserNumberSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSNodesUserNumberSetCPtr

!!==================================================================================================================================
!!
!! PROBLEM_ROUTINES
!!
!!==================================================================================================================================

  !>Finishes the process of a problem identified by user number for C.
  FUNCTION CMISSProblemCreateFinishCNum(ProblemUserNumber) BIND(C, NAME="CMISSProblemCreateFinishNum")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to finish creating for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemCreateFinishCNum !<Error Code.
    !Local variables

    CALL CMISSProblemCreateFinish(ProblemUserNumber, CMISSProblemCreateFinishCNum)

    RETURN

  END FUNCTION CMISSProblemCreateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a problem identified by an object for C.
  FUNCTION CMISSProblemCreateFinishCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemCreateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to finish the creation of.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemCreateFinishCPtr !<Error Code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemCreateFinish(Problem, CMISSProblemCreateFinishCPtr)
      ELSE
        CMISSProblemCreateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts the process of a problem identified by user number for C.
  FUNCTION CMISSProblemCreateStartCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to start creating for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemCreateStartCNum !<Error code.
    !Local variables

    CALL CMISSProblemCreateStart(ProblemUserNumber, CMISSProblemCreateStartCNum)

    RETURN

  END FUNCTION CMISSProblemCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts the creation of a problem identified by an object for C.
  FUNCTION CMISSProblemCreateStartCPtr(ProblemUserNumber,ProblemPtr) BIND(C, NAME = "CMISSProblemCreateStart")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to start creating for C.
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the created problem.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemCreateStartCPtr !<Error code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemCreateStart(ProblemUserNumber, Problem, CMISSProblemCreateStartCPtr)
      ELSE
        CMISSProblemCreateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemCreateStartCPtr

  !
  !================================================================================================================================
  !

  !>Destroys a problem identified by an user number for C.
  FUNCTION CMISSProblemDestroyCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemDestroyNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to destroy for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemDestroyCNum !<Error code.
    !Local variables

    CALL CMISSProblemDestroy(ProblemUserNumber, CMISSProblemDestroyCNum)

    RETURN

  END FUNCTION CMISSProblemDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys a problem identified by an object for C.
  FUNCTION CMISSProblemDestroyCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemDestroy")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the destroyed problem.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemDestroyCPtr !<Error code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemDestroy(Problem, CMISSProblemDestroyCPtr)
      ELSE
        CMISSProblemDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a control loop for a problem identified by user number for C.
  FUNCTION CMISSProblemControlLoopCreateFinishCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemControlLoopCreateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to finish creating the control loop for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemControlLoopCreateFinishCNum !<Error code.
    !Local variables

    CALL CMISSProblemControlLoopCreateFinish(ProblemUserNumber, CMISSProblemControlLoopCreateFinishCNum)

    RETURN

  END FUNCTION CMISSProblemControlLoopCreateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a control loop on a problem identified by an object for C.
  FUNCTION CMISSProblemControlLoopCreateFinishCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemControlLoopCreateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to finish creating the control loop for.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemControlLoopCreateFinishCPtr !<Error code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemControlLoopCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemControlLoopCreateFinish(Problem, CMISSProblemControlLoopCreateFinishCPtr)
      ELSE
        CMISSProblemControlLoopCreateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemControlLoopCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemControlLoopCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a control loop for a problem identified by user number for C.
  FUNCTION CMISSProblemControlLoopCreateStartCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemControlLoopCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to start creating the control loop for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemControlLoopCreateStartCNum !<Error code.
    !Local variables

    CALL CMISSProblemControlLoopCreateStart(ProblemUserNumber, CMISSProblemControlLoopCreateStartCNum)

    RETURN

  END FUNCTION CMISSProblemControlLoopCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts the creation of a control loop on a problem identified by an object for C.
  FUNCTION CMISSProblemControlLoopCreateStartCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemControlLoopCreateStart")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to start creating the control loop for.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemControlLoopCreateStartCPtr !<Error code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemControlLoopCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemControlLoopCreateStart(Problem, CMISSProblemControlLoopCreateStartCPtr)
      ELSE
        CMISSProblemControlLoopCreateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemControlLoopCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemControlLoopCreateStartCPtr

  !
  !================================================================================================================================
  !

  !>Destroys the control loops for a problem identified by user number for C.
  FUNCTION CMISSProblemControlLoopDestroyCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemControlLoopDestroyNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to destroy the control loops for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemControlLoopDestroyCNum !<Error code.
    !Local variables

    CALL CMISSProblemControlLoopDestroy(ProblemUserNumber, CMISSProblemControlLoopDestroyCNum)

    RETURN

  END FUNCTION CMISSProblemControlLoopDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys the control loops on a problem identified by an object for C.
  FUNCTION CMISSProblemControlLoopDestroyCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemControlLoopDestroy")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to destroy the control loops for.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemControlLoopDestroyCPtr !<Error code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemControlLoopDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemControlLoopDestroy(Problem, CMISSProblemControlLoopDestroyCPtr)
      ELSE
        CMISSProblemControlLoopDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemControlLoopDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemControlLoopDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Returns a control loop from a problem identified by an user number for C.
  FUNCTION CMISSProblemControlLoopGetCNum(ProblemUserNumber,ControlLoopIdentifiersSize,ControlLoopIdentifiersPtr,ControlLoopPtr) &
    & BIND(C, NAME = "CMISSProblemControlLoopGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the control loop for, for C.
    INTEGER(C_INT), INTENT(IN) :: ControlLoopIdentifiersSize(1) !<Size of the control loop identifiers for C.
    TYPE(C_PTR), INTENT(IN) :: ControlLoopIdentifiersPtr !<C pointer to the location of the control loop identifier, for C.
    TYPE(C_PTR), INTENT(OUT) :: ControlLoopPtr !<C pointer to the specified problem control loop.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemControlLoopGetCNum !<Error Code.
    !Local variable
    INTEGER(C_INT), POINTER :: ControlLoopIdentifiers(:)
    TYPE(CMISSControlLoopType), POINTER :: ControlLoop

    CMISSProblemControlLoopGetCNum = CMISSNoError
    IF(C_ASSOCIATED(ControlLoopIdentifiersPtr)) THEN
      CALL C_F_POINTER(ControlLoopIdentifiersPtr, ControlLoopIdentifiers, ControlLoopIdentifiersSize)
      IF(ASSOCIATED(ControlLoopIdentifiers)) THEN
        CALL CMISSProblemControlLoopGet(ProblemUserNumber, ControlLoopIdentifiers, ControlLoop, &
          & CMISSProblemControlLoopGetCNum)
        IF(ASSOCIATED(ControlLoop)) THEN
          ControlLoopPtr = C_LOC(ControlLoop)
        ELSE
          CMISSProblemControlLoopGetCNum = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSProblemControlLoopGetCNum = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemControlLoopGetCNum = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemControlLoopGetCNum

  !
  !================================================================================================================================
  !

  !>Returns a control loop from a problem identified by an object for C.
  FUNCTION CMISSProblemControlLoopGetCPtr(ProblemPtr,ControlLoopIdentifiersSize,ControlLoopIdentifiersPtr,ControlLoopPtr) BIND(C, &
    & NAME = "CMISSProblemControlLoopGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<The user number of the problem to get the control loop for, for C.
    INTEGER(C_INT), INTENT(IN) :: ControlLoopIdentifiersSize(1) !<Size of the control loop identifiers for C.
    TYPE(C_PTR), INTENT(IN) :: ControlLoopIdentifiersPtr !<C pointer to the location of the control loop identifier, for C.
    TYPE(C_PTR), INTENT(OUT) :: ControlLoopPtr !<C pointer to the specified problem control loop.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemControlLoopGetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSProblemType), POINTER :: Problem
    INTEGER(C_INT), POINTER :: ControlLoopIdentifiers(:)
    TYPE(CMISSControlLoopType), POINTER :: ControlLoop

    CMISSProblemControlLoopGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(ControlLoopPtr)) THEN
      CMISSProblemControlLoopGetCPtr = CMISSPointerNotNULL
    ELSE
      IF(C_ASSOCIATED(ProblemPtr)) THEN
        CALL C_F_POINTER(ProblemPtr, Problem)
        IF(ASSOCIATED(Problem)) THEN
          IF(C_ASSOCIATED(ControlLoopIdentifiersPtr)) THEN
            CALL C_F_POINTER(ControlLoopIdentifiersPtr, ControlLoopIdentifiers, ControlLoopIdentifiersSize)
            IF(ASSOCIATED(ControlLoopIdentifiers)) THEN
              CALL CMISSProblemControlLoopGet(Problem, ControlLoopIdentifiers, ControlLoop, &
                & CMISSProblemControlLoopGetCPtr)
              IF(ASSOCIATED(ControlLoop)) THEN
                ControlLoopPtr = C_LOC(ControlLoop)
              ELSE
                CMISSProblemControlLoopGetCPtr = CMISSPointerIsNULL
              ENDIF
            ELSE
              CMISSProblemControlLoopGetCPtr = CMISSErrorConvertingPointer
            ENDIF
          ELSE
            CMISSProblemControlLoopGetCPtr = CMISSPointerIsNULL
          ENDIF
        ELSE
          CMISSProblemControlLoopGetCPtr = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSProblemControlLoopGetCPtr = CMISSPointerIsNULL
      ENDIF
    ENDIF

    RETURN

  END FUNCTION CMISSProblemControlLoopGetCPtr

  !
  !================================================================================================================================
  !

  !>Solves a problem identified by an user number for C.
  FUNCTION CMISSProblemSolveCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemSolveNum")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to solve, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolveCNum !<Error Code.
    !Local variables

    CALL CMISSProblemSolve(ProblemUserNumber, CMISSProblemSolveCNum)

    RETURN

  END FUNCTION CMISSProblemSolveCNum

  !
  !================================================================================================================================
  !

  !>Solves a problem identified by an object for C.
  FUNCTION CMISSProblemSolveCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemSolve")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to solve, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolveCPtr !<Error Code.
    !Local variable
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemSolveCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemSolve(Problem, CMISSProblemSolveCPtr)
      ELSE
        CMISSProblemSolveCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSolveCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolveCPtr

  !
  !================================================================================================================================
  !

  !>Returns a solver from a problem identified by an user number for C.
  FUNCTION CMISSProblemSolverGetCNum(ProblemUserNumber,ControlLoopIdentifiersSize, ControlLoopIdentifiersPtr,SolverIndex,SolverPtr)&
    & BIND(C,NAME = "CMISSProblemSolverGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to get the solver for, for C.
    INTEGER(C_INT), INTENT(IN) :: ControlLoopIdentifiersSize(1) !<Size of the control loop identifiers for C.
    TYPE(C_PTR), INTENT(IN) :: ControlLoopIdentifiersPtr !<C pointer to the location of the control loop identifiers, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: SolverIndex !<The solver index to get the solver for, for C.
    TYPE(C_PTR), INTENT(OUT) :: SolverPtr !<C pointer to the specified solver.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolverGetCNum !<Error Code.
    !Local variable
    INTEGER(C_INT), POINTER :: ControlLoopIdentifiers(:)
    TYPE(CMISSSolverType), POINTER :: Solver

    CMISSProblemSolverGetCNum = CMISSNoError
    IF(C_ASSOCIATED(SolverPtr)) THEN
      CMISSProblemSolverGetCNum = CMISSPointerNotNULL
      IF(C_ASSOCIATED(ControlLoopIdentifiersPtr)) THEN
        CALL C_F_POINTER(ControlLoopIdentifiersPtr, ControlLoopIdentifiers, ControlLoopIdentifiersSize)
        IF(ASSOCIATED(ControlLoopIdentifiers)) THEN
          CALL CMISSProblemSolverGet(ProblemUserNumber,ControlLoopIdentifiers,SolverIndex, &
            & Solver, CMISSProblemSolverGetCNum)
          IF(ASSOCIATED(Solver)) THEN
            SolverPtr = C_LOC(Solver)
          ELSE
            CMISSProblemSolverGetCNum = CMISSPointerIsNULL
          ENDIF
        ELSE
          CMISSProblemSolverGetCNum = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSProblemSolverGetCNum = CMISSPointerIsNULL
      ENDIF
    ELSE
      CMISSProblemSolverGetCNum = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolverGetCNum

  !
  !================================================================================================================================
  !

  !>Returns a solver from a problem identified by an object for C.
  FUNCTION CMISSProblemSolverGetCPtr(ProblemPtr,ControlLoopIdentifiersSize,ControlLoopIdentifiersPtr,SolverIndex,SolverPtr) BIND(C,&
    & NAME = "CMISSProblemSolverGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to get the solver for.
    INTEGER(C_INT), INTENT(IN) :: ControlLoopIdentifiersSize(1) !<Size of the control loop identifiers for C.
    TYPE(C_PTR), INTENT(IN) :: ControlLoopIdentifiersPtr !<C pointer to the location of the control loop identifiers, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: SolverIndex !<The solver index to get the solver for, for C.
    TYPE(C_PTR), INTENT(OUT) :: SolverPtr !<C pointer to the specified solver.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolverGetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSProblemType), POINTER :: Problem
    INTEGER(C_INT), POINTER :: ControlLoopIdentifiers(:)
    TYPE(CMISSSolverType), POINTER :: Solver

    CMISSProblemSolverGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        IF(C_ASSOCIATED(SolverPtr)) THEN
          CMISSProblemSolverGetCPtr = CMISSPointerNotNULL
          IF(C_ASSOCIATED(ControlLoopIdentifiersPtr)) THEN
            CALL C_F_POINTER(ControlLoopIdentifiersPtr, ControlLoopIdentifiers, ControlLoopIdentifiersSize)
            IF(ASSOCIATED(ControlLoopIdentifiers)) THEN
              CALL CMISSProblemSolverGet(Problem,ControlLoopIdentifiers,SolverIndex, &
                & Solver, CMISSProblemSolverGetCPtr)
              IF(ASSOCIATED(Solver)) THEN
                SolverPtr = C_LOC(Solver)
              ELSE
                CMISSProblemSolverGetCPtr = CMISSPointerIsNULL
              ENDIF
            ELSE
              CMISSProblemSolverGetCPtr = CMISSErrorConvertingPointer
            ENDIF
          ELSE
            CMISSProblemSolverGetCPtr = CMISSPointerIsNULL
          ENDIF
        ELSE
          CMISSProblemSolverGetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSProblemSolverGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSolverGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolverGetCPtr

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating solver equations for a problem identified by user number for C.
  FUNCTION CMISSProblemSolverEquationsCreateFinishCNum(ProblemUserNumber) BIND(C, NAME = &
    & "CMISSProblemSolverEquationsCreateFinishNum")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber  !<The user number of the problem to finish the creation of solver equations for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolverEquationsCreateFinishCNum !<Error Code.
    !Local variable

    CALL CMISSProblemSolverEquationsCreateFinish(ProblemUserNumber, CMISSProblemSolverEquationsCreateFinishCNum)

    RETURN

  END FUNCTION CMISSProblemSolverEquationsCreateFinishCNum


  !
  !================================================================================================================================
  !

  !>Finishes the creation of solver equations for problem identified by an object for C.
  FUNCTION CMISSProblemSolverEquationsCreateFinishCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemSolverEquationsCreateFinish")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to finish the creation of solver equations for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolverEquationsCreateFinishCPtr !<Error Code.
    !Local variable
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemSolverEquationsCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemSolverEquationsCreateFinish(Problem, CMISSProblemSolverEquationsCreateFinishCPtr)
      ELSE
        CMISSProblemSolverEquationsCreateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSolverEquationsCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolverEquationsCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts the process of creating solver equations for a problem identified by user number for C.
  FUNCTION CMISSProblemSolverEquationsCreateStartCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemSolverEquationsCreateStartNum")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to start the creation of solver equations for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolverEquationsCreateStartCNum !<Error Code.
    !Local variable

    CALL CMISSProblemSolverEquationsCreateStart(ProblemUserNumber, CMISSProblemSolverEquationsCreateStartCNum)

    RETURN

  END FUNCTION CMISSProblemSolverEquationsCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts the creation of solver equations for problem identified by an object for C.
  FUNCTION CMISSProblemSolverEquationsCreateStartCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemSolverEquationsCreateStart")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolverEquationsCreateStartCPtr
    !Local variable
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemSolverEquationsCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr,Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemSolverEquationsCreateStart(Problem, CMISSProblemSolverEquationsCreateStartCPtr)
      ELSE
        CMISSProblemSolverEquationsCreateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSolverEquationsCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolverEquationsCreateStartCPtr

  !
  !================================================================================================================================
  !

  !>Destroys the solver equations for a problem identified by an user number for C.
  FUNCTION CMISSProblemSolverEquationsDestroyCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemSolverEquationsDestroyNum")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to destroy the solver equations for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolverEquationsDestroyCNum
    !Local variable

    CALL CMISSProblemSolverEquationsDestroy(ProblemUserNumber, CMISSProblemSolverEquationsDestroyCNum)

    RETURN

  END FUNCTION CMISSProblemSolverEquationsDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys the solver equations for problem identified by an object for C.
  FUNCTION CMISSProblemSolverEquationsDestroyCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemSolverEquationsDestroy")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolverEquationsDestroyCPtr
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemSolverEquationsDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemSolverEquationsDestroy(Problem, CMISSProblemSolverEquationsDestroyCPtr)
      ELSE
        CMISSProblemSolverEquationsDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSolverEquationsDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolverEquationsDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Returns the solver equations from a problem identified by an user number.
  FUNCTION CMISSProblemSolverEquationsGetCNum(ProblemUserNumber,ControlLoopIdentifiersSize,ControlLoopIdentifiersPtr,SolverIndex, &
    & SolverEquationsPtr) BIND(C, NAME = "CMISSProblemSolverEquationsGetNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<C pointer to the problem to get the solver for.
    INTEGER(C_INT), INTENT(IN) :: ControlLoopIdentifiersSize(1) !<Size of the control loop identifiers for C.
    TYPE(C_PTR), INTENT(IN) :: ControlLoopIdentifiersPtr !<C pointer to the location of the control loop identifiers, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: SolverIndex !<The solver index to get the solver for, for C.
    TYPE(C_PTR), INTENT(OUT) :: SolverEquationsPtr !<C pointer to the specified solver equations.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolverEquationsGetCNum !<Error Code.
    !Local variable
    INTEGER(C_INT), POINTER :: ControlLoopIdentifiers(:)
    TYPE(CMISSSolverEquationsType), POINTER :: SolverEquations

    CMISSProblemSolverEquationsGetCNum = CMISSNoError
    IF(C_ASSOCIATED(SolverEquationsPtr)) THEN
      CMISSProblemSolverEquationsGetCNum = CMISSPointerNotNULL
      IF(C_ASSOCIATED(ControlLoopIdentifiersPtr)) THEN
        CALL C_F_POINTER(ControlLoopIdentifiersPtr, ControlLoopIdentifiers, ControlLoopIdentifiersSize)
        IF(ASSOCIATED(ControlLoopIdentifiers)) THEN
          CALL CMISSProblemSolverEquationsGet(ProblemUserNumber,ControlLoopIdentifiers,SolverIndex, &
            & SolverEquations, CMISSProblemSolverEquationsGetCNum)
          IF(ASSOCIATED(SolverEquations)) THEN
            SolverEquationsPtr = C_LOC(SolverEquations)
          ELSE
            CMISSProblemSolverEquationsGetCNum = CMISSPointerIsNULL
          ENDIF
        ELSE
          CMISSProblemSolverEquationsGetCNum = CMISSErrorConvertingPointer
        ENDIF
      ELSE
        CMISSProblemSolverEquationsGetCNum = CMISSPointerIsNULL
      ENDIF
    ELSE
      CMISSProblemSolverEquationsGetCNum = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolverEquationsGetCNum

  !
  !================================================================================================================================
  !

  !>Returns the solver equations from a problem identified by an object.
  FUNCTION CMISSProblemSolverEquationsGetCPtr(ProblemPtr,ControlLoopIdentifiersSize,ControlLoopIdentifiersPtr,SolverIndex, &
    & SolverEquationsPtr)BIND(C, NAME = "CMISSProblemSolverEquationsGet")

    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to get the solver for.
    INTEGER(C_INT), INTENT(IN) :: ControlLoopIdentifiersSize(1) !<Size of the control loop identifiers for C.
    TYPE(C_PTR), INTENT(IN) :: ControlLoopIdentifiersPtr !<C pointer to the location of the control loop identifiers, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: SolverIndex !<The solver index to get the solver for, for C.
    TYPE(C_PTR), INTENT(OUT) :: SolverEquationsPtr !<C pointer to the specified solver equations.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolverEquationsGetCPtr !<Error Code.
    !Local variable
    TYPE(CMISSProblemType), POINTER :: Problem
    INTEGER(C_INT), POINTER :: ControlLoopIdentifiers(:)
    TYPE(CMISSSolverEquationsType), POINTER :: SolverEquations

    CMISSProblemSolverEquationsGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        IF(C_ASSOCIATED(SolverEquationsPtr)) THEN
          CMISSProblemSolverEquationsGetCPtr = CMISSPointerNotNULL
          IF(C_ASSOCIATED(ControlLoopIdentifiersPtr)) THEN
            CALL C_F_POINTER(ControlLoopIdentifiersPtr, ControlLoopIdentifiers, ControlLoopIdentifiersSize)
            IF(ASSOCIATED(ControlLoopIdentifiers)) THEN
              CALL CMISSProblemSolverEquationsGet(Problem,ControlLoopIdentifiers,SolverIndex, &
                & SolverEquations, CMISSProblemSolverEquationsGetCPtr)
              IF(ASSOCIATED(SolverEquations)) THEN
                SolverEquationsPtr = C_LOC(SolverEquations)
              ELSE
                CMISSProblemSolverEquationsGetCPtr = CMISSPointerIsNULL
              ENDIF
            ELSE
              CMISSProblemSolverEquationsGetCPtr = CMISSErrorConvertingPointer
            ENDIF
          ELSE
            CMISSProblemSolverEquationsGetCPtr = CMISSPointerIsNULL
          ENDIF
        ELSE
          CMISSProblemSolverEquationsGetCPtr = CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSProblemSolverEquationsGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSolverEquationsGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolverEquationsGetCPtr

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating solvers for a problem identified by user numbe for Cr.
  FUNCTION CMISSProblemSolversCreateFinishCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemSolversCreateFinishNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to finish the creation of solvers for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolversCreateFinishCNum
    !Local variable

    CALL CMISSProblemSolversCreateFinish(ProblemUserNumber, CMISSProblemSolversCreateFinishCNum)

    RETURN

  END FUNCTION CMISSProblemSolversCreateFinishCNum

  !
  !================================================================================================================================
  !

  !>Finishes the creation of solvers for problem identified by an object for C.
  FUNCTION CMISSProblemSolversCreateFinishCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemSolversCreateFinish")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to finish the creation of solvers for.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolversCreateFinishCPtr !<Error Code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemSolversCreateFinishCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemSolversCreateFinish(Problem, CMISSProblemSolversCreateFinishCPtr)
      ELSE
        CMISSProblemSolversCreateFinishCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSolversCreateFinishCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolversCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Starts the process of creating solvers for a problem identified by user number for C.
  FUNCTION CMISSProblemSolversCreateStartCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemSolversCreateStartNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to start the creation of solvers for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolversCreateStartCNum !<Error Code.
    !Local variables

    CALL CMISSProblemSolversCreateStart(ProblemUserNumber, CMISSProblemSolversCreateStartCNum)

    RETURN

  END FUNCTION CMISSProblemSolversCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Starts the creation of solvers for problem identified by an object for C.
  FUNCTION CMISSProblemSolversCreateStartCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemSolversCreateStart")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to start the creation of solvers for.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolversCreateStartCPtr !<Error Code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemSolversCreateStartCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemSolversCreateStart(Problem, CMISSProblemSolversCreateStartCPtr)
      ELSE
        CMISSProblemSolversCreateStartCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSolversCreateStartCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolversCreateStartCPtr

  !
  !================================================================================================================================
  !

  !>Destroys the solvers for a problem identified by an user number for C.
  FUNCTION CMISSProblemSolversDestroyCNum(ProblemUserNumber) BIND(C, NAME = "CMISSProblemSolversDestroyNum")

    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number of the problem to destroy the solvers for, for C.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolversDestroyCNum !<Error Code.
    !Local variables

    CALL CMISSProblemSolversDestroy(ProblemUserNumber, CMISSProblemSolversDestroyCNum)

    RETURN

  END FUNCTION CMISSProblemSolversDestroyCNum

  !
  !================================================================================================================================
  !

  !>Destroys the solvers for problem identified by an object.
  FUNCTION CMISSProblemSolversDestroyCPtr(ProblemPtr) BIND(C, NAME = "CMISSProblemSolversDestroy")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to destroy the solvers for.
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSolversDestroyCPtr !<Error Code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemSolversDestroyCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemSolversDestroy(Problem, CMISSProblemSolversDestroyCPtr)
      ELSE
        CMISSProblemSolversDestroyCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSolversDestroyCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSolversDestroyCPtr

  !
  !================================================================================================================================
  !

  !>Returns the specification i.e., problem class, type and subtype for a problem identified by an user number for C.
  FUNCTION CMISSProblemSpecificationGetCNum(ProblemUserNumber,ProblemClass,ProblemType,ProblemSubtype) BIND(C, NAME = &
    & "CMISSProblemSpecificationGetNum")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number for the problem to get the specification for, for C.
    INTEGER(C_INT), INTENT(OUT) :: ProblemClass !<The problem class for C. \see OPENCMISS_ProblemClasses
    INTEGER(C_INT), INTENT(OUT) :: ProblemType !<The problem type for C. \see OPENCMISS_ProblemTypes
    INTEGER(C_INT), INTENT(OUT) :: ProblemSubtype !<The problem subtype for C. \see OPENCMISS_ProblemSubTypes
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSpecificationGetCNum !<Error Code.
    !Local variables

    CALL CMISSProblemSpecificationGet(ProblemUserNumber,ProblemClass,ProblemType,ProblemSubtype,CMISSProblemSpecificationGetCNum)

    RETURN

  END FUNCTION CMISSProblemSpecificationGetCNum

    !
  !================================================================================================================================
  !

  !>Returns the specification i.e., problem class, type and subtype for a problem identified by an object for C..
  FUNCTION CMISSProblemSpecificationGetCPtr(ProblemPtr,ProblemClass,ProblemType,ProblemSubtype) BIND(C, NAME = &
    & "CMISSProblemSpecificationGet")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to get the specification for.
    INTEGER(C_INT), INTENT(OUT) :: ProblemClass !<The problem class for C. \see OPENCMISS_ProblemClasses
    INTEGER(C_INT), INTENT(OUT) :: ProblemType !<The problem type for C. \see OPENCMISS_ProblemTypes
    INTEGER(C_INT), INTENT(OUT) :: ProblemSubtype !<The problem subtype for C. \see OPENCMISS_ProblemSubTypes
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSpecificationGetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemSpecificationGetCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemSpecificationGet(Problem, ProblemClass,ProblemType,ProblemSubtype,CMISSProblemSpecificationGetCPtr)
      ELSE
        CMISSProblemSpecificationGetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSpecificationGetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSpecificationGetCPtr

    !
  !================================================================================================================================
  !

  !>Sets/changes the specification i.e., problem class, type and subtype for a problem identified by an user number for C.
  FUNCTION CMISSProblemSpecificationSetCNum(ProblemUserNumber,ProblemClass,ProblemType,ProblemSubtype) BIND(C, NAME = &
    & "CMISSProblemSpecificationSetNum")

    !Argument variable
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemUserNumber !<The user number for the problem to set the specification for, for C.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemClass !<The problem class to set for C. \see OPENCMISS_ProblemClasses
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemType !<The problem type to set for C. \see OPENCMISS_ProblemTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemSubtype !<The problem subtype to set for C. \see OPENCMISS_ProblemSubTypes
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSpecificationSetCNum !<Error Code.
    !Local variables

    CALL CMISSProblemSpecificationSet(ProblemUserNumber,ProblemClass,ProblemType,ProblemSubtype,CMISSProblemSpecificationSetCNum)

    RETURN

  END FUNCTION CMISSProblemSpecificationSetCNum

  !
  !================================================================================================================================
  !

  !>Sets/changes the specification i.e., problem class, type and subtype for a problem identified by an object for C.
  FUNCTION CMISSProblemSpecificationSetCPtr(ProblemPtr,ProblemClass,ProblemType,ProblemSubtype) BIND(C, NAME = &
    & "CMISSProblemSpecificationSet")

    !Argument variable
    TYPE(C_PTR), VALUE, INTENT(IN) :: ProblemPtr !<C pointer to the problem to set the specification for.
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemClass !<The problem class to set for C. \see OPENCMISS_ProblemClasses
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemType !<The problem type to set for C. \see OPENCMISS_ProblemTypes
    INTEGER(C_INT), VALUE, INTENT(IN) :: ProblemSubtype !<The problem subtype to set for C. \see OPENCMISS_ProblemSubTypes
    !Function variable
    INTEGER(C_INT) :: CMISSProblemSpecificationSetCPtr !<Error Code.
    !Local variables
    TYPE(CMISSProblemType), POINTER :: Problem

    CMISSProblemSpecificationSetCPtr = CMISSNoError
    IF(C_ASSOCIATED(ProblemPtr)) THEN
      CALL C_F_POINTER(ProblemPtr, Problem)
      IF(ASSOCIATED(Problem)) THEN
        CALL CMISSProblemSpecificationSet(Problem, ProblemClass,ProblemType,ProblemSubtype,CMISSProblemSpecificationSetCPtr)
      ELSE
        CMISSProblemSpecificationSetCPtr = CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSProblemSpecificationSetCPtr = CMISSPointerIsNULL
    ENDIF

    RETURN

  END FUNCTION CMISSProblemSpecificationSetCPtr

!!==================================================================================================================================
!!
!! REGION_ROUTINES
!!
!!==================================================================================================================================

  !>Finishes the process of creating a region identified by a pointer for C.
  FUNCTION CMISSRegionCreateFinishCPtr(RegionPtr) BIND(C,NAME="CMISSRegionCreateFinish")

    !Argument variablesCMISSAnalyticAnalysisOutputNumber
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr
    !Function variable
    INTEGER(C_INT) :: CMISSRegionCreateFinishCPtr !<Error Code.
    !Local variables
    TYPE(CMISSRegionType), POINTER :: Region

    CMISSRegionCreateFinishCPtr=CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr,Region)
      IF(ASSOCIATED(Region)) THEN        
        CALL CMISSRegionCreateFinish(Region,CMISSRegionCreateFinishCPtr)
      ELSE
        CMISSRegionCreateFinishCPtr=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionCreateFinishCPtr=CMISSPointerIsNULL
    ENDIF

    RETURN
    
  END FUNCTION CMISSRegionCreateFinishCPtr

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a region identified by an user number for C.
  FUNCTION CMISSRegionCreateFinishCNum(RegionUserNumber) BIND(C,NAME="CMISSRegionCreateFinishNum")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    !Function variable
    INTEGER(C_INT) :: CMISSRegionCreateFinishCNum
    !Local variables

    CALL CMISSRegionCreateFinish(RegionUserNumber,CMISSRegionCreateFinishCNum)

    RETURN
    
  END FUNCTION CMISSRegionCreateFinishCNum

  !
  !================================================================================================================================
  !
  
  !>Starts the process of creating a region identified by a pointer for C.
  FUNCTION CMISSRegionCreateStartCPtr(RegionUserNumber,ParentRegionPtr,RegionPtr) BIND(C,NAME="CMISSRegionCreateStart")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    TYPE(C_PTR), VALUE, INTENT(IN) :: ParentRegionPtr
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr
    !Function variable
    INTEGER(C_INT) :: CMISSRegionCreateStartCPtr
    !Local variables
    TYPE(CMISSRegionType), POINTER :: Region,ParentRegion

    CMISSRegionCreateStartCPtr=CMISSNoError
    IF(C_ASSOCIATED(ParentRegionPtr)) THEN
      CALL C_F_POINTER(ParentRegionPtr,ParentRegion)
      IF(ASSOCIATED(ParentRegion)) THEN        
        IF(C_ASSOCIATED(RegionPtr)) THEN
          CALL C_F_POINTER(RegionPtr,Region)
          IF(ASSOCIATED(Region)) THEN        
            CALL CMISSRegionCreateStart(RegionUserNumber,ParentRegion,Region,CMISSRegionCreateStartCPtr)
          ELSE
            CMISSRegionCreateStartCPtr=CMISSErrorConvertingPointer
          ENDIF
        ELSE
          CMISSRegionCreateStartCPtr=CMISSPointerIsNULL
        ENDIF
      ELSE
        CMISSRegionCreateStartCPtr=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionCreateStartCPtr=CMISSPointerIsNULL
    ENDIF
    
    RETURN
    
  END FUNCTION CMISSRegionCreateStartCPtr

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a region identified by an user number for C.
  FUNCTION CMISSRegionCreateStartCNum(RegionUserNumber,ParentRegionUserNumber) BIND(C,NAME="CMISSRegionCreateStartNum")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: ParentRegionUserNumber
    !Function variable
    INTEGER(C_INT) :: CMISSRegionCreateStartCNum
    !Local variables

    CALL CMISSRegionCreateStart(RegionUserNumber,ParentRegionUserNumber,CMISSRegionCreateStartCNum)

    RETURN
    
  END FUNCTION CMISSRegionCreateStartCNum

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a region identified by a pointer for C.
  FUNCTION CMISSRegionLabelGetCPtr(RegionPtr,LabelSize,Label) BIND(C,NAME="CMISSRegionLabelGet")
  
    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize)
    !Function variable
    INTEGER(C_INT) :: CMISSRegionLabelGetCPtr
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel
    TYPE(CMISSRegionType), POINTER :: Region    

    CMISSRegionLabelGetCPtr=CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr,Region)
      IF(ASSOCIATED(Region)) THEN        
        CALL CMISSRegionLabelGet(Region,FLabel,CMISSRegionLabelGetCPtr)
        CALL CMISSF2CString(Flabel,Label)
      ELSE
        CMISSRegionLabelGetCPtr=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionLabelGetCPtr=CMISSPointerIsNULL
    ENDIF

    RETURN
    
  END FUNCTION CMISSRegionLabelGetCPtr

  !
  !================================================================================================================================
  !

  !>Returns the character string label for a region identified by an user number for C.
  FUNCTION CMISSRegionLabelGetCNum(RegionUserNumber,LabelSize,Label) BIND(C,NAME="CMISSRegionLabelGetNum")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Label(LabelSize)
    !Function variable
    INTEGER(C_INT) :: CMISSRegionLabelGetCNum
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel
 
    CALL CMISSRegionLabelGet(RegionUserNumber,FLabel,CMISSRegionLabelGetCNum)
    CALL CMISSF2CString(Flabel,Label)
 
    RETURN
    
  END FUNCTION CMISSRegionLabelGetCNum

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the label for a region identified by a pointer for C.
  FUNCTION CMISSRegionLabelSetCPtr(RegionPtr,LabelSize,Label) BIND(C,NAME="CMISSRegionLabelSet")
  
    !Argument variables
    TYPE(C_PTR), VALUE, INTENT(IN) :: RegionPtr
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Label(LabelSize)
    !Function variable
    INTEGER(C_INT) :: CMISSRegionLabelSetCPtr
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel
    TYPE(CMISSRegionType), POINTER :: Region    

    CMISSRegionLabelSetCPtr=CMISSNoError
    IF(C_ASSOCIATED(RegionPtr)) THEN
      CALL C_F_POINTER(RegionPtr,Region)
      IF(ASSOCIATED(Region)) THEN        
        CALL CMISSC2FString(Label,Flabel)
        CALL CMISSRegionLabelSet(Region,FLabel,CMISSRegionLabelSetCPtr)
      ELSE
        CMISSRegionLabelSetCPtr=CMISSErrorConvertingPointer
      ENDIF
    ELSE
      CMISSRegionLabelSetCPtr=CMISSPointerIsNULL
    ENDIF

    RETURN
    
  END FUNCTION CMISSRegionLabelSetCPtr

  !
  !================================================================================================================================
  !

  !>Sets/changes the label for a region identified by an user number for C.
  FUNCTION CMISSRegionLabelSetCNum(RegionUserNumber,LabelSize,Label) BIND(C,NAME="CMISSRegionLabelSetNum")
  
    !Argument variables
    INTEGER(C_INT), VALUE, INTENT(IN) :: RegionUserNumber
    INTEGER(C_INT), VALUE, INTENT(IN) :: LabelSize
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Label(LabelSize)
    !Function variable
    INTEGER(C_INT) :: CMISSRegionLabelSetCNum
    !Local variables
    CHARACTER(LEN=LabelSize-1) :: FLabel
 
    CALL CMISSC2FString(Label,Flabel)
    CALL CMISSRegionLabelSet(RegionUserNumber,FLabel,CMISSRegionLabelSetCNum)
    
    RETURN
    
  END FUNCTION CMISSRegionLabelSetCNum

  !
  !================================================================================================================================
  !

END MODULE OPENCMISS_C
