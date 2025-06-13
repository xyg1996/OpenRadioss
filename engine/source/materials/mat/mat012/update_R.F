Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Altair Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
      !||====================================================================
      !||    m14gtf   ../engine/source/materials/mat/mat014/m14gtf.F
      !||--- called by ------------------------------------------------------
      !||    m12law   ../engine/source/materials/mat/mat012/m12law.F
      !||    m14law   ../engine/source/materials/mat/mat014/m14law.F
      !||====================================================================
      SUBROUTINE UPDATE_R( 
     1   AX, AY, AZ, BX, BY,  
     2   BZ, CX, CY, CZ, NEL,
     3   R_OLD, ORTH_DIR)
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   G l o b a l   P a r a m e t e r s
C-----------------------------------------------
#include      "mvsiz_p.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com01_c.inc"
#include      "param_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER NEL
C     REAL
      my_real
     .   R_OLD(NEL,9), ORTH_DIR(NEL,9),
     .   AX(*),AY(*),AZ(*),BX(*),BY(*),BZ(*),CX(*),CY(*),CZ(*)
C     REAL
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER I, J, K
      my_real
     .   TEMPQ(3,3), R(3,3), Q1(3), Q2(3), Q3(3)
c     .   AX_temp,AY_temp,AZ_temp,BX_temp,BY_temp,BZ_temp,CX_temp,CY_temp,CZ_temp,
c     .   RR, TX(MVSIZ), TY(MVSIZ),TZ(MVSIZ)
C-----------------------------------------------
      DO I = 1, NEL
C--- 还原 R_OLD(I,:) 成为 3x3 矩阵 R(J,K)
        DO J = 1, 3
          DO K = 1, 3
            R(J,K) = R_OLD(I, (J-1)*3 + K)
          END DO
        END DO

c        DO J = 1, 9
c          PRINT *, 'ORTH_DIR is', ORTH_DIR(I, J)
c        ENDDO 

c        TX(I)=RY(I)*SZ(I)-RZ(I)*SY(I)
c        TY(I)=RZ(I)*SX(I)-RX(I)*SZ(I)
c        TZ(I)=RX(I)*SY(I)-RY(I)*SX(I)
c        RR=ONE/SQRT(TX(I)**2+TY(I)**2+TZ(I)**2)
c        TX(I)=TX(I)*RR
c        TY(I)=TY(I)*RR
c        TZ(I)=TZ(I)*RR

c        AX_temp = A(I,1)*RX(I)+A(I,2)*SX(I)+A(I,3)*TX(I)
c        AY_temp = A(I,1)*RY(I)+A(I,2)*SY(I)+A(I,3)*TY(I)
c       AZ_temp = A(I,1)*RZ(I)+A(I,2)*SZ(I)+A(I,3)*TZ(I)
C
c        BX_temp = A(I,4)*RX(I)+A(I,5)*SX(I)+A(I,6)*TX(I)
c        BY_temp = A(I,4)*RY(I)+A(I,5)*SY(I)+A(I,6)*TY(I)
c        BZ_temp = A(I,4)*RZ(I)+A(I,5)*SZ(I)+A(I,6)*TZ(I)

c        CX_temp = AY(I)*BZ(I)-AZ(I)*BY(I)
c        CY_temp = AZ(I)*BX(I)-AX(I)*BZ(I)
c        CZ_temp = AX(I)*BY(I)-AY(I)*BX(I)

C--- 构造原 Q 矩阵列为 [AX AY AZ; BX BY BZ; CX CY CZ]
        Q1(1) = ORTH_DIR(I, 1)
        Q1(2) = ORTH_DIR(I, 4)
        Q1(3) = ORTH_DIR(I, 7)

        Q2(1) = ORTH_DIR(I, 2)
        Q2(2) = ORTH_DIR(I, 5)
        Q2(3) = ORTH_DIR(I, 8)

        Q3(1) = ORTH_DIR(I, 3)
        Q3(2) = ORTH_DIR(I, 6)
        Q3(3) = ORTH_DIR(I, 9)
c        PRINT *, 'Q1 is', Q1
c        PRINT *, 'Q2 is', Q2
c        PRINT *, 'Q3 is', Q3
C--- 计算新 Q 矩阵列：TEMPQ(:,1) = R * Q1, etc.
        DO J = 1, 3
          TEMPQ(J,1) = R(J,1)*Q1(1) + R(J,2)*Q1(2) + R(J,3)*Q1(3)
          TEMPQ(J,2) = R(J,1)*Q2(1) + R(J,2)*Q2(2) + R(J,3)*Q2(3)
          TEMPQ(J,3) = R(J,1)*Q3(1) + R(J,2)*Q3(2) + R(J,3)*Q3(3)
        END DO

C--- 更新新列向量 AX...CZ
        AX(I) = TEMPQ(1,1)
        AY(I) = TEMPQ(2,1)
        AZ(I) = TEMPQ(3,1)

        BX(I) = TEMPQ(1,2)
        BY(I) = TEMPQ(2,2)
        BZ(I) = TEMPQ(3,2)

        CX(I) = TEMPQ(1,3)
        CY(I) = TEMPQ(2,3)
        CZ(I) = TEMPQ(3,3)
      END DO
C
      RETURN
      END
