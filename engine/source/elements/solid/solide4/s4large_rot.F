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
      !||    s4defot3   ../engine/source/elements/solid/solide4/s4defot3.F
      !||--- called by ------------------------------------------------------
      !||    s4forc3    ../engine/source/elements/solid/solide4/s4forc3.F
      !||====================================================================
      SUBROUTINE S4LARGE_ROT(
     1   PX1,     PX2,     PX3,     PX4,
     2   PY1,     PY2,     PY3,     PY4,
     3   PZ1,     PZ2,     PZ3,     PZ4,
     4   VX1,     VX2,     VX3,     VX4,
     5   VY1,     VY2,     VY3,     VY4,
     6   VZ1,     VZ2,     VZ3,     VZ4,
     7   DXX,     DXY,     DXZ,     DYX,
     8   DYY,     DYZ,     DZX,     DZY,
     9   DZZ,     NEL,     DT1,     R_OLD)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MESSAGE_MOD
      USE DT_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
#include      "mvsiz_p.inc"
#include      "comlock.inc"
#include      "com01_c.inc"
#include      "com04_c.inc"
#include      "scr19_c.inc"
#include      "task_c.inc"
#include      "userlib.inc"
#include      "timeri_c.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "units_c.inc"
#include      "param_c.inc"
C-----------------------------------------------
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN) :: NEL
      my_real
     .   VX1(*), VX2(*), VX3(*), VX4(*),
     .   VY1(*), VY2(*), VY3(*), VY4(*),
     .   VZ1(*), VZ2(*), VZ3(*), VZ4(*),
     .   PX1(*), PX2(*), PX3(*), PX4(*),
     .   PY1(*), PY2(*), PY3(*), PY4(*),
     .   PZ1(*), PZ2(*), PZ3(*), PZ4(*), 
     .   DXX(*), DXY(*), DXZ(*),
     .   DYX(*), DYY(*), DYZ(*),
     .   DZX(*), DZY(*), DZZ(*),
     .   TIMESTEP, DT1, R_OLD(NEL,9)
     
C      TYPE(DT_), INTENT(INOUT) :: DT
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  I, J, K, L, M
C      INTEGER NODADT
      my_real
c     .   F_inc11, F_inc12, F_inc13,
c     .   F_inc21, F_inc22, F_inc23,
c     .   F_inc31, F_inc32, F_inc33,
c     .   F_inc_inv11, F_inc_inv12, F_inc_inv13,
c     .   F_inc_inv21, F_inc_inv22, F_inc_inv23,
c     .   F_inc_inv31, F_inc_inv32, F_inc_inv33, det_F_inc_inv,
     .   F_inc(3,3), F_inc_inv(3,3), det_F,
     .   E(3,3), A(3,3), alpha(3),
     .   Q, P, tr, cos2, cos, term21, term31, eijk(3,3,3), term32(3,3), 
     .   term1(3,3), term2(3,3), term3(3,3), R_inc(3,3), temp(9), sum
c      .   nodadt, MIN_ASPECT, MIN_DEFV, IDELM

C-----------------------------------------------
C
      DO I=1,NEL
        DXX(I)=PX1(I)*VX1(I)+PX2(I)*VX2(I)+PX3(I)*VX3(I)+PX4(I)*VX4(I)
        DYY(I)=PY1(I)*VY1(I)+PY2(I)*VY2(I)+PY3(I)*VY3(I)+PY4(I)*VY4(I)
        DZZ(I)=PZ1(I)*VZ1(I)+PZ2(I)*VZ2(I)+PZ3(I)*VZ3(I)+PZ4(I)*VZ4(I)
        DXY(I)=PY1(I)*VX1(I)+PY2(I)*VX2(I)+PY3(I)*VX3(I)+PY4(I)*VX4(I)
        DXZ(I)=PZ1(I)*VX1(I)+PZ2(I)*VX2(I)+PZ3(I)*VX3(I)+PZ4(I)*VX4(I)
        DYX(I)=PX1(I)*VY1(I)+PX2(I)*VY2(I)+PX3(I)*VY3(I)+PX4(I)*VY4(I)
        DYZ(I)=PZ1(I)*VY1(I)+PZ2(I)*VY2(I)+PZ3(I)*VY3(I)+PZ4(I)*VY4(I)
        DZX(I)=PX1(I)*VZ1(I)+PX2(I)*VZ2(I)+PX3(I)*VZ3(I)+PX4(I)*VZ4(I)
        DZY(I)=PY1(I)*VZ1(I)+PY2(I)*VZ2(I)+PY3(I)*VZ3(I)+PY4(I)*VZ4(I)
        
c        DO M = 1, 9
c          PRINT *, 'R_OLD is', R_OLD(I, M)
c        ENDDO

c        R_OLD(I, 1) = ONE
c        R_OLD(I, 5) = ONE
c        R_OLD(I, 9) = ONE

c        DO M = 1, 6
c          PRINT *, 'sig is', SIG(I, M)
c        ENDDO
c        F_inc11 = ONE + (PX1(I)*VX1(I)+PX2(I)*VX2(I)+PX3(I)*VX3(I)+PX4(I)*VX4(I))*DT1
c        F_inc22 = ONE + (PY1(I)*VY1(I)+PY2(I)*VY2(I)+PY3(I)*VY3(I)+PY4(I)*VY4(I))*DT1
c        F_inc33 = ONE + (PZ1(I)*VZ1(I)+PZ2(I)*VZ2(I)+PZ3(I)*VZ3(I)+PZ4(I)*VZ4(I))*DT1

c        F_inc12 = (PY1(I)*VX1(I)+PY2(I)*VX2(I)+PY3(I)*VX3(I)+PY4(I)*VX4(I))*DT1
c        F_inc13 = (PZ1(I)*VX1(I)+PZ2(I)*VX2(I)+PZ3(I)*VX3(I)+PZ4(I)*VX4(I))*DT1
c        F_inc21 = (PX1(I)*VY1(I)+PX2(I)*VY2(I)+PX3(I)*VY3(I)+PX4(I)*VY4(I))*DT1

c        F_inc23 = (PZ1(I)*VY1(I)+PZ2(I)*VY2(I)+PZ3(I)*VY3(I)+PZ4(I)*VY4(I))*DT1
c        F_inc31 = (PX1(I)*VZ1(I)+PX2(I)*VZ2(I)+PX3(I)*VZ3(I)+PX4(I)*VZ4(I))*DT1
c        F_inc32 = (PY1(I)*VZ1(I)+PY2(I)*VZ2(I)+PY3(I)*VZ3(I)+PY4(I)*VZ4(I))*DT1

c        det_F_inc_inv = F_inc11*(F_inc22*F_inc33 - F_inc23*F_inc32)
c     1                - F_inc12*(F_inc21*F_inc33 - F_inc23*F_inc31)
c     1                + F_inc13*(F_inc21*F_inc32 - F_inc22*F_inc31)


c        if (det_F_inc_inv == 0.0d0) then
c              print *, 'F_inc 不可逆'
c              stop
c        end if

        ! 计算逆矩阵的各个分量（使用代数余子式并除以 det）
c        F_inc_inv11 =  (F_inc22*F_inc33 - F_inc23*F_inc32) / det_F_inc_inv
c        F_inc_inv12 = -(F_inc12*F_inc33 - F_inc13*F_inc32) / det_F_inc_inv
c        F_inc_inv13 =  (F_inc12*F_inc23 - F_inc13*F_inc22) / det_F_inc_inv

c        F_inc_inv21 = -(F_inc21*F_inc33 - F_inc23*F_inc31) / det_F_inc_inv
c        F_inc_inv22 =  (F_inc11*F_inc33 - F_inc13*F_inc31) / det_F_inc_inv
c        F_inc_inv23 = -(F_inc11*F_inc23 - F_inc13*F_inc21) / det_F_inc_inv

c        F_inc_inv31 =  (F_inc21*F_inc32 - F_inc22*F_inc31) / det_F_inc_inv
c        F_inc_inv32 = -(F_inc11*F_inc32 - F_inc12*F_inc31) / det_F_inc_inv
c        F_inc_inv33 =  (F_inc11*F_inc22 - F_inc12*F_inc21) / det_F_inc_inv

        F_inc(1,1) = ONE + DT1*(PX1(I)*VX1(I) + PX2(I)*VX2(I) + PX3(I)*VX3(I) + PX4(I)*VX4(I))
        F_inc(2,2) = ONE + DT1*(PY1(I)*VY1(I) + PY2(I)*VY2(I) + PY3(I)*VY3(I) + PY4(I)*VY4(I))
        F_inc(3,3) = ONE + DT1*(PZ1(I)*VZ1(I) + PZ2(I)*VZ2(I) + PZ3(I)*VZ3(I) + PZ4(I)*VZ4(I))

        F_inc(1,2) = DT1*(PY1(I)*VX1(I) + PY2(I)*VX2(I) + PY3(I)*VX3(I) + PY4(I)*VX4(I))
        F_inc(1,3) = DT1*(PZ1(I)*VX1(I) + PZ2(I)*VX2(I) + PZ3(I)*VX3(I) + PZ4(I)*VX4(I))
        F_inc(2,1) = DT1*(PX1(I)*VY1(I) + PX2(I)*VY2(I) + PX3(I)*VY3(I) + PX4(I)*VY4(I))

        F_inc(2,3) = DT1*(PZ1(I)*VY1(I) + PZ2(I)*VY2(I) + PZ3(I)*VY3(I) + PZ4(I)*VY4(I))
        F_inc(3,1) = DT1*(PX1(I)*VZ1(I) + PX2(I)*VZ2(I) + PX3(I)*VZ3(I) + PX4(I)*VZ4(I))
        F_inc(3,2) = DT1*(PY1(I)*VZ1(I) + PY2(I)*VZ2(I) + PY3(I)*VZ3(I) + PY4(I)*VZ4(I))

c        PRINT *, 'F_inc is', F_inc
c     ——— 计算行列式 det_F ————————————————————————————————
        det_F = F_inc(1,1)*(F_inc(2,2)*F_inc(3,3) - F_inc(2,3)*F_inc(3,2))
     1        - F_inc(1,2)*(F_inc(2,1)*F_inc(3,3) - F_inc(2,3)*F_inc(3,1))
     1        + F_inc(1,3)*(F_inc(2,1)*F_inc(3,2) - F_inc(2,2)*F_inc(3,1))

        if (det_F .eq. ZERO) then
          print *, 'F_inc 不可逆 at I =', I
          stop
        endif

c     ——— 用代数余子式法求逆矩阵 F_inc_inv ————————————————————
        F_inc_inv(1,1) =  (F_inc(2,2)*F_inc(3,3) - F_inc(2,3)*F_inc(3,2))/det_F
        F_inc_inv(1,2) = -(F_inc(1,2)*F_inc(3,3) - F_inc(1,3)*F_inc(3,2))/det_F
        F_inc_inv(1,3) =  (F_inc(1,2)*F_inc(2,3) - F_inc(1,3)*F_inc(2,2))/det_F

        F_inc_inv(2,1) = -(F_inc(2,1)*F_inc(3,3) - F_inc(2,3)*F_inc(3,1))/det_F
        F_inc_inv(2,2) =  (F_inc(1,1)*F_inc(3,3) - F_inc(1,3)*F_inc(3,1))/det_F
        F_inc_inv(2,3) = -(F_inc(1,1)*F_inc(2,3) - F_inc(1,3)*F_inc(2,1))/det_F

        F_inc_inv(3,1) =  (F_inc(2,1)*F_inc(3,2) - F_inc(2,2)*F_inc(3,1))/det_F
        F_inc_inv(3,2) = -(F_inc(1,1)*F_inc(3,2) - F_inc(1,2)*F_inc(3,1))/det_F
        F_inc_inv(3,3) =  (F_inc(1,1)*F_inc(2,2) - F_inc(1,2)*F_inc(2,1))/det_F

c     ——— 1. 构造单位矩阵 I        
        do 10 J = 1,3
          do 20 K = 1,3
            if (J .eq. K) then
              E(J,K) = ONE
            else
              E(J,K) = ZERO
            endif
   20     continue
   10 continue
c     ——— 2. A = I - F_inc_inv
        do 30 J = 1,3
          do 40 K = 1,3
            A(J,K) = E(J,K) - F_inc_inv(J,K)
   40 continue
   30 continue
c     ——— 3. alpha_k = -ε_{i j k} A_{j k}
        alpha(1) = - ( A(2,3) - A(3,2) )
        alpha(2) = - ( A(3,1) - A(1,3) )
        alpha(3) = - ( A(1,2) - A(2,1) )
c     ——— 4. Q = 0.25 * α·α
        Q = ONE/FOUR * ( alpha(1)**2 + alpha(2)**2 + alpha(3)**2 )    
c     ——— 5. P = 0.25*(tr(F_inv)-1)^2
        tr = F_inc_inv(1,1) + F_inc_inv(2,2) + F_inc_inv(3,3)
        P = ONE/FOUR * ( tr - ONE )**2
c    6. cos2, cos
        cos2 = P + THREE*P**2*(ONE-(P+Q))/(P+Q)**2 - TWO*P**3*(ONE-(P+Q))/(P+Q)**3
        cos = SIGN(SQRT(cos2), tr - ONE)
c    7. term1 = I * cos
        do 15 J = 1,3
          do 16 K = 1,3
            term1(J,K) = E(J,K) * cos
   16 continue
   15 continue
c    8. term2 =  (Q small?) or else
        if (Q .le. EM02) then
          term21 = ONE/EIGHT
     1            + Q*( (P**2-TWELVE*(P-ONE)) / (THIRTY2*P**2) )
     1            + Q**2*(P-TWO)*(P**2-TEN*P+THIRTY2)/(SIXTY4*P**3)
     1            + Q**3*(1104.0-992.0*P+376.0*P**2
     1                      -72.0*P**3+FIVE*P**4)/(FIVEHUNDRED12*P**4)
          do 17 J = 1,3
            do 18 K = 1,3
              term2(J,K) = term21 * alpha(J) * alpha(K)
   18 continue
   17 continue
        else
          do 29 J = 1,3
            do 49 K = 1,3
              term2(J,K) = (ONE - cos) / (FOUR*Q) * alpha(J) * alpha(K)
   49 continue
   29 continue
        endif
        term31 = HALF * SQRT( (P*Q*(THREE-Q) + P**3 + Q**2) / (P+Q)**3 )
C     初始化 Levi-Civita 张量
c        call init_levi(eijk)
C     初始化为全 0
        do 11 J = 1,3
        do 12 K = 1,3
        do 13 L = 1,3
          eijk(J,K,L) = ZERO
  13  continue
  12  continue
  11  continue

C     设置非零项
        eijk(1,2,3) =  ONE
        eijk(2,3,1) =  ONE
        eijk(3,1,2) =  ONE
        eijk(1,3,2) = -ONE
        eijk(2,1,3) = -ONE
        eijk(3,2,1) = -ONE
C     term32(i,j) = sum_k e(i,j,k) * alpha(k)
      do 59 J = 1,3
      do 69 K = 1,3
        term32(J,K) = ZERO
        do 79 L = 1,3
          term32(J,K) = term32(J,K) + eijk(J,K,L) * alpha(L)
 79   continue
 69   continue
 59   continue

C     term3(i,j) = term31 * term32(i,j)
      do 89 J = 1,3
      do 99 K = 1,3
        term3(J,K) = term31 * term32(J,K)
 99   continue
 89   continue
c    14. R_inc = term1 + term2 - term3
        do 240 J = 1,3
          do 250 K = 1,3
            R_inc(J,K) = term1(J,K) + term2(J,K) - term3(J,K)
250   continue
240   continue
C        PRINT *, 'R_inc is', R_inc
C     计算 temp = R_inc * R_OLD，保持行优先存储
        DO J = 1, 3
          DO K = 1, 3
            sum = 0.0
            DO L = 1, 3
              sum = sum + R_inc(J, L) * R_OLD(I, (L - 1) * 3 + K)
c            PRINT *, 'sum is', sum
            END DO
            temp((J - 1) * 3 + K) = sum
          END DO
        END DO

C     写回到 R_OLD
        DO J = 1, 9
          R_OLD(I, J) = temp(J)
        END DO
C     
C      PRINT *, 'TIMESTEP is', TIMESTEP
C        DO J = 1, 9
C          PRINT *, 'R_OLD is', R_OLD(I,J)
C        END DO
      ENDDO
c      PRINT *, 'DT1 is', DT1
      RETURN
      END
c       my_real
c      .   intent(out) eijk(3,3,3)

c       integer i, j, k

c C     初始化为全 0
c       do 10 i = 1,3
c       do 20 j = 1,3
c       do 30 k = 1,3
c         eijk(i,j,k) = ZERO
c  30   continue
c  20   continue
c  10   continue

c C     设置非零项
c       eijk(1,2,3) =  ONE
c       eijk(2,3,1) =  ONE
c       eijk(3,1,2) =  ONE
c       eijk(1,3,2) = -ONE
c       eijk(2,1,3) = -ONE
c       eijk(3,2,1) = -ONE

c       return
c       end

