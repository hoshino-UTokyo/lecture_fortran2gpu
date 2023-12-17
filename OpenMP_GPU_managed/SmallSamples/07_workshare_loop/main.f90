!********************************************
! 単純な一重ループとしてDAXPY (Y = a*X + Y)を並列化する例
! DAXPY計算をNT回繰り返すプログラム
!********************************************

module test
  implicit none
contains
  
  subroutine daxpy(x, y, a, n)
    real(KIND=8),dimension(:),intent(out) :: y ! intent属性はできるだけつけた方がコンパイラの最適化が効きやすい
    real(KIND=8),dimension(:),intent(in) :: x
    real(KIND=8),intent(in) :: a
    integer,intent(in) :: n
    integer :: i

#if 0
    !***********************************************
    ! OpenACCと同じノリで以下のように書くと、正しい結果は出るが非常に遅い。
    ! 並列化されずに１スレッドで動く。
    !***********************************************

    !$omp target 
    y(:) = a * x(:) + y(:)
    !$omp end target
#endif

#if 0
    !***********************************************
    ! 以下もコンパイルメッセージでは並列化されたようにも見えるが非常に遅い。
    !***********************************************

    !$omp target
    !$omp workshare
    y(:) = a * x(:) + y(:)
    !$omp end workshare
    !$omp end target workshare
#endif

    print *, "Not available"
    stop

  end subroutine daxpy
    
end module test

program main
  use omp_lib
  use test
  implicit none
  integer :: nt = 100
  integer :: n = 1000000000
  real(KIND=8),allocatable,dimension(:) :: y,x
  real(KIND=8) :: a
  real(KIND=8) :: t1,t2,t3,t4,t5,t6
  real(KIND=8) :: ans, dummy(1)
  integer :: i,times
  
  allocate(y(n),x(n))

  y(:) = 0.0d0
  x(:) = 1.0d0
  a = 2.0d0
  
  !******************************
  ! 実行時間計測時の注意点。
  ! プログラム中でGPUに最初に触れる際、GPU起動に多少時間がかかる。
  ! Unified Memoryを用いないOpenMPプログラムでは、target指示文を触った時に表面化する。
  ! 次のループはGPUの起動時間を除外するためだけのもので、プログラム上の意味はない。
  !******************************

  t1 = omp_get_wtime()
  !$omp target 
  dummy(1) = 0
  !$omp end target
  t2 = omp_get_wtime()

  print *, "GPU boot [s]: ", t2-t1

  !******************************
  ! DAXPY() NT回の繰り返しに対して一回だけCPU-GPU間のデータ転送を行う。
  !******************************

  t1 = omp_get_wtime()
  !******************************
  ! Unified Memoryの場合、data指示文はいらない。
  ! ページフォルト（扱おうとしたGPU/CPUのデータがCPU/GPUより古い）が発生した時に
  ! 自動で転送される仕組みなので、データ転送の時間を正確に図るのは難しい。
  ! データ転送はdaxpyのカーネル関数内で発生する。
  ! データ転送時間を除外するなら、以下のように最初の一回目を計測から外す。
  !******************************
  call daxpy(x,y,a,n)
  t2 = omp_get_wtime()
  !$omp target
  !$omp loop 
  do i = 1, n
     y(i) = 0.0d0
  end do
  !$omp end target
  t3 = omp_get_wtime()
  do times = 1, nt
     call daxpy(x,y,a,n)
  end do
  t4 = omp_get_wtime()

  
  !******************************
  ! 答えはyの平均値を出力。2*ntになるはず。
  ! 配列yはGPUで更新されておりCPU側が古くなっているので、
  ! sum(y)の計算時にGPU->CPUの通信が発生する。
  !******************************
  ans = sum(y)/n
  print *, "ans :", ans
  print *, "First run [s]: ", (t2-t1)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(x,y)

end program main

