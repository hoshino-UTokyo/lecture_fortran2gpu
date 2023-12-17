!********************************************
! 単純な一重ループとしてDAXPY (Y = a*X + Y)を並列化する例
! DAXPY計算をNT回繰り返すプログラム
!********************************************

module test
  implicit none
contains

  real(KIND=8) function madd(a,x,y) 
    real(KIND=8),intent(in) :: a,x,y
    !***********************************************
    ! この関数はpure function（引数が全てintent(in)で戻り値がスカラ一つ。グローバル変数の書き換えなども無い）
    ! だと思うので、do concurrentの仕様上は呼び出せるはずだが、NVIDIAコンパイラ22.7にはnot pureだと怒られる。
    !***********************************************

    madd = a * x + y
  end function madd
  
  subroutine daxpy(x, y, a, n)
    real(KIND=8),dimension(:),intent(out) :: y ! intent属性はできるだけつけた方がコンパイラの最適化が効きやすい
    real(KIND=8),dimension(:),intent(in) :: x
    real(KIND=8),intent(in) :: a
    integer,intent(in) :: n
    integer :: i

#if 0
    !***********************************************
    ! 下のコードはpure fucntionの呼び出しなので、仕様上は合ってると思う。
    !***********************************************
    
    do concurrent (i = 1: n)
       y(i) = madd(a,x(i),y(i))
    end do
#else
    print *, "Not available"
    stop
#endif
    
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
  ! do concurrentはUnified Memoryを用いるので、
  ! プログラムの開始時にすでに起動処理がかかり、
  ! 最初のGPU実行で起動時間が表面化しない。
  !******************************

  t1 = omp_get_wtime()
  do concurrent (i=1:1)
     dummy(1) = 0.0d0
  end do
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
  do concurrent (i=1:n)
     y(i) = 0.0d0
  end do
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
