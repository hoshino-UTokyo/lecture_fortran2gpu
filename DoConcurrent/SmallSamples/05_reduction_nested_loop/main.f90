!********************************************
! 総和など、reduction計算を並列化する例
!********************************************

module test
  implicit none
contains
  
  subroutine reduction(y,sum,n)
    real(KIND=8),dimension(:,:,:),intent(in) :: y ! intent属性はできるだけつけた方がコンパイラの最適化が効きやすい
    real(KIND=8),intent(inout) :: sum
    integer,intent(in) :: n
    integer :: i, j, k
    
    !***********************************************
    ! 多重ループの場合は、reductionが必要なループ全てに付ける。
    ! collapseを使って一重化しても良い。
    !***********************************************
    
    do concurrent (k = 1: n) local(i,j)
       do concurrent (j = 1: n) local(i)
          do concurrent (i = 1: n)
             sum = sum + y(i,j,k)
          end do
       end do
    end do
    
  end subroutine reduction
    
end module test

program main
  use omp_lib
  use test
  implicit none
  integer :: nt = 100
  integer :: n = 1000
  real(KIND=8),allocatable,dimension(:,:,:) :: y
  real(KIND=8) :: sum
  real(KIND=8) :: t1,t2,t3,t4
  real(KIND=8) :: ans, dummy(1)
  integer :: i,times
  
  allocate(y(n,n,n))

  y(:,:,:) = 1.0d0
  sum = 0.0d0
  
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
  ! NT回の繰り返しに対して一回だけCPU-GPU間のデータ転送を行う。
  !******************************

  !******************************
  ! Unified Memoryの場合、data指示文はいらない。
  ! ページフォルト（扱おうとしたGPU/CPUのデータがCPU/GPUより古い）が発生した時に
  ! 自動で転送される仕組みなので、データ転送の時間を正確に図るのは難しい。
  ! データ転送はカーネル関数内で発生する。
  ! データ転送時間を除外するなら、以下のように最初の一回目を計測から外す。
  !******************************
  t1 = omp_get_wtime()
  call reduction(y,sum,n)
  t2 = omp_get_wtime()
  sum = 0.0d0
  t3 = omp_get_wtime()
  do times = 1, nt
     call reduction(y,sum,n)
  end do
  t4 = omp_get_wtime()
  
  !******************************
  ! 答えはxとyの総和nt回足し合わせたもの。
  !******************************
  print *, "ans :", sum
  print *, "First run [s]: ", (t2-t1)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(y)

end program main

