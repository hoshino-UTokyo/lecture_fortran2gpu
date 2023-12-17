!********************************************
! atomic命令を利用してhistgramを並列化する例
!********************************************

module test
  implicit none
contains
  
  subroutine histgram(x, y)
    integer,dimension(:),intent(in) :: x
    integer,dimension(:),intent(out) :: y 
    integer :: i,j
    
    !***********************************************
    ! atomic指示文は並列化したループの内側に現れる。
    ! NVIDIAのハイエンドGPUでは、atomic addに対して
    ! ハードウェアサポートがあるので、それなりに速い。
    !***********************************************
    
    !$acc kernels
    !$acc loop independent
    do i = 1, size(x)
       j = x(i)
       !$acc atomic 
       y(j) = y(j) + 1
       !$acc end atomic
    end do
    !$acc end kernels
    
  end subroutine histgram
    
end module test

program main
  use omp_lib
  use test
  implicit none
  integer :: nt = 100
  integer :: n = 10000000
  integer :: ny = 100
  integer,allocatable,dimension(:) :: y,x
  real(KIND=8) :: a
  real(KIND=8) :: t1,t2,t3,t4
  real(KIND=8) :: ans, dummy(1)
  integer :: i,times
  
  allocate(x(n),y(ny))

  y(:) = 0
  do i = 1, n
     call random_number(a)
     x(i) = int(a*100)+1
  end do
  
  !******************************
  ! 実行時間計測時の注意点。
  ! プログラム中でGPUに最初に触れる際、GPU起動に多少時間がかかる。
  ! Unified Memoryを用いる場合、プログラムの開始時にすでに起動処理がかかるので、
  ! プログラム中の最初の指示文では表面化しない。
  ! 下のkernelsにかかる時間は通常のOpenACCより小さい
  !******************************

  t1 = omp_get_wtime()
  !$acc kernels
  dummy(:) = 0
  !$acc end kernels
  t2 = omp_get_wtime()

  print *, "GPU boot [s]: ", t2-t1

  !******************************
  ! NT回の繰り返しに対して一回だけCPU-GPU間のデータ転送を行う。
  !******************************

  t1 = omp_get_wtime()
  !******************************
  ! Unified Memoryの場合、data指示文はいらない。
  ! ページフォルト（扱おうとしたGPU/CPUのデータがCPU/GPUより古い）が発生した時に
  ! 自動で転送される仕組みなので、データ転送の時間を正確に図るのは難しい。
  ! データ転送はカーネル関数内で発生する。
  ! データ転送時間を除外するなら、以下のように最初の一回目を計測から外す。
  !******************************
  call histgram(x,y)
  t2 = omp_get_wtime()
  !$acc kernels
  y(:) = 0
  !$acc end kernels
  t3 = omp_get_wtime()
  do times = 1, nt
     call histgram(x,y)
  end do
  t4 = omp_get_wtime()

  !******************************
  ! 答えはyの総和を出力。n*ntになるはず。
  !******************************
  print *, "ans :", sum(y)
  print *, "First run [s]: ", (t2-t1)
  print *, "Time / iteration [s]: ", (t4-t3)/nt

  deallocate(x,y)

end program main

