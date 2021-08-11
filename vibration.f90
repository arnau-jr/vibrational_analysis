module vibration
      use molecule
      use force_field
      implicit none

      real*8,allocatable :: normal_base(:,:),normal_eigenvalues(:),normal_frequencies(:)
      
      contains

      subroutine comp_normal_modes()
            implicit none
            real*8  :: H(3*Natoms,3*Natoms),Hm(3*Natoms,3*Natoms),d(3*Natoms)
            real*8  :: work(100)
            integer :: a,b,p,q,i,j

            !Build mass weighted matrix of force constants
            H = build_hessian(Natoms,xyz_eq)
            do a=1,Natoms
                  do b=1,Natoms
                        do p=0,2
                              do q=0,2
                                    i = 3*a - p
                                    j = 3*b - q
      
                                    Hm(i,j) = H(i,j)/sqrt(M_mol(a)*M_mol(b))
                              end do
                        end do
                  end do
            end do

            allocate(normal_base(3*Natoms,3*Natoms),normal_eigenvalues(3*Natoms),normal_frequencies(3*Natoms))
            !Diagonalize using LAPACK
            call dsyev("V","U",3*Natoms,Hm,3*Natoms,d,work,100,i)

            ! normal_eigenvalues = d

            !Reorganize normal modes to put the first 6 (0 eigenvalue) in the last 6
            do i=1,6
                  normal_base(:,3*Natoms-i+1) = Hm(:,i)
                  normal_eigenvalues(3*Natoms-i+1) = d(i)
            end do
            normal_base(:,:3*Natoms-6) = Hm(:,7:)
            normal_eigenvalues(:3*Natoms-6) = d(7:)

            !Convert eigenvalues to cm-1 frequencies
            normal_frequencies(:3*Natoms-6) = sqrt(normal_eigenvalues(:3*Natoms-6))*hbar_cm_dps
            normal_frequencies(3*Natoms-5:) = sqrt(abs(normal_eigenvalues(3*Natoms-5:)))*hbar_cm_dps
      end subroutine comp_normal_modes

end module vibration