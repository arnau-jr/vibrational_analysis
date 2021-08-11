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

      function unpack_coords(packed) result (unpacked)
            implicit none
            real*8 :: packed(:,:)
            real*8 :: unpacked(3*size(packed,2))
            integer :: i
            do i=1,size(packed,2)
                  unpacked(3*i-2) = packed(1,i)
                  unpacked(3*i-1) = packed(2,i)
                  unpacked(3*i  ) = packed(3,i)
            enddo
      end function unpack_coords

      function pack_coords(unpacked) result(packed)
            implicit none
            real*8 :: unpacked(:)
            real*8 :: packed(3,int(size(unpacked)/3))
            integer :: i
            do i=1,int(size(unpacked)/3)
                  packed(1,i) = unpacked(3*i-2)
                  packed(2,i) = unpacked(3*i-1)
                  packed(3,i) = unpacked(3*i  )
            enddo
      end function pack_coords

      function cart_to_normal(xyz) result(xyz_nm)
            implicit none
            real*8              :: xyz(3,Natoms),xyz_nm(3*Natoms)
            real*8              :: xyz_mw(3,Natoms),xyz_mw_unpacked(3*Natoms)

            xyz_mw(1,:) = xyz(1,:)*sqrt(M_mol)
            xyz_mw(2,:) = xyz(2,:)*sqrt(M_mol)
            xyz_mw(3,:) = xyz(3,:)*sqrt(M_mol)

            xyz_mw_unpacked = unpack_coords(xyz_mw)

            xyz_nm = matmul(transpose(normal_base),xyz_mw_unpacked)
      end function cart_to_normal

            function normal_to_cart(xyz_nm) result(xyz)
            implicit none
            real*8              :: xyz_nm(3*Natoms),xyz(3,Natoms)
            real*8              :: xyz_mw(3,Natoms),xyz_mw_unpacked(3*Natoms)

            xyz_mw_unpacked = matmul(normal_base,xyz_nm)

            xyz_mw = pack_coords(xyz_mw_unpacked)

            xyz(1,:) = xyz_mw(1,:)/sqrt(M_mol)
            xyz(2,:) = xyz_mw(2,:)/sqrt(M_mol)
            xyz(3,:) = xyz_mw(3,:)/sqrt(M_mol)
      end function normal_to_cart

      function build_eckart_matrix() result(EM)
            implicit none
            real*8  :: xpa,xma,ypa,yma,zpa,zma
            real*8  :: EM(4,4)
            integer :: a
            EM = 0.d0
            do a=1,Natoms
                  xpa = xyz_eq(1,a) + xyz_cm(1,a)
                  xma = xyz_eq(1,a) - xyz_cm(1,a)

                  ypa = xyz_eq(2,a) + xyz_cm(2,a)
                  yma = xyz_eq(2,a) - xyz_cm(2,a)

                  zpa = xyz_eq(3,a) + xyz_cm(3,a)
                  zma = xyz_eq(3,a) - xyz_cm(3,a)

                  !Diagonal
                  EM(1,1) = EM(1,1) + M_mol(a)*(xma**2 + yma**2 + zma**2)
                  EM(2,2) = EM(2,2) + M_mol(a)*(xma**2 + ypa**2 + zpa**2)
                  EM(3,3) = EM(3,3) + M_mol(a)*(xpa**2 + yma**2 + zpa**2)
                  EM(4,4) = EM(4,4) + M_mol(a)*(xpa**2 + ypa**2 + zma**2)

                  !Rest of 1st row
                  EM(1,2) = EM(1,2) + M_mol(a)*(ypa*zma - yma*zpa)
                  EM(1,3) = EM(1,3) + M_mol(a)*(xma*zpa - xpa*zma)
                  EM(1,4) = EM(1,4) + M_mol(a)*(xpa*yma - xma*ypa)

                  !Rest of 2nd row
                  EM(2,3) = EM(2,3) + M_mol(a)*(xma*yma - xpa*ypa)
                  EM(2,4) = EM(2,4) + M_mol(a)*(xma*zma - xpa*zpa)

                  !Rest of 3rd row
                  EM(3,4) = EM(3,4) + M_mol(a)*(yma*zma - ypa*zpa)
            end do

            !Symmetrize
            EM(2,1) = EM(1,2)
            EM(3,1) = EM(1,3)
            EM(4,1) = EM(1,4)

            EM(3,2) = EM(2,3)
            EM(4,2) = EM(2,4)

            EM(4,3) = EM(3,4)
      end function build_eckart_matrix

      function build_direction_cosine_matrix(V) result (U)
            implicit none
            real*8 :: V(4)
            real*8 :: U(3,3)
            real*8 :: q0,q1,q2,q3
            q0 = V(1)
            q1 = V(2)
            q2 = V(3)
            q3 = V(4)

            !1st row
            U(1,1) = q0**2 + q1**2 - q2**2 - q3**2
            U(1,2) = 2.d0*(q1*q2 + q0*q3)
            U(1,3) = 2.d0*(q1*q3 - q0*q2)
            
            !2nd row
            U(2,1) = 2.d0*(q1*q2 - q0*q3)
            U(2,2) = q0**2 - q1**2 + q2**2 - q3**2
            U(2,3) = 2.d0*(q2*q3 + q0*q1)

            !3rd row
            U(3,1) = 2.d0*(q1*q3 + q0*q2)
            U(3,2) = 2.d0*(q2*q3 - q0*q1)
            U(3,3) = q0**2 - q1**2 - q2**2 + q3**2
      end function build_direction_cosine_matrix

      subroutine get_eckart_frame(U)
            implicit none
            real*8             :: EM(4,4),U(3,3)
            real*8             :: work(100),d(4)
            integer            :: i,ierror

            call update_cm_coords()

            EM = build_eckart_matrix()

            call dsyev("V","U",4,EM,4,d,work,100,ierror)

            U = build_direction_cosine_matrix(EM(:,1))

            do i=1,Natoms
                  xyz_eckart(:,i) = matmul(transpose(U),xyz_eq(:,i))
            end do
            call check_eckart_conditions()
      end subroutine get_eckart_frame

      subroutine check_eckart_conditions()
            implicit none
            real*8           :: tra_cond(3),rot_cond(3),comb_cond(3),disp(3,Natoms)
            real*8,parameter :: eps=1.d-10
            integer          :: i

            disp = xyz_cm - xyz_eckart
            tra_cond = 0.d0
            rot_cond = 0.d0
            comb_cond = 0.d0
            do i=1,Natoms
                  tra_cond  = tra_cond  + M_mol(i)*disp(:,i)
                  rot_cond  = rot_cond  + M_mol(i)*cross_product(xyz_eckart(:,i),disp(:,i))
                  comb_cond = comb_cond + M_mol(i)*cross_product(xyz_eckart(:,i),xyz_cm(:,i))
            end do

            if(any(tra_cond > eps) .or. any(rot_cond > eps) .or. any(comb_cond > eps)) then
                  print*,"Eckart conditions not satisfied"
                  print*,"Eckart Conditons:"
                  print"(A,2X,3(E16.8,2X))","Translational:",tra_cond
                  print"(A,2X,3(E16.8,2X))","Rotational:   ",rot_cond
                  print"(A,2X,3(E16.8,2X))","Combined:     ",comb_cond
                  print*,""
            end if
      end subroutine check_eckart_conditions

      function build_pseudo_inertia_moment() result(I)
            implicit none
            real*8  :: I(3,3)
            real*8  :: delta_term
            integer :: j,a,b

            delta_term = 0.d0
            do j=1,Natoms
                  delta_term = delta_term + M_mol(j)*sum(xyz_eckart(:,j)*xyz_cm(:,j))
            end do

            I = 0.d0
            do a=1,3
                  do b=1,3
                        if(a==b) I(a,b) = I(a,b) + delta_term
                        do j=1,Natoms
                              I(a,b) = I(a,b) - M_mol(j)*xyz_cm(a,j)*xyz_eckart(b,j)
                        end do
                  end do
            end do
      end function build_pseudo_inertia_moment

      function invert_matrix(A) result(Ainv)
            implicit none
            real*8 :: A(3,3),Ainv(3,3)
            real*8 :: det
            det = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                  - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                  + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
            det = 1.d0/det

            Ainv(1,1) = +det * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
            Ainv(2,1) = -det * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
            Ainv(3,1) = +det * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
            Ainv(1,2) = -det * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
            Ainv(2,2) = +det * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
            Ainv(3,2) = -det * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
            Ainv(1,3) = +det * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
            Ainv(2,3) = -det * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
            Ainv(3,3) = +det * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
      end function invert_matrix

      function comp_kinetic_energy(vel) result(K)
            implicit none
            real*8 :: vel(:,:),K
            integer :: i
            K = sum(0.5d0*M_mol*sum(vel**2,1))
      end function comp_kinetic_energy

      function comp_cm_vel(vel) result(vel_cm)
            implicit none
            real*8 :: vel(3,Natoms),vel_cm(3)
            integer :: i
            vel_cm = 0.d0
            do i=1,Natoms
                  vel_cm = vel_cm + M_mol(i)*vel(:,i)
            end do
            vel_cm = vel_cm/sum(M_mol)
      end function comp_cm_vel

      function comp_cm_energy(vel) result(Ecm)
            implicit none
            real*8 :: vel(3,Natoms),vel_cm(3)
            real*8 :: Ecm
            vel_cm = comp_cm_vel(vel)
            Ecm = 0.5d0*sum(M_mol)*sum(vel_cm**2)
      end function comp_cm_energy

      function comp_pseudo_angular_moment(vel) result(l)
            implicit none
            real*8  :: vel(3,Natoms)
            real*8  :: l(3)
            integer :: i
            l = 0.d0
            do i=1,Natoms
                  ! l = l + cross_product(xyz_eckart(:,i),vel(:,i))
                  l = l + M_mol(i)*cross_product(xyz_eckart(:,i),vel(:,i))
            end do
      end function comp_pseudo_angular_moment

      function comp_angular_vel(vel) result(omega)
            implicit none
            real*8 :: vel(3,Natoms)
            real*8 :: omega(3)
            real*8 :: Iab(3,3),Iab_inv(3,3),l(3)

            Iab = build_pseudo_inertia_moment()
            Iab_inv = invert_matrix(Iab)

            l = comp_pseudo_angular_moment(vel)

            omega = matmul(Iab_inv,l)
      end function comp_angular_vel

      function comp_rotational_energy(omega) result(Erot)
            implicit none
            real*8  :: omega(3),xyz_cm(3,Natoms)
            real*8  :: Erot
            integer :: i
            Erot = 0.d0
            do i=1,Natoms
                  Erot = Erot + M_mol(i)*sum(cross_product(omega,xyz_cm(:,i))**2)
            end do
            Erot = 0.5d0*Erot
      end function comp_rotational_energy

end module vibration