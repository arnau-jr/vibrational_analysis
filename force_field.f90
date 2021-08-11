module force_field
      implicit none 

      !Constants and parameters
      real*8,parameter :: PI = 4.d0*atan(1.d0)

      !Number of bonds, angles, etc.
      integer :: Nbonds,Nangles,Ntorsions,Nimpropers

      !Pairs and values of each type
      integer,allocatable :: bond_pairs(:,:),angle_pairs(:,:),torsion_pairs(:,:),improper_pairs(:,:)
      real*8,allocatable  :: bond_vals(:),angle_vals(:),torsion_vals(:),improper_vals(:)
      !Units are standarized to be Angstroms for bonds and degrees for the rest
      !Energy units are internally kJ/mol, so the time unit is the dps (0.1 ps)
      !(Using g/mol as mass unit, as it also standard).

      !Selection of potential type, from the 4 terms we decide which type of potential
      !it is following a one character string
      character           :: pot_types(4)*1

      !Coefficients for every type of potential supported, not all will be used at runtime
      real*8,allocatable  :: harmonic_stretching(:,:) !Harmonic stretching (k,r)
      real*8,allocatable  :: morse_stretching(:,:) !Morse stretching (De,beta,r)
      real*8,allocatable  :: harmonic_bending(:,:) !Harmonic angles (k,a)
      real*8,allocatable  :: cosine_torsions(:,:) !Cosine torsions (An,delta,n)
      real*8,allocatable  :: cosine_impropers(:,:) !Cosine impropers (An,delta,n)

      !Unit conversions and other constants
      real*8,parameter :: kcal_to_kj = 4.184d0
      real*8,parameter :: R_kJ_mol_K = 8.31446261815324d-3
      real*8,parameter :: c_cm_dps = 2.99792458d-3 
      real*8,parameter :: h_cm_dps = 8065.6d0*4.135667696d-2
      ! real*8,parameter :: hbar_cm_dps = 8065.6d0*6.582119569d-3
      real*8,parameter :: hbar_cm_dps = 1.d0/(2.d0*pi*c_cm_dps)

      contains

      subroutine init_forcefield(unit,Natoms,xyz)
            implicit none
            integer,intent(in) :: unit,Natoms
            real*8,intent(in)  :: xyz(3,Natoms)
            character          :: dummy*90,units*90,pot_type*1
            integer            :: i,a,b,c,d,dummy2

            if(Natoms == 0) then
                  print*,"Molecule has 0 atoms, not initialized correctly, aborting..."
                  stop
            end if

            !Initialize
            Nbonds = 0
            Nangles = 0
            Ntorsions = 0
            Nimpropers = 0

            !Read stretching info
            read(unit,*)Nbonds
            read(unit,*)pot_type
            read(unit,*)units

            allocate(bond_pairs(2,Nbonds),bond_vals(Nbonds))
            if(pot_type=="H") then
                  allocate(harmonic_stretching(2,Nbonds))
            elseif(pot_type=="M") then
                  allocate(morse_stretching(3,Nbonds))
            else
                  print*,"Stretching potential not supported, aborting..."
                  stop
            end if
            pot_types(1) = pot_type

            !Read streching coefs
            do i=1,Nbonds
                  if(pot_type=="H") then
                        read(unit,*)a,b,dummy,harmonic_stretching(1,i),harmonic_stretching(2,i)
                        !Correct units if necessary
                        if(units=="KC") then
                              harmonic_stretching(1,i) = harmonic_stretching(1,i)*kcal_to_kj
                        end if
                  elseif(pot_type=="M") then
                        read(unit,*)a,b,dummy,morse_stretching(1,i),morse_stretching(2,i),&
                        morse_stretching(3,i)
                        if(units=="KC") then
                              morse_stretching(1,i) = morse_stretching(1,i)*kcal_to_kj
                        end if
                  end if
                  bond_pairs(:,i) = (/a,b/)
            end do
            read(unit,*,end=10) !Check if end file, blank line should be there if continuing

            !Read bending info
            read(unit,*)Nangles
            read(unit,*)pot_type
            read(unit,*)units

            allocate(angle_pairs(3,Nangles),angle_vals(Nangles))
            if(pot_type=="H") then
                  allocate(harmonic_bending(2,Nangles))
            else
                  print*,"Bending potential not supported, aborting..."
                  stop
            end if
            pot_types(2) = pot_type

            !Read bending coefs
            do i=1,Nangles
                  if(pot_type=="H") then
                        read(unit,*)a,b,c,dummy,harmonic_bending(1,i),harmonic_bending(2,i)
                        if(units=="KC") then
                              harmonic_bending(1,i) = harmonic_bending(1,i)*kcal_to_kj
                        endif
                  end if
                  angle_pairs(:,i) = (/a,b,c/)
            end do
            read(unit,*,end=10)

            !Read torsion info
            read(unit,*)Ntorsions
            read(unit,*)pot_type
            read(unit,*)units

            allocate(torsion_pairs(4,Ntorsions),torsion_vals(Ntorsions))
            if(pot_type=="C") then
                  allocate(cosine_torsions(3,Ntorsions))
            else
                  print*,"Torsion potential not supported, aborting..."
                  stop
            end if
            pot_types(3) = pot_type

            do i=1,Ntorsions
                  if(pot_type=="C") then
                        read(unit,*)a,b,c,d,dummy,dummy2,cosine_torsions(1,i),cosine_torsions(2,i),&
                        cosine_torsions(3,i)
                        if(units=="KC") then
                              cosine_torsions(1,i) = cosine_torsions(1,i)*kcal_to_kj
                        endif
                  end if
                  torsion_pairs(:,i) = (/a,b,c,d/)
            end do
            read(unit,*,end=10)

            !Read impropers info
            read(unit,*)Nimpropers
            read(unit,*)pot_type
            read(unit,*)units

            allocate(improper_pairs(4,Nimpropers),improper_vals(Nimpropers))
            if(pot_type=="C") then
                  allocate(cosine_impropers(3,Nimpropers))
            else
                  print*,"Improper potential not supported, aborting..."
                  stop
            end if
            pot_types(4) = pot_type

            do i=1,Nimpropers
                  if(pot_type=="C") then
                        read(unit,*)a,b,c,d,dummy,cosine_impropers(1,i),cosine_impropers(2,i),&
                        cosine_impropers(3,i)
                        if(units=="KC") then
                              cosine_impropers(1,i) = cosine_impropers(1,i)*kcal_to_kj
                        endif
                  end if
                  improper_pairs(:,i) = (/a,b,c,d/)
                  ! improper_vals(i) = comp_improper(xyz(:,a),xyz(:,b),&
                  ! xyz(:,c),xyz(:,d))
            enddo
            10 continue
            call comp_all(Natoms,xyz)
      end subroutine init_forcefield

      function unit_cross(u1,u2) result(u3)
            implicit none
            real*8 :: u1(3),u2(3)
            real*8 :: u3(3)

            u3(1) = u1(2)*u2(3) - u1(3)*u2(2)
            u3(2) = - u1(1)*u2(3) + u1(3)*u2(1)
            u3(3) = u1(1)*u2(2) - u1(2)*u2(1)

            u3 = u3/sqrt(sum(u3**2))
      end function unit_cross

      function cross_product(v1,v2) result(v3)
            implicit none
            real*8 :: v1(3),v2(3)
            real*8 :: v3(3)

            v3(1) =   v1(2)*v2(3) - v1(3)*v2(2)
            v3(2) = - v1(1)*v2(3) + v1(3)*v2(1)
            v3(3) =   v1(1)*v2(2) - v1(2)*v2(1)
      end function cross_product

      function comp_norm(u) result(a)
            implicit none
            real*8 :: u(:),a
            a = sqrt(sum(u**2))
      end function comp_norm

      function comp_dist(c1,c2) result (d)
            implicit none
            real*8 :: c1(3),c2(3)
            real*8 :: d

            d = sqrt(sum((c1-c2)**2))
      end function comp_dist

      function comp_angle(c1,c2,c3) result(A)
            implicit none
            real*8 :: c1(3),c2(3),c3(3)
            real*8 :: u1(3),u2(3),proj
            real*8 :: A

            u1 = c2-c1
            u2 = c3-c2

            u1 = u1/comp_norm(u1)
            u2 = u2/comp_norm(u2)

            proj = sum(u1*u2)
            if(proj>=1.d0) then
                  A = 180.d0
            elseif(proj<=-1.d0) then
                  A = 0.d0
            else
                  A = 180.d0 - (180.d0/pi)*acos(proj)
            endif
      end function comp_angle

      function comp_torsion(c1,c2,c3,c4) result(T)
            implicit none
            real*8 :: c1(3),c2(3),c3(3),c4(3)
            real*8 :: u12(3),u23(3),u32(3),u43(3)
            real*8 :: u1232(3),u2343(3)
            real*8 :: proj,proj2
            real*8 :: T

            u12 = c2-c1
            u23 = c3-c2
            u32 = c2-c3
            u43 = c3-c4

            u12 = u12/sqrt(sum(u12**2))
            u23 = u23/sqrt(sum(u23**2))
            u32 = u32/sqrt(sum(u32**2))
            u43 = u43/sqrt(sum(u43**2))

            u1232 = unit_cross(u12,u32)
            u2343 = unit_cross(u23,u43)
            
            proj = sum(u1232*u2343)
            proj2 = sum(u1232*u43)

            ! if(proj>=1.d0) then
            !       T = 180.d0
            ! elseif(proj<=-1.d0) then
            !       T = 180.d0 - 180.d0*sign(1.d0,proj2)
            ! else
            !       T = 180.d0 - (180.d0/pi)*sign(1.d0,proj2)*acos(proj)
            ! endif
            if(proj>=1.d0) then
                  T = 0.d0*sign(1.d0,-proj2)
            elseif(proj<=-1.d0) then
                  T = 180.d0*sign(1.d0,-proj2)
            else
                  T = (180.d0/pi)*sign(1.d0,-proj2)*acos(proj)
            endif
      end function comp_torsion

      function comp_improper(c4,c1,c3,c2) result(T)
            implicit none
            real*8 :: c1(3),c2(3),c3(3),c4(3)
            real*8 :: u12(3),u23(3),u32(3),u43(3)
            real*8 :: u1232(3),u2343(3)
            real*8 :: proj,proj2
            real*8 :: T

            u12 = c2-c1
            u23 = c3-c2
            u32 = c2-c3
            u43 = c3-c4

            u12 = u12/sqrt(sum(u12**2))
            u23 = u23/sqrt(sum(u23**2))
            u32 = u32/sqrt(sum(u32**2))
            u43 = u43/sqrt(sum(u43**2))

            u1232 = unit_cross(u12,u32)
            u2343 = unit_cross(u23,u43)
            
            proj = sum(u1232*u2343)
            proj2 = sum(u1232*u43)
            ! proj2 = -1.d0

            if(proj>=1.d0) then
                  T = 0.d0*sign(1.d0,-proj2)
            elseif(proj<=-1.d0) then
                  T = 180.d0*sign(1.d0,-proj2)
            else
                  T = (180.d0/pi)*sign(1.d0,-proj2)*acos(proj)
            endif
      end function comp_improper

      subroutine comp_bonds(Natoms,xyz)
            implicit none
            integer,intent(in) :: Natoms
            real*8,intent(in)  :: xyz(3,Natoms)
            integer            :: bond

            do bond=1,Nbonds
                  bond_vals(bond) = comp_dist(xyz(:,bond_pairs(1,bond)),&
                  xyz(:,bond_pairs(2,bond)))
            enddo
      end subroutine comp_bonds

       subroutine comp_angles(Natoms,xyz)
            implicit none
            integer,intent(in) :: Natoms
            real*8,intent(in)  :: xyz(3,Natoms)
            integer            :: angle

            do angle=1,Nangles
                  angle_vals(angle) = comp_angle(xyz(:,angle_pairs(1,angle)),&
                  xyz(:,angle_pairs(2,angle)),xyz(:,angle_pairs(3,angle)))
            enddo
      end subroutine comp_angles

      subroutine comp_torsions(Natoms,xyz)
            implicit none
            integer,intent(in) :: Natoms
            real*8,intent(in)  :: xyz(3,Natoms)
            integer            :: torsion

            do torsion=1,Ntorsions
                  torsion_vals(torsion) = comp_torsion(xyz(:,torsion_pairs(1,torsion)),&
                  xyz(:,torsion_pairs(2,torsion)),&
                  xyz(:,torsion_pairs(3,torsion)),&
                  xyz(:,torsion_pairs(4,torsion)))
            enddo
      end subroutine comp_torsions

      subroutine comp_impropers(Natoms,xyz)
            implicit none
            integer,intent(in) :: Natoms
            real*8,intent(in)  :: xyz(3,Natoms)
            integer            :: improper

            do improper=1,Nimpropers
                  improper_vals(improper) = comp_improper(xyz(:,improper_pairs(1,improper)),&
                  xyz(:,improper_pairs(2,improper)),&
                  xyz(:,improper_pairs(3,improper)),&
                  xyz(:,improper_pairs(4,improper)))
            enddo
      end subroutine comp_impropers

      subroutine comp_all(Natoms,xyz)
            implicit none
            integer,intent(in) :: Natoms
            real*8,intent(in)  :: xyz(3,Natoms)
            call comp_bonds    (Natoms,xyz)
            call comp_angles   (Natoms,xyz)
            call comp_torsions (Natoms,xyz)
            call comp_impropers(Natoms,xyz)
      end subroutine comp_all

      real*8 function comp_bonds_energy() result(E)
            implicit none
            real*8  :: k,req
            real*8  :: De,beta
            integer :: i
            E = 0.
            if(pot_types(1)=="H") then
                  do i=1,Nbonds
                        k = harmonic_stretching(1,i)
                        req = harmonic_stretching(2,i)
                        E = E + k*(bond_vals(i)-req)**2
                  enddo
            elseif(pot_types(1)=="M") then
                  do i=1,Nbonds
                        De = morse_stretching(1,i)
                        beta = morse_stretching(2,i)
                        req = morse_stretching(3,i)
                        ! E = E + De*((1.d0-exp(-beta*(bond_vals(i)-req)))**2-1.d0)
                        E = E + De*((1.d0-exp(-beta*(bond_vals(i)-req)))**2)
                  enddo
            endif
      end function comp_bonds_energy


      real*8 function comp_angles_energy() result(E)
            implicit none
            real*8  :: k,aeq
            integer :: i
            E = 0.
            if(pot_types(2)=="H") then
                  do i=1,Nangles
                        k = harmonic_bending(1,i)
                        aeq = harmonic_bending(2,i)
                        E = E + k*((pi/180.d0)*(angle_vals(i)-aeq))**2
                  enddo
            end if
      end function comp_angles_energy


      real*8 function comp_torsions_energy() result(E)
            implicit none
            real*8  :: An,n,delta
            integer :: i
            E = 0.
            if(pot_types(3)=="C") then
                  do i=1,Ntorsions
                        An = cosine_torsions(1,i)
                        delta = cosine_torsions(2,i)
                        n = cosine_torsions(3,i)
                        E = E + An*(1. + cos((pi/180.d0)*(n*torsion_vals(i) - delta)))
                  enddo
            end if
      end function comp_torsions_energy

      real*8 function comp_impropers_energy() result(E)
            implicit none 
            real*8  :: An,n,delta
            integer :: i
            E = 0.
            if(pot_types(4)=="C") then
                  do i=1,Nimpropers
                        An = cosine_impropers(1,i)
                        delta = cosine_impropers(2,i)
                        n = cosine_impropers(3,i)
                        E = E + An*(1. + cos((pi/180.d0)*(n*improper_vals(i) - delta)))
                  enddo
            end if
      end function comp_impropers_energy

      real*8 function comp_energy(Natoms,xyz,recomp_flag) result (E)
            implicit none
            integer,intent(in) :: Natoms
            real*8,intent(in)  :: xyz(3,Natoms)
            logical,optional   :: recomp_flag
            if(present(recomp_flag))then
                  if(recomp_flag) then
                        call comp_all(Natoms,xyz)
                  end if
            else
                  call comp_all(Natoms,xyz)
            endif

            E = 0.
            E = E + comp_bonds_energy()
            E = E + comp_angles_energy()
            E = E + comp_torsions_energy()
            E = E + comp_impropers_energy()
      end function comp_energy

      function build_gradient(Natoms,xyz) result(G)
                  implicit none
                  integer :: Natoms
                  real*8 :: xyz(3,Natoms),xyz_mod(3,Natoms)
                  real*8 :: G(3,Natoms)
                  real*8,parameter :: hi = 5.d-6
                  real*8 :: Vp,Vm
                  integer :: a,p

                  G = 0.d0
                  do a=1,Natoms
                        do p=1,3
                              xyz_mod = xyz
                              xyz_mod(p,a) = xyz(p,a) + hi

                              Vp = comp_energy(Natoms,xyz_mod)

                              xyz_mod = xyz
                              xyz_mod(p,a) = xyz(p,a) - hi

                              Vm = comp_energy(Natoms,xyz_mod)

                              G(p,a) = (Vp-Vm)/(2.d0*hi)
                        end do

                  end do
                  call comp_all(Natoms,xyz)
      end function build_gradient

      function build_hessian(Natoms,xyz) result(H)
                  implicit none
                  integer :: Natoms
                  real*8 :: xyz(3,Natoms),xyz_mod(3,Natoms)
                  real*8 :: H(3*Natoms,3*Natoms)
                  real*8,parameter :: hi = 5.d-6,hj = 5.d-6
                  ! real*8,parameter :: hi = 1.d-5,hj = 1.d-5
                  real*8 :: V,Vpp,Vmm,Vpm,Vmp
                  integer :: a,b,p,q,i,j

                  H = 0.d0
                  do a=1,Natoms
                        do b=1,Natoms
                              do p=0,2
                              do q=0,2
                              i = 3*a - p
                              j = 3*b - q
                              
                              if(i/=j) then
                              xyz_mod = xyz
                              xyz_mod(3-p,a) = xyz(3-p,a) + hi
                              xyz_mod(3-q,b) = xyz(3-q,b) + hj
                              Vpp = comp_energy(Natoms,xyz_mod)


                              xyz_mod = xyz
                              xyz_mod(3-p,a) = xyz(3-p,a) - hi
                              xyz_mod(3-q,b) = xyz(3-q,b) - hj
                              Vmm = comp_energy(Natoms,xyz_mod)


                              xyz_mod = xyz
                              xyz_mod(3-p,a) = xyz(3-p,a) + hi
                              xyz_mod(3-q,b) = xyz(3-q,b) - hj
                              Vpm = comp_energy(Natoms,xyz_mod)


                              xyz_mod = xyz
                              xyz_mod(3-p,a) = xyz(3-p,a) - hi
                              xyz_mod(3-q,b) = xyz(3-q,b) + hj
                              Vmp = comp_energy(Natoms,xyz_mod)

                              H(i,j) = (Vpp - Vpm - Vmp + Vmm)/(4.d0*hi*hj)
                              else

                              xyz_mod = xyz
                              xyz_mod(3-p,a) = xyz(3-p,a) + hi
                              Vpp = comp_energy(Natoms,xyz_mod)


                              xyz_mod = xyz
                              xyz_mod(3-p,a) = xyz(3-p,a) - hi
                              Vmm = comp_energy(Natoms,xyz_mod)


                              xyz_mod = xyz
                              V = comp_energy(Natoms,xyz_mod)

                              H(i,j) = (Vpp - 2.d0*V + Vmm)/(hi*hj)
                              
                              end if
                              end do
                              end do
                        end do
                  end do
                  call comp_all(Natoms,xyz)
      end function build_hessian



end module force_field

