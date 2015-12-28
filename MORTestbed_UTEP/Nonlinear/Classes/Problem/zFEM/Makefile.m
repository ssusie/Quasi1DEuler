function  [] = Makefile(flag)

if nargin < 1 || isempty(flag)
    flag = 2:7;
end

for i = flag
    switch i
        case 1
            mex zFEM.cpp array.hpp convertML2Cstruct.cpp matops.cpp zfem_forces.cpp zfem_elems.cpp zfem_mass.cpp zfem_precomp_rom.cpp
        case 2
            mex zFEM_mxMass.cpp array.hpp convertML2Cstruct.cpp matops.cpp zfem_elems.cpp zfem_mass.cpp
        case 3
            mex zFEM_mxMassVec.cpp array.hpp convertML2Cstruct.cpp matops.cpp zfem_elems.cpp zfem_mass.cpp
        case 4
            mex zFEM_mxIntForces.cpp array.hpp convertML2Cstruct.cpp matops.cpp zfem_elems.cpp zfem_mass.cpp zfem_forces.cpp
        case 5
            mex zFEM_mxStresses.cpp array.hpp convertML2Cstruct.cpp matops.cpp zfem_elems.cpp zfem_mass.cpp zfem_forces.cpp
        case 6
            mex zFEM_mxBodyForces.cpp array.hpp convertML2Cstruct.cpp matops.cpp zfem_elems.cpp zfem_mass.cpp zfem_forces.cpp
        case 7
            mex zFEM_mxPrecompRom.cpp array.hpp convertML2Cstruct.cpp matops.cpp zfem_elems.cpp zfem_precomp_rom.cpp
    end
end

