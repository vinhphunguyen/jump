# ----------------------------------------------------------------------
#
#                    ***       JUMP       ***
#                Material Point Method in Julia
#
# Copyright (2020) Vinh Phu Nguyen, phu.nguyen@monash.edu
# Civil Engineering, Monash University
# Clayton VIC 3800, Australia
# This software is distributed under the GNU General Public License.
#
# -----------------------------------------------------------------------

module Algorithm
	abstract type AlgorithmType end

	struct USL  <: AlgorithmType
		tolerance::Float64
	end

	struct FEM  <: AlgorithmType
		tolerance::Float64
	end

	struct TLFEM  <: AlgorithmType  
		tolerance::Float64
		alpha::Float64
    end

    struct TLFEM_MUSL  <: AlgorithmType  
		tolerance::Float64
    end

    struct TLFEMFull  <: AlgorithmType  
		tolerance::Float64
    end

	struct MUSL <: AlgorithmType
		alpha::Float64   # FLIP/PIC mixing param
	end
	
	struct APIC <: AlgorithmType end


    export AlgorithmType, USL, MUSL, APIC, FEM, TLFEM, TLFEMFull, TLFEM_MUSL
end
