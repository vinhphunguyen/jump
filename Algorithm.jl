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
