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

push!(LOAD_PATH,"/Users/vingu/my-codes/julia-codes/juMP")
#
import PyPlot
using Printf
using LinearAlgebra
using StaticArrays   # if not yet installed, in REPL, do import Pkg and Pkd.add("StaticArrays")
using TimerOutputs
using Tensors
using BenchmarkTools

struct Solid3D_1
        F :: Vector{SMatrix{3,3,Float64,9}}  # F, 2x2 matrix
        ϵ :: Vector{SMatrix{3,3,Float64,9}}  # F, 2x2 matrix
        σ :: Vector{SMatrix{3,3,Float64,9}}  # F, 2x2 matrix
        parcount::Int64
        function Solid3D_1(count)
            strain    = fill(zeros(3,3),count)
            defor     = fill(zeros(3,3),count)
            stress    = fill(zeros(3,3),count)
            return new(defor,strain,stress,count)
        end
end

struct Solid3D_2
        F :: Vector{Tensor{2,3,Float64,9}}           # Tensor{order,dim,T<:Real} F, 2x2 matrix
        ϵ :: Vector{SymmetricTensor{2,3,Float64,6}}  # Tensor{order,dim,T<:Real} F, 2x2 matrix
        σ :: Vector{SymmetricTensor{2,3,Float64,6}}  # F, 2x2 matrix
        parcount::Int64
        function Solid3D_2(count)
            defor     = zeros(Tensor{2,3,Float64,9}, count) #fill(zero(Tensor{2, 3}),         count)
            strain    = zeros(SymmetricTensor{2,3,Float64,6},count)
            stress    = zeros(SymmetricTensor{2,3,Float64,6},count)
            return new(defor,strain,stress,count)
        end
end

function main_static_arrays(solid)
        defor = solid.F
        epsi  = solid.ϵ
        stres = solid.σ
        for ip = 1:solid.parcount
            defor[ip]     = SMatrix{3,3}(1,1,1,1,1,1,1,1,1)
            J             = det(defor[ip])
            Finv          = inv(defor[ip])    # no memory alloc
            P             = J*stres[ip]*Finv'  # convert to Piola Kirchoof stress, no memory alloc
            L             = Finv*Finv 
            D             = 0.5  * (L + L')  # no memory alloc
        end
end

function main_tensor(solid)
        defor = solid.F
        epsi  = solid.ϵ
        stres = solid.σ
        for ip = 1:solid.parcount
            defor[ip]     = ones(Tensor{2,3})
            J             = det(defor[ip])
            Finv          = inv(defor[ip])    # no memory alloc
            P             = J*stres[ip]⋅Finv'  # convert to Piola Kirchoof stress, no memory alloc
            L             = Finv⋅Finv 
            D             = symmetric(L)  # no memory alloc
        end
end

function main1()
    parCount = 100000
    s1       = Solid3D_1(parCount)
    @btime main_static_arrays($s1)
end

function main2()
    parCount = 100000
    s2       = Solid3D_2(parCount)
    @btime main_tensor($s2)
end

main1()
main2()
