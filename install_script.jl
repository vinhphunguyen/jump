#!/usr/bin/env julia
import Pkg

Pkg.add(["OhMyREPL",
         "Revise",
         "LinearAlgebra",
         "Gmsh", 
         "Glob", 
         "Images",
	     "TimerOutputs", 
         "StaticArrays",
         "WriteVTK",
         "PyPlot",
         "BenchmarkTools",
         "DelimitedFiles",
         "RecursiveArrayTools"])
