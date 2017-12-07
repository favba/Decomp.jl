__precompile__()
module Decomp

using StaticArrays: MVector, SVector
export eigen!, eigen, inphdecomp, propdecomp, inphmodel, propmodel, mixmodel

 function eigen!(e11::T,e22::T,e33::T,e12::T,e13::T,e23::T,eig::MVector{3,T},eigv1::MVector{3,T},eigv2::MVector{3,T},eigv3::MVector{3,T}) where T<:AbstractFloat

  p1 = e12^2 + e13^2 + e23^2
  if p1 == 0.0
    # E is diagonal.
    # eig[1] = e11
    # eig[2] = e22
    # eig[3] = e33

    if abs(e11) >= abs(e22)
      if abs(e11) >= abs(e33)
        eig[1]=e11

        eigv1[1] = oneunit(T)
        eigv1[2] = zero(T)
        eigv1[3] = zero(T)

        if abs(e22) >= abs(e33)
          eig[2] = e22
          eig[3] = e33

          eigv2[1] = zero(T)
          eigv2[2] = oneunit(T)
          eigv2[3] = zero(T)
        else
          eig[3] = e22
          eig[2] = e33

          eigv2[1] = zero(T)
          eigv2[2] = zero(T)
          eigv2[3] = oneunit(T)
        end
      else
        eig[1] = e33
        eig[2] = e11
        eig[3] = e22

        eigv1[1] = zero(T)
        eigv1[2] = zero(T)
        eigv1[3] = oneunit(T)

        eigv2[1] = oneunit(T)
        eigv2[2] = zero(T)
        eigv2[3] = zero(T)
      end
    else
      if abs(e22) >= abs(e33)
        eig[1] = e22

        eigv1[1] = zero(T)
        eigv1[2] = oneunit(T)
        eigv1[3] = zero(T)

      elseif abs(e33) >= abs(e11)
        eig[2] = e33
        eig[3] = e11

        eigv2[1] = zero(T)
        eigv2[2] = zero(T)
        eigv2[3] = oneunit(T)

      else
        eig[3] = e33
        eig[2] = e11

        eigv2[1] = oneunit(T)
        eigv2[2] = zero(T)
        eigv2[3] = zero(T)
      end
    end

    eigv3[1] = eigv1[2]*eigv2[3] - (eigv1[3]*eigv2[2])
    eigv3[2] = eigv1[3]*eigv2[1] - (eigv1[1]*eigv2[3])
    eigv3[3] = eigv1[1]*eigv2[2] - (eigv1[2]*eigv2[1])
  else
    q = (e11 + e22 + e33)/3
    p2 = (e11-q)^2 + (e22-q)^2 + (e33-q)^2 + 2*p1
    p = sqrt(p2/6)
    r = ((e11-q)*(e22-q)*(e33-q) - (e11-q)*(e23^2) - (e12^2)*(e33-q) + 2*(e12*e13*e23) - (e13^2)*(e22-q))/(2*p*p*p)

    # In exact arithmetic for a symmetric matrix  -1 <= r <= 1
    # but computation error can leave it slightly outside this range.
    if r <= -1
      ϕ = T(π/3)
    elseif r >= 1
      ϕ = zero(T)
    else
      ϕ = acos(r)/3
    end

    # the eigenvalues satisfy eig[3] <= eig[2] <= eig[1]
    eig[1] = q + 2*p*cos(ϕ)
    eig[3] = q + 2*p*cos(ϕ+(2*π/3))
    eig[2] = 3*q - eig[1] - eig[3]     # since trace(E) = eig[1] + eig[2] + eig[3] = 3q

    sort!(eig,by=abs,rev=true)

    bla = ((e22 - eig[1])*(e33 - eig[1]) - e23*e23)
    if bla != 0
      eigv1[1] = oneunit(T)
      eigv1[2] = (e23*e13 - (e33-eig[1])*e12)/bla
      eigv1[3] = (-e13 -e23*eigv1[2])/(e33-eig[1])
      aux = sqrt(1 + eigv1[2]^2 + eigv1[3]^2)
      eigv1[1] = 1/aux
      eigv1[2] = eigv1[2]/aux
      eigv1[3] = eigv1[3]/aux
    else
      bla = ((e11 - eig[1])*(e22 - eig[1]) - e12*e12)
      eigv1[3] = oneunit(T)
      eigv1[1] = (e23*e12 - (e22-eig[1])*e13)/bla
      eigv1[2] = (-e23 -e12*eigv1[1])/(e22-eig[1])
      aux = sqrt(1 + eigv1[2]^2 + eigv1[1]^2)
      eigv1[1] = eigv1[1]/aux
      eigv1[2] = eigv1[2]/aux
      eigv1[3] = 1/aux
    end
    bla = ((e22 - eig[2])*(e33 - eig[2]) - e23*e23)
    if bla != 0
      eigv2[1] = oneunit(T)
      eigv2[2] = (e23*e13 - (e33-eig[2])*e12)/bla
      eigv2[3] = (-e13 -e23*eigv2[2])/(e33-eig[2])
      aux = sqrt(1 + eigv2[2]^2 + eigv2[3]^2)
      eigv2[1] = 1/aux
      eigv2[2] = eigv2[2]/aux
      eigv2[3] = eigv2[3]/aux
    else
      bla = ((e11 - eig[2])*(e22 - eig[2]) - e12*e12)
      eigv2[3] = oneunit(T)
      eigv2[1] = (e23*e12 - (e22-eig[2])*e13)/bla
      eigv2[2] = (-e23 -e12*eigv2[1])/(e22-eig[2])
      aux = sqrt(1 + eigv2[2]^2 + eigv2[1]^2)
      eigv2[1] = eigv2[1]/aux
      eigv2[2] = eigv2[2]/aux
      eigv2[3] = 1.0/aux
    end

    eigv3[1] = eigv1[2]*eigv2[3] - (eigv1[3]*eigv2[2])
    eigv3[2] = eigv1[3]*eigv2[1] - (eigv1[1]*eigv2[3])
    eigv3[3] = eigv1[1]*eigv2[2] - (eigv1[2]*eigv2[1])

  end
  return 0
end


function eigen(e11::T,e22::T,e33::T,e12::T,e13::T,e23::T) where T<:AbstractFloat
  eig = MVector{3,T}()
  eigv1 = MVector{3,T}()
  eigv2 = MVector{3,T}()
  eigv3 = MVector{3,T}()
  eigen!(e11,e22,e33,e12,e13,e23,eig,eigv1,eigv2,eigv3)
  return SVector(eig), SVector(eigv1), SVector(eigv2), SVector(eigv3)
end

function eigen(e11::Real,e22::Real,e33::Real,e12::Real,e13::Real,e23::Real)
  return eigen(Float64(e11),Float64(e22),Float64(e33),Float64(e12),Float64(e13),Float64(e23))
end

function eigen(A::AbstractArray{<:Real,2})
  size(A) == (3,3) || throw(ArgumentError("Input must be a 3x3 Matrix"))
  (A[1,2] == A[2,1] && A[1,3] == A[3,1] && A[2,3] == A[3,2]) || throw(ArgumentError("Input must be symmetric matrix"))
  return eigen(A[1,1],A[2,2],A[3,3],A[1,2],A[1,3],A[2,3])
end

 function inphdecomp!(t11::T,t22::T,t33::T,t12::T,t13::T,t23::T,
                e11::T,e22::T,e33::T,e12::T,e13::T,e23::T,
                eig::MVector{3,T},eigv1::MVector{3,T},eigv2::MVector{3,T},eigv3::MVector{3,T}) where T<:AbstractFloat

  eigen!(e11,e22,e33,e12,e13,e23,eig, eigv1, eigv2, eigv3)

  m1 = t11*eigv1[1]*eigv1[1] + t22*eigv1[2]*eigv1[2] + t33*eigv1[3]*eigv1[3] + 2*(
       t12*eigv1[1]*eigv1[2] + t13*eigv1[1]*eigv1[3] + t23*eigv1[2]*eigv1[3])
  m2 = t11*eigv2[1]*eigv2[1] + t22*eigv2[2]*eigv2[2] + t33*eigv2[3]*eigv2[3] + 2*(
       t12*eigv2[1]*eigv2[2] + t13*eigv2[1]*eigv2[3] + t23*eigv2[2]*eigv2[3])
  m3 = t11 + t22 + t33 - m1 - m2

  lam2 = eig[1]^2 + eig[2]^2 + eig[3]^2

  a = 1 - (3*(eig[1]^2)/lam2)
  b = 1 - (3*(eig[2]^2)/lam2)
  c = b*eig[1] - a*eig[2]

  alpha1 = (b*m1 - a*m2)/c
  alpha0 = (m1 - eig[1]*alpha1)/a
  alpha2 = -3*alpha0/lam2

  tm11 = m1*eigv1[1]*eigv1[1] + m2*eigv2[1]*eigv2[1] + m3*eigv3[1]*eigv3[1]
  tm22 = m1*eigv1[2]*eigv1[2] + m2*eigv2[2]*eigv2[2] + m3*eigv3[2]*eigv3[2]
  tm33 = m1*eigv1[3]*eigv1[3] + m2*eigv2[3]*eigv2[3] + m3*eigv3[3]*eigv3[3]
  tm12 = m1*eigv1[1]*eigv1[2] + m2*eigv2[1]*eigv2[2] + m3*eigv3[1]*eigv3[2]
  tm13 = m1*eigv1[1]*eigv1[3] + m2*eigv2[1]*eigv2[3] + m3*eigv3[1]*eigv3[3]
  tm23 = m1*eigv1[2]*eigv1[3] + m2*eigv2[2]*eigv2[3] + m3*eigv3[2]*eigv3[3]

  return tm11, tm22, tm33, tm12, tm13, tm23, alpha0, alpha1, alpha2
end

function inphdecomp(t11::T,t22::T,t33::T,t12::T,t13::T,t23::T,
                e11::T,e22::T,e33::T,e12::T,e13::T,e23::T) where T<:AbstractFloat

  eig = MVector{3,T}()
  eigv1 = MVector{3,T}()
  eigv2 = MVector{3,T}()
  eigv3 = MVector{3,T}()

  return inphdecomp!(t11,t22,t33,t12,t13,t23,e11,e22,e33,e12,e13,e23,eig,eigv1,eigv2,eigv3)
end

inphdecomp(var...) = inphdecomp((Float64(v) for v in var)...)

 function propdecomp(t11::T,t22::T,t33::T,t12::T,t13::T,t23::T,
                    e11::T,e22::T,e33::T,e12::T,e13::T,e23::T) where T<:AbstractFloat

  mode2 = (e11^2 + 2*e12^2 + 2*e13^2 + e22^2 + 2*e23^2 + e33^2)
  alphat = (e11*t11 + 2*e12*t12 + 2*e13*t13 + e22*t22 + 2*e23*t23 + e33*t33)/mode2

  m12 = alphat*e12
  m13 = alphat*e13
  m23 = alphat*e23
  m11 = alphat*e11
  m22 = alphat*e22
  m33 = alphat*e33

  return m11, m22, m33, m12, m13, m23, alphat
end

propdecomp(var...) = propdecomp((Float64(v) for v in var)...)

 function inphmodel(t11::Vector{T},t22::Vector{T},t33::Vector{T},t12::Vector{T},t13::Vector{T},t23::Vector{T}, #input
                   e11::Vector{T},e22::Vector{T},e33::Vector{T},e12::Vector{T},e13::Vector{T},e23::Vector{T}, #input
                   tm11::Vector{T},tm22::Vector{T},tm33::Vector{T},tm12::Vector{T},tm13::Vector{T},tm23::Vector{T}, #output
                   alpha0::Vector{T},alpha1::Vector{T},alpha2::Vector{T},ratio::Vector{T}) where T<:Real

  eigg  = [MVector{3,T}() for i in 1:Threads.nthreads()]
  eigv1 = [MVector{3,T}() for i in 1:Threads.nthreads()]
  eigv2 = [MVector{3,T}() for i in 1:Threads.nthreads()]
  eigv3 = [MVector{3,T}() for i in 1:Threads.nthreads()]

  Threads.@threads for i in 1:length(e11)
    j = Threads.threadid()
    @inbounds tm11[i], tm22[i], tm33[i], tm12[i], tm13[i], tm23[i], alpha0[i], alpha1[i], alpha2[i] = inphdecomp!(t11[i],t22[i],t33[i],t12[i],t13[i],t23[i],
                            e11[i],e22[i],e33[i],e12[i],e13[i],e23[i],
                            eigg[j],eigv1[j],eigv2[j],eigv3[j])

    @inbounds ratio[i] = (tm11[i]^2 + 2*tm12[i]^2 + 2*tm13[i]^2 + tm22[i]^2 + 2*tm23[i]^2 + tm33[i]^2)/(t11[i]^2 + 2*t12[i]^2 + 2*t13[i]^2 + t22[i]^2 + 2*t23[i]^2 + t33[i]^2)
  end
  return 0
end

 function propmodel(t11::A,t22::A,t33::A,t12::A,t13::A,t23::A, #input
                   e11::A,e22::A,e33::A,e12::A,e13::A,e23::A, #input
                   tm11::A,tm22::A,tm33::A,tm12::A,tm13::A,tm23::A, #output
                   alpha::A,ratio::A) where {A<:Vector{<:Real}}

  Threads.@threads for i in 1:length(e11)
    @inbounds tm11[i], tm22[i], tm33[i], tm12[i], tm13[i], tm23[i], alpha[i] = propdecomp(t11[i],t22[i],t33[i],t12[i],t13[i],t23[i],e11[i],e22[i],e33[i],e12[i],e13[i],e23[i])

    @inbounds ratio[i] = (tm11[i]^2 + 2*tm12[i]^2 + 2*tm13[i]^2 + tm22[i]^2 + 2*tm23[i]^2 + tm33[i]^2)/(t11[i]^2 + 2*t12[i]^2 + 2*t13[i]^2 + t22[i]^2 + 2*t23[i]^2 + t33[i]^2)
  end
  return 0
end

 function mixmodel(t11::A,t22::A,t33::A,t12::A,t13::A,t23::A, #input
                  e11::A,e22::A,e33::A,e12::A,e13::A,e23::A, #input
                  c11::A,c22::A,c33::A,c12::A,c13::A,c23::A, #input
                  tm11::A,tm22::A,tm33::A,tm12::A,tm13::A,tm23::A, #output
                  ratio::A) where {A<:Vector{<:Real}}

  Threads.@threads for i in 1:length(e11)
    @inbounds tm11[i] = e11[i] + c11[i]
    @inbounds tm22[i] = e22[i] + c22[i]
    @inbounds tm33[i] = e33[i] + c33[i]
    @inbounds tm12[i] = e12[i] + c12[i]
    @inbounds tm13[i] = e13[i] + c13[i]
    @inbounds tm23[i] = e23[i] + c23[i]

    @inbounds ratio[i] = (tm11[i]^2 + 2*tm12[i]^2 + 2*tm13[i]^2 + tm22[i]^2 + 2*tm23[i]^2 + tm33[i]^2)/(t11[i]^2 + 2*t12[i]^2 + 2*t13[i]^2 + t22[i]^2 + 2*t23[i]^2 + t33[i]^2)
  end
  return 0
end

end # module
