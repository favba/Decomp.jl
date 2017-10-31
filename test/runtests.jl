using Decomp
using Base.Test

a = zeros(3,3)
for i=1:100
  a[1,1] = rand(-2.:10e-8:2.)
  a[2,2] = rand(-2.:10e-8:2.)
  a[3,3] = rand(-2.:10e-8:2.)
  a[1,2] = rand(-2.:10e-8:2.)
  a[1,3] = rand(-2.:10e-8:2.)
  a[2,3] = rand(-2.:10e-8:2.)
  a[2,1] = a[1,2]
  a[3,1] = a[1,3]
  a[3,2] = a[2,3]
  eigv,eigvec1,eigvec2,eigvec3 = eigen(a)
  eigvc,eigvec123 = eig(a)
  @test eigv ≈ sort!(eigvc,by=abs,rev=true)
  @test eigv[1] * eigvec1*eigvec1' .+ eigv[2] * eigvec2*eigvec2' .+ eigv[3] * eigvec3*eigvec3' ≈ a
end
