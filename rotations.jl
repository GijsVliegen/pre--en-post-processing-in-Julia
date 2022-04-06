using LinearAlgebra
using Einsum

P = -1

"""TODO
    -voorwaarden checken in sommige omzetfuncties
    -documentatie schrijven
    -from_basis() en as_basis() schrijven
    -extra argumenten toevoegen zoals degrees, pair...
    -efficiëntie vergelijken met python-implementatie

    ambiguiteit bij quaternion: als angle/omega 0 is,
"""

#om documentatie te schrijven
#https://julia-doc.readthedocs.io/en/latest/manual/documentation/

#python damask documentatie
#https://damask.mpie.de/documentation/processing_tools/pre-processing.html#damask.Rotation.from_quaternion

struct Rotation
    #eigenlijk gewoon een quaternion
    ω
    i
    j
    k
end

function tupleSize(tuple)
    total = 1
    for i in tuple
        total *= i
    end
    return total
end

function normalize(rot ::Rotation)
    abs_val = sqrt(rot.ω^2 + rot.i^2 + rot.j^2 + rot.k^2)
    #return multiply(rot, abs_val)
    return Rotation(rot.ω/abs_val, rot.i/abs_val, rot.j/abs_val, rot.k/abs_val)
end

function getComponents(rot ::Rotation)
    return [rot.ω, rot.i, rot.j, rot.k]
end

function Base.:copy(rot ::Rotation)
    return Rotation(rot.ω, rot.i, rot.j, rot.k)

end
#in tegenstelling tot in damask wordt er niet rekening gehouden met NaN values
#dit doen we omdat het niet is ingebouwd in isapprox()

function isClose(first ::Rotation, other ::Rotation, rtol, atol, nanEquals = true)
    fcom = getComponents(first)
    scom = getComponents(other)
    for i in 1:4
        if !isapprox(fcom[i], scom[i], rtol = rtol , atol = atol, nans = nanEquals)
            return false
        end
    end
    return true
end

function isClose(first ::Array{Rotation, 1}, other ::Array{Rotation, 1}, rtol, atol, nanEquals = true)
    if length(first) != length(other)
        return false
    end
    result = true
    for i in 1:length(first)
        if !isClose(first[i], other[i], rtol, atol, nanEquals)
            result = false
            println("oke het is fout: de twee quaternionen zijn:")
            println(first[i])
            println(other[i])
            break
        end
    end
    return result
end

#wat volgt zijn algoritmes voor het omzetten van rotaties naar een bepaalde representatie
#enkel de directe omzettingen zijn gegeven, andere omzettingen kunnen opgebouwd worden
#volgens de tabel pagina 11 van de paper "Consistent representations of and conversions
#between 3D rotations"

#voor de gebruiksvriendelijkheid kan het idee van bovenstaande code ook gebruikt worden.

#Euler angle heeft ZXZ structuur

function as_euler_angle!(rotations ::Array{Rotation})
    return as!(qu2eu!, (3,), rotations)
end

function fill!(result, low, high)
    result[low:high] = [(low:high)...]
end

function testInPlace()
    result = zeros(50)
    fill!(result, 2, 48)
    return result
end

function as!(qu2x::Function, xSize ::Tuple, rotations ::Array{Rotation})
    l = tupleSize(xSize)
    sizes = size(rotations)
    flat_array = vec(rotations)
    result = Array{Float64}(undef, tupleSize(sizes)*l)
    for i in 1:length(flat_array)
        qu2x(result, i*l-l+1, i*l, flat_array[i])
    end
    return reshape(result, (xSize..., sizes...))
end

#input = [phi1, PHI, phi2]
function qu2eu!(result, low, high, rot ::Rotation)
    #moet unit quaternion zijn
    q₀ = rot.ω
    q₁ = rot.i
    q₂ = rot.j
    q₃ = rot.k
    q₀₃ = q₀^2 + q₃^2
    q₁₂ = q₁^2 + q₂^2
    χ = sqrt(q₀₃ * q₁₂)
    if (χ == 0 && q₁₂ == 0)
        result[low:high] = [atan(-2*P*q₀*q₃, q₀^2 - q₃^2), 0.0, 0.0]
    elseif (χ == 0 && q₀₃ == 0)
        result[low:high] = [atan(2*q₁*q₂, q₁^2 - q₂^2), π, 0.0]
    else #als χ != 0
        result[low:high] = [atan((q₁*q₃ - P*q₀*q₂)/χ,(-P*q₀*q₁ - q₂*q₃)/χ),
            atan(2*χ, q₀₃-q₁₂), atan((P*q₀*q₂ + q₁*q₃)/χ, (q₂*q₃ - P*q₀*q₁)/χ)]
    end
end

#input = [phi1, PHI, phi2]
function eu2qu(rot, P = -1) ::Rotation
    ϕ1 = rot[1]
    Φ = rot[2]
    ϕ2 = rot[3]
    σ = (ϕ1 + ϕ2)/2.0
    δ = (ϕ1 - ϕ2)/2.0
    c = cos(Φ/2)
    s = sin(Φ/2)
    q₀ = c*cos(σ)
    if q₀ > 0
        return Rotation(q₀, -P*s*cos(δ), -P*s*sin(δ), -P*c*sin(σ))
    else
        return Rotation(-q₀, P*s*cos(δ), P*s*sin(δ), P*c*sin(σ))
    end
end

#input = matrix
function om2qu(α, P = -1) ::Rotation
    q₀ = 1/2f0*sqrt(1 + α[1] + α[5] + α[9])
    q₁ = P/2f0*sqrt(1 + α[1] - α[5] - α[9])
    q₂ = P/2f0*sqrt(1 - α[1] + α[5] - α[9])
    q₃ = P/2f0*sqrt(1 - α[1] - α[5] + α[9])
    if (α[6] < α[8])
        q₁ = -q₁
    end
    if (α[7] < α[3])
        q₂ = -q₂
    end
    if (α[2] < α[4])
        q₃ = -q₃
    end
    return normalize(Rotation(q₀, q₁, q₂, q₃))
end

#input = [n_0, n_1, n_2, w]
function ax2qu(rot, P = -1) ::Rotation
    n̂ = rot[1:3]
    ω = rot[4]
    n = n̂*sin(ω/2)
    return Rotation(cos(ω/2), n[1], n[2], n[3])
end

#input = [n_0, n_1, n_2, tan(ω/2)]
#om de een of andere rede klopt dit niet, god mag weten waarom
function ro2ax(rotation)
    n = rotation[1:3]
    ρ = norm(n)
    #als ik de paper gebruik zou het 2*atan(rho) moeten zijn,
    #dit werkt echter niet, 2*atan(rotation[4]) werkt wel.
    return vcat(n/ρ, [2*atan(rotation[4])])
end

function ro2qu(rotation, P = -1)
    ax = ro2ax(rotation)
    return ax2qu(ax)
end

#deze functie is er enkel voor debugging
function ax2ro(rotation)
    #wat als omega = π???
    n̂ = rotation[1:3]
    f = tan(rotation[4]/2)
    return vcat(n̂, [f])
end

#input = [n_0, n_1, n_2]
function ho2ax(rotation)
    γ = [+1.0000000000018852,      -0.5000000002194847,
    -0.024999992127593126,    -0.003928701544781374,
    -0.0008152701535450438,   -0.0002009500426119712,
    -0.00002397986776071756,  -0.00008202868926605841,
    +0.00012448715042090092,  -0.0001749114214822577,
    +0.0001703481934140054,   -0.00012062065004116828,
    +0.000059719705868660826, -0.00001980756723965647,
    +0.000003953714684212874, -0.00000036555001439719544]
    small_h = norm(rotation)
    small_h²  = small_h^2
    sh²_copy = 1
    if small_h == 0
        return [0.0, 0.0, 1.0, 0.0]
    end
    h′ = rotation/small_h
    s = 0
    for i in 1:16
        s += γ[i]*sh²_copy
        sh²_copy *= small_h²
    end
    return vcat(h′, [2*acos(s)])
end

#input = 
function ho2qu(rotation, P = -1) #::rotation
    ax = ho2ax(rotation)
    return ax2qu(ax) #als ik rotation laat staan, komt er een error lol
end

function qu2eu(rot ::Rotation)
    #moet unit quaternion zijn
    q₀ = rot.ω
    q₁ = rot.i
    q₂ = rot.j
    q₃ = rot.k
    q₀₃ = q₀^2 + q₃^2
    q₁₂ = q₁^2 + q₂^2
    χ = sqrt(q₀₃ * q₁₂)
    if (χ == 0 && q₁₂ == 0)
        return [atan(-2*P*q₀*q₃, q₀^2 - q₃^2), 0.0, 0.0]
    elseif (χ == 0 && q₀₃ == 0)
        return [atan(2*q₁*q₂, q₁^2 - q₂^2), π, 0.0]
    else #als χ != 0
        return [atan((q₁*q₃ - P*q₀*q₂)/χ,(-P*q₀*q₁ - q₂*q₃)/χ),
            atan(2*χ, q₀₃-q₁₂), atan((P*q₀*q₂ + q₁*q₃)/χ, (q₂*q₃ - P*q₀*q₁)/χ)]
    end
end

function qu2om(rot ::Rotation)
    #moet unit quaternion zijn
    #geeft passieve interpretatie
    q₀ = rot.ω
    q₁ = rot.i
    q₂ = rot.j
    q₃ = rot.k
    q̄ = q₀^2 - (q₁^2 + q₂^2 + q₃^2)
    α = [q̄+2*q₁^2 2*(q₁*q₂-P*q₀*q₃) 2*(q₁*q₃+P*q₀*q₂);
        2(q₁*q₂+P*q₀*q₃) q̄+2*q₂^2 2*(q₂*q₃-P*q₀*q₁);
        2*(q₁*q₃-P*q₀*q₂) 2*(q₂*q₃+P*q₀*q₁) q̄+2*q₃^2]
    return α
end

function qu2ax(rot ::Rotation)
    q₀ = rot.ω
    q₁ = rot.i
    q₂ = rot.j
    q₃ = rot.k
    ω = 2*acos(q₀)
    if (ω == 0)
        #[n_0, n_1, n_2, w]
        return [0.0, 0.0, 1.0, 0.0]
    end
    if (q₀ == 0)
        return [q₁, q₂, q₃, π]
    end
    s = sign(q₀)/sqrt(q₁^2 + q₂^2 + q₃^2)
    return [s*q₁, s*q₂, s*q₃, ω]
end

function qu2ro(rotation ::Rotation)
    q₀ = rotation.ω
    q₁ = rotation.i
    q₂ = rotation.j
    q₃ = rotation.k
    s = sqrt(q₁^2 + q₂^2 + q₃^2)
    t = tan(acos(q₀))
    #oppassen als s klein wordt!!!
    #rodriguesFrank opslaan als vector van 4 elementen
    #[n_0, n_1, n_2, t]
    return [q₁/s, q₂/s, q₃/s, t]
end

function qu2ho(rotation ::Rotation)
    q₀ = rotation.ω
    q₁ = rotation.i
    q₂ = rotation.j
    q₃ = rotation.k
    ω = 2*acos(q₀)
    if ω == 0
        return [0.0, 0.0, 0.0]
    end
    s = 1/sqrt(q₁^2 + q₂^2 + q₃^2)
    n̂ = [s*q₁, s*q₂, s*q₃]
    f = 3(ω-sin(ω))/4
    return n̂*f^(1/3)
end

#omega, i, j, k, omega > 0
function qu2rot(rot, P) ::Rotation
    return Rotation(rot[1], rot[2]*-P, rot[3]*-P, rot[4]*-P)
end

function rot2qu(rot ::Rotation)
    return [rot.ω, rot.i, rot.j, rot.k]
end

#Rotate a vector
"""
    apply(rot, vector) -> vector

Rotate vector, second order tensor, or fourth order tensor.

# Arguments

   -`rot ::Rotation`

   -`array ::Array{<:Number}`: shape (3), (3,3), or (3,3,3,3)

        Vector or tensor on which to apply the rotation.

# Returns

    rotated vector or 2nd/4th order tensor, shape (3), (3,3), or (3,3,3,3)

        Rotated vector or tensor, i.e. transformed to frame defined by rotation.

"""
function apply(rot ::Rotation, vector ::Array{<:Number, 1})
    qvector = Rotation(vector[1], vector[2], vector[3], 0)
    result = rot*qvector*(inv(rot))
    return [result.i, result.j, result.k]
end

#Rotate a matrix
function apply(rot ::Rotation, vector ::Array{<:Number, 2})
    R = as_matrix(rot)
    return R*vector
end

# Rotate a fourth order tensor
function apply(rot ::Rotation, vector ::Array{<:Number, 4})
    R = as_matrix(rot)
    result = zeros(3,3,3,3)
    @einsum result[i,j,k,l] += R[i,m]*R[j,n]*R[k,o]*R[l,p]*vector[m,n,o,p]
    return result
end

function Base.:*(r ::Rotation, nr ::Number) ::Rotation
    return Rotation(r.ω * nr, r.i * nr, r.j * nr, r.k * nr)
end

function Base.:*(nr ::Number, r ::Rotation) ::Rotation
    return r*nr
end

#niet commutatief
function Base.:*(f ::Rotation, s ::Rotation) ::Rotation
    ω = f.ω*s.ω - (f.i * s.i + f.j * s.j + f.k * s.k)
    i = (f.j * s.k - f.k * s.j) + f.ω * s.i + s.ω * f.i
    j = (f.k * s.i - f.i * s.k) + f.ω * s.j + s.ω * f.j
    k = (f.i * s.j - f.j * s.i) + f.ω * s.k + s.ω * f.k
    return Rotation(ω, i, j, k)
end

function Base.:inv(rot ::Rotation)
    return Rotation(rot.ω, -rot.i, -rot.j, -rot.k)
end


function from_random()
    return from_random(1)[1]
end
"""
    from_random(n, sizes = n) -> Array{Rotation}

Construeert een array van random rotaties, met gegeven dimensies
"""
function from_random(n, sizes = n) ::Array{Rotation}
    if (n != tupleSize(sizes))
        print("invalid dimensions")
    end
    rotations = Rotation[]
    for i = 1:n
        u1 = rand()
        u2 = rand()
        u3 = rand()
        h = Rotation(sqrt(1-u1)*sin(pi*u2*2), sqrt(1-u1)*cos(pi*u2*2), sqrt(u1)*sin(pi*u3*2), sqrt(u1)*cos(pi*u3*2))
        if (h.ω < 0)
            h = h*-1
        end
        push!(rotations, copy(h))
    end
    return reshape(rotations, sizes)
end


#algemene from functie
#x2qu is een v.d. functies: qu2qu, eu2qu, ax2qu, om2qu, ho2qu of ro2qu
#nrOfDimension is aantal dimensies waarin resultaat wordt voorgesteld (1, behalve bij om is het 2)
#l is aantal elementen per resultaat (3, 4 of 9)
function from(x2qu::Function, nrOfDimension ::Int, l ::Int, array, P = -1)
    sizes = size(array) #is een tupel, in de vorm van (4, ...)
    result = Rotation[]
    for i in 1:l:length(array)
        qu = x2qu(array[i:i+l-1], P)
        push!(result, qu)
    end
    reshape(result, sizes[nrOfDimension+1:length(sizes)])
end

function from_quaternion(array ::Array{<:Number, 1}, accept_homomorph = false, P = -1)
    if accept_homomorph
        for i in 1:4:length(array)
            if array[i] < 0
                array[i:i+3] = -array[i:i+3]
            end
        end
    end
    return qu2qu(array, P)
end

#moeten deze in positive real hemisphere zijn? en wat is dan de logica hierachter?
#voorlopig wordt er nergens gecheckt of het unit quaternionen zijn
#misschien hiervoor een voorwaarde schrijven in de struct zelf?

#in de python damask is het shape (..., 4), wij doen (4,...), zou dit veel uitmaken?
"""
    from_quaternion(array, accept_homomorph = false, P = -1) -> Array{Rotation}

Initialize from quaternion.

# Arguments
   - `array ::Array`: shape (4,…)

        Unit quaternion (q₀, q₁, q₂, q₃) in positive real hemisphere, i.e. ǀqǀ = 1, q_0 ≥ 0.

   - `accept_homomorph ::Bool`

        Allow homomorphic variants, i.e. q_0 < 0 (negative real hemisphere). Defaults to false.

   - `P ::Int`, -1 or 1, optional

        Sign convention. Defaults to -1.

# Examples
```julia-repl
julia> from_quaternion(reshape([(1:16)...], 4, 2, 2))
2×2 Array{Rotation,2}:
 Rotation(1, 2, 3, 4)  Rotation(9, 10, 11, 12)
 Rotation(5, 6, 7, 8)  Rotation(13, 14, 15, 16)
```
"""
function from_quaternion(array ::Array{<:Number}, accept_homomorph = false, P = - 1)
    if accept_homomorph
        for i in 1:4:length(array)
            if array[i] < 0
                array[i:i+3] = array[i:i+3]*-1
            end
        end
    end
    return from(qu2qu, 1, 4, array, P)
end

function from_Euler_angles(array ::Array{<:Number, 1}, degrees = false)
    if degrees
        return eu2qu(array/180*pi, -1)
    end
    return eu2qu(array, -1)
end

#deze voorwaarden worden opnieuw nergens gechecked, misschien dit checken in de functie eu2qu?
"""
    from_Euler_angles(array, degrees = false) -> Array{Rotation}

Initialize from Bungle Euler angles.

# Arguments
   - `array ::Array`: shape (3,…)

        Euler angles (φ1 ∈ [0,2π], ϕ ∈ [0,π], φ2 ∈ [0,2π]) or (φ1 ∈ [0,360], ϕ ∈ [0,180], φ2 ∈ [0,360]) if degrees == true.

   - `degrees ::Bool`, optional

        Euler angles are given in degrees. Defaults to false.
"""
function from_Euler_angles(array ::Array{<:Number}, degrees = false) ::Array{Rotation}
    if degrees
        array = array/180*pi
    end
    return from(eu2qu, 1, 3, array)
end

function from_axis_angle(array ::Array{<:Number, 1}, degrees = false, normalize = false, P = -1) ::Rotation
    if degrees
        array[4] = array[4]/180*pi
    end
    if normalize
        LinearAlgebra.normalize!(array[1:3])
    end
    return ax2qu(array)
end

#deze voorwaarden worden opnieuw nergens gechecked, misschien dit checken in de functie ax2qu?
"""
    from_axis_angle(array, degrees = false, normalize = false, P = -1) -> Array{Rotation}

Initialize from Axis angle pair.

# Arguments
   - `array ::Array{<:Number}`: shape (4,…)

        Axis and angle (n_1, n_2, n_3, ω) with ǀnǀ = 1 and ω ∈ [0,π] or ω ∈ [0,180] if degrees == true.

   - `degrees ::Bool`, optional

        Angle ω is given in degrees. Defaults to false.

   - `normalize ::Bool`, optional

        Allow ǀnǀ ≠ 1. Defaults to false.
    
   - `P ::Int`, -1 or 1, optional

        Sign convention. Defaults to -1.
"""

function from_axis_angle(array ::Array{<:Number}, degrees ::Bool = false, normalize ::Bool = false, P ::Int = -1) ::Array{Rotation}
    if degrees
        for i in 1:4:length(array)
            array[i+3] = array[i+3]/180*pi
        end
    end
    if normalize
        for i in 1:4:length(array)
            array[i:i+2] = LinearAlgebra.normalize(array[i:i+2])
        end
    end
    return from(ax2qu, 1, 4, array)
end

#deze voorwaarden worden nergens gechecked, misschien dit checken in de functie om2qu?
"""
    from_basis(basis, orthonormal = true, reciprocal = false) -> Array{Rotation}

Initialize from lattice basis vectors.

# Arguments
   - `basis ::Array{<:Number}`: shape (3, 3,…)

        Three three-dimensional lattice basis vectors.

   - `orthonormal ::Bool`, optional

        Basis is strictly orthonormal, i.e. is free of stretch components. Defaults to true.

   - `reciprocal ::Bool`, optional

        Basis vectors are given in reciprocal (instead of real) space. Defaults to false.
"""
function from_basis(basis ::Array{<:Number}, orthonormal ::Bool = true, reciprocal ::Bool = false) ::Array{rotation}
    # TODO check dimensions
    om = copy(basis)
    if reciprocal
        om = LinearAlgebra.inv!(LinearAlgebra.transpose!(om)/π)
        orthonormal = false
    end
    if !orthonormal
        svd = LinearAlgebra.svd!(om)
        om = svd.U*svd.Vt
    end
    # TODO check determinant == 1
    # TODO check orthogonality
    return from_matrix(om)
end

"""
    from_parallel(a,b) -> Array{Rotation}

Initialize from pairs of two orthogonal lattice basis vectors.

# Arguments

    - `a ::Array`: shape (3, 2,…)

        Two three-dimensional lattice vectors of first orthogonal basis.

    - `b ::Array`: shape (3, 2,…)

        Corresponding three-dimensional lattice vectors of second basis.
"""
function from_parallel(a, b) ::Array{rotation}

end


function from_matrix(array ::Array{<:Number, 2}) ::Rotation
    return om2qu(array)
end

#deze voorwaarden worden nergens gechecked, misschien dit checken in de functie om2qu?
"""
    from_matrix(array ::Array{<:Number}) -> Array{Rotation}

Initialize from rotation matrix.

# Arguments
   - `array ::Array{<:Number}`: shape (3, 3,…)

        Rotation matrix with det(R) = 1, R.T ∙ R = I.
"""
function from_matrix(array ::Array{<:Number})
    return from(om2qu, 2, 9, array)
end

function from_Rodrigues_vector(array ::Array{<:Number, 1}, normalize ::Bool = False, P ::Int = -1)
    if normalize
        LinearAlgebra.normalize!(array[1:3])
    end
    return ro2qu(array, P)
end

"""
    from_Rodrigues_vector(array, normalize = false, P = -1) -> Array{Rotation}

Initialize from Rodrigues–Frank vector (angle separated from axis).

# Arguments
   - `array ::Array{<:Number}`: shape (4,…)

        Rodrigues–Frank vector (n_1, n_2, n_3, tan(ω/2)) with ǀnǀ = 1 and ω ∈ [0,π].
   
   - `normalize ::Bool`, optional

        Allow ǀnǀ ≠ 1. Defaults to false.
   
   - `P ::Int`, -1 or 1, optional

        Sign convention. Defaults to -1.

"""
function from_Rodrigues_vector(array ::Array{<:Number}, normalize ::Bool = false, P ::Int = -1)
    if normalize
        for i in 1:4:length(array)
            array[i:i+2] = LinearAlgebra.normalize(array[i:i+2])
        end
    end
    from(ro2qu, 1, 4, array, P)
end

function from_homochoric(array ::Array{<:Number, 1}, P ::Int = -1)
    return ho2qu(array, P)
end

"""
    from_homochoric(array, P = -1) -> Array{Rotation}

Initialize from homochoric vector.

# Arguments
   - `array ::Array{<:Number}`: shape (3,…)

        Homochoric vector (h_1, h_2, h_3) with ǀhǀ < (3/4*π)^(1/3).
   
   - `P ::Int`, -1 or 1, optional

        Sign convention. Defaults to -1.

"""
function from_homochoric(array ::Array{<:Number}, P ::Int = -1)
    from(ho2qu, 1, 3, array, P)
end

function as(qu2x::Function, xSize ::Tuple, rotations ::Array{Rotation})
    l = tupleSize(xSize)
    sizes = size(rotations)
    result = Array{Float64}(undef, tupleSize(sizes)*l)
    for i in 1:length(rotations)
        x = qu2x(rotations[i])
        result[i*l-l+1:i*l] = x
    end
    return reshape(result, (xSize..., sizes...))
end

function as_quaternion(rotation ::Rotation)
    return rot2qu(rotation)
end

"""
    as_quaternion(rotations)

Represent as unit quaternion.

# Arguments
   - `rotations ::Array{Rotation}`: shape(…)
    
        Rotations that need to be converted

# Returns
   - `q ::Array{<:Number}`: shape (4,…)

        Unit quaternion (q_0, q_1, q_2, q_3) in positive real hemisphere, i.e. ǀqǀ = 1, q_0 ≥ 0.
"""

function as_quaternion(rotations ::Array{Rotation})
    return as(rot2qu, (4,), rotations)
end

function as_Euler_angles(rotation ::Rotation, degrees = false)#nodig als er maar een element is
    return qu2eu(rotation)
end

"""
    as_Euler_angles(rotations, degrees = false)

Represent as Bunge Euler angles.

# Arguments
   - `rotations ::Array{Rotation}`: shape(…)
    
        Rotations that need to be converted

   - `degrees ::Bool`, optional

        Return angles in degrees. Defaults to false.

# Returns
   - `phi ::Array{<:Number}`: shape (3,…)

        Bunge Euler angles (φ_1 ∈ [0,2π], ϕ ∈ [0,π], φ_2 ∈ [0,2π]) or (φ_1 ∈ [0,360], ϕ ∈ [0,180], φ_2 ∈ [0,360]) if degrees == true.

"""
function as_Euler_angles(rotations ::Array{Rotation}, degrees = false)
    if degrees
        return as(qu2eu, (3,), rotations)*180/pi
    end
    return as(qu2eu, (3,), rotations)
end

function as_axis_angle(rotation ::Rotation, degrees = false, pair = false)#nodig als er maar een element is
    result = qu2ax(rotation)
    if degrees
        result = result*180/pi
    end
    if pair
        return (result[1:3], result[4])
    end
    return result
end

"""
    as_axis_angle(rotations, degrees = false, pair = false)

Represent as axis–angle pair.

# Arguments
   - `rotations ::Array{Rotation}`: shape(…)
    
        Rotations that need to be converted

   - `degrees ::Bool`, optional

        Return rotation angle in degrees. Defaults to false.
    
   - `pair ::Bool`
        Return tuple of axis and angle. Defaults to false.

# Returns
   - `axis_angle ::Array{<:Number}`: shape (4,…) or tuple ((3,…), (…)) if pair == true

   Axis and angle [n_1, n_2, n_3, ω] with ǀnǀ = 1 and ω ∈ [0,π] or ω ∈ [0,180] if degrees == true.
"""
function as_axis_angle(rotations ::Array{Rotation}, degrees = false, pair = false)
    result = as(qu2ax, (4,), rotations)
    if degrees
        for i in 1:4:length(result)
            result[i+3] = result[i+3]*180/pi
        end
    end
    if pair
        return (result[1:3, :], result[4, :])
    end
    return result
end

function as_matrix(rotation ::Rotation) #nodig als er maar een element is
    return qu2om(rotation)
end

"""
    as_matrix(rotations)

Represent as rotation matrix.

# Arguments
   - `rotations ::Array{Rotation}`: shape(…)
    
        Rotations that need to be converted

# Returns
   - `R ::Array{<:Number}`: shape (3, 3,…)

        Rotation matrix R with det(R) = 1, R.T ∙ R = I.
"""
function as_matrix(rotations ::Array{Rotation})
    return as(qu2om, (3,3), rotations)
end

function as_homochoric(rotation ::Rotation)#nodig als er maar een element is
    return qu2ho(rotation)
end

"""
    as_homochoric(rotations)

Represent as homochoric vector.

# Arguments
   - `rotations ::Array{Rotation}`: shape(…)
    
        Rotations that need to be converted

   - `rotations ::Array{Rotations}`

        Rotations

# Returns
   - `h ::Array{<:Number}`: shape (3,…)

        Homochoric vector (h_1, h_2, h_3) with ǀhǀ < (3/4*π)^(1/3).
"""
function as_homochoric(rotations ::Array{Rotation})
    return as(qu2ho, (3,), rotations)
end

function as_Rodrigues_vector(rotation ::Rotation, compact = false)#nodig als er maar een element is
    if compact
        result = qu2ro(rotation)
        return result[1:3]*result[4]
    end
    return qu2ro(rotation)
end

"""
    as_Rodrigues_vector(rotations, compact = false)

Represent as Rodrigues–Frank vector with separate axis and angle argument.

# Arguments
   - `rotations ::Array{Rotation}`: shape(…)
    
        Rotations that need to be converted

   - `compact ::Bool`, optional

        Return three-component Rodrigues–Frank vector, i.e. axis and angle argument are not separated. Defaults to false

# Returns
   - `rho ::Array{<:Number}`: shape (4,…) or (3,…) if compact == true

        Rodrigues–Frank vector [n_1, n_2, n_3, tan(ω/2)] with ǀnǀ = 1 and ω ∈ [0,π] or [n_1, n_2, n_3] with ǀnǀ = tan(ω/2) and ω ∈ [0,π] if compact == true.
"""
function as_Rodrigues_vector(rotations ::Array{Rotation}, compact = false)
    if compact
        result = as(qu2ro, (4, ), rotations)
        sizes = size(result)
        new_result = zeros(((3,)..., sizes[2: end]...))
        for i in 1:length(rotations)
            new_result[i*3-2:i*3] = result[i*4-3:i*4-1]*result[i*4]
        end
        return new_result
    end
    return as(qu2ro, (4, ), rotations)
end

"""
function eu2om(rotation ::eulerAngle)
    c1 = cos(rotation.phi1)
    s1 = sin(rotation.phi1)
    c2 = cos(rotation.phi2)
    s2 = sin(rotation.phi2)
    C = cos(rotation.PHI)
    S = sin(rotation.PHI)
    matrix = [c1*c2-s1*C*s2 s1*c2+c1*C*s2 S*s2; -c1*s2-s1*C*c2 -s1*s2+c1*C*c2 S*c2; s1*S -c1*S C]
    return rotationMatrix(matrix)
end

function eu2ax(rotation ::eulerAngle)::axisAnglePair
    t = tan(rotation.PHI/2)
    sigma = (rotation.phi1 + rotation.phi2)/2.0
    delta = (rotation.phi1 - rotation.phi2)/2.0
    tau = sqrt(t^2 + sin(sigma)^2)
    alpha = 2*atan(tau/(cos(sigma)))
    if (alpha > π)
        alpha = 2*π - alpha
    end
    return axisAnglePair([P/t*cos(delta) P/t*sin(delta) P/t*sin(sigma)], a)
end

function eu2ro(rotation ::eulerAngle)
    #EA naar AA en dan AA naar RF

    #nog eens nakijken
    axisAngle = toAxisAngle(rotation)
    return rodriguesFrank(axisAngle.n, tan(axisAngle.w/2))
end

function om2eu(rotation ::rotationMatrix)
    a_33 = rotation.matrix[3, 3]
    if abs(a_33) != 1
        zeta = 1/sqrt(1-a_33^2)
        phi1 = atan(rotation.matrix[3, 1]*zeta, -rotation.matrix[3, 2]*zeta)
        PHI = acos(a_33)
        phi2 = atan(rotation.matrix[1, 3]*zeta, rotation.matrix[2, 3]*zeta)
        return eulerAngle(phi1, PHI, phi2)
    else
        phi1 = atan(rotation.matrix[1, 2], rotation.matrix[1, 1])
        PHI = π/2*(1 - a_33)
        return eulerAngle(phi1, PHI, 0)
    end
end

function om2ax(rotation ::rotationMatrix)
    omega = acos((tr(rotation.matrix) - 1)/2)
    #inspiratie voor de code: https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/index.htm
    #dit algo komt niet uit de paper, die onduidelijk was, opletten voor omega = 0 of omega = π
    #omega = 0 -> x y en z zijn arbitrair
    #omega = π -> x y en z maken uit dus moeten berekend worden
    #handelen van singularities: https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/david.htm
end

function ax2om(rotation ::axisAnglePair)
    n = rotation.n
    c = cos(rotation.omega)
    s = sin(rotation.omega)
    a = [c+(1-c)*n[1]^2 (1-c)*n[1]*n[2]+s*n[3] (1-c)*n[1]*n[3]-s*n[2];
    (1-c)*n[1]*n[2]-s*n[3] c+(1-c)*n[2]^2 (1-c)*n[2]*n[3]+s*n[1];
    (1-c)*n[1]*n[3]+s*n[2] (1-c)*n[2]*n[3]-sn[1] c+(1-c)*n[3]^2]
    if (P == 1)
        a = transpose(a)
    end
    return rotationMatrix(a)
end

function ax2ho(rotation ::axisAnglePair)
    f = (3/4(rotation.omega - sin(rotation.omega)))^(1/3)
    return homochoric(rotation.n*f)
end


function ro2ho(rotation ::rodriguesFrank)
    rho = abs(rotation.rho)
    if (rho == 0)
        return homochoric([0, 0, 0])
    end
    #is het nodig om rho == oneindig te behandelen
    w = 2*atan(rho)
    f = 3(w - sin(w))/4
    return homochoric(rotation.n*f^(1/3))
end


function ho2cu()
end

function cu2ho()
end"""

#yeet
