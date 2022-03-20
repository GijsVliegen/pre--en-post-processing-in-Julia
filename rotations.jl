using LinearAlgebra
using Einsum

P = 1

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

struct rotation
    #eigenlijk gewoon een quaternion
    ω #change to omega
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

function normalize(rot ::rotation)
    abs_val = sqrt(rot.ω^2 + rot.i^2 + rot.j^2 + rot.k^2)
    #return multiply(rot, abs_val)
    return rotation(rot.ω/abs_val, rot.i/abs_val, rot.j/abs_val, rot.k/abs_val)
end

function getComponents(rot ::rotation)
    return [rot.ω, rot.i, rot.j, rot.k]
end

function Base.:copy(rot ::rotation)
    return rotation(rot.ω, rot.i, rot.j, rot.k)

end
#in tegenstelling tot in damask wordt er niet rekening gehouden met NaN values
#dit doen we omdat het niet is ingebouwd in isapprox()

function isClose(first ::rotation, other ::rotation, rtol, atol, nanEquals = true)
    fcom = getComponents(first)
    scom = getComponents(other)
    for i in 1:4
        if !isapprox(fcom[i], scom[i], rtol = rtol , atol = atol, nans = nanEquals)
            return false
        end
    end
    return true
end

function isClose(first ::Array{rotation, 1}, other ::Array{rotation, 1}, rtol, atol, nanEquals = true)
    println("dit is slechts een test")
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
#input: [phi1, PHI, phi2]
function eu2qu(rot) ::rotation
    ϕ1 = rot[1]
    Φ = rot[2]
    ϕ2 = rot[3]
    σ = (ϕ1 + ϕ2)/2.0
    δ = (ϕ1 - ϕ2)/2.0
    c = cos(Φ/2)
    s = sin(Φ/2)
    q₀ = c*cos(σ)
    if q₀ > 0
        return rotation(q₀, -P*s*cos(δ), -P*s*sin(δ), -P*c*sin(σ))
    else
        return rotation(-q₀, P*s*cos(δ), P*s*sin(δ), P*c*sin(σ))
    end
end

function om2qu(α) ::rotation
    q₀ = 1/2f0*sqrt(1 + α[1, 1] + α[2, 2] + α[3,3])
    q₁ = P/2f0*sqrt(1 + α[1, 1] - α[2, 2] - α[3,3])
    q₂ = P/2f0*sqrt(1 - α[1, 1] + α[2, 2] - α[3,3])
    q₃ = P/2f0*sqrt(1 - α[1, 1] - α[2, 2] + α[3,3])
    if (α[3,2] < α[2, 3])
        q₁ = -q₁
    end
    if (α[1, 3] < α[3, 1])
        q₂ = -q₂
    end
    if (α[2, 1] < α[1, 2])
        q₃ = -q₃
    end
    return normalize(rotation(q₀, q₁, q₂, q₃))
end

#input = [n_0, n_1, n_2, w]
function ax2qu(rot) ::rotation
    n̂ = rot[1:3]
    ω = rot[4]
    n = n̂*sin(ω/2)
    return rotation(cos(ω/2), n[1], n[2], n[3])
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

function ro2qu(rotation)
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

γ = [+1.0000000000018852,      -0.5000000002194847,
-0.024999992127593126,    -0.003928701544781374,
-0.0008152701535450438,   -0.0002009500426119712,
-0.00002397986776071756,  -0.00008202868926605841,
+0.00012448715042090092,  -0.0001749114214822577,
+0.0001703481934140054,   -0.00012062065004116828,
+0.000059719705868660826, -0.00001980756723965647,
+0.000003953714684212874, -0.00000036555001439719544]
#[n_0, n_1, n_2]
function ho2ax(rotation)
    small_h = norm(rotation)
    small_h²  = small_h^2
    sh²_copy = 1
    if small_h == 0
        return [0, 0, 1, 0]
    end
    h′ = rotation/small_h
    s = 0
    for i in 1:16
        s += γ[i]*sh²_copy
        sh²_copy *= small_h²
    end
    return vcat(h′, [2*acos(s)])
end

function ho2qu(rotation) #::rotation
    ax = ho2ax(rotation)
    return ax2qu(ax) #als ik rotation laat staan, komt er een error lol
end

function qu2eu(rot ::rotation)
    #moet unit quaternion zijn
    q₀ = rot.ω
    q₁ = rot.i
    q₂ = rot.j
    q₃ = rot.k
    q₀₃ = q₀^2 + q₃^2
    q₁₂ = q₁^2 + q₂^2
    χ = sqrt(q₀₃ * q₁₂)
    if (χ == 0 && q₁₂ == 0)
        return [atan(-2*P*q₀*q₃, q₀^2 - q₃^2), 0, 0]
    elseif (χ == 0 && q₀₃ == 0)
        return [atan(2*q₁*q₂, q₁^2 - q₂^2), π, 0]
    else #als χ != 0
        return [atan((q₁*q₃ - P*q₀*q₂)/χ,(-P*q₀*q₁ - q₂*q₃)/χ),
            atan(2*χ, q₀₃-q₁₂), atan((P*q₀*q₂ + q₁*q₃)/χ, (q₂*q₃ - P*q₀*q₁)/χ)]
    end
end

function qu2om(rot ::rotation)
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

function qu2ax(rot ::rotation)
    q₀ = rot.ω
    q₁ = rot.i
    q₂ = rot.j
    q₃ = rot.k
    ω = 2*acos(q₀)
    if (ω == 0)
        #[n_0, n_1, n_2, w]
        return [0, 0, 1, 0]
    end
    if (q₀ == 0)
        return [q₁, q₂, q₃, π]
    end
    s = sign(q₀)/sqrt(q₁^2 + q₂^2 + q₃^2)
    return [s*q₁, s*q₂, s*q₃, ω]
end

function qu2ro(rotation ::rotation)
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

function qu2ho(rotation ::rotation)
    q₀ = rotation.ω
    q₁ = rotation.i
    q₂ = rotation.j
    q₃ = rotation.k
    ω = 2*acos(q₀)
    if ω == 0
        return [0, 0, 0]
    end
    s = 1/sqrt(q₁^2 + q₂^2 + q₃^2)
    n̂ = [s*q₁, s*q₂, s*q₃]
    f = 3(ω-sin(ω))/4
    return n̂*f^(1/3)
end

#moeten deze in positive real hemisphere zijn? en wat is dan de logica hierachter?
#voorlopig wordt er nergens gecheckt of het unit quaternionen zijn
#misschien hiervoor een voorwaarde schrijven in de struct zelf?

#in de python damask is het shape (..., 4), wij doen (4,...), zou dit veel uitmaken?
"""
    from_quaternion(array) ::Array{rotation}

Initialize from quaternion.

# Arguments
   - `array ::Array`: shape (4, ...)

        Unit quaternion (q₀, q₁, q₂, q₃) in positive real hemisphere, i.e. ǀqǀ = 1, q_0 ≥ 0.

# Examples
```julia-repl
julia> from_quaternion(reshape([(1:16)...], 4, 2, 2))
2×2 Array{rotation,2}:
 rotation(1, 2, 3, 4)  rotation(9, 10, 11, 12)
 rotation(5, 6, 7, 8)  rotation(13, 14, 15, 16)
```
"""
function from_quaternion(array) #greek letters?
    sizes = size(array) #is een tupel, in de vorm van (4, ...)
    flat_array = vec(array)
    rotations = rotation[]
    if (length(flat_array) == 4) #als er maar een rotatie is, geen array teruggeven
        return rotation(flat_array[1], flat_array[2], flat_array[3], flat_array[4])
    end
    for i in 1:4:length(flat_array)
        push!(rotations, rotation(flat_array[i], flat_array[i+1], flat_array[i+2], flat_array[i+3]))
    end
    return reshape(rotations, sizes[2:length(sizes)])
end

#deze voorwaarden worden opnieuw nergens gechecked, misschien dit checken in de functie eu2qu?
"""
    from_eulerAngle(array, degrees = false) ::Array{rotation}

Initialize from Bungle Euler angles.

# Arguments
   - `array ::Array`: shape (3, ...)

        Euler angles (φ1 ∈ [0,2π], ϕ ∈ [0,π], φ2 ∈ [0,2π]) or (φ1 ∈ [0,360], ϕ ∈ [0,180], φ2 ∈ [0,360]) if degrees == True.

   - `degrees ::bool`, optional

        Euler angles are given in degrees. Defaults to False.
"""
function from_euler_angle(array, degrees = false)
    sizes = size(array)
    flat_array = vec(array)
    rotations = rotation[]
    if degrees
        flat_array = flat_array/(2*π)
    end
    if (length(flat_array) == 3) #als er maar een rotatie is, geen array teruggeven
        return eu2qu(flat_array)
    end
    for i in 1:3:length(flat_array)
        q = eu2qu(flat_array[i:i+2])
        push!(rotations, q)
    end
    return reshape(rotations, sizes[2:length(sizes)])
end



#deze voorwaarden worden opnieuw nergens gechecked, misschien dit checken in de functie ax2qu?
"""
    from_axisAngle(array, degrees = false) ::Array{rotation}

Initialize from Axis angle pair.

# Arguments
   - `array ::Array`: shape (4, ...)

        Axis and angle (n_1, n_2, n_3, ω) with ǀnǀ = 1 and ω ∈ [0,π] or ω ∈ [0,180] if degrees == True.

   - `degrees ::bool`, optional

        Euler angles are given in degrees. Defaults to False.
"""
function from_axis_angle(array, degrees = false)
    sizes = size(array)
    flat_array = vec(array)
    rotations = rotation[]
    if degrees
        flat_array = flat_array/(2*π)
    end
    if (length(flat_array) == 4) #als er maar een rotatie is, geen array teruggeven
        return ax2qu(flat_array)
    end
    for i in 1:4:length(flat_array)
        q = ax2qu(flat_array[i:i+3])
        push!(rotations, q)
    end
    return reshape(rotations, sizes[2:length(sizes)])
end

#deze voorwaarden worden nergens gechecked, misschien dit checken in de functie om2qu?
"""
    from_basis(array, degrees = false) ::Array{rotation}

Initialize from lattice basis vectors.

# Arguments
   - `array ::Array`: shape (3, 3, ...)

        Three three-dimensional lattice basis vectors.

   - `orthonormal ::bool`, optional

        Basis is strictly orthonormal, i.e. is free of stretch components. Defaults to True.

   - `reciprocal ::bool`, optional

        Basis vectors are given in reciprocal (instead of real) space. Defaults to False.
"""
function from_basis(basis, orthonormal = true, reciprocal = false)
    """# TODO check dimensions
    om = copy(basis)
    if reciprocal
        om = LinearAlgebra.inv!(LinearAlgebra.transpose!(om)/π)
        orthonormal = false
    end
    if !orthonormal
        svd = LinearAlgebra.svd!(om)
        @einsum om[i,j] = svd.U[i,j]*svd.Vt[j,l]
    end"""
end


#deze voorwaarden worden nergens gechecked, misschien dit checken in de functie om2qu?
"""
    from_basis(array, degrees = false) ::Array{rotation}

Initialize from rotation matrix.

# Arguments
   - `array ::Array`: shape (3, 3, ...)

        Rotation matrix with det(R) = 1, R.T ∙ R = I.
"""

function from_matrix(array, degrees = false)
    sizes = size(array)
    flat_array = vec(array)
    rotations = rotation[]
    if degrees
        flat_array = flat_array/(2*π)
    end
    if (length(flat_array) == 9) #als er maar een rotatie is, geen array teruggeven
        return om2qu(array)
    end
    for i in 1:9:length(flat_array)
        #zou deze reshape veel tijd in beslag nemen? anders gewoon een om2qu maken die op een vector werkt?
        #of zou ge kunnen flattenen naar een 2dimensionale array ipv naar een vector
        q = om2qu(reshape(flat_array[i:i+8], (3, 3)))
        push!(rotations, q)
    end
    return reshape(rotations, sizes[3:length(sizes)])
end

function from_rodriguesfrank(array ::Array{<:Number, 1})
    return ro2qu(array)
end

function from_rodriguesfrank(array ::Array{<:Number, 2})
    sizes = size(array)
    flat_array = vec(array)
    rotations = rotation[]
    if (length(flat_array) == 3) #als er maar een rotatie is, geen array teruggeven
        return ro2qu(flat_array)
    end
    for i in 1:4:length(flat_array)
        q = ro2qu(flat_array[i:i+3])
        push!(rotations, q)
    end
    return reshape(rotations, sizes[2:length(sizes)])
end

function from_homochoric(array)
    sizes = size(array)
    flat_array = vec(array)
    rotations = rotation[]
    if (length(flat_array) == 3) #als er maar een rotatie is, geen array teruggeven
        return ho2qu(flat_array)
    end
    for i in 1:3:length(flat_array)
        q = ho2qu(flat_array[i:i+2])
        push!(rotations, q)
    end
    return reshape(rotations, sizes[2:length(sizes)])
end

# misschien voor later
function from_cubochoric(array)
end

#Rotate a vector
function apply(rot ::rotation, vector ::Array{<:Number, 1})
    qvector = rotation(vector[1], vector[2], vector[3], 0)
    result = rot*qvector*(inv(rot))
    return [result.i, result.j, result.k]
end

#Rotate a matrix
function apply(rot ::rotation, vector ::Array{<:Number, 2})
    R = as_matrix(rot)
    return R*vector
end

function apply(rot ::rotation, vector ::Array{<:Number, 4})
    R = as_matrix(rot)
    @einsum result[i,j,k,l] := R[i,m]*R[j,n]*R[k,o]*R[l,p]*vector[m,n,o,p]
    return result
end

function Base.:*(r ::rotation, nr ::Number) ::rotation
    return rotation(r.ω * nr, r.i * nr, r.j * nr, r.k * nr)
end

function Base.:*(nr ::Number, r ::rotation) ::rotation
    return r*nr
end

#niet commutatief
function Base.:*(f ::rotation, s ::rotation) ::rotation
    angle = f.ω*s.ω - (f.i * s.i + f.j * s.j + f.k * s.k)
    i = (f.j * s.k - f.k * s.j) + f.ω * s.i + s.ω * f.i
    j = (f.k * s.i - f.i * s.k) + f.ω * s.j + s.ω * f.j
    k = (f.i * s.j - f.j * s.i) + f.ω * s.k + s.ω * f.k
    return rotation(angle, i, j, k)
end

function Base.:inv(rot ::rotation)
    return rotation(rot.ω, -rot.i, -rot.j, -rot.k)
end

function from_random()
    return from_random(1)[1]
end
"""
    from_random(n, sizes = n, type) ::Array{rotation}

    type = Float32, Float64...

Construeert een array van random rotaties, met gegeven dimensies
"""
function from_random(n, sizes = n, type = Float64) ::Array{rotation}
    if (n != tupleSize(sizes))
        print("invalid dimensions")
    end
    rotations = rotation[]
    for i = 1:n
        u1 = rand(type)
        u2 = rand(type)
        u3 = rand(type)
        h = rotation(sqrt(1-u1)*sin(π*u2*2), sqrt(1-u1)*cos(π*u2*2), sqrt(u1)*sin(π*u3*2), sqrt(u1)*cos(π*u3*2))
        if (h.ω < 0)
            h = h*-1
        end
        push!(rotations, copy(h))
    end
    return reshape(rotations, sizes)
end

function as(qu2x::Function, xSize ::Tuple, rotations ::Array{rotation})
    l = tupleSize(xSize)
    sizes = size(rotations)
    flat_array = vec(rotations)
    result = Array{Float64}(undef, tupleSize(sizes)*l)
    for i in 1:length(flat_array)
        x = qu2x(flat_array[i])
        result[i*l-l+1:i*l] = x
    end
    return reshape(result, (xSize..., sizes...))
end

function as_axis_angle(rotation ::rotation)
    return as_axis_angle([rotation])
end

function as_axis_angle(rotations ::Array{rotation}, degrees = false, pair = false)
    return as(qu2ax, (4,), rotations)
end

function as_matrix(rotation ::rotation) #nodig als er maar een element is
    return as_matrix([rotation])
end

function as_matrix(rotations ::Array{rotation})
    return as(qu2om, (3,3), rotations)
end

function as_euler_angle(rotation ::rotation)
    return as_euler_angle([rotation])
end

function as_euler_angle(rotations ::Array{rotation})
    return as(qu2eu, (3,), rotations)
end

function as_homochoric(rotation ::rotation)
    return as_homochoric([rotation])
end

function as_homochoric(rotations ::Array{rotation})
    return as(qu2ho, (3,), rotations)
end

function as_rodriguesfrank(rotations ::rotation)
    return as_rodriguesfrank([rotations])
end

function as_rodriguesfrank(rotations ::Array{rotation})
    return as(qu2ro, (4, ), rotations)
end

function as_rodriguesfrank(rotation ::rotation)
    return qu2ro(rotation)
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
