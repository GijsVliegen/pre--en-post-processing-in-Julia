using LinearAlgebra
using Einsum

P = 1

"""TODO
    -voorwaarden checken in sommige omzetfuncties
    -documentatie schrijven
    -tests schrijven
    -from_basis() en as_basis() schrijven
    -extra argumenten toevoegen zoals degrees, pair...
    -efficiëntie vergelijken met python-implementatie

    row major and column major

    ambiguiteit bij quaternion: als angle/omega 0 is, 
"""

#om documentatie te schrijven
#https://julia-doc.readthedocs.io/en/latest/manual/documentation/

#python damask documentatie
#https://damask.mpie.de/documentation/processing_tools/pre-processing.html#damask.Rotation.from_quaternion

struct rotation
    #eigenlijk gewoon een quaternion
    angle #change to omega
    i
    j
    k
end

function normalize(rot ::rotation)
    abs_val = sqrt(rot.angle^2 + rot.i^2 + rot.j^2 + rot.k^2)
    #return multiply(rot, abs_val)
    return rotation(rot.i/abs_val, rot.j/abs_val, rot.k/abs_val, rot.angle/abs_val)
end

function getComponents(rot ::rotation)
    return [rot.angle, rot.i, rot.j, rot.k]
end

function Base.:copy(rot ::rotation)
    return rotation(rot.angle, rot.i, rot.j, rot.k)

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

#wat volgt zijn algoritmes voor het omzetten van rotaties naar een bepaalde representatie
#enkel de directe omzettingen zijn gegeven, andere omzettingen kunnen opgebouwd worden
#volgens de tabel pagina 11 van de paper "Consistent representations of and conversions
#between 3D rotations"

#voor de gebruiksvriendelijkheid kan het idee van bovenstaande code ook gebruikt worden.

#Euler angle heeft ZXZ structuur
#input: [phi1, PHI, phi2]
function eu2qu(rot) ::rotation
    phi1 = rot[1]
    PHI = rot[2]
    phi2 = rot[3]
    sigma = (phi1 + phi2)/2.0
    delta = (phi1 - phi2)/2.0
    c = cos(PHI/2)
    s = sin(PHI/2)
    q0 = c*cos(sigma)
    if q0 < 0
        return rotation(q0, P*s*cos(delta), P*s*sin(delta), P*c*sin(sigma))
    else
        return rotation(q0, -P*s*cos(delta), -P*s*sin(delta), -P*c*sin(sigma))
    end
end

function om2qu(a) ::rotation
    if (1 + a[1] + a[5] + a[9] < 0 || 1 + a[1] - a[5] - a[9] < 0
        || 1 - a[1] + a[5] - a[9] < 0 || 1 - a[1] - a[5] + a[9] < 0)
        print(a)
        return rotation(1, 1, 0, 0)
    end

    q0 = 1/2f0*sqrt(1 + a[1] + a[5] + a[9])
    q1 = P/2f0*sqrt(1 + a[1] - a[5] - a[9])
    q2 = P/2f0*sqrt(1 - a[1] + a[5] - a[9])
    q3 = P/2f0*sqrt(1 - a[1] - a[5] + a[9])
    if (a[6] < a[8])
        q1 = -q1
    end
    if (a[7] < a[3])
        q2 = -q2
    end
    if (a[2] < a[4])
        q3 = -q3
    end
    #weet niet waarom deze volgorde werkt, maar het werkt
    return normalize(rotation(q3, q0, q1, q2))
end

#input = [n_0, n_1, n_2, w]
function ax2qu(rot) ::rotation
    n = rot[1:3]
    n = n*sin(rot[4]/2)
    return rotation(cos(rot[4]/2), n[1], n[2], n[3])
end

#input = [n_0, n_1, n_2, tan(ω/2)]
#om de een of andere rede klopt dit niet, god mag weten waarom
function ro2ax(rotation)
    n = rotation[1:3]
    rho = norm(n)
    #als ik de paper gebruik zou het 2*atan(rho) moeten zijn,
    #dit werkt echter niet, 2*atan(rotation[4]) werkt wel.
    return vcat(n/rho, [2*atan(rotation[4])])
end

function ro2qu(rotation)
    ax = ro2ax(rotation)
    return ax2qu(ax)
end

#deze functie is er enkel voor debugging
function ax2ro(rotation)
    #wat als omega = pi???
    n = rotation[1:3]
    f = tan(rotation[4]/2)
    return vcat(n, [f])
end

gamma = [1.0000000000018805, -0.500000000219485, -0.024999992127593, -0.003928701544781,
    -0.000815270153545, -0.000200950042612, -0.000023979867761, -0.000082028689266,
    0.000124487150421, -0.000174911421482, 0.000170348193414, -0.000120620650041,
    0.000059719705869, -0.000019807567240, 0.000003953714684, -0.000000365550014]
#[n_0, n_1, n_2]
function ho2ax(rotation)
    small_h = norm(rotation)
    small_h_squared = small_h^2
    if small_h == 0
        return [0, 0, 1, 0]
    end
    h_prime = rotation/small_h
    s = 0
    for i in 1:16
        s += gamma[i]*small_h_squared^(i-1)
    end
    return vcat(h_prime, [2*acos(s)])
end

function ho2qu(rotation) #::rotation
    ax = ho2ax(rotation)
    return ax2qu(ax) #als ik rotation laat staan, komt er een error lol
end

function qu2eu(rotation ::rotation)
    #moet unit quaternion zijn
    q0 = rotation.angle
    q1 = rotation.i
    q2 = rotation.j
    q3 = rotation.k
    q03 = q0^2 + q3^2
    q12 = q1^2 + q2^2
    x = sqrt(q03 * q12)
    if (x == 0 && q12 == 0)
        return [atan(-2*P*q0*q3, q0^2 - q3^2), 0, 0]
    elseif (x == 0 && q03 == 0)
        return [atan(2*q1*q2, q1^2 - q2^2), pi, 0]
    else #als x != 0
        return [atan((q1*q3 - P*q0*q2)/x,(-P*q0*q1 - q2*q3)/x),
            atan(2*x, q03-q12), atan((P*q0*q2 + q1*q3)/x, (q2*q3 - P*q0*q1)/x)]
    end
end

function qu2om(rot ::rotation)
    #moet unit quaternion zijn
    #geeft passieve interpretatie
    q0 = convert(Float64, rot.angle)
    q1 = convert(Float64, rot.i)
    q2 = convert(Float64, rot.j)
    q3 = convert(Float64, rot.k)
    q = q0^2 - (q1^2 + q2^2 + q3^2)
    a = [q+2*q1^2, 2*(q1*q2-P*q0*q3), 2*(q1*q3+P*q0*q2),
        2(q1*q2+P*q0*q3), q+2*q2^2, 2*(q2*q3-P*q0*q1),
        2*(q1*q3-P*q0*q2), 2*(q2*q3+P*q0*q1), q+2*q3^2]
    return a
end

function qu2ax(rotation ::rotation)
    q0 = rotation.angle
    q1 = rotation.i
    q2 = rotation.j
    q3 = rotation.k
    w = 2*acos(q0)
    if (w == 0)
        #[n_0, n_1, n_2, w]
        return [0, 0, 1, 0]
    end
    if (q0 == 0)
        return [q1, q2, q3, pi]
    end
    s = sign(q0)/sqrt(q1^2 + q2^2 + q3^2)
    return [s*q1, s*q2, s*q3, w]
end

function qu2ro(rotation ::rotation)
    q0 = rotation.angle
    q1 = rotation.i
    q2 = rotation.j
    q3 = rotation.k
    s = sqrt(q1^2 + q2^2 + q3^2)
    t = tan(acos(q0))
    #oppassen als s klein wordt!!!
    #rodriguesFrank opslaan als vector van 4 elementen
    #[n_0, n_1, n_2, t]
    return [q1/s, q2/s, q3/s, t]
end

function qu2ho(rotation ::rotation)
    q0 = rotation.angle
    q1 = rotation.i
    q2 = rotation.j
    q3 = rotation.k
    w = 2*acos(q0)
    if w == 0
        return [0, 0, 0]
    end
    s = 1/sqrt(q1^2 + q2^2 + q3^2)
    n = [s*q1, s*q2, s*q3]
    f = 3(w-sin(w))/4
    return n*f^(1/3)
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

        Unit quaternion (q0, q1, q2, q3) in positive real hemisphere, i.e. ǀqǀ = 1, q_0 ≥ 0.

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
function from_Euler_angles(array, degrees = false)
    sizes = size(array)
    flat_array = vec(array)
    rotations = rotation[]
    if degrees
        flat_array = flat_array/(2*pi)
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
        flat_array = flat_array/(2*pi)
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
function from_basis(array, orthonormal = true, reciprocal = false)
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
        flat_array = flat_array/(2*pi)
    end
    if (length(flat_array) == 9) #als er maar een rotatie is, geen array teruggeven
        return om2qu(array)
    end
    for i in 1:9:length(flat_array)
        #zou deze reshape veel tijd in beslag nemen? anders gewoon een om2qu maken die op een vector werkt?
        #of zou ge kunnen flattenen naar een 2dimensionale array ipv naar een vector
        q = om2qu(flat_array[i:i+8])
        push!(rotations, q)
    end
    return reshape(rotations, sizes[3:length(sizes)])
end

function from_RodriguesFrank(array)
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
function apply(rotation ::rotation, vector ::Array{Number, 1})
    qvector = rotation(vector[1], vector[2], vector[3], 0)
    result = multiply(multiply(rotation, qvector), inv(rotation))
    return [result.i, result.j, result.k]
end

#Rotate a matrix
function apply(rotation ::rotation, vector ::Array{Number, 2})
    R = as_matrix(rotation)
    return R*vector
end

#explanation of ij,j notation:
#v'[1] = R[1,1]*v[1] + R[1,2]*v[2] + R[1,3]*v[3]

#Rotate a fourth order tensor TODO
function apply(rotation ::rotation, vector ::Array{Number, 4})
    R = as_matrix(rotation)
end

function Base.:*(r ::rotation, nr ::Number) ::rotation
    return rotation(r.angle * nr, r.i * nr, r.j * nr, r.k * nr)
end

function Base.:*(nr ::Number, r ::rotation) ::rotation
    return r*nr
end

#niet commutatief
function Base.:*(f ::rotation, s ::rotation) ::rotation
    angle = f.angle*s.angle - (f.i * s.i + f.j * s.j + f.k * s.k)
    i = (f.j * s.k - f.k * s.j) + f.angle * s.i + s.angle * f.i
    j = (f.k * s.i - f.i * s.k) + f.angle * s.j + s.angle * f.j
    k = (f.i * s.j - f.j * s.i) + f.angle * s.k + s.angle * f.k
    return rotation(angle, i, j, k)
end

function Base.:-(rot ::rotation)
    return rotation(rot.angle, -rot.i, -rot.j, -rot.k)
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
    total = 1
    for i in sizes
        total *= i
    end
    if (n != total)
        print("invalid dimensions")
    end
    rotations = rotation[]
    for i = 1:n
        u1 = rand(type)
        u2 = rand(type)
        u3 = rand(type)
        h = rotation(sqrt(1-u1)*sin(pi*u2*2), sqrt(1-u1)*cos(pi*u2*2), sqrt(u1)*sin(pi*u3*2), sqrt(u1)*cos(pi*u3*2))
        push!(rotations, copy(h))
    end
    return reshape(rotations, sizes)
end

function as_axis_angle(rotation ::rotation)
    return qu2ax(rotation)
end

function as_axis_angle(rotations ::Array{rotation}, degrees = false, pair = false)
    sizes = size(rotations)
    flat_array = vec(rotations)
    rotationmatrices = Float32[]
    for i in 1:length(flat_array)
        ax = qu2ax(flat_array[i])
        rotationmatrices = vcat(rotationmatrices, ax)
    end
    return reshape(rotationmatrices, ((4)..., sizes...))
end

function as_matrix(rotation ::rotation) #nodig als er maar een element is
    return reshape(qu2om(rotation), (3,3))
end

function as_matrix(rotations ::Array{rotation})
    sizes = size(rotations)
    flat_array = vec(rotations)
    rotationmatrices = Float32[]
    for i in 1:length(flat_array)
        om = qu2om(flat_array[i])
        rotationmatrices = vcat(rotationmatrices, om)
    end
    return reshape(rotationmatrices, ((3,3)..., sizes...))
end

function as_euler_angle(rotation ::rotation)
    return qu2eu(rotation)
end

function as_euler_angle(rotations ::Array{rotation})
    sizes = size(rotations)
    flat_array = vec(rotations)
    rotationmatrices = Float32[]
    for i in 1:length(flat_array)
        eu = qu2eu(flat_array[i])
        rotationmatrices = vcat(rotationmatrices, eu)
    end
    return reshape(rotationmatrices, ((3)..., sizes...))
end

function as_homochoric(rotations ::Array{rotation})
    sizes = size(rotations)
    flat_array = vec(rotations)
    rotationmatrices = Float32[]
    for i in 1:length(flat_array)
        ho = qu2ho(flat_array[i])
        rotationmatrices = vcat(rotationmatrices, ho)
    end
    return reshape(rotationmatrices, ((3)..., sizes...))
end

function as_RodriguezFrank(rotations ::Array{rotation})
    sizes = size(rotations)
    flat_array = vec(rotations)
    rotationmatrices = Float32[]
    for i in 1:length(flat_array)
        ro = qu2ro(flat_array[i])
        rotationmatrices = vcat(rotationmatrices, ro)
    end
    return reshape(rotationmatrices, ((4)..., sizes...))
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
    if (alpha > pi)
        alpha = 2*pi - alpha
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
        PHI = pi/2*(1 - a_33)
        return eulerAngle(phi1, PHI, 0)
    end
end

function om2ax(rotation ::rotationMatrix)
    omega = acos((tr(rotation.matrix) - 1)/2)
    #inspiratie voor de code: https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/index.htm
    #dit algo komt niet uit de paper, die onduidelijk was, opletten voor omega = 0 of omega = pi
    #omega = 0 -> x y en z zijn arbitrair
    #omega = pi -> x y en z maken uit dus moeten berekend worden
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
