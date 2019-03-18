import Foundation

// MARK: Perlin Noise

var PERLIN_N:Int = (1 << 8)

func permuteAxis(_ axis: inout [Int], _ i: Int) {
    let tar = Int(drand48f()*Float(i + 1))
    let tmp = axis[i]
    axis[i] = axis[tar]
    axis[tar] = tmp
}

// NOTE: (Kapsy) This is an incomplete implementation!
class Perlin {
    var permX: [Int]
    var permY: [Int]
    var permZ: [Int]

    var randFloat: [Float]

    init() {
        let N = PERLIN_N

        self.permX = [Int](repeating: 0, count: PERLIN_N)
        self.permY = [Int](repeating: 0, count: PERLIN_N)
        self.permZ = [Int](repeating: 0, count: PERLIN_N)

        self.randFloat = [Float](repeating: 0.0, count: PERLIN_N)

        for i in 0..<N {
            self.randFloat[i] = drand48f()
            self.permX[i] = i
            self.permY[i] = i
            self.permZ[i] = i
        }

        for i in (0..<N).reversed() {
            permuteAxis(&self.permX, i)
            permuteAxis(&self.permY, i)
            permuteAxis(&self.permZ, i)
        }
    }
}

func getNoise(_ per: Perlin, _ p0: V3) -> Float {

    let p1 = p0*20.0

    let u = Float(p1.x - floor(p1.x))
    let v = Float(p1.y - floor(p1.y))
    let w = Float(p1.z - floor(p1.z))

    let i = Int(floor(p1.x))
    let j = Int(floor(p1.y))
    let k = Int(floor(p1.z))

    // Crazy Swift multi dimensional array syntax.
    var c:[[[Float]]] = Array(repeating:
                        Array(repeating:
                        Array(repeating: 0.0,
                            count: 2),
                            count: 2),
                            count: 2)

    for di in 0..<2 {
    for dj in 0..<2 {
    for dk in 0..<2 {
        c[di][dj][dk] =
            per.randFloat[per.permX[(i+di) & (PERLIN_N - 1)] ^
                          per.permY[(j+dj) & (PERLIN_N - 1)] ^
                          per.permZ[(k+dk) & (PERLIN_N - 1)]]
    }}}

    var accum = Float(0.0)

    for i in 0..<2 {
    for j in 0..<2 {
    for k in 0..<2 {
        let I = (Float(i)*u + (1.0 - Float(i))*(1.0 - u))
        let J = (Float(j)*v + (1.0 - Float(j))*(1.0 - v))
        let K = (Float(k)*w + (1.0 - Float(k))*(1.0 - w))

        accum += I*J*K*Float(c[i][j][k])
    }}}

    assert(accum <= 1.0)

    return accum
}

// MARK: Texture

enum TextureType {

    case plain
    case checker
    case perlin
}

class Texture {

    var type: TextureType = .plain
    var albedo = V3(0)
    var perlin: Perlin?
}

// MARK: Material

enum MaterialType {

    case lambertian
    case metal
    case dielectric
}

class Material {

    var type: MaterialType
    var texture: Texture
    var fuzz = Float(0) // For metal only, leaving in same class for now.
    var refIndex = Float(1.1) // For dielectrics only, leaving here for now.

    init(type: MaterialType, texture: Texture) {
        self.type = type
        self.texture = texture
    }
}

func getAlbedo(_ texture: Texture, _ u: Float, _ v: Float, _ p: V3) -> V3 {

    var res = V3(0)

    switch texture.type {

        case .checker:
            let selector = Float(sin(10.0*p.x)*sin(10.0*p.z))
            if selector > 0.0 {
                res = V3(0,0,0)
            } else {
                res = V3(1.0,1.0,1.0)
            }

        case .plain:
            res = texture.albedo

        case .perlin:
            if let perlin = texture.perlin {
                res = V3(1.0)*getNoise(perlin, p)
            }
    }

    return res
}

func schlick(_ cos: Float, _ refIndex: Float) -> Float {
    var r0 = (1.0 - refIndex)/(1.0 + refIndex);
    r0 = r0*r0;
    r0 = r0 + (1.0 - r0)*pow((1.0 - cos), 5.0);

    return r0
}

func reflect(_ v: V3, _ N: V3) -> V3 {
    return (v - 2*dot(v, N)*N)
}

// MARK: Sphere

class Sphere
{
    var center: V3
    var rad: Float
    var material: Material

    init(center: V3, rad: Float, material: Material) {
        self.center = center
        self.rad = rad
        self.material = material
    }
}

var globalSpheres: [Sphere] = []

enum PrimativeType {
    case sphere
    case triangle
}

// MARK: Camera

struct Camera {

    var origin = V3(0)
    var lowerLeft = V3(0)
    var horiz = V3(0)
    var vert = V3(0)

    var w = V3(0)
    var u = V3(0)
    var v = V3(0)

    var lensRad = Float(0)
}

// MARK: Ray

class HitRecord {
    var dist = Float(0)
    var primRef:Sphere? = nil
    var primType:PrimativeType = .sphere
}

struct Ray {

    var A: V3
    var B: V3

    // NOTE: (Kapsy) For C union like behavior.
    var orig: V3 { get { return A } }
    var dir: V3 { get { return B } }

    init(_ A: V3, _ B: V3) {
        self.A = A
        self.B = B
    }

    init() {
        self.A = V3(0)
        self.B = V3(0)
    }
}

func getRayTemp(_ c: inout Camera, _ s: Float, _ t: Float) -> Ray {

    let res = Ray(c.origin, c.lowerLeft + s*c.horiz + t*c.vert - c.origin)
    return res
}

func getRay(_ c: inout Camera, _ s: Float, _ t: Float) -> Ray {

    // NOTE: (Kapsy) Random in unit disk.
    var rand = V3(0)
    repeat  {
        rand = 2.0*V3(drand48f(), drand48f(), 0) - V3(1,1,0)
    } while dot(rand, rand) >= 1.0

    let rd = c.lensRad*rand;
    let offset = c.u*rd.x + c.v*rd.y

    let res = Ray(c.origin + offset, c.lowerLeft + s*c.horiz + t*c.vert - c.origin - offset)
    return res
}

func pointAt(_ ray: inout Ray, _ t: Float) -> V3 {
    let res = ray.A + t*ray.B
    return res
}

func traverseSpheres(_ ray: inout Ray, _ hit: HitRecord) {

    let tnear = Float(0.001)
    var tfar = Float.greatestFiniteMagnitude

    for sphere in globalSpheres {

        let rad = sphere.rad
        let center = sphere.center
        let oc = ray.orig - center

        let a = dot(ray.dir, ray.dir)
        let b = dot(oc, ray.dir)
        let c = dot(oc, oc) - rad*rad
        let discriminant = b*b - a*c;

        if discriminant > 0.0 {
            var t = (-b - sqrt(discriminant))/a
            if (tnear < t && t < tfar)
            {
                tfar = t

                hit.dist = t;
                hit.primRef = sphere
                hit.primType = .sphere;
            }

            t = (-b + sqrt(discriminant))/a
            if (tnear < t && t < tfar)
            {
                tfar = t

                hit.dist = t;
                hit.primRef = sphere
                hit.primType = .sphere;
            }
        }
    }
}

var MAX_DEPTH = Int(10)

func getColorForRay(_ ray: inout Ray, _ depth: Int) -> V3 {

    var res = V3()

    let hit = HitRecord()
    hit.dist = Float.greatestFiniteMagnitude

    traverseSpheres(&ray, hit)

    if hit.dist < Float.greatestFiniteMagnitude {

        var p = V3(0)
        var N = V3(0)
        var mat: Material? = nil;

        switch hit.primType {

        case .sphere:

            if let sphere = hit.primRef {

                let rad = sphere.rad
                let center = sphere.center

                p = pointAt(&ray, hit.dist)
                N = (p - center)/rad
                mat = sphere.material
            }

        case .triangle:
            break
        }

        if let _mat = mat {

            switch _mat.type {

            case .lambertian:

                let rand = randomInUnitSphere()
                let target = p + N + rand
                let albedo = getAlbedo(_mat.texture, 0, 0, p)
                var scattered = Ray(p, target - p)

                if depth < MAX_DEPTH {
                    res = albedo*getColorForRay(&scattered, depth+1)
                } else {
                    res = V3(0)
                }

            case .metal:

                let v = unit(ray.dir)
                let reflected = v - 2*dot(v, N)*N
                let bias = N*1e-4

                var scattered = Ray(p + bias, reflected + _mat.fuzz*randomInUnitSphere())

                let albedo = getAlbedo(_mat.texture, 0, 0, p)

                // NOTE: (Kapsy) Checking if < 90 deg.
                let result = (dot(scattered.dir, N) > 0.0)
                if (depth < MAX_DEPTH && result) {
                    res = albedo*getColorForRay(&scattered, depth+1)
                } else {
                    res = V3(0.0)
                }

            case .dielectric:

                    var scattered = Ray()

                    var outwardNormal = V3(0)
                    var niOverNt = Float(0)
                    let reflected = reflect(ray.dir, N)

                    var reflectProb = Float(0)
                    var cos = Float(0)

                    if dot(ray.dir, N) > 0.0 {
                        outwardNormal = -N
                        niOverNt = _mat.refIndex
                        cos = _mat.refIndex*dot(ray.dir, N)/length(ray.dir)
                    }
                    else
                    {
                        outwardNormal = N
                        niOverNt = 1.0/_mat.refIndex
                        cos = -dot(ray.dir, N)/length(ray.dir)
                    }

                    var refracted = V3(0.0)

                    let bias = outwardNormal*1e-2
                    p = p - bias

                    let uv = unit(ray.dir)
                    let dt = dot(uv, outwardNormal)
                    let discriminant = 1.0 - niOverNt*niOverNt*(1.0 - dt*dt)

                    // NOTE: (Kapsy) Approximate reflection/refraction probability.
                    if discriminant > 0.0 {
                        refracted = niOverNt*(uv - outwardNormal*dt) - outwardNormal*sqrt(discriminant)
                        reflectProb = schlick(cos, _mat.refIndex)
                    } else {
                        scattered = Ray(p, reflected)
                        reflectProb = 1.0
                    }

                    if drand48f() < reflectProb {
                        scattered = Ray(p, reflected)
                    } else {
                        scattered = Ray(p, refracted)
                    }

                    if depth < MAX_DEPTH {
                        res = getColorForRay(&scattered, depth+1)
                    } else {
                        res = V3 (0.0)
                    }
            }
        }

    } else {

        // NOTE: (Kapsy) Draw our psuedo sky background.
        let rdir = ray.dir

        let t = (unit(rdir).y + 1.0)*0.5
        let cola = V3(1.0)
        let colb = (1.0/255.0)*V3(255.0, 128.0, 0.0)

        res = (1.0 - t)*cola + t*colb
    }

    return res
}

func main() {

    // MARK: Init spheres

    let perlinTexture = Texture()
    perlinTexture.albedo = V3(1,1,1)
    perlinTexture.perlin = Perlin()
    perlinTexture.type = .perlin
    let sphere0Mat = Material(type: .lambertian, texture: perlinTexture)
    let sphere0 = Sphere(center: V3(0, 0.32, 0), rad: 0.34, material: sphere0Mat)
    globalSpheres.append(sphere0)

    let glassTexture = Texture()
    glassTexture.albedo = V3(1)
    let sphere1Mat = Material(type: .dielectric, texture: glassTexture)
    let sphere1 = Sphere(center: V3(0.53, 0.3, -0.33), rad: -0.23, material: sphere1Mat)
    globalSpheres.append(sphere1)

    let whiteTexture = Texture()
    whiteTexture.albedo = V3(1,0.97,0.97)
    let sphere2Mat = Material(type: .metal, texture: whiteTexture)
    sphere2Mat.fuzz = 0.24
    let sphere2 = Sphere(center: V3(-0.7, 0.3, 0), rad: 0.24, material: sphere2Mat)
    globalSpheres.append(sphere2)

    let groundTexture = Texture()
    groundTexture.albedo = V3(0.2,0.5,0.3)
    groundTexture.type = .checker
    let sphere3Mat = Material(type: .lambertian, texture: groundTexture)
    sphere3Mat.refIndex = 1.3
    let sphere3 = Sphere(center: V3(0, -99.99, 0), rad: 100.0, material: sphere3Mat)
    globalSpheres.append(sphere3)

    let greenTexture = Texture()
    greenTexture.albedo = V3(0,1.3,0)
    let sphere4Mat = Material(type: .lambertian, texture: greenTexture)
    let sphere4 = Sphere(center: V3(0.0, 0.3, 0.5), rad: 0.13, material: sphere4Mat)
    globalSpheres.append(sphere4)

    let redTexture = Texture()
    redTexture.albedo = V3(2,0.3,0.3)
    let sphere5Mat = Material(type: .lambertian, texture: redTexture)
    let sphere5 = Sphere(center: V3(0.1, 0.3, -0.6), rad: 0.16, material: sphere5Mat)
    globalSpheres.append(sphere5)

    let purpleTexture = Texture()
    purpleTexture.albedo = V3(1,0,1)
    let sphere6Mat = Material(type: .metal, texture: purpleTexture)
    sphere6Mat.fuzz = 0.2
    let sphere6 = Sphere(center: V3(0.68, 0.33, 0.79), rad: 0.33, material: sphere6Mat)
    globalSpheres.append(sphere6)

    let blueTexture = Texture()
    blueTexture.albedo = V3(0.2,0.2,3)
    let sphere7Mat = Material(type: .lambertian, texture: blueTexture)
    let sphere7 = Sphere(center: V3(-0.5, 0.3, -0.9), rad: 0.13, material: sphere7Mat)
    globalSpheres.append(sphere7)

    let purple2Texture = Texture()
    purple2Texture.albedo = V3(1,1,1)
    let sphere8Mat = Material(type: .dielectric, texture: purple2Texture)
    let sphere8 = Sphere(center: V3(-0.6, 0.24, 0.6), rad: 0.18, material: sphere8Mat)
    globalSpheres.append(sphere8)

    let metalTexture = Texture()
    metalTexture.albedo = V3(0,1,1)
    let sphere9Mat = Material(type: .metal, texture: metalTexture)
    sphere9Mat.fuzz = 0.3
    let sphere9 = Sphere(center: V3(0.5, 0.3, -0.9), rad: 0.10, material: sphere9Mat)
    globalSpheres.append(sphere9)

    let frameRate = Float(25)
    // NOTE: (Kapsy) Uncomment for animation preview renders.
    //let frameRate = Float(0.5)

    let durationSeconds = Float(12)
    let frameCount = frameRate*durationSeconds

    // NOTE: (Kapsy) Camera rotation step
    let omega = (2*Float.pi)/frameCount
    let k = V3 (0,1,0)
    var ellipsephase = Float(0)

    // NOTE: (Kapsy) Image plane dimensions
    let nx = Int(600)
    let ny = Int(300)

    //// let nx = Int(200)
    //// let ny = Int(100)

    // NOTE: (Kapsy) Primary rays per pixel
    let ns = Int(30)

    var lookFrom = V3(0.001,0.39,-1.0)
    let lookAt = V3(0.0, 0.3, 0.0)

    for f in 0..<Int(frameCount) {

        var data = Data()

        let header = "P3\n\(nx) \(ny)\n255\n".data(using: .ascii)!
        data.append(header)

        // NOTE: (Kapsy) Camera setup stuff.

        var cam = Camera()

        var lookFromRes = lookFrom
        lookFromRes = lookFromRes*((-cos(ellipsephase) + 1.0)*0.07 + 1.3);

        let vup = V3(0.18, 1, 0)
        let vfov = Float(60)
        let aspect = Float(nx)/Float(ny)
        let aperture = Float(0.09)
        let focusDist = length(lookFromRes - lookAt)

        cam.lensRad = aperture/2.0;

        let theta = vfov*Float.pi/180
        let halfHeight = tan(theta/2)
        let halfWidth = Float(aspect*halfHeight)

        cam.origin = lookFromRes
        cam.w = unit(lookFromRes - lookAt)
        cam.u = unit(cross(vup, cam.w))
        cam.v = cross(cam.w, cam.u)

        cam.lowerLeft = cam.origin - halfWidth*focusDist*cam.u - halfHeight*focusDist*cam.v - focusDist*cam.w
        cam.horiz = 2*halfWidth*focusDist*cam.u
        cam.vert = 2*halfHeight*focusDist*cam.v

        for j in (0..<ny).reversed() {

            for i in 0..<nx {

                var col = V3(0)

                for _ in 0..<ns {

                    let u = (Float(i) + drand48f())/Float(nx)
                    let v = (Float(j) + drand48f())/Float(ny)

                    var r = getRay(&cam, u, v)

                    col += getColorForRay(&r, 0)

                    assertNaN(col.r)
                    assertNaN(col.g)
                    assertNaN(col.b)
                }

                col /= Float(ns)

                col.r = clamp01(col.r)
                col.g = clamp01(col.g)
                col.b = clamp01(col.b)

                let ir = Int(255.99*col.r)
                let ig = Int(255.99*col.g)
                let ib = Int(255.99*col.b)

                assert((ir >= 0) && (ir <= 255))
                assert((ig >= 0) && (ig <= 255))
                assert((ib >= 0) && (ib <= 255))

                let pixel = "\(ir) \(ig) \(ib)\n".data(using: .ascii)!
                data.append(pixel)
            }
        }

        let outPath = "temp/out_\(String(format: "%03d", f + 1)).ppm"

        FileManager.default.createFile(atPath: outPath, contents: data)

        // NOTE: (Kapsy) Rodrigues Rotation formula
        var v = lookFrom
        v = v*cos(omega) + cross(k, v)*sin(omega) + k*dot(k, v)*(1.0 - cos(omega))
        lookFrom = v

        ellipsephase += omega
    }
}

main()
