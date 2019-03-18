import Foundation

enum TextureType {

    case plain
    case checker
    case perlin
}

class Texture {

    var type: TextureType = .plain
    var albedo = V3(0)

}

enum MaterialType {

    case lambertian
}

class Material {

    var type: MaterialType
    var texture: Texture

    init(type: MaterialType, texture: Texture) {
        self.type = type
        self.texture = texture
    }
}

func getAlbedo(_ texture: Texture, _ u: Float, _ v: Float, _ p: V3) -> V3 {

    var res = V3(0)

    switch texture.type {

        case .checker:
            let selector = Float(sin(10.0*p.x)*sin(10.0*p.y)*sin(10.0*p.z))
            if selector > 0.0 {
                res = V3(0,0,0)
            } else {
                res = V3(1.0,1.0,1.0)
            }

        case .plain:
            res = texture.albedo

        case .perlin:
            break
            //res = V3(1.0)*GetNoise(tex->perlin, p);
    }

    return res
}

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

    var hit = HitRecord()
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
                let attenuation = getAlbedo(_mat.texture, 0, 0, p)
                var scattered = Ray(p, target - p)

                if depth < MAX_DEPTH {
                    res = attenuation*getColorForRay(&scattered, depth+1)
                } else {
                    res = V3(0)
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

    //// var s0 = sphereCount
    //// spheres[s0].center = V3(0, -99.99, 0)
    //// spheres[s0].rad = 100.f
    //// spheres[s0].mat = Mat(type: MatLambertian, t3)
    //// sphereCount += 1

    //// var s3 = sphereCount
    //// spheres[s3].center = V3(0, 0, -1)
    //// spheres[s3].rad = 0.5
    //// sphereCount += 1

    //let perlinTexture = PerlinTexture()

    let perlinTexture = Texture()
    perlinTexture.albedo = V3(1,1,0)
    let sphere0Mat = Material(type: .lambertian, texture: perlinTexture)
    let sphere0 = Sphere(center: V3(0, 0.3, 0), rad: 0.3, material: sphere0Mat)
    globalSpheres.append(sphere0)

    let sphere1Mat = Material(type: .lambertian, texture: perlinTexture)
    let sphere1 = Sphere(center: V3(0.5, 0.3, -0.3), rad: 0.2, material: sphere1Mat)
    globalSpheres.append(sphere1)

    let sphere2Mat = Material(type: .lambertian, texture: perlinTexture)
    let sphere2 = Sphere(center: V3(-0.7, 0.3, 0), rad: 0.2, material: sphere2Mat)
    globalSpheres.append(sphere2)

    let groundTexture = Texture()
    groundTexture.albedo = V3(0.2,0.5,0.3)
    let sphere3Mat = Material(type: .lambertian, texture: groundTexture)
    let sphere3 = Sphere(center: V3(0, -99.99, 0), rad: 100.0, material: sphere3Mat)
    globalSpheres.append(sphere3)

    //// sphere *s3 = &spheres[spherecount++];
    //// s3->center = V3 (0.6f, 0.8f, -0.3f);
    //// s3->rad = 0.2f;
    //// s3->mat = (mat_t) { MAT_METAL, t3 };

    //// sphere *s4 = &spheres[spherecount++];
    //// s4->center = V3 (0.7f, 0.32f, -0.4f);
    //// s4->rad = 0.23f;
    //// s4->mat = (mat_t) { MAT_DIELECTRIC, t3 };
    //// s4->mat.refindex = 1.1f;

    //let frameRate = Float(25)
    let frameRate = Float(1)
    let durationS = Float(12)
    let frameCount = frameRate*durationS

    let omega = (2*Float.pi)/frameCount
    let k = V3 (0,1,0)

    var ellipsephase = Float(0)

    let nx = Int(200)
    let ny = Int(100)
    let ns = Int(100)

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
        //let T = M44Trans(0.f, 0.4, 0.0);
        //lookFromRes = lookFromRes*T

        let vup = V3(0.18, 1, 0)
        let vfov = Float(60)
        let aspect = Float(nx)/Float(ny)
        let aperture = Float(0) // (0.04)
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


        // Rodrigues Rotation formula
        var v = lookFrom
        v = v*cos(omega) + cross(k, v)*sin(omega) + k*dot(k, v)*(1.0 - cos(omega))
        lookFrom = v

        ellipsephase += omega
    }
}

main()
