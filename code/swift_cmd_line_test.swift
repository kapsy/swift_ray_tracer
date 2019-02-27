

import Foundation

// MARK: Number functions

func assertNaN(_ value: Float) {
    assert(!value.isNaN)
}

func drand48f() -> Float {
    return Float(drand48())
}

// MARK: V3

struct V3
{
    var x = Float(0)
    var y = Float(0)
    var z = Float(0)

    // NOTE: (Kapsy) For C union like behavior.
    var r: Float { get { return x } }
    var g: Float { get { return y } }
    var b: Float { get { return z } }

    init() {}

    init(_ val: Float) {
        self.x = val
        self.y = val
        self.z = val
    }

    init(_ x: Float, _ y: Float, _ z: Float) {
        self.x = x
        self.y = y
        self.z = z
    }
}


func +(a: V3, b: V3) -> V3 {
    let res = V3(a.x + b.x, a.y + b.y, a.z + b.z)
    return res
}

func +=(a: inout V3, b: V3) {
    a = a + b
}

func -(a: V3, b: V3) -> V3 {
    let res = V3(a.x - b.x, a.y - b.y, a.z - b.z)
    return res
}

func -=(a: inout V3, b: V3) {
    a = a - b
}

func *(a: V3, b: V3) -> V3 {
    let res = V3(a.x*b.x, a.y*b.y, a.z*b.z)
    return res
}

func *=(a: inout V3, b: V3) {
    a = a * b
}

func /(a: V3, b: Float) -> V3 {
    let res = V3(a.x/b, a.y/b, a.z/b)
    return res
}

func /=(a: inout V3, b: Float) {
    a = a/b
}

func *(a: Float, b: V3) -> V3 {
    let res = V3(a*b.x, a*b.y, a*b.z)
    return res
}

func *(a: V3, b: Float) -> V3 {
    let res = b*a
    return res
}

func dot(_ a: V3, _ b: V3) -> Float {
    let res = a.x*b.x + a.y*b.y + a.z*b.z
    return res
}

func cross(_ a: V3, _ b: V3) -> V3 {
    let res = V3(a.y*b.z - a.z*b.y,
                 a.z*b.x - a.x*b.z,
                 a.x*b.y - a.y*b.x);
    return res
}


func length(_ a: V3) -> Float {
    let res = sqrt(dot(a, a))
    return res
}

func unit(_ a: V3) -> V3 {
    let res = a/length(a)
    return res
}

// MARK: Static arrays


struct Sphere {
    var center = V3(0)
    var rad = Float(0)
    //material mat;
}

var MAX_SPHERES = Int(1 << 6)
var spheres = [Sphere](repeating: Sphere(), count: MAX_SPHERES)
var sphereCount = 0

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

// NOTE: (Kapsy) Ray-sphere intersection:
//
// x^2 + y^2 + z^2 = R^2
//
// A.A = R^2
// (A - C).(A - C) = R^2
//
// since A = A + tB
//
// (A + tB - C).(A + tB - C) = R^2
//
// v =def A - C
//
// starting another way:
// in vector notation:
//
// (x - c)^2 = r^2
// sub x = A + tB
// (A + tB - C)^2 = r^2
//
// def v = A - C
//
// |v + tB|^2 = r^2
//
// (v + tB)*(v + tB)
//
// v^2 + vtB + vtB + t^2B^2
// v^2 + 2v.tB + t^2B^2
//
// solve for:
//
// t^2(B^2) + t(2v.B) + (v^2 - r^2) = 0
//

//    a.b = |a||b|cos(?)
//
//
//    a.a = |a||a|
//    a.a = |a|^2

// |v + tB|^2 = r^2
// (v + tB).(v + tB)
// v.(v + tB) + tB.(v + tB)
// v.v + v.tB + tB.v + tB.tB
// v.v + 2v.tB + t^2B.B
//
//
// t*t*B.B + 2*t*B.v + v.v - R*R = 0
//
// |v|^2 + 2v.tB + t^2B^2

func getColor(_ ray: inout Ray, _ depth: Int) -> V3 {

    var res = V3()
    let tmin = Float(0.001)
    var tmax = Float.greatestFiniteMagnitude

    var collision = false
    var N = V3()
    var p = V3()

    for i in 0..<sphereCount {

        let rad = spheres[i].rad
        let center = spheres[i].center

        let oc = ray.orig - center

        let a = dot(ray.dir, ray.dir)
        let b = dot(oc, ray.dir)
        let c = dot(oc, oc) - rad*rad
        let discriminant = b*b - a*c;

        if discriminant > 0.0 {
            var t = (-b - sqrt(discriminant))/a
            if (tmin < t && t < tmax)
            {
                tmax = t
                p = pointAt(&ray, t)
                N = (p - center)/rad
                // matIndex = i
                collision = true
            }

            t = (-b + sqrt(discriminant))/a
            if (tmin < t && t < tmax)
            {
                tmax = t
                p = pointAt(&ray, t)
                N = (p - center)/rad
                // matIndex = i
                collision = true
            }

        }
    }

    if (collision) {
        res = 0.5*V3 (N.x+1, N.y+1, N.z+1);
    } else {
        res = V3(0,0,1)
    }

    return res
}


func main() {

    var s0 = sphereCount
    spheres[s0].center = V3(0, 0, -1)
    spheres[s0].rad = 0.5
    sphereCount += 1


    let frameRate = Float(25)
    let durationS = Float(5)
    let frameCount = 1//frameRate*durationS

    let testStart = Float(-2)
    let testEnd = Float(2)
    let testDelta = (testEnd - testStart)/Float(frameCount)

    let lookFromTest = testStart

    for f in 0..<frameCount {

        var data = Data()

        let nx = Int(600)
        let ny = Int(300)
        let ns = Int(100)

        // NOTE: (Kapsy) Camera setup stuff.

        let lookFrom = V3(0,1,2)
        let lookAt = V3(0,0.25,-0.2)
        let vup = V3(0,1,0)
        let vfov = Float(50)
        let aspect = Float(nx)/Float(ny)
        let aperture = Float(0)//(0.04)
        let focusDist = length(lookFrom - lookAt)

        let theta = vfov*Float.pi/180
        let halfHeight = tan(theta/2)
        let halfWidth = Float(aspect*halfHeight)

        var cam = Camera()

        cam.origin = lookFrom
        cam.w = unit(lookFrom - lookAt)
        cam.u = unit(cross(vup, cam.w))
        cam.v = cross(cam.w, cam.u)

        cam.lowerLeft = cam.origin - halfWidth*focusDist*cam.u - halfHeight*focusDist*cam.v - focusDist*cam.w
        cam.horiz = 2*halfWidth*focusDist*cam.u
        cam.vert = 2*halfHeight*focusDist*cam.v


        var lensRad: Float

        let header = "P3\n\(nx) \(ny)\n255\n".data(using: .ascii)!
        data.append(header)

        for j in (0..<ny).reversed() {

            for i in 0..<nx {

                var col = V3(0)

                for s in 0..<ns {

                    let u = (Float(i) + drand48f())/Float(nx)
                    let v = (Float(j) + drand48f())/Float(ny)

                    var r = getRay(&cam, u, v)

                    var depth = 0

                    col += getColor(&r, depth)

                    assertNaN(col.r)
                    assertNaN(col.g)
                    assertNaN(col.b)
                }

                col /= Float(ns)

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
    }
}

main()
